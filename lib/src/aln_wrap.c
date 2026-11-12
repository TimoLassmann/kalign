#include <unistd.h>
#include "tldevel.h"
#include "tlmisc.h"
#include <strings.h>
#include "esl_stopwatch.h"
#include "task.h"
#include "msa_struct.h"

#include "kalign_simd.h"
#include "msa_op.h"
#include "msa_alloc.h"
#include "msa_check.h"
#include "msa_sort.h"
#include "msa_io.h"

#include "alphabet.h"
#include "bisectingKmeans.h"
#include "sl_tree.h"
#include "guide_dist.h"

#include "aln_param.h"
#include "aln_run.h"
#include "aln_refine.h"
#include "aln_apair_dist.h"
#include "anchor_consistency.h"
#include "kalign/kalign.h"
#include "kalign/kalign_config.h"
#include "ensemble.h"
#include "threadpool.h"

#define ALN_WRAP_IMPORT
#include "aln_wrap.h"

/* Resolve KALIGN_MATRIX_AUTO into a concrete matrix constant.
   For protein: use length-ratio heuristic (PFASUM43 vs PFASUM60).
   For DNA/RNA: pick the standard matrix for that biotype. */
static int resolve_matrix_auto(struct msa *msa, int *type)
{
        int i;
        int min_len, max_len;
        float len_ratio;

        if(*type != KALIGN_MATRIX_AUTO){
                return OK;
        }

        if(msa->biotype == ALN_BIOTYPE_DNA){
                /* DNA biotype covers both DNA and RNA sequences.
                   Default to RNA matrix (broader scoring range). */
                *type = KALIGN_MATRIX_RNA;
                return OK;
        }

        /* Protein: use sequence length ratio (max/min) to select matrix.
           Similar-length sequences (ratio < 1.5) prefer PFASUM43;
           high length variation (insertions/extensions) prefers PFASUM60. */
        min_len = msa->sequences[0]->len;
        max_len = msa->sequences[0]->len;
        for(i = 1; i < msa->numseq; i++){
                int l = msa->sequences[i]->len;
                if(l < min_len) min_len = l;
                if(l > max_len) max_len = l;
        }
        len_ratio = (min_len > 0) ? (float)max_len / (float)min_len : 1.0f;

        if(len_ratio < 1.5f){
                *type = KALIGN_MATRIX_PFASUM43;
        }else{
                *type = KALIGN_MATRIX_PFASUM60;
        }
        if(!msa->quiet){
                LOG_MSG("Auto matrix: len_ratio=%.2f -> %s",
                        len_ratio,
                        *type == KALIGN_MATRIX_PFASUM60 ? "PFASUM60" : "PFASUM43");
        }
        return OK;
}

static int compute_tree_weights(struct msa* msa, struct aln_tasks* tasks)
{
        float* nw = NULL;
        int i;

        MMALLOC(nw, sizeof(float) * msa->num_profiles);
        for(i = 0; i < msa->num_profiles; i++){
                nw[i] = 0.0f;
        }

        /* Root gets total weight = numseq */
        nw[tasks->list[tasks->n_tasks - 1]->c] = (float)msa->numseq;

        /* Walk root → leaves: at each split, distribute parent's weight
           to children in proportion to the OTHER child's size */
        for(i = tasks->n_tasks - 1; i >= 0; i--){
                int a = tasks->list[i]->a;
                int b = tasks->list[i]->b;
                int c = tasks->list[i]->c;
                float total = (float)(msa->nsip[a] + msa->nsip[b]);
                nw[a] = nw[c] * (float)msa->nsip[b] / total;
                nw[b] = nw[c] * (float)msa->nsip[a] / total;
        }

        /* Copy leaf weights to seq_weights */
        if(msa->seq_weights){
                MFREE(msa->seq_weights);
        }
        MMALLOC(msa->seq_weights, sizeof(float) * msa->numseq);
        for(i = 0; i < msa->numseq; i++){
                msa->seq_weights[i] = nw[i];
        }

        MFREE(nw);
        return OK;
ERROR:
        if(nw) MFREE(nw);
        return FAIL;
}

/* ======================================================================== */
/* kalign_single_run — unified internal entry point for one alignment run.   */
/*                                                                           */
/* Handles the full pipeline: input check, sort, alphabet conversion, tree   */
/* building, alignment, refinement, and optional realign iterations.          */
/*                                                                           */
/* The caller owns msa->pool (NULL for one thread).  This function does NOT */
/* create or destroy the threadpool — it reads msa->pool for parallel work.  */
/* ======================================================================== */

int kalign_single_run(struct msa *msa,
                      const struct kalign_run_config *cfg,
                      int n_threads)
{ KPERF_REGION(KP_RUN);
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;
        struct guide_dist_ctx* gctx = NULL;
        int type;
        int iter;

        ASSERT(cfg->guide_dist == KALIGN_GUIDE_DIST_BPM ||
               cfg->guide_dist == KALIGN_GUIDE_DIST_PAIRALN,
               "Unknown guide_dist %d", cfg->guide_dist);
        ASSERT(cfg->gap_model == 0 || cfg->gap_model == 1, "Unknown gap_model %d", cfg->gap_model);

        /* Guide-tree options of this run */
        RUN(sl_opts_resolve(&msa->tree, cfg));

        /* Input check (also sets ranks) */
        RUN(kalign_essential_input_check(msa, 0));

        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }
        RUN(msa_sort_len_name(msa));

        /* Convert to reduced alphabet for tree building */
        if(msa->biotype == ALN_BIOTYPE_DNA){
                msa->L = ALPHA_defDNA;
                RUN(convert_msa_to_internal(msa, ALPHA_defDNA));
        }else if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                msa->L = ALPHA_redPROTEIN;
                RUN(convert_msa_to_internal(msa, ALPHA_redPROTEIN));
        }else{
                ERROR_MSG("Unable to determine what alphabet to use.");
        }

        RUN(alloc_tasks(&tasks, msa->numseq));

        /* Resolve auto matrix selection (sequence-length heuristic) and
           init alignment parameters.  Done before the guide tree because
           guide_dist = PAIRALN aligns sequence pairs with these parameters. */
        type = cfg->matrix;
        RUN(resolve_matrix_auto(msa, &type));

        RUN(aln_param_init(&ap, msa->biotype, n_threads, type,
                           cfg->gpo, cfg->gpe, cfg->tgpe));
        ap->gap_model = cfg->gap_model;
        ap->pool = msa->pool;

        if(cfg->guide_dist == KALIGN_GUIDE_DIST_PAIRALN){
                RUN(guide_dist_ctx_create(msa, ap, &gctx));
        }

        /* Build guide tree — noisy variant for ensemble diversity */
        if(cfg->tree_seed != 0 && cfg->tree_noise > 0.0f){
                RUN(build_tree_kmeans_noisy(msa, &tasks, cfg->tree_seed, cfg->tree_noise, gctx));
        }else{
                RUN(build_tree_kmeans(msa, &tasks, gctx));
        }
        guide_dist_ctx_free(gctx);
        gctx = NULL;

        /* Convert to full alphabet for alignment */
        if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }

        /* Apply config overrides (config values are concrete, not sentinels) */
        ap->adaptive_budget = cfg->adaptive_budget;
        ap->use_seq_weights = cfg->seq_weights;
        ap->dist_scale = cfg->dist_scale;
        ap->vsm_amax = cfg->vsm_amax;
        ap->vsm_mode = cfg->vsm_mode;
        ap->vsm_kappa = cfg->vsm_kappa;
        ap->vsm_dzero = cfg->vsm_dzero;

        if(ap->use_seq_weights > 0.0f){
                RUN(compute_tree_weights(msa, tasks));
        }

        /* Build anchor consistency table if requested */
        if(cfg->consistency_anchors > 0){
                ap->consistency_anchors = cfg->consistency_anchors;
                ap->consistency_weight = cfg->consistency_weight;
                RUN(anchor_consistency_build(msa, ap, cfg->consistency_anchors,
                                             cfg->consistency_weight,
                                             (struct consistency_table**)&msa->consistency_table));
        }

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Aligning");
        }
        START_TIMER(t1);

        /* Initial alignment */
        if(cfg->refine == KALIGN_REFINE_INLINE){
                RUN(create_msa_tree_inline_refine(msa, ap, tasks, 3));
        }else{
                RUN(create_msa_tree(msa, ap, tasks));
        }
        msa->aligned = ALN_STATUS_ALIGNED;

        /* Iterative realignment: align → N×K anchor distances → bisecting
           k-means tree → dealign → re-encode → re-align.
           Sequences stay aligned until after tree build so the identity
           distance callback can read gap characters from seq->seq. */
        msa->tree.in_realign = 1;
        for(iter = 0; iter < cfg->realign; iter++){
                float** dm = NULL;
                int* realign_anchors = NULL;
                int n_realign_anchors = 0;
                int si;

                RUN(finalise_alignment(msa));

                /* Select K diverse anchors and compute N×K distances
                   (sequences still aligned — gap chars in seq->seq). */
                RUN(pick_anchor_from_alignment(msa, REALIGN_NUM_ANCHORS,
                                               &realign_anchors,
                                               &n_realign_anchors));
                RUN(compute_aln_anchor_dist(msa, realign_anchors,
                                             n_realign_anchors, &dm));
                MFREE(realign_anchors);

                /* Rebuild guide tree via bisecting k-means.
                   The leaf-cluster callback uses identity distances from
                   the aligned sequences (seq->seq still has gap chars). */
                free_tasks(tasks);
                tasks = NULL;
                {
                        /* dm: rows from kalign_aligned_malloc, outer array
                           from MMALLOC; freed here on success and failure. */
                        int fi;
                        int rc = FAIL;
                        if(alloc_tasks(&tasks, msa->numseq) == OK){
                                rc = build_tree_kmeans_from_dm(msa, &tasks, dm,
                                                               n_realign_anchors,
                                                               aln_identity_pair_dist);
                        }
                        for(fi = 0; fi < msa->numseq; fi++){
                                kalign_aligned_free(dm[fi]);
                        }
                        MFREE(dm);
                        dm = NULL;
                        if(rc != OK){
                                ERROR_MSG("Realign guide tree failed");
                        }
                }

                /* NOW dealign: strip gap characters from seq->seq. */
                RUN(dealign_msa(msa));
                for(si = 0; si < msa->numseq; si++){
                        struct msa_seq* seq = msa->sequences[si];
                        int r, w = 0;
                        for(r = 0; seq->seq[r] != '\0'; r++){
                                if(seq->seq[r] != '-'){
                                        seq->seq[w++] = seq->seq[r];
                                }
                        }
                        seq->seq[w] = '\0';
                        seq->len = w;
                }

                /* Re-encode for alignment */
                if(msa->biotype == ALN_BIOTYPE_DNA){
                        RUN(convert_msa_to_internal(msa, ALPHA_defDNA));
                }else if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                        RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
                }

                RUN(set_sip_nsip(msa));

                if(ap->use_seq_weights > 0.0f){
                        RUN(compute_tree_weights(msa, tasks));
                }

                /* Re-align with new tree */
                if(cfg->refine == KALIGN_REFINE_INLINE){
                        RUN(create_msa_tree_inline_refine(msa, ap, tasks, 3));
                }else{
                        RUN(create_msa_tree(msa, ap, tasks));
                }
                msa->aligned = ALN_STATUS_ALIGNED;
        }
        msa->tree.in_realign = 0;

        /* Two-pass refinement after the final alignment pass */
        if(cfg->refine != KALIGN_REFINE_NONE && cfg->refine != KALIGN_REFINE_INLINE){
                RUN(refine_alignment(msa, ap, tasks, cfg->refine));
        }

        /* Free consistency table AFTER refinement */
        if(msa->consistency_table){
                anchor_consistency_free((struct consistency_table*)msa->consistency_table);
                msa->consistency_table = NULL;
        }

        if(!msa->defer_finalise){
                RUN(finalise_alignment(msa));
        }
        RUN(msa_sort_rank(msa));

        STOP_TIMER(t1);
        if(!msa->quiet){
                GET_TIMING(t1);
        }
        DESTROY_TIMER(t1);

        aln_param_free(ap);
        free_tasks(tasks);
        return OK;
ERROR:
        if(msa && msa->consistency_table){
                anchor_consistency_free((struct consistency_table*)msa->consistency_table);
                msa->consistency_table = NULL;
        }
        guide_dist_ctx_free(gctx);
        aln_param_free(ap);
        if(tasks) free_tasks(tasks);
        return FAIL;
}

/* Run a mode on an msa exactly as the command line does: apply the sequence
   type, check the size limit, take the fitted preset for the sequence type and
   run it.  gpo, gpe and tgpe < 0 keep the fitted values, >= 0 override them in
   every run.  Shared by kalign() and kalign_align_mode(). */
static int run_mode_on_msa(struct msa *msa, const char *mode, int seq_type,
                           float gpo, float gpe, float tgpe, int n_threads)
{
        struct kalign_run_config runs[KALIGN_MAX_PRESET_RUNS];
        struct kalign_ensemble_config ens;
        int n_runs = 0;
        int k;
        const char *m = mode ? mode : "fast";

        ASSERT(msa != NULL, "No MSA");

        if(kalign_set_seq_type(msa, seq_type, NULL) != 0){
                ERROR_MSG("Invalid sequence type %d.", seq_type);
        }
        if(msa->biotype == ALN_BIOTYPE_UNDEF){
                RUN(detect_alphabet(msa));
        }
        if(kalign_get_mode_preset(m, msa->biotype, runs, &n_runs, &ens) != 0){
                ERROR_MSG("Unknown mode: '%s'. Use: fast or accurate.", m);
        }
        if(!kalign_mode_size_ok(m, msa->biotype, msa->numseq)){
                ERROR_MSG("Mode '%s' is not feasible above %d sequences (input: %d).",
                          m, KALIGN_ACCURATE_MAX_SEQ, msa->numseq);
        }
        for(k = 0; k < n_runs; k++){
                if(gpo >= 0.0f)  runs[k].gpo  = gpo;
                if(gpe >= 0.0f)  runs[k].gpe  = gpe;
                if(tgpe >= 0.0f) runs[k].tgpe = tgpe;
        }
        RUN(kalign_align_full(msa, runs, n_runs, &ens, n_threads));
        return OK;
ERROR:
        return FAIL;
}

int kalign_align_mode(struct msa *msa, const char *mode, int seq_type, int n_threads)
{
        return run_mode_on_msa(msa, mode, seq_type, -1.0f, -1.0f, -1.0f, n_threads);
}

/* Array entry point: the "fast" mode, as on the command line and in Python.
   type selects the sequence type like --type; the fitted fast setting for that
   type is used, and gpo, gpe and tgpe >= 0 override its gap penalties. */
int kalign(char **seq, int *len, int numseq, int n_threads, int type,
           float gpo, float gpe, float tgpe, char ***aligned, int *out_aln_len)
{
        struct msa *msa = NULL;
        int seq_type;

        switch(type){
        case KALIGN_TYPE_UNDEFINED:
                seq_type = KALIGN_SEQ_AUTO;
                break;
        case KALIGN_TYPE_PROTEIN:
        case KALIGN_TYPE_PROTEIN_PFASUM60:
        case KALIGN_TYPE_PROTEIN_CORBLOSUM66:
        case KALIGN_TYPE_PROTEIN_DIVERGENT:
                seq_type = KALIGN_SEQ_PROTEIN;
                break;
        case KALIGN_TYPE_DNA:
        case KALIGN_TYPE_DNA_INTERNAL:
        case KALIGN_MATRIX_NUC_1PAM:
        case KALIGN_MATRIX_NUC_20PAM:
        case KALIGN_MATRIX_NUC_200PAM:
                seq_type = KALIGN_SEQ_DNA;
                break;
        case KALIGN_TYPE_RNA:
                seq_type = KALIGN_SEQ_RNA;
                break;
        default:
                ERROR_MSG("Unknown sequence type %d.", type);
        }

        RUN(kalign_arr_to_msa(seq, len, numseq, &msa));
        msa->quiet = 1;
        if(n_threads < 1){
                n_threads = 1;
        }

        RUN(run_mode_on_msa(msa, "fast", seq_type, gpo, gpe, tgpe, n_threads));

        RUN(kalign_msa_to_arr(msa, aligned, out_aln_len));
        kalign_free_msa(msa);
        return OK;
ERROR:
        if(msa) kalign_free_msa(msa);
        return FAIL;
}

/* ======================================================================== */
/* Config defaults and unified entry point                                   */
/* ======================================================================== */

struct kalign_run_config kalign_run_config_defaults(void)
{
        struct kalign_run_config cfg;
        cfg.matrix = KALIGN_MATRIX_PFASUM43;
        cfg.gpo = 7.0f;          /* PFASUM43 penalty */
        cfg.gpe = 1.25f;
        cfg.tgpe = 1.0f;
        cfg.vsm_amax = 2.0f;     /* VSM amplitude */
        cfg.seq_weights = 0.0f;
        cfg.dist_scale = 0.0f;
        cfg.refine = KALIGN_REFINE_NONE;
        cfg.adaptive_budget = 0;
        cfg.realign = 0;
        cfg.tree_seed = 0;
        cfg.tree_noise = 0.0f;
        cfg.consistency_anchors = 0;
        cfg.consistency_weight = 2.0f;
        cfg.vsm_mode = 0;        /* VSM with vsm_amax */
        cfg.vsm_kappa = 1.0f;
        cfg.vsm_dzero = 1.0f;
        cfg.guide_dist = KALIGN_GUIDE_DIST_BPM;  /* BPM guide distance (default) */
        cfg.gap_model = 0;       /* Kalign gap rule */
        cfg.tree_method = -1;    /* -1: default tree (sl_opts_resolve) */
        cfg.tree_dist = -1;
        cfg.tree_alpha = -1;
        cfg.tree_maxlen = -1;
        cfg.tree_leaf = -1;
        cfg.tree_kmeans_dist = -1;
        cfg.realign_tree = -1;
        cfg.tree_sl_noise = 0.0f;
        return cfg;
}

struct kalign_ensemble_config kalign_ensemble_config_defaults(void)
{
        struct kalign_ensemble_config ens;
        ens.min_support = 0;
        ens.consistency_merge = 0;
        ens.consistency_merge_weight = 2.0f;
        ens.fill_in = 0;
        ens.fill_anchor_frac = 0.0f;
        ens.fill_support = 1;
        return ens;
}

int kalign_align_full(struct msa* msa,
                      const struct kalign_run_config* runs,
                      int n_runs,
                      const struct kalign_ensemble_config* ens,
                      int n_threads)
{
        threadpool_t *pool = NULL;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(runs != NULL, "No run configs");
        ASSERT(n_runs >= 1, "n_runs must be >= 1");

        if(n_threads < 1) n_threads = 1;
        {
                /* More threads than cores only slows the run down: idle workers spin before
                   they park and then compete with the busy ones for the cores. The output
                   does not depend on the number of threads. */
                const long cores = sysconf(_SC_NPROCESSORS_ONLN);
                if(cores >= 1 && n_threads > cores){
                        if(!msa->quiet){
                                WARNING_MSG("%d threads requested; using %ld, the number of cores.", n_threads, cores);
                        }
                        n_threads = (int)cores;
                }
        }

        RUN(msa_report_foreign_chars(msa, msa->biotype));

        /* n_threads is the total number of threads that run tasks.
           tp_create(n) starts n workers besides the calling thread, which
           also executes tasks, so the pool gets n_threads - 1 workers; with
           n_threads = 1 no pool is created and all pool calls run serially
           on this thread (threadpool.h, NULL pool). */
        if(n_threads > 1){
                pool = tp_create(n_threads - 1);
                ASSERT(pool != NULL, "Could not create thread pool");
        }
        msa->pool = pool;

        if(n_runs > 1){
                /* The ensemble works on finalised member alignments */
                msa->defer_finalise = 0;
                RUN(kalign_ensemble_from_configs(msa, runs, n_runs, ens, n_threads));
        }else{
                RUN(kalign_single_run(msa, &runs[0], n_threads));
        }

        msa->pool = NULL;
        tp_destroy(pool);
        return OK;
ERROR:
        msa->pool = NULL;
        tp_destroy(pool);
        return FAIL;
}

/* ======================================================================== */
/* Mode presets (Kalign 4)                                                   */
/* ======================================================================== */

/* The table is generated from presets/results/v2/presets_v2.json by
   presets/release/emit_c_presets.py: two modes per sequence type, chosen at
   the knees of the pooled nested cross-validation Pareto front
   (presets/results/v2/RESULTS_v2.md); protein fast is the fast-fit knee
   setting (presets/results/fastfit/final_pick.json: member gap model and
   hybrid single-linkage guide tree); protein accurate is the confidence re-fit's
   pick (presets/results/ensconf/final_pick.json: the ensemble fit's 5
   members, member 1 = protein fast, all on member 1's guide tree with
   per-member single-linkage tree noise, consensus min_support 2). Each run is built exactly as
   python-kalign/_core.cpp ensemble_custom_file_to_file builds it from the
   same parameters (seed 42 -> tree_seed 42 + k; per-member tree_noise), so
   the preset path reproduces the evaluated alignments byte for byte. */
#include "mode_presets.h"

static void preset_apply(const struct kalign_mode_preset *p,
                         struct kalign_run_config *runs,
                         int *n_runs,
                         struct kalign_ensemble_config *ens)
{
        *n_runs = p->n_runs;
        for(int k = 0; k < p->n_runs; k++){
                struct kalign_run_config *r = &runs[k];
                *r = kalign_run_config_defaults();
                r->matrix = p->matrix[k];
                r->gpo = p->gpo[k];
                r->gpe = p->gpe[k];
                r->tgpe = p->tgpe[k];
                r->tree_seed = 42 + (uint64_t)k;
                r->tree_noise = p->tree_noise[k];
                r->vsm_amax = p->vsm_amax[k];
                r->dist_scale = 0.0f;
                r->seq_weights = p->seq_weights[k];
                r->refine = p->refine[k];
                r->realign = p->realign[k];
                r->consistency_anchors = p->consistency_anchors[k];
                r->consistency_weight = p->consistency_weight[k];
                r->vsm_mode = p->vsm_mode[k];
                r->vsm_kappa = p->vsm_kappa[k];
                r->vsm_dzero = p->vsm_dzero[k];
                r->guide_dist = p->guide_dist[k];
                r->gap_model = p->gap_model;
                r->tree_method = p->tree_method[k];
                r->tree_dist = p->tree_dist[k];
                r->tree_alpha = p->tree_alpha[k];
                r->tree_maxlen = p->tree_maxlen[k];
                r->tree_leaf = p->tree_leaf[k];
                r->tree_kmeans_dist = p->tree_kmeans_dist[k];
                r->realign_tree = p->realign_tree[k];
                r->tree_sl_noise = p->tree_sl_noise[k];
        }
        ens->min_support = p->min_support;
        ens->consistency_merge = p->consistency_merge;
        ens->consistency_merge_weight = p->consistency_merge_weight;
        ens->fill_in = p->fill_in;
        ens->fill_anchor_frac = p->fill_anchor_frac;
        ens->fill_support = p->fill_support;
}

int kalign_get_mode_preset(const char *mode,
                            int biotype,
                            struct kalign_run_config *runs,
                            int *n_runs,
                            struct kalign_ensemble_config *ens)
{
        const char *m = mode ? mode : "fast";
        int nuc = (biotype == ALN_BIOTYPE_DNA);
        int idx;

        *ens = kalign_ensemble_config_defaults();

        if(strcasecmp(m, "fast") == 0 ||
           strcasecmp(m, "default") == 0){
                /* "default" is an alias of "fast" (the default mode) */
                idx = nuc ? KALIGN_PRESET_NUCLEOTIDE_FAST : KALIGN_PRESET_PROTEIN_FAST;
        }else if(strcasecmp(m, "accurate") == 0 ||
                 strcasecmp(m, "recall") == 0){
                /* "recall" is an alias of "accurate" */
                idx = nuc ? KALIGN_PRESET_NUCLEOTIDE_ACCURATE : KALIGN_PRESET_PROTEIN_ACCURATE;
        }else{
                return -1;
        }
        preset_apply(&kalign_mode_presets[idx], runs, n_runs, ens);
        return 0;
}

int kalign_mode_size_ok(const char *mode, int biotype, int numseq)
{
        const char *m = mode ? mode : "fast";

        (void) biotype;  /* same limit for protein and nucleotide */
        if(strcasecmp(m, "accurate") == 0 || strcasecmp(m, "recall") == 0){
                return numseq <= KALIGN_ACCURATE_MAX_SEQ;
        }
        return 1;
}

int kalign_parse_seq_type(const char *name, int *seq_type)
{
        if(name == NULL || strcasecmp(name, "auto") == 0){
                *seq_type = KALIGN_SEQ_AUTO;
        }else if(strcasecmp(name, "protein") == 0){
                *seq_type = KALIGN_SEQ_PROTEIN;
        }else if(strcasecmp(name, "dna") == 0){
                *seq_type = KALIGN_SEQ_DNA;
        }else if(strcasecmp(name, "rna") == 0){
                *seq_type = KALIGN_SEQ_RNA;
        }else{
                return -1;
        }
        return 0;
}

int kalign_set_seq_type(struct msa *msa, int seq_type, int *overridden)
{
        int want;

        if(overridden) *overridden = 0;
        if(msa == NULL) return -1;

        switch(seq_type){
        case KALIGN_SEQ_AUTO:
                return 0;
        case KALIGN_SEQ_PROTEIN:
                want = ALN_BIOTYPE_PROTEIN;
                break;
        case KALIGN_SEQ_DNA:
        case KALIGN_SEQ_RNA:
                want = ALN_BIOTYPE_DNA;  /* DNA and RNA share one biotype */
                break;
        default:
                return -1;
        }
        if(overridden && msa->biotype != ALN_BIOTYPE_UNDEF && msa->biotype != want){
                *overridden = 1;
        }
        msa->biotype = want;
        return 0;
}

const char *kalign_seq_type_name(struct msa *msa)
{
        if(msa == NULL) return "unknown";
        if(msa->biotype == ALN_BIOTYPE_PROTEIN) return "protein";
        if(msa->biotype == ALN_BIOTYPE_DNA) return "DNA/RNA";
        return "unknown";
}
