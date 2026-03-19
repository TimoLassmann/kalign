#include "tldevel.h"
#include "tlmisc.h"
#include <strings.h>
#include "esl_stopwatch.h"
#include "task.h"
#include "msa_struct.h"
#include "msa_op.h"
#include "msa_alloc.h"
#include "msa_check.h"
#include "msa_sort.h"
#include "msa_io.h"

#include "alphabet.h"
#include "bisectingKmeans.h"

#include "aln_param.h"
#include "aln_run.h"
#include "aln_refine.h"
#include "aln_apair_dist.h"
#include "anchor_consistency.h"
#include "kalign/kalign.h"
#include "kalign/kalign_config.h"
#include "ensemble.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

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

int kalign(char **seq, int *len, int numseq,int n_threads, int type, float gpo, float gpe, float tgpe, char ***aligned, int *out_aln_len)
{
        struct msa *msa = NULL;
        RUN(kalign_arr_to_msa(seq, len,numseq, &msa));

        msa->quiet = 1;
        if(n_threads < 1){
                n_threads = 1;
        }
        RUN(kalign_run(msa,n_threads, type,  gpo, gpe, tgpe, KALIGN_REFINE_NONE, 0));

        RUN(kalign_msa_to_arr(msa, aligned, out_aln_len));

        kalign_free_msa(msa);

        return OK;
ERROR:
        if(msa){
                kalign_free_msa(msa);
        }
        return FAIL;
}

int kalign_run_seeded(struct msa *msa, int n_threads, int type,
                      float gpo, float gpe, float tgpe,
                      int refine, int adaptive_budget,
                      uint64_t tree_seed, float tree_noise,
                      float dist_scale, float vsm_amax,
                      float use_seq_weights,
                      int consistency_anchors, float consistency_weight)
{
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;
        /* This also adds the ranks of the sequences !  */
        RUN(kalign_essential_input_check(msa, 0));

        /* If already aligned unalign ! */
        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }
        /* Make sure sequences are in order  */
        RUN(msa_sort_len_name(msa));

        /* Convert into internal representation  */
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

#ifdef HAVE_OPENMP
        omp_set_num_threads(n_threads);
#endif
        /* Build guide tree - noisy variant if seed != 0 */
        if(tree_seed != 0 && tree_noise > 0.0f){
                RUN(build_tree_kmeans_noisy(msa, &tasks, tree_seed, tree_noise));
        }else{
                RUN(build_tree_kmeans(msa, &tasks));
        }

        /* Convert to full alphabet after having converted to reduced alphabet for tree building above  */
        if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }

        /* Resolve auto matrix selection using BPM distances */
        RUN(resolve_matrix_auto(msa, &type));

        /* align  */
        RUN(aln_param_init(&ap,
                           msa->biotype,
                           n_threads,
                           type,
                           gpo,
                           gpe,
                           tgpe));
        ap->adaptive_budget = adaptive_budget;
        if(use_seq_weights >= 0.0f){
                ap->use_seq_weights = use_seq_weights;
        }
        if(dist_scale > 0.0f){
                ap->dist_scale = dist_scale;
        }
        if(vsm_amax >= 0.0f){
                ap->vsm_amax = vsm_amax;
        }

        if(ap->use_seq_weights > 0.0f){
                RUN(compute_tree_weights(msa, tasks));
        }

        /* Build anchor consistency table if requested */
        if(consistency_anchors > 0){
                ap->consistency_anchors = consistency_anchors;
                ap->consistency_weight = consistency_weight;
                RUN(anchor_consistency_build(msa, ap, consistency_anchors,
                                             consistency_weight,
                                             (struct consistency_table**)&msa->consistency_table));
        }

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Aligning");
        }
        START_TIMER(t1);

        if(refine == KALIGN_REFINE_INLINE){
                RUN(create_msa_tree_inline_refine(msa, ap, tasks, 3));
        }else{
                RUN(create_msa_tree(msa, ap, tasks));
        }
        msa->aligned = ALN_STATUS_ALIGNED;

        /* Optional iterative refinement (two-pass approach) */
        if(refine != KALIGN_REFINE_NONE && refine != KALIGN_REFINE_INLINE){
                RUN(refine_alignment(msa, ap, tasks, refine));
        }

        /* Free consistency table AFTER refinement */
        if(msa->consistency_table){
                anchor_consistency_free((struct consistency_table*)msa->consistency_table);
                msa->consistency_table = NULL;
        }

        RUN(finalise_alignment(msa));

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
        if(msa->consistency_table){
                anchor_consistency_free((struct consistency_table*)msa->consistency_table);
                msa->consistency_table = NULL;
        }
        aln_param_free(ap);
        free_tasks(tasks);
        return FAIL;
}

int kalign_run(struct msa *msa, int n_threads, int type, float gpo, float gpe, float tgpe, int refine, int adaptive_budget)
{
        return kalign_run_seeded(msa, n_threads, type, gpo, gpe, tgpe, refine, adaptive_budget, 0, 0.0f, 0.0f, -1.0f, -1.0f, 0, 2.0f);
}

int kalign_run_dist_scale(struct msa *msa, int n_threads, int type,
                          float gpo, float gpe, float tgpe,
                          int refine, int adaptive_budget,
                          float dist_scale, float vsm_amax,
                          float use_seq_weights)
{
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;
        RUN(kalign_essential_input_check(msa, 0));

        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }
        RUN(msa_sort_len_name(msa));

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

#ifdef HAVE_OPENMP
        omp_set_num_threads(n_threads);
#endif
        RUN(build_tree_kmeans(msa, &tasks));

        if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }

        RUN(resolve_matrix_auto(msa, &type));

        RUN(aln_param_init(&ap,
                           msa->biotype,
                           n_threads,
                           type,
                           gpo,
                           gpe,
                           tgpe));
        ap->adaptive_budget = adaptive_budget;
        if(use_seq_weights >= 0.0f){
                ap->use_seq_weights = use_seq_weights;
        }
        ap->dist_scale = dist_scale;
        if(vsm_amax >= 0.0f){
                ap->vsm_amax = vsm_amax;
        }

        if(ap->use_seq_weights > 0.0f){
                RUN(compute_tree_weights(msa, tasks));
        }

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Aligning (dist_scale=%.2f, vsm_amax=%.2f)", dist_scale, vsm_amax);
        }
        START_TIMER(t1);

        if(refine == KALIGN_REFINE_INLINE){
                RUN(create_msa_tree_inline_refine(msa, ap, tasks, 3));
        }else{
                RUN(create_msa_tree(msa, ap, tasks));
        }
        msa->aligned = ALN_STATUS_ALIGNED;

        if(refine != KALIGN_REFINE_NONE && refine != KALIGN_REFINE_INLINE){
                RUN(refine_alignment(msa, ap, tasks, refine));
        }

        RUN(finalise_alignment(msa));
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
        aln_param_free(ap);
        free_tasks(tasks);
        return FAIL;
}

int kalign_run_realign(struct msa *msa, int n_threads, int type,
                       float gpo, float gpe, float tgpe,
                       int refine, int adaptive_budget,
                       float dist_scale, float vsm_amax,
                       int realign_iterations,
                       float use_seq_weights,
                       int consistency_anchors, float consistency_weight)
{
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;
        int iter;

        RUN(kalign_essential_input_check(msa, 0));

        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }
        RUN(msa_sort_len_name(msa));

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

#ifdef HAVE_OPENMP
        omp_set_num_threads(n_threads);
#endif
        /* Initial guide tree from BPM anchor distances */
        RUN(build_tree_kmeans(msa, &tasks));

        if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }

        RUN(resolve_matrix_auto(msa, &type));

        RUN(aln_param_init(&ap,
                           msa->biotype,
                           n_threads,
                           type,
                           gpo,
                           gpe,
                           tgpe));
        ap->adaptive_budget = adaptive_budget;
        if(use_seq_weights >= 0.0f){
                ap->use_seq_weights = use_seq_weights;
        }
        ap->dist_scale = dist_scale;
        if(vsm_amax >= 0.0f){
                ap->vsm_amax = vsm_amax;
        }

        if(ap->use_seq_weights > 0.0f){
                RUN(compute_tree_weights(msa, tasks));
        }

        /* Build anchor consistency table if requested */
        if(consistency_anchors > 0){
                ap->consistency_anchors = consistency_anchors;
                ap->consistency_weight = consistency_weight;
                RUN(anchor_consistency_build(msa, ap, consistency_anchors,
                                             consistency_weight,
                                             (struct consistency_table**)&msa->consistency_table));
        }

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Aligning (realign=%d, dist_scale=%.2f, vsm_amax=%.2f)",
                        realign_iterations, dist_scale, vsm_amax);
        }
        START_TIMER(t1);

        /* First alignment with BPM-based guide tree */
        if(refine == KALIGN_REFINE_INLINE){
                RUN(create_msa_tree_inline_refine(msa, ap, tasks, 3));
        }else{
                RUN(create_msa_tree(msa, ap, tasks));
        }
        msa->aligned = ALN_STATUS_ALIGNED;

        /* Iterative realignment: align -> compute distances -> new tree -> re-align */
        for(iter = 0; iter < realign_iterations; iter++){
                float** dm = NULL;
                int si;

                /* Finalize to get character sequences with gap characters */
                RUN(finalise_alignment(msa));

                /* Compute NxN pairwise identity distances from alignment */
                RUN(compute_aln_pairwise_dist(msa, &dm));

                /* Remove gaps, reset alignment status.
                   dealign_msa zeroes the gaps[] array but does NOT strip '-'
                   from seq->seq (which was linearized by finalise_alignment).
                   We must rebuild seq->seq without gap characters. */
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

                /* Re-encode internal representation for alignment */
                if(msa->biotype == ALN_BIOTYPE_DNA){
                        RUN(convert_msa_to_internal(msa, ALPHA_defDNA));
                }else if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                        RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
                }

                /* Reset profile tracking */
                RUN(set_sip_nsip(msa));

                /* Rebuild guide tree from alignment-derived distances */
                free_tasks(tasks);
                tasks = NULL;
                RUN(alloc_tasks(&tasks, msa->numseq));
                RUN(build_tree_from_pairwise(msa, &tasks, dm));
                free_aln_dm(dm, msa->numseq);

                if(ap->use_seq_weights > 0.0f){
                        RUN(compute_tree_weights(msa, tasks));
                }

                /* Re-align with new tree */
                if(refine == KALIGN_REFINE_INLINE){
                        RUN(create_msa_tree_inline_refine(msa, ap, tasks, 3));
                }else{
                        RUN(create_msa_tree(msa, ap, tasks));
                }
                msa->aligned = ALN_STATUS_ALIGNED;
        }

        /* Refinement after all realign iterations (two-pass, skip for inline) */
        if(refine != KALIGN_REFINE_NONE && refine != KALIGN_REFINE_INLINE){
                RUN(refine_alignment(msa, ap, tasks, refine));
        }

        /* Free consistency table AFTER refinement */
        if(msa->consistency_table){
                anchor_consistency_free((struct consistency_table*)msa->consistency_table);
                msa->consistency_table = NULL;
        }

        RUN(finalise_alignment(msa));
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
        if(msa->consistency_table){
                anchor_consistency_free((struct consistency_table*)msa->consistency_table);
                msa->consistency_table = NULL;
        }
        aln_param_free(ap);
        free_tasks(tasks);
        return FAIL;
}

int kalign_post_realign(struct msa *msa, int n_threads, int type,
                        float gpo, float gpe, float tgpe,
                        int refine, int adaptive_budget,
                        float dist_scale, float vsm_amax,
                        int realign_iterations,
                        float use_seq_weights)
{
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;
        int iter;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(realign_iterations > 0, "Need at least 1 realign iteration");

        /* Detect biotype if not set */
        if(msa->biotype == ALN_BIOTYPE_UNDEF){
                RUN(detect_alphabet(msa));
        }

        /* seq_distances available from prior alignment */
        RUN(resolve_matrix_auto(msa, &type));

        RUN(aln_param_init(&ap,
                           msa->biotype,
                           n_threads,
                           type,
                           gpo,
                           gpe,
                           tgpe));
        ap->adaptive_budget = adaptive_budget;
        if(use_seq_weights >= 0.0f){
                ap->use_seq_weights = use_seq_weights;
        }
        ap->dist_scale = dist_scale;
        if(vsm_amax >= 0.0f){
                ap->vsm_amax = vsm_amax;
        }

#ifdef HAVE_OPENMP
        omp_set_num_threads(n_threads);
#endif

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Post-realign (%d iterations, vsm_amax=%.2f)",
                        realign_iterations, ap->vsm_amax);
        }
        START_TIMER(t1);

        for(iter = 0; iter < realign_iterations; iter++){
                float** dm = NULL;
                int si;

                /* Finalize if not already (first iter may already be FINAL from ensemble) */
                if(msa->aligned != ALN_STATUS_FINAL){
                        RUN(finalise_alignment(msa));
                }

                /* Compute NxN pairwise identity distances from alignment */
                RUN(compute_aln_pairwise_dist(msa, &dm));

                /* Strip gap characters from seq->seq and fix seq->len.
                   Consensus alignment may have set len to alignment length,
                   so we recompute from the ungapped sequence.
                   We also zero gaps[] and reset alignment status manually
                   (rather than calling dealign_msa which uses the possibly
                   wrong len to bound the gaps[] loop). */
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
                        /* Zero gaps array (len+1 entries) */
                        for(r = 0; r <= w; r++){
                                seq->gaps[r] = 0;
                        }
                }
                msa->aligned = ALN_STATUS_UNALIGNED;

                /* Re-encode to internal representation */
                if(msa->biotype == ALN_BIOTYPE_DNA){
                        RUN(convert_msa_to_internal(msa, ALPHA_defDNA));
                }else if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                        RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
                }

                /* Reset profile tracking */
                RUN(set_sip_nsip(msa));

                /* Build UPGMA tree from alignment-derived distances */
                if(tasks){ free_tasks(tasks); tasks = NULL; }
                RUN(alloc_tasks(&tasks, msa->numseq));
                RUN(build_tree_from_pairwise(msa, &tasks, dm));
                free_aln_dm(dm, msa->numseq);

                if(ap->use_seq_weights > 0.0f){
                        RUN(compute_tree_weights(msa, tasks));
                }

                /* Re-align with new tree */
                if(refine == KALIGN_REFINE_INLINE){
                        RUN(create_msa_tree_inline_refine(msa, ap, tasks, 3));
                }else{
                        RUN(create_msa_tree(msa, ap, tasks));
                }
                msa->aligned = ALN_STATUS_ALIGNED;
        }

        /* Refinement after all realign iterations (two-pass, skip for inline) */
        if(refine != KALIGN_REFINE_NONE && refine != KALIGN_REFINE_INLINE){
                RUN(refine_alignment(msa, ap, tasks, refine));
        }

        RUN(finalise_alignment(msa));
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
        aln_param_free(ap);
        if(tasks) free_tasks(tasks);
        return FAIL;
}

/* ======================================================================== */
/* Config defaults and unified entry point                                   */
/* ======================================================================== */

struct kalign_run_config kalign_run_config_defaults(void)
{
        struct kalign_run_config cfg;
        cfg.matrix = KALIGN_MATRIX_PFASUM43;
        cfg.gpo = 7.0f;          /* PFASUM43 default */
        cfg.gpe = 1.25f;
        cfg.tgpe = 1.0f;
        cfg.vsm_amax = 2.0f;     /* protein default */
        cfg.seq_weights = 0.0f;
        cfg.dist_scale = 0.0f;
        cfg.refine = KALIGN_REFINE_NONE;
        cfg.adaptive_budget = 0;
        cfg.realign = 0;
        cfg.tree_seed = 0;
        cfg.tree_noise = 0.0f;
        cfg.consistency_anchors = 0;
        cfg.consistency_weight = 2.0f;
        return cfg;
}

struct kalign_ensemble_config kalign_ensemble_config_defaults(void)
{
        struct kalign_ensemble_config ens;
        ens.min_support = 0;
        ens.consistency_merge = 0;
        ens.consistency_merge_weight = 2.0f;
        return ens;
}

int kalign_align_full(struct msa* msa,
                      const struct kalign_run_config* runs,
                      int n_runs,
                      const struct kalign_ensemble_config* ens,
                      int n_threads)
{
        ASSERT(msa != NULL, "No MSA");
        ASSERT(runs != NULL, "No run configs");
        ASSERT(n_runs >= 1, "n_runs must be >= 1");

        if(n_runs > 1){
                /* Ensemble path */
                RUN(kalign_ensemble_from_configs(msa, runs, n_runs, ens, n_threads));
        }else{
                /* Single-run path */
                const struct kalign_run_config* r = &runs[0];
                if(r->realign > 0){
                        RUN(kalign_run_realign(msa, n_threads, r->matrix,
                                              r->gpo, r->gpe, r->tgpe,
                                              r->refine, r->adaptive_budget,
                                              r->dist_scale, r->vsm_amax,
                                              r->realign, r->seq_weights,
                                              r->consistency_anchors,
                                              r->consistency_weight));
                }else{
                        RUN(kalign_run_seeded(msa, n_threads, r->matrix,
                                             r->gpo, r->gpe, r->tgpe,
                                             r->refine, r->adaptive_budget,
                                             r->tree_seed, r->tree_noise,
                                             r->dist_scale, r->vsm_amax,
                                             r->seq_weights,
                                             r->consistency_anchors,
                                             r->consistency_weight));
                }
        }

        return OK;
ERROR:
        return FAIL;
}

/* ======================================================================== */
/* NSGA-III optimized protein mode presets                                   */
/* ======================================================================== */

/* Helper: fill one run config with preset values */
static void preset_run(struct kalign_run_config *r,
                        int matrix, float gpo, float gpe, float tgpe,
                        float vsm_amax, float seq_weights,
                        int realign, int refine,
                        uint64_t seed, float noise)
{
        *r = kalign_run_config_defaults();
        r->matrix = matrix;
        r->gpo = gpo;
        r->gpe = gpe;
        r->tgpe = tgpe;
        r->vsm_amax = vsm_amax;
        r->seq_weights = seq_weights;
        r->realign = realign;
        r->refine = refine;
        r->tree_seed = seed;
        r->tree_noise = noise;
}

/* Helper: set shared consistency params on all runs */
static void preset_consistency(struct kalign_run_config *runs, int n_runs,
                                int anchors, float weight)
{
        for(int i = 0; i < n_runs; i++){
                runs[i].consistency_anchors = anchors;
                runs[i].consistency_weight = weight;
        }
}

/* ---- Protein presets (NSGA-III optimized on BAliBASE v4, gen 41) ---- */

static int preset_protein(const char *m,
                           struct kalign_run_config *runs,
                           int *n_runs,
                           struct kalign_ensemble_config *ens)
{
        /* fast: single run, ~34s on BAliBASE. F1=0.737 R=0.818 P=0.670 */
        if(strcasecmp(m, "fast") == 0){
                *n_runs = 1;
                preset_run(&runs[0],
                           KALIGN_MATRIX_PFASUM60,
                           8.626f, 0.843f, 0.433f,
                           0.592f, 1.534f,
                           0, KALIGN_REFINE_NONE,
                           0, 0.0f);
                return 0;
        }

        /* default: single run with inline refinement, ~130s.
           F1=0.750 R=0.779 P=0.724 */
        if(strcasecmp(m, "default") == 0){
                *n_runs = 1;
                preset_run(&runs[0],
                           KALIGN_MATRIX_PFASUM60,
                           8.121f, 0.684f, 0.560f,
                           1.661f, 2.356f,
                           0, KALIGN_REFINE_INLINE,
                           0, 0.0f);
                preset_consistency(runs, 1, 1, 1.640f);
                return 0;
        }

        /* recall: 5-run ensemble optimized for recall, ~1175s.
           F1=0.777 R=0.837 P=0.726 */
        if(strcasecmp(m, "recall") == 0){
                *n_runs = 5;
                preset_run(&runs[0], KALIGN_MATRIX_CORBLOSUM66,
                           5.416f, 1.071f, 1.620f,
                           2.536f, 0.282f,
                           1, KALIGN_REFINE_ALL, 42, 0.0f);
                preset_run(&runs[1], KALIGN_MATRIX_CORBLOSUM66,
                           12.091f, 1.024f, 1.284f,
                           2.078f, 0.282f,
                           1, KALIGN_REFINE_INLINE, 43, 0.0f);
                preset_run(&runs[2], KALIGN_MATRIX_PFASUM60,
                           6.970f, 2.919f, 0.632f,
                           0.877f, 0.282f,
                           1, KALIGN_REFINE_NONE, 44, 0.0f);
                preset_run(&runs[3], KALIGN_MATRIX_PFASUM43,
                           6.173f, 1.102f, 0.510f,
                           1.276f, 0.282f,
                           1, KALIGN_REFINE_ALL, 45, 0.0f);
                preset_run(&runs[4], KALIGN_MATRIX_PFASUM43,
                           5.278f, 1.764f, 1.088f,
                           2.315f, 0.282f,
                           1, KALIGN_REFINE_CONFIDENT, 46, 0.0f);
                preset_consistency(runs, 5, 3, 2.120f);
                ens->min_support = 0;
                return 0;
        }

        /* accurate: 5-run ensemble optimized for F1, ~2005s.
           F1=0.814 R=0.782 P=0.848 */
        if(strcasecmp(m, "accurate") == 0){
                *n_runs = 5;
                preset_run(&runs[0], KALIGN_MATRIX_PFASUM43,
                           8.682f, 0.650f, 1.465f,
                           2.532f, 1.468f,
                           1, KALIGN_REFINE_INLINE, 42, 0.0f);
                preset_run(&runs[1], KALIGN_MATRIX_PFASUM43,
                           6.148f, 1.174f, 1.297f,
                           2.011f, 1.468f,
                           1, KALIGN_REFINE_ALL, 43, 0.0f);
                preset_run(&runs[2], KALIGN_MATRIX_CORBLOSUM66,
                           3.622f, 1.004f, 0.631f,
                           0.872f, 1.468f,
                           1, KALIGN_REFINE_CONFIDENT, 44, 0.0f);
                preset_run(&runs[3], KALIGN_MATRIX_CORBLOSUM66,
                           5.988f, 0.552f, 0.495f,
                           1.914f, 1.468f,
                           1, KALIGN_REFINE_CONFIDENT, 45, 0.0f);
                preset_run(&runs[4], KALIGN_MATRIX_PFASUM43,
                           13.939f, 2.629f, 1.406f,
                           1.830f, 1.468f,
                           1, KALIGN_REFINE_NONE, 46, 0.0f);
                preset_consistency(runs, 5, 8, 2.457f);
                ens->min_support = 3;
                return 0;
        }

        return -1;
}

/* ---- Nucleotide presets (NSGA-III optimized on BRAliBASE + MDSA, gen 100) ----
 *
 * Unified presets for both DNA and RNA. Optimized on combined
 * BRAliBASE (599 RNA cases) + MDSA (325 DNA cases).
 * All presets use Kimura two-parameter matrices.
 */

static int preset_nucleotide(const char *m,
                              struct kalign_run_config *runs,
                              int *n_runs,
                              struct kalign_ensemble_config *ens)
{
        /* fast: single run with realign, ~5s. F1=0.790 R=0.792 P=0.788 */
        if(strcasecmp(m, "fast") == 0){
                *n_runs = 1;
                preset_run(&runs[0],
                           KALIGN_MATRIX_NUC_200PAM,
                           19.873f, 0.642f, 2.582f,
                           0.874f, 0.199f,
                           1, KALIGN_REFINE_NONE,
                           0, 0.0f);
                return 0;
        }

        /* default: 3-run ensemble with realign, ~26s.
           F1=0.806 R=0.773 P=0.842 */
        if(strcasecmp(m, "default") == 0){
                *n_runs = 3;
                preset_run(&runs[0], KALIGN_MATRIX_NUC_20PAM,
                           19.275f, 0.628f, 2.854f,
                           0.445f, 0.448f,
                           1, KALIGN_REFINE_NONE, 42, 0.0f);
                preset_run(&runs[1], KALIGN_MATRIX_NUC_200PAM,
                           16.581f, 0.274f, 2.491f,
                           0.706f, 0.448f,
                           1, KALIGN_REFINE_NONE, 43, 0.0f);
                preset_run(&runs[2], KALIGN_MATRIX_NUC_200PAM,
                           16.499f, 3.754f, 2.896f,
                           0.491f, 0.448f,
                           1, KALIGN_REFINE_CONFIDENT, 44, 0.0f);
                ens->min_support = 2;
                return 0;
        }

        /* recall: single run with realign, ~17s.
           F1=0.798 R=0.800 P=0.796 */
        if(strcasecmp(m, "recall") == 0){
                *n_runs = 1;
                preset_run(&runs[0],
                           KALIGN_MATRIX_NUC_200PAM,
                           19.882f, 0.646f, 2.854f,
                           0.874f, 0.449f,
                           1, KALIGN_REFINE_CONFIDENT,
                           0, 0.0f);
                preset_consistency(runs, 1, 1, 1.616f);
                return 0;
        }

        /* accurate: 5-run ensemble with realign, ~100s.
           F1=0.810 R=0.760 P=0.867 */
        if(strcasecmp(m, "accurate") == 0){
                *n_runs = 5;
                preset_run(&runs[0], KALIGN_MATRIX_NUC_200PAM,
                           19.475f, 0.728f, 2.856f,
                           0.492f, 0.590f,
                           1, KALIGN_REFINE_CONFIDENT, 42, 0.0f);
                preset_run(&runs[1], KALIGN_MATRIX_NUC_1PAM,
                           19.377f, 1.379f, 2.486f,
                           1.211f, 0.590f,
                           1, KALIGN_REFINE_CONFIDENT, 43, 0.0f);
                preset_run(&runs[2], KALIGN_MATRIX_NUC_200PAM,
                           15.869f, 0.353f, 2.430f,
                           0.721f, 0.590f,
                           1, KALIGN_REFINE_NONE, 44, 0.0f);
                preset_run(&runs[3], KALIGN_MATRIX_NUC_200PAM,
                           15.989f, 2.331f, 2.764f,
                           0.188f, 0.590f,
                           1, KALIGN_REFINE_CONFIDENT, 45, 0.0f);
                preset_run(&runs[4], KALIGN_MATRIX_NUC_1PAM,
                           5.346f, 3.483f, 0.610f,
                           1.406f, 0.590f,
                           1, KALIGN_REFINE_NONE, 46, 0.0f);
                preset_consistency(runs, 5, 2, 0.800f);
                ens->min_support = 3;
                return 0;
        }

        return -1;
}

int kalign_get_mode_preset(const char *mode,
                            int biotype,
                            struct kalign_run_config *runs,
                            int *n_runs,
                            struct kalign_ensemble_config *ens)
{
        const char *m = mode ? mode : "default";

        *ens = kalign_ensemble_config_defaults();

        if(biotype == ALN_BIOTYPE_DNA){
                return preset_nucleotide(m, runs, n_runs, ens);
        }

        /* Default: protein presets */
        return preset_protein(m, runs, n_runs, ens);
}
