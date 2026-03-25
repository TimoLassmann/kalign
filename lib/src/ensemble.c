#include "tldevel.h"
#include "esl_stopwatch.h"

#include "msa_struct.h"
#include "msa_op.h"
#include "msa_alloc.h"
#include "msa_sort.h"
#include "msa_check.h"

#include "aln_param.h"
#include "aln_wrap.h"
#include "poar.h"
#include "consensus_msa.h"
#include "anchor_consistency.h"
#include "msa_consistency.h"

#include "kalign/kalign.h"
#include "kalign/kalign_config.h"

#define ENSEMBLE_IMPORT
#include "ensemble.h"

/* Each ensemble run varies gap-open and gap-extend independently,
   with optional iterative refinement on select runs.
   This explores "few long gaps" vs "many short gaps" rather than
   just uniformly scaling all penalties together.
   Entry 0 is unused (run 0 always uses defaults). */
struct ensemble_params {
        float gpo_scale;    /* gap open multiplier */
        float gpe_scale;    /* gap extend multiplier */
        float tgpe_scale;   /* terminal gap extend multiplier */
        float noise;        /* tree noise sigma */
};

static struct ensemble_params run_params[] = {
        {1.0f, 1.0f, 1.0f, 0.0f},    /* 0: default (unused, handled specially) */
        {0.5f, 1.5f, 0.8f, 0.20f},    /* 1: fewer gap-opens, longer extensions */
        {1.5f, 0.5f, 1.2f, 0.20f},    /* 2: more gap-opens, shorter extensions */
        {0.7f, 0.7f, 0.5f, 0.25f},    /* 3: globally relaxed gaps */
        {1.4f, 1.4f, 1.5f, 0.25f},    /* 4: globally strict gaps */
        {0.8f, 1.2f, 1.0f, 0.30f},    /* 5: slight open-relax, extend-strict */
        {1.3f, 0.8f, 0.7f, 0.30f},    /* 6: strict open, relaxed extend+terminal */
        {0.6f, 1.0f, 1.3f, 0.15f},    /* 7: relaxed open only */
        {1.0f, 0.6f, 0.6f, 0.15f},    /* 8: relaxed extend only */
        {1.8f, 1.0f, 1.0f, 0.35f},    /* 9: very strict open, big tree perturbation */
        {1.0f, 1.8f, 1.8f, 0.35f},    /* 10: very strict extend */
        {0.4f, 0.4f, 0.3f, 0.20f},    /* 11: very relaxed all */
};
#define N_RUN_PARAMS 12

/* ---------------------------------------------------------------------------
 * Helper: score_alignments
 *
 * Score all N alignments against the POAR table.  Returns the scores
 * array (caller must MFREE) and the index of the best run (best_k).
 * Run 0 is preferred unless another run exceeds its score by >5%.
 * ------------------------------------------------------------------------- */
static int score_alignments(struct msa** alignments,
                            struct poar_table* poar,
                            int numseq, int n_runs, int quiet,
                            double** out_scores, int* out_best_k
#ifdef USE_THREADPOOL
                            , threadpool_t* pool
#endif
                            )
{
        struct pos_matrix* pm = NULL;
        double* scores = NULL;
        int k;

        MMALLOC(scores, sizeof(double) * n_runs);

        for(k = 0; k < n_runs; k++){
                char** aln_seqs = NULL;

                MMALLOC(aln_seqs, sizeof(char*) * numseq);
                for(int i = 0; i < numseq; i++){
                        aln_seqs[i] = alignments[k]->sequences[i]->seq;
                }

                RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, alignments[k]->alnlen));
#ifdef USE_THREADPOOL
                RUN(score_alignment_poar(poar, pm, numseq, n_runs, &scores[k], pool));
#else
                RUN(score_alignment_poar(poar, pm, numseq, n_runs, &scores[k]));
#endif

                if(!quiet){
                        LOG_MSG("  Run %d score: %.1f", k + 1, scores[k]);
                }

                pos_matrix_free(pm);
                pm = NULL;
                MFREE(aln_seqs);
        }

        /* Select: prefer run 0 (default params) unless another run
           has a meaningfully higher score (>5% improvement). */
        int best_k = 0;
        double baseline = scores[0];
        for(k = 1; k < n_runs; k++){
                if(scores[k] > scores[best_k] && scores[k] > baseline * 1.05){
                        best_k = k;
                }
        }

        *out_scores = scores;
        *out_best_k = best_k;
        return OK;
ERROR:
        if(pm) pos_matrix_free(pm);
        if(scores) MFREE(scores);
        return FAIL;
}

/* ---------------------------------------------------------------------------
 * Helper: build_consensus_from_poar
 *
 * Build a consensus MSA from a POAR table at the given min_support
 * threshold.  Allocates and returns consensus_msa via *out_consensus.
 * The source MSA is used as the template (deep-copied internally).
 * ------------------------------------------------------------------------- */
static int build_consensus_from_poar(struct poar_table* poar,
                                     struct msa* msa,
                                     int numseq, int min_support,
                                     struct msa** out_consensus
#ifdef USE_THREADPOOL
                                     , threadpool_t* pool
#endif
                                     )
{
        struct msa* consensus_msa = NULL;
        int* seq_lens = NULL;

        RUN(msa_cpy(&consensus_msa, msa));

        MMALLOC(seq_lens, sizeof(int) * numseq);
        for(int i = 0; i < numseq; i++){
                seq_lens[i] = msa->sequences[i]->len;
        }

#ifdef USE_THREADPOOL
        RUN(build_consensus(poar, seq_lens, numseq, min_support, consensus_msa, pool));
#else
        RUN(build_consensus(poar, seq_lens, numseq, min_support, consensus_msa));
#endif
        MFREE(seq_lens);

        *out_consensus = consensus_msa;
        return OK;
ERROR:
        if(seq_lens) MFREE(seq_lens);
        if(consensus_msa) kalign_free_msa(consensus_msa);
        return FAIL;
}

/* ---------------------------------------------------------------------------
 * Helper: copy_alignment_to_msa
 *
 * Move alignment sequences from src into dst, transferring ownership.
 * After this call, src's seq pointers are replaced with empty strings
 * so it can be safely freed.
 * ------------------------------------------------------------------------- */
static int copy_alignment_to_msa(struct msa* dst, struct msa* src, int numseq)
{
        for(int i = 0; i < numseq; i++){
                MFREE(dst->sequences[i]->seq);
                dst->sequences[i]->seq = src->sequences[i]->seq;
                dst->sequences[i]->len = src->sequences[i]->len;
                src->sequences[i]->seq = NULL;
                MMALLOC(src->sequences[i]->seq, 1);
                src->sequences[i]->seq[0] = '\0';
        }
        dst->alnlen = src->alnlen;
        dst->aligned = src->aligned;
        return OK;
ERROR:
        return FAIL;
}

/* ---------------------------------------------------------------------------
 * Helper: score_single_msa
 *
 * Score a single MSA against a POAR table.  Returns the score via
 * out_score.  This avoids repeating the aln_seqs + pos_matrix pattern.
 * ------------------------------------------------------------------------- */
static int score_single_msa(struct msa* aln, struct poar_table* poar,
                            int numseq, int n_runs, double* out_score
#ifdef USE_THREADPOOL
                            , threadpool_t* pool
#endif
                            )
{
        struct pos_matrix* pm = NULL;
        char** aln_seqs = NULL;

        MMALLOC(aln_seqs, sizeof(char*) * numseq);
        for(int i = 0; i < numseq; i++){
                aln_seqs[i] = aln->sequences[i]->seq;
        }

        RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, aln->alnlen));
#ifdef USE_THREADPOOL
        RUN(score_alignment_poar(poar, pm, numseq, n_runs, out_score, pool));
#else
        RUN(score_alignment_poar(poar, pm, numseq, n_runs, out_score));
#endif

        pos_matrix_free(pm);
        MFREE(aln_seqs);
        return OK;
ERROR:
        if(pm) pos_matrix_free(pm);
        if(aln_seqs) MFREE(aln_seqs);
        return FAIL;
}

/* kalign_ensemble and kalign_ensemble_custom removed —
   use kalign_align_full with per-run configs instead. */

/* ---- Parallel ensemble run support ---- */
#ifdef USE_THREADPOOL
#include "threadpool/threadpool.h"

struct ensemble_run_arg {
        struct msa* copy;
        const struct kalign_run_config* cfg;
        int n_threads;
        int error;
};

static void ensemble_run_task_fn(void* arg)
{
        struct ensemble_run_arg* ra = (struct ensemble_run_arg*)arg;
        if(kalign_single_run(ra->copy, ra->cfg, ra->n_threads) != 0){
                ra->error = 1;
        }
}
#endif

/* ======================================================================== */

/* kalign_generate_ensemble_runs: expand base config into N diversified runs.
 *
 * IMPORTANT: base.gpo/gpe/tgpe must be resolved (non-sentinel) values.
 * If they are -1.0 (sentinel), the scale factors will produce garbage.
 * The caller should resolve sentinels via aln_param_init before calling this.
 */
int kalign_generate_ensemble_runs(const struct kalign_run_config* base,
                                  int n_runs, uint64_t seed,
                                  struct kalign_run_config* out)
{
        int k;

        ASSERT(base != NULL, "base config is NULL");
        ASSERT(out != NULL, "output array is NULL");
        ASSERT(n_runs >= 1, "n_runs must be >= 1");

        for(k = 0; k < n_runs; k++){
                /* Start with a copy of the base config */
                out[k] = *base;

                if(k == 0){
                        /* Run 0: base params, deterministic tree */
                        out[k].tree_seed = 0;
                        out[k].tree_noise = 0.0f;
                }else{
                        /* Apply diversity table scale factors */
                        struct ensemble_params ep = run_params[k % N_RUN_PARAMS];
                        out[k].gpo = base->gpo * ep.gpo_scale;
                        out[k].gpe = base->gpe * ep.gpe_scale;
                        out[k].tgpe = base->tgpe * ep.tgpe_scale;
                        out[k].tree_seed = seed + (uint64_t)k;
                        out[k].tree_noise = ep.noise;
                }
        }

        return OK;
ERROR:
        return FAIL;
}

/* ======================================================================== */
/* kalign_ensemble_from_configs: run ensemble alignment with per-run configs.
 *
 * This is the core ensemble implementation used by kalign_align_full.
 * Each runs[k] is a fully-specified run configuration.
 */
int kalign_ensemble_from_configs(struct msa* msa,
                                 const struct kalign_run_config* runs,
                                 int n_runs,
                                 const struct kalign_ensemble_config* ens,
                                 int n_threads)
{
        struct msa* copy = NULL;
        struct msa* consensus_msa = NULL;
        struct msa** alignments = NULL;
        struct poar_table* poar = NULL;
        struct pos_matrix* pm = NULL;
        double* scores = NULL;
        int numseq;
        int k;
        int best_k = 0;
        int use_consensus = 0;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(runs != NULL, "No run configs");
        ASSERT(n_runs >= 1, "n_runs must be >= 1");

        RUN(kalign_essential_input_check(msa, 0));

        numseq = msa->numseq;

        DECLARE_TIMER(t_ensemble);
        if(!msa->quiet){
                LOG_MSG("Ensemble alignment with %d runs", n_runs);
        }
        START_TIMER(t_ensemble);

        if(msa->biotype == ALN_BIOTYPE_UNDEF){
                RUN(detect_alphabet(msa));
        }

        RUN(poar_table_alloc(&poar, numseq));
        MMALLOC(alignments, sizeof(struct msa*) * n_runs);
        for(k = 0; k < n_runs; k++){
                alignments[k] = NULL;
        }

        /* Phase timing instrumentation */
        DECLARE_TIMER(t_phase);

        /* Run N alignments concurrently — all runs share the one global
           threadpool, giving the pool N× more tasks to keep workers busy.
           POAR extraction is sequential (sorted insert not thread-safe). */
        START_TIMER(t_phase);
#ifdef USE_THREADPOOL
        if(msa->pool != NULL && n_runs > 1){
                struct ensemble_run_arg* run_args = NULL;
                MMALLOC(run_args, sizeof(struct ensemble_run_arg) * n_runs);

                for(k = 0; k < n_runs; k++){
                        run_args[k].copy = NULL;
                        RUN(msa_cpy(&run_args[k].copy, msa));
                        run_args[k].copy->quiet = 1;
                        run_args[k].copy->pool = msa->pool;
                        run_args[k].cfg = &runs[k];
                        run_args[k].n_threads = n_threads;
                        run_args[k].error = 0;

                        if(!msa->quiet){
                                LOG_MSG("  Run %d/%d (gpo=%.1f gpe=%.1f tgpe=%.1f noise=%.2f)",
                                        k + 1, n_runs,
                                        runs[k].gpo, runs[k].gpe, runs[k].tgpe,
                                        runs[k].tree_noise);
                        }
                }

                /* Fork all runs into the shared pool */
                {
                        tp_group_t *g = tp_group_create(msa->pool);
                        for(k = 0; k < n_runs; k++){
                                tp_group_submit(g, ensemble_run_task_fn, &run_args[k]);
                        }
                        tp_group_wait(g);
                        tp_group_destroy(g);
                }

                /* Check for errors and collect results */
                for(k = 0; k < n_runs; k++){
                        if(run_args[k].error){
                                for(int j = 0; j < n_runs; j++){
                                        if(run_args[j].copy) kalign_free_msa(run_args[j].copy);
                                }
                                MFREE(run_args);
                                ERROR_MSG("Ensemble run %d failed", k + 1);
                        }
                        alignments[k] = run_args[k].copy;
                        run_args[k].copy = NULL;
                }
                MFREE(run_args);

                /* Extract POARs sequentially */
                for(k = 0; k < n_runs; k++){
                        char** aln_seqs = NULL;
                        MMALLOC(aln_seqs, sizeof(char*) * numseq);
                        for(int i = 0; i < numseq; i++){
                                aln_seqs[i] = alignments[k]->sequences[i]->seq;
                        }
                        RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, alignments[k]->alnlen));
                        {
#ifdef USE_THREADPOOL
                                int _ep_ret = extract_poars(poar, pm, k, msa->pool);
#else
                                int _ep_ret = extract_poars(poar, pm, k);
#endif
                                if(_ep_ret != OK) goto ERROR;
                        }
                        pos_matrix_free(pm);
                        pm = NULL;
                        MFREE(aln_seqs);
                }
        }else
#endif
        {
                /* Sequential fallback (no threadpool or single run) */
                for(k = 0; k < n_runs; k++){
                        copy = NULL;
                        RUN(msa_cpy(&copy, msa));
                        copy->quiet = 1;
#ifdef USE_THREADPOOL
                        copy->pool = msa->pool;
#endif
                        if(!msa->quiet){
                                LOG_MSG("  Run %d/%d (gpo=%.1f gpe=%.1f tgpe=%.1f noise=%.2f)",
                                        k + 1, n_runs,
                                        runs[k].gpo, runs[k].gpe, runs[k].tgpe,
                                        runs[k].tree_noise);
                        }
                        RUN(kalign_single_run(copy, &runs[k], n_threads));

                        char** aln_seqs = NULL;
                        MMALLOC(aln_seqs, sizeof(char*) * numseq);
                        for(int i = 0; i < numseq; i++){
                                aln_seqs[i] = copy->sequences[i]->seq;
                        }
                        RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, copy->alnlen));
                        {
#ifdef USE_THREADPOOL
                                int _ep_ret = extract_poars(poar, pm, k, msa->pool);
#else
                                int _ep_ret = extract_poars(poar, pm, k);
#endif
                                if(_ep_ret != OK) goto ERROR;
                        }
                        pos_matrix_free(pm);
                        pm = NULL;
                        MFREE(aln_seqs);

                        alignments[k] = copy;
                        copy = NULL;
                }
        }

        STOP_TIMER(t_phase);
        if(!msa->quiet){ LOG_MSG("  [time] alignment runs + POAR extraction: "); GET_TIMING(t_phase); }

        /* Score all alignments and select the best */
        START_TIMER(t_phase);
        #ifdef USE_THREADPOOL
        RUN(score_alignments(alignments, poar, numseq, n_runs, msa->quiet,
                             &scores, &best_k, msa->pool));
#else
        RUN(score_alignments(alignments, poar, numseq, n_runs, msa->quiet,
                             &scores, &best_k));
#endif

        if(!msa->quiet){
                LOG_MSG("  Selected run %d (score=%.1f)", best_k + 1, scores[best_k]);
        }

        STOP_TIMER(t_phase);
        if(!msa->quiet){ LOG_MSG("  [time] scoring:                         "); GET_TIMING(t_phase); }

        /* Determine merge strategy */
        START_TIMER(t_phase);
        int min_support = (ens != NULL) ? ens->min_support : 0;
        int use_consistency_merge = (ens != NULL) ? ens->consistency_merge : 0;

        if(use_consistency_merge){
                /* ---- POAR consistency merge path ----
                 * Use the already-built POAR table as a source of pairwise
                 * residue consistency scores for a fresh progressive alignment.
                 * Uses best_k's gap penalties and matrix. */
                float cm_weight = (ens != NULL) ? ens->consistency_merge_weight : 2.0f;
                struct poar_consistency_ctx poar_ctx;
                poar_ctx.poar = poar;
                poar_ctx.n_runs = n_runs;
                poar_ctx.weight = cm_weight;

                copy = NULL;
                RUN(msa_cpy(&copy, msa));
                copy->quiet = msa->quiet ? 1 : 0;

                /* Attach POAR consistency context — the progressive alignment
                   will pick it up via msa->poar_consistency in aln_run.c */
                copy->poar_consistency = &poar_ctx;

                if(!msa->quiet){
                        LOG_MSG("  Consistency merge (weight=%.1f) using run %d params",
                                cm_weight, best_k + 1);
                }

                /* Run a fresh progressive alignment with best_k's params.
                   No additional anchor consistency or realign — the POAR
                   consistency signal is the main guide. */
                {
                        struct kalign_run_config cm_cfg = runs[best_k];
                        cm_cfg.refine = KALIGN_REFINE_NONE;
                        cm_cfg.tree_seed = 0;
                        cm_cfg.tree_noise = 0.0f;
                        cm_cfg.seq_weights = 0.0f;
                        cm_cfg.consistency_anchors = 0;
                        cm_cfg.realign = 0;
#ifdef USE_THREADPOOL
                        copy->pool = msa->pool;
#endif
                        RUN(kalign_single_run(copy, &cm_cfg, n_threads));
                }

                /* Clear the non-owning pointer before freeing the copy */
                copy->poar_consistency = NULL;

                RUN(copy_alignment_to_msa(msa, copy, numseq));
                kalign_free_msa(copy);
                copy = NULL;

        }else{
                /* ---- POAR consensus / selection path (existing) ---- */

                if(min_support > 0){
                        #ifdef USE_THREADPOOL
                        RUN(build_consensus_from_poar(poar, msa, numseq, min_support,
                                                      &consensus_msa, msa->pool));
                        #else
                        RUN(build_consensus_from_poar(poar, msa, numseq, min_support,
                                                      &consensus_msa));
                        #endif
                        use_consensus = 1;
                        if(!msa->quiet){
                                LOG_MSG("  Using consensus alignment (min_support=%d)", min_support);
                        }
                }else{
                        double consensus_score = 0.0;
                        int min_sup = (n_runs + 2) / 3;
                        if(min_sup < 2) min_sup = 2;

                        #ifdef USE_THREADPOOL
                        RUN(build_consensus_from_poar(poar, msa, numseq, min_sup,
                                                      &consensus_msa, msa->pool));
                        #else
                        RUN(build_consensus_from_poar(poar, msa, numseq, min_sup,
                                                      &consensus_msa));
                        #endif

                        #ifdef USE_THREADPOOL
                        RUN(score_single_msa(consensus_msa, poar, numseq, n_runs,
                                             &consensus_score, msa->pool));
                        #else
                        RUN(score_single_msa(consensus_msa, poar, numseq, n_runs,
                                             &consensus_score));
                        #endif

                        if(!msa->quiet){
                                LOG_MSG("  Consensus score: %.1f (selection: %.1f)",
                                        consensus_score, scores[best_k]);
                        }

                        if(consensus_score > scores[best_k]){
                                use_consensus = 1;
                                if(!msa->quiet){
                                        LOG_MSG("  Using consensus alignment");
                                }
                        }else{
                                kalign_free_msa(consensus_msa);
                                consensus_msa = NULL;
                                if(!msa->quiet){
                                        LOG_MSG("  Keeping selection winner");
                                }
                        }
                }

                /* Post-selection refinement: re-run the winner with REFINE_CONFIDENT */
                if(!use_consensus){
                        struct kalign_run_config ref_cfg = runs[best_k];
                        ref_cfg.refine = KALIGN_REFINE_CONFIDENT;

                        copy = NULL;
                        RUN(msa_cpy(&copy, msa));
                        copy->quiet = 1;
#ifdef USE_THREADPOOL
                        copy->pool = msa->pool;
#endif

                        if(!msa->quiet){
                                LOG_MSG("  Refining run %d...", best_k + 1);
                        }

                        RUN(kalign_single_run(copy, &ref_cfg, n_threads));

                        double refined_score = 0.0;
                        #ifdef USE_THREADPOOL
                        RUN(score_single_msa(copy, poar, numseq, n_runs,
                                             &refined_score, msa->pool));
                        #else
                        RUN(score_single_msa(copy, poar, numseq, n_runs,
                                             &refined_score));
                        #endif

                        if(!msa->quiet){
                                LOG_MSG("  Refined score: %.1f (was %.1f)",
                                        refined_score, scores[best_k]);
                        }

                        if(refined_score > scores[best_k]){
                                kalign_free_msa(alignments[best_k]);
                                alignments[best_k] = copy;
                                copy = NULL;
                                if(!msa->quiet){
                                        LOG_MSG("  Using refined alignment");
                                }
                        }else{
                                kalign_free_msa(copy);
                                copy = NULL;
                                if(!msa->quiet){
                                        LOG_MSG("  Keeping original alignment");
                                }
                        }
                }

                MFREE(scores);
                scores = NULL;

                if(use_consensus){
                        RUN(copy_alignment_to_msa(msa, consensus_msa, numseq));
                        kalign_free_msa(consensus_msa);
                        consensus_msa = NULL;
                }else{
                        RUN(copy_alignment_to_msa(msa, alignments[best_k], numseq));
                }
        }

        if(scores){
                MFREE(scores);
                scores = NULL;
        }

        STOP_TIMER(t_phase);
        if(!msa->quiet){ LOG_MSG("  [time] consensus/selection:             "); GET_TIMING(t_phase); }

        START_TIMER(t_phase);
        #ifdef USE_THREADPOOL
        RUN(compute_residue_confidence(poar, msa, msa->pool));
#else
        RUN(compute_residue_confidence(poar, msa));
#endif
        STOP_TIMER(t_phase);
        if(!msa->quiet){ LOG_MSG("  [time] confidence:                      "); GET_TIMING(t_phase); }

        RUN(msa_sort_rank(msa));

        STOP_TIMER(t_ensemble);
        if(!msa->quiet){
                GET_TIMING(t_ensemble);
        }
        DESTROY_TIMER(t_ensemble);
        DESTROY_TIMER(t_phase);

        for(k = 0; k < n_runs; k++){
                if(alignments[k]) kalign_free_msa(alignments[k]);
        }
        MFREE(alignments);
        poar_table_free(poar);
        return OK;
ERROR:
        if(copy) kalign_free_msa(copy);
        if(consensus_msa) kalign_free_msa(consensus_msa);
        if(pm) pos_matrix_free(pm);
        if(alignments){
                for(k = 0; k < n_runs; k++){
                        if(alignments[k]) kalign_free_msa(alignments[k]);
                }
                MFREE(alignments);
        }
        poar_table_free(poar);
        if(scores) MFREE(scores);
        return FAIL;
}

/* ======================================================================== */

int kalign_consensus_from_poar(struct msa* msa,
                               const char* poar_path,
                               int min_support)
{
        struct msa* consensus_msa = NULL;
        struct poar_table* poar = NULL;
        int numseq;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(poar_path != NULL, "No POAR file path");
        ASSERT(min_support >= 1, "min_support must be >= 1");

        RUN(kalign_essential_input_check(msa, 0));
        numseq = msa->numseq;

        /* Read POAR table from file */
        RUN(poar_table_read(&poar, poar_path));

        if(poar->numseq != numseq){
                ERROR_MSG("POAR file has %d sequences, input has %d",
                          poar->numseq, numseq);
        }

        /* Build consensus at given min_support threshold */
        #ifdef USE_THREADPOOL
                        RUN(build_consensus_from_poar(poar, msa, numseq, min_support,
                                                      &consensus_msa, msa->pool));
#else
                        #ifdef USE_THREADPOOL
                        RUN(build_consensus_from_poar(poar, msa, numseq, min_support,
                                                      &consensus_msa, msa->pool));
                        #else
                        RUN(build_consensus_from_poar(poar, msa, numseq, min_support,
                                                      &consensus_msa));
                        #endif
#endif

        /* Copy consensus alignment back into original MSA */
        RUN(copy_alignment_to_msa(msa, consensus_msa, numseq));
        kalign_free_msa(consensus_msa);
        consensus_msa = NULL;

        /* Compute per-residue and per-column confidence */
        #ifdef USE_THREADPOOL
        RUN(compute_residue_confidence(poar, msa, msa->pool));
#else
        RUN(compute_residue_confidence(poar, msa));
#endif

        RUN(msa_sort_rank(msa));

        poar_table_free(poar);
        return OK;
ERROR:
        if(consensus_msa) kalign_free_msa(consensus_msa);
        if(poar) poar_table_free(poar);
        return FAIL;
}
