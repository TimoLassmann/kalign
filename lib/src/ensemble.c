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

#include "kalign/kalign.h"

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
 * Helper: resolve_run_params
 *
 * Given base gap penalties and a run index, compute the run-specific
 * gap-open, gap-extend, terminal-gap-extend, seed, and noise values.
 * Run 0 always uses defaults (deterministic, no noise).
 * ------------------------------------------------------------------------- */
static void resolve_run_params(float base_gpo, float base_gpe, float base_tgpe,
                               int k, uint64_t seed,
                               float* out_gpo, float* out_gpe, float* out_tgpe,
                               uint64_t* out_seed, float* out_noise)
{
        if(k == 0){
                /* Run 0: default params, deterministic */
                *out_gpo = base_gpo;
                *out_gpe = base_gpe;
                *out_tgpe = base_tgpe;
                *out_seed = 0;
                *out_noise = 0.0f;
        }else{
                /* Subsequent runs: independent per-penalty scaling + tree noise */
                struct ensemble_params ep = run_params[k % N_RUN_PARAMS];
                *out_gpo = base_gpo * ep.gpo_scale;
                *out_gpe = base_gpe * ep.gpe_scale;
                *out_tgpe = base_tgpe * ep.tgpe_scale;
                *out_seed = seed + (uint64_t)k;
                *out_noise = ep.noise;
        }
}

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
                            double** out_scores, int* out_best_k)
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
                RUN(score_alignment_poar(poar, pm, numseq, n_runs, &scores[k]));

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
                                     struct msa** out_consensus)
{
        struct msa* consensus_msa = NULL;
        int* seq_lens = NULL;

        RUN(msa_cpy(&consensus_msa, msa));

        MMALLOC(seq_lens, sizeof(int) * numseq);
        for(int i = 0; i < numseq; i++){
                seq_lens[i] = msa->sequences[i]->len;
        }

        RUN(build_consensus(poar, seq_lens, numseq, min_support, consensus_msa));
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
                            int numseq, int n_runs, double* out_score)
{
        struct pos_matrix* pm = NULL;
        char** aln_seqs = NULL;

        MMALLOC(aln_seqs, sizeof(char*) * numseq);
        for(int i = 0; i < numseq; i++){
                aln_seqs[i] = aln->sequences[i]->seq;
        }

        RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, aln->alnlen));
        RUN(score_alignment_poar(poar, pm, numseq, n_runs, out_score));

        pos_matrix_free(pm);
        MFREE(aln_seqs);
        return OK;
ERROR:
        if(pm) pos_matrix_free(pm);
        if(aln_seqs) MFREE(aln_seqs);
        return FAIL;
}

/* ======================================================================== */

int kalign_ensemble(struct msa* msa, int n_threads, int type,
                    int n_runs, float gpo, float gpe, float tgpe,
                    uint64_t seed, int min_support,
                    const char* save_poar_path,
                    int refine, float dist_scale, float vsm_amax,
                    int realign, float use_seq_weights,
                    int consistency_anchors, float consistency_weight)
{
        struct msa* copy = NULL;
        struct msa* consensus_msa = NULL;
        struct msa** alignments = NULL;
        struct poar_table* poar = NULL;
        struct pos_matrix* pm = NULL;
        struct aln_param* ap = NULL;
        double* scores = NULL;
        int numseq;
        int k;
        int best_k = 0;
        int use_consensus = 0;
        float base_gpo, base_gpe, base_tgpe;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(n_runs >= 1, "n_runs must be >= 1");

        /* Seq_weights hurts ensemble performance (POAR consensus already
           handles profile imbalance).  Default to OFF in ensemble mode. */
        if(use_seq_weights < 0.0f){
                use_seq_weights = 0.0f;
        }

        /* Essential input check + detect alphabet */
        RUN(kalign_essential_input_check(msa, 0));

        numseq = msa->numseq;

        DECLARE_TIMER(t_ensemble);
        if(!msa->quiet){
                LOG_MSG("Ensemble alignment with %d runs", n_runs);
        }
        START_TIMER(t_ensemble);

        /* Resolve default gap penalties using aln_param_init.
           We need to detect biotype first. */
        if(msa->biotype == ALN_BIOTYPE_UNDEF){
                RUN(detect_alphabet(msa));
        }

        /* Use aln_param_init to resolve defaults */
        RUN(aln_param_init(&ap, msa->biotype, n_threads, type, gpo, gpe, tgpe));
        base_gpo = ap->gpo;
        base_gpe = ap->gpe;
        base_tgpe = ap->tgpe;
        aln_param_free(ap);
        ap = NULL;

        /* Allocate POAR table and array to store completed alignments */
        RUN(poar_table_alloc(&poar, numseq));
        MMALLOC(alignments, sizeof(struct msa*) * n_runs);
        for(k = 0; k < n_runs; k++){
                alignments[k] = NULL;
        }

        /* Run N alignments, extract POARs, and keep each alignment */
        for(k = 0; k < n_runs; k++){
                float run_gpo, run_gpe, run_tgpe, run_noise;
                uint64_t run_seed;

                resolve_run_params(base_gpo, base_gpe, base_tgpe, k, seed,
                                   &run_gpo, &run_gpe, &run_tgpe,
                                   &run_seed, &run_noise);

                /* Deep-copy MSA */
                copy = NULL;
                RUN(msa_cpy(&copy, msa));
                copy->quiet = 1;

                if(!msa->quiet){
                        LOG_MSG("  Run %d/%d (gpo=%.1f gpe=%.1f tgpe=%.1f noise=%.2f)",
                                k + 1, n_runs, run_gpo, run_gpe, run_tgpe, run_noise);
                }

                /* Run alignment */
                if(realign > 0){
                        RUN(kalign_run_realign(copy, n_threads, type,
                                              run_gpo, run_gpe, run_tgpe,
                                              refine, 0,
                                              dist_scale, vsm_amax,
                                              realign, use_seq_weights,
                                              consistency_anchors, consistency_weight));
                }else{
                        RUN(kalign_run_seeded(copy, n_threads, type,
                                              run_gpo, run_gpe, run_tgpe,
                                              refine, 0,
                                              run_seed, run_noise,
                                              dist_scale, vsm_amax,
                                              use_seq_weights,
                                              consistency_anchors, consistency_weight));
                }

                /* Extract POARs from the finalized alignment */
                char** aln_seqs = NULL;
                MMALLOC(aln_seqs, sizeof(char*) * numseq);
                for(int i = 0; i < numseq; i++){
                        aln_seqs[i] = copy->sequences[i]->seq;
                }

                RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, copy->alnlen));
                RUN(extract_poars(poar, pm, k));

                pos_matrix_free(pm);
                pm = NULL;
                MFREE(aln_seqs);

                /* Keep this alignment for scoring later */
                alignments[k] = copy;
                copy = NULL;
        }

        /* Score all alignments and select the best */
        RUN(score_alignments(alignments, poar, numseq, n_runs, msa->quiet,
                             &scores, &best_k));

        if(!msa->quiet){
                LOG_MSG("  Selected run %d (score=%.1f)", best_k + 1, scores[best_k]);
        }

        /* Save POAR table if requested */
        if(save_poar_path != NULL){
                RUN(poar_table_write(poar, save_poar_path));
                if(!msa->quiet){
                        LOG_MSG("  Saved POAR table to %s", save_poar_path);
                }
        }

        /* When min_support > 0: explicit consensus threshold, skip selection.
           When min_support == 0: auto behavior (selection vs consensus). */
        if(min_support > 0){
                /* Explicit consensus: force consensus path */
                RUN(build_consensus_from_poar(poar, msa, numseq, min_support,
                                              &consensus_msa));
                use_consensus = 1;
                if(!msa->quiet){
                        LOG_MSG("  Using consensus alignment (min_support=%d)", min_support);
                }
        }else{
                /* Try consensus approach: build a new alignment from POAR table.
                   This can combine correct pairs from multiple runs, potentially
                   outperforming any single run. */
                double consensus_score = 0.0;
                int min_sup = (n_runs + 2) / 3;
                if(min_sup < 2) min_sup = 2;

                RUN(build_consensus_from_poar(poar, msa, numseq, min_sup,
                                              &consensus_msa));

                /* Score consensus against POAR table */
                RUN(score_single_msa(consensus_msa, poar, numseq, n_runs,
                                     &consensus_score));

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

        /* Post-selection refinement: only when using selection (not consensus),
           re-run the winner with REFINE_CONFIDENT and keep if it scores higher. */
        if(!use_consensus){
                float ref_gpo, ref_gpe, ref_tgpe, ref_noise;
                uint64_t ref_seed;

                resolve_run_params(base_gpo, base_gpe, base_tgpe, best_k, seed,
                                   &ref_gpo, &ref_gpe, &ref_tgpe,
                                   &ref_seed, &ref_noise);

                copy = NULL;
                RUN(msa_cpy(&copy, msa));
                copy->quiet = 1;

                if(!msa->quiet){
                        LOG_MSG("  Refining run %d...", best_k + 1);
                }

                RUN(kalign_run_seeded(copy, n_threads, type,
                                      ref_gpo, ref_gpe, ref_tgpe,
                                      KALIGN_REFINE_CONFIDENT, 0,
                                      ref_seed, ref_noise,
                                      dist_scale, vsm_amax,
                                      use_seq_weights,
                                      consistency_anchors, consistency_weight));

                /* Score the refined alignment against the same POAR table */
                double refined_score = 0.0;
                RUN(score_single_msa(copy, poar, numseq, n_runs,
                                     &refined_score));

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

        /* Copy the winning alignment back into the original MSA */
        if(use_consensus){
                RUN(copy_alignment_to_msa(msa, consensus_msa, numseq));
                kalign_free_msa(consensus_msa);
                consensus_msa = NULL;
        }else{
                RUN(copy_alignment_to_msa(msa, alignments[best_k], numseq));
        }

        /* Compute per-residue and per-column confidence from POAR table */
        RUN(compute_residue_confidence(poar, msa));

        /* Sort back to original rank order */
        RUN(msa_sort_rank(msa));

        STOP_TIMER(t_ensemble);
        if(!msa->quiet){
                GET_TIMING(t_ensemble);
        }
        DESTROY_TIMER(t_ensemble);

        /* Free all alignments */
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
        if(ap) aln_param_free(ap);
        return FAIL;
}

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
        RUN(build_consensus_from_poar(poar, msa, numseq, min_support,
                                       &consensus_msa));

        /* Copy consensus alignment back into original MSA */
        RUN(copy_alignment_to_msa(msa, consensus_msa, numseq));
        kalign_free_msa(consensus_msa);
        consensus_msa = NULL;

        /* Compute per-residue and per-column confidence */
        RUN(compute_residue_confidence(poar, msa));

        RUN(msa_sort_rank(msa));

        poar_table_free(poar);
        return OK;
ERROR:
        if(consensus_msa) kalign_free_msa(consensus_msa);
        if(poar) poar_table_free(poar);
        return FAIL;
}
