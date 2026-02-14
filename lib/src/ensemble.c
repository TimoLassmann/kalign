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

int kalign_ensemble(struct msa* msa, int n_threads, int type,
                    int n_runs, float gpo, float gpe, float tgpe,
                    uint64_t seed)
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
        int use_consensus = 0;
        float base_gpo, base_gpe, base_tgpe;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(n_runs >= 1, "n_runs must be >= 1");

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
                struct ensemble_params ep = run_params[k % N_RUN_PARAMS];
                float run_gpo, run_gpe, run_tgpe;
                uint64_t run_seed;
                float run_noise;

                if(k == 0){
                        /* Run 0: default params, deterministic */
                        run_gpo = base_gpo;
                        run_gpe = base_gpe;
                        run_tgpe = base_tgpe;
                        run_seed = 0;
                        run_noise = 0.0f;
                }else{
                        /* Subsequent runs: independent per-penalty scaling + tree noise */
                        run_gpo = base_gpo * ep.gpo_scale;
                        run_gpe = base_gpe * ep.gpe_scale;
                        run_tgpe = base_tgpe * ep.tgpe_scale;
                        run_seed = seed + (uint64_t)k;
                        run_noise = ep.noise;
                }

                /* Deep-copy MSA */
                copy = NULL;
                RUN(msa_cpy(&copy, msa));
                copy->quiet = 1;

                if(!msa->quiet){
                        LOG_MSG("  Run %d/%d (gpo=%.1f gpe=%.1f tgpe=%.1f noise=%.2f)",
                                k + 1, n_runs, run_gpo, run_gpe, run_tgpe, run_noise);
                }

                /* Run alignment */
                RUN(kalign_run_seeded(copy, n_threads, type,
                                      run_gpo, run_gpe, run_tgpe,
                                      KALIGN_REFINE_NONE, 0,
                                      run_seed, run_noise));

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

        /* Score each alignment against the POAR table.
           Use run 0 (default params) as baseline; only switch to
           another run if it scores strictly higher. */
        MMALLOC(scores, sizeof(double) * n_runs);

        for(k = 0; k < n_runs; k++){
                char** aln_seqs = NULL;

                MMALLOC(aln_seqs, sizeof(char*) * numseq);
                for(int i = 0; i < numseq; i++){
                        aln_seqs[i] = alignments[k]->sequences[i]->seq;
                }

                RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, alignments[k]->alnlen));
                RUN(score_alignment_poar(poar, pm, numseq, n_runs, &scores[k]));

                if(!msa->quiet){
                        LOG_MSG("  Run %d score: %.1f", k + 1, scores[k]);
                }

                pos_matrix_free(pm);
                pm = NULL;
                MFREE(aln_seqs);
        }

        /* Select: prefer run 0 (default params) unless another run
           has a meaningfully higher score (>2% improvement). */
        int best_k = 0;
        double baseline = scores[0];
        for(k = 1; k < n_runs; k++){
                if(scores[k] > scores[best_k] && scores[k] > baseline * 1.05){
                        best_k = k;
                }
        }

        if(!msa->quiet){
                LOG_MSG("  Selected run %d (score=%.1f)", best_k + 1, scores[best_k]);
        }

        /* Try consensus approach: build a new alignment from POAR table.
           This can combine correct pairs from multiple runs, potentially
           outperforming any single run. */
        {
                int* seq_lens = NULL;
                double consensus_score = 0.0;
                int min_sup = (n_runs + 2) / 3;
                if(min_sup < 2) min_sup = 2;

                RUN(msa_cpy(&consensus_msa, msa));

                MMALLOC(seq_lens, sizeof(int) * numseq);
                for(int i = 0; i < numseq; i++){
                        seq_lens[i] = msa->sequences[i]->len;
                }

                RUN(build_consensus(poar, seq_lens, numseq, min_sup,
                                    consensus_msa));
                MFREE(seq_lens);

                /* Score consensus against POAR table */
                {
                        char** aln_seqs = NULL;
                        MMALLOC(aln_seqs, sizeof(char*) * numseq);
                        for(int i = 0; i < numseq; i++){
                                aln_seqs[i] = consensus_msa->sequences[i]->seq;
                        }
                        RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq,
                                                consensus_msa->alnlen));
                        RUN(score_alignment_poar(poar, pm, numseq, n_runs,
                                                 &consensus_score));
                        pos_matrix_free(pm);
                        pm = NULL;
                        MFREE(aln_seqs);
                }

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
                struct ensemble_params ep = run_params[best_k % N_RUN_PARAMS];
                float ref_gpo, ref_gpe, ref_tgpe;
                uint64_t ref_seed;
                float ref_noise;

                if(best_k == 0){
                        ref_gpo = base_gpo;
                        ref_gpe = base_gpe;
                        ref_tgpe = base_tgpe;
                        ref_seed = 0;
                        ref_noise = 0.0f;
                }else{
                        ref_gpo = base_gpo * ep.gpo_scale;
                        ref_gpe = base_gpe * ep.gpe_scale;
                        ref_tgpe = base_tgpe * ep.tgpe_scale;
                        ref_seed = seed + (uint64_t)best_k;
                        ref_noise = ep.noise;
                }

                copy = NULL;
                RUN(msa_cpy(&copy, msa));
                copy->quiet = 1;

                if(!msa->quiet){
                        LOG_MSG("  Refining run %d...", best_k + 1);
                }

                RUN(kalign_run_seeded(copy, n_threads, type,
                                      ref_gpo, ref_gpe, ref_tgpe,
                                      KALIGN_REFINE_CONFIDENT, 0,
                                      ref_seed, ref_noise));

                /* Score the refined alignment against the same POAR table */
                double refined_score = 0.0;
                char** aln_seqs = NULL;
                MMALLOC(aln_seqs, sizeof(char*) * numseq);
                for(int i = 0; i < numseq; i++){
                        aln_seqs[i] = copy->sequences[i]->seq;
                }
                RUN(pos_matrix_from_msa(&pm, aln_seqs, numseq, copy->alnlen));
                RUN(score_alignment_poar(poar, pm, numseq, n_runs,
                                         &refined_score));
                pos_matrix_free(pm);
                pm = NULL;
                MFREE(aln_seqs);

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

        /* Copy the winning alignment back into the original MSA */
        if(use_consensus){
                for(int i = 0; i < numseq; i++){
                        MFREE(msa->sequences[i]->seq);
                        msa->sequences[i]->seq = consensus_msa->sequences[i]->seq;
                        msa->sequences[i]->len = consensus_msa->sequences[i]->len;
                        consensus_msa->sequences[i]->seq = NULL;
                        MMALLOC(consensus_msa->sequences[i]->seq, 1);
                        consensus_msa->sequences[i]->seq[0] = '\0';
                }
                msa->alnlen = consensus_msa->alnlen;
                msa->aligned = consensus_msa->aligned;
                kalign_free_msa(consensus_msa);
                consensus_msa = NULL;
        }else{
                struct msa* best = alignments[best_k];
                for(int i = 0; i < numseq; i++){
                        MFREE(msa->sequences[i]->seq);
                        msa->sequences[i]->seq = best->sequences[i]->seq;
                        msa->sequences[i]->len = best->sequences[i]->len;
                        best->sequences[i]->seq = NULL;
                        MMALLOC(best->sequences[i]->seq, 1);
                        best->sequences[i]->seq[0] = '\0';
                }
                msa->alnlen = best->alnlen;
                msa->aligned = best->aligned;
        }

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

int kalign_ensemble_consensus(struct msa* msa, int n_threads, int type,
                              int n_runs, float gpo, float gpe, float tgpe,
                              uint64_t seed, int min_support)
{
        struct msa* copy = NULL;
        struct msa* consensus_msa = NULL;
        struct msa** alignments = NULL;
        struct poar_table* poar = NULL;
        struct pos_matrix* pm = NULL;
        struct aln_param* ap = NULL;
        int numseq;
        int k;
        float base_gpo, base_gpe, base_tgpe;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(n_runs >= 1, "n_runs must be >= 1");
        ASSERT(min_support >= 1, "min_support must be >= 1");

        RUN(kalign_essential_input_check(msa, 0));
        numseq = msa->numseq;

        if(msa->biotype == ALN_BIOTYPE_UNDEF){
                RUN(detect_alphabet(msa));
        }

        RUN(aln_param_init(&ap, msa->biotype, n_threads, type, gpo, gpe, tgpe));
        base_gpo = ap->gpo;
        base_gpe = ap->gpe;
        base_tgpe = ap->tgpe;
        aln_param_free(ap);
        ap = NULL;

        RUN(poar_table_alloc(&poar, numseq));
        MMALLOC(alignments, sizeof(struct msa*) * n_runs);
        for(k = 0; k < n_runs; k++){
                alignments[k] = NULL;
        }

        /* Run N alignments and collect POARs */
        for(k = 0; k < n_runs; k++){
                struct ensemble_params ep = run_params[k % N_RUN_PARAMS];
                float run_gpo, run_gpe, run_tgpe;
                uint64_t run_seed;
                float run_noise;

                if(k == 0){
                        run_gpo = base_gpo;
                        run_gpe = base_gpe;
                        run_tgpe = base_tgpe;
                        run_seed = 0;
                        run_noise = 0.0f;
                }else{
                        run_gpo = base_gpo * ep.gpo_scale;
                        run_gpe = base_gpe * ep.gpe_scale;
                        run_tgpe = base_tgpe * ep.tgpe_scale;
                        run_seed = seed + (uint64_t)k;
                        run_noise = ep.noise;
                }

                copy = NULL;
                RUN(msa_cpy(&copy, msa));
                copy->quiet = 1;

                RUN(kalign_run_seeded(copy, n_threads, type,
                                      run_gpo, run_gpe, run_tgpe,
                                      KALIGN_REFINE_NONE, 0,
                                      run_seed, run_noise));

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

                alignments[k] = copy;
                copy = NULL;
        }

        /* Build consensus at given min_support threshold */
        RUN(msa_cpy(&consensus_msa, msa));

        int* seq_lens = NULL;
        MMALLOC(seq_lens, sizeof(int) * numseq);
        for(int i = 0; i < numseq; i++){
                seq_lens[i] = msa->sequences[i]->len;
        }

        RUN(build_consensus(poar, seq_lens, numseq, min_support, consensus_msa));
        MFREE(seq_lens);

        /* Copy consensus alignment back into original MSA */
        for(int i = 0; i < numseq; i++){
                MFREE(msa->sequences[i]->seq);
                msa->sequences[i]->seq = consensus_msa->sequences[i]->seq;
                msa->sequences[i]->len = consensus_msa->sequences[i]->len;
                consensus_msa->sequences[i]->seq = NULL;
                MMALLOC(consensus_msa->sequences[i]->seq, 1);
                consensus_msa->sequences[i]->seq[0] = '\0';
        }
        msa->alnlen = consensus_msa->alnlen;
        msa->aligned = consensus_msa->aligned;
        kalign_free_msa(consensus_msa);
        consensus_msa = NULL;

        RUN(msa_sort_rank(msa));

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
        if(ap) aln_param_free(ap);
        return FAIL;
}
