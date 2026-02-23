#include "tldevel.h"
#include "tlmisc.h"
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
#include "probmsa.h"
#include "kalign/kalign.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#define ALN_WRAP_IMPORT
#include "aln_wrap.h"

/* Resolve PFASUM_AUTO: pick PFASUM43 for divergent, PFASUM60 for closer.
   Must be called after build_tree_kmeans which fills msa->seq_distances[]. */
static int resolve_pfasum_auto(struct msa *msa, int *type)
{
        int i;
        int min_len, max_len;
        float len_ratio;

        if(*type != KALIGN_TYPE_PROTEIN_PFASUM_AUTO){
                return OK;
        }
        if(msa->biotype != ALN_BIOTYPE_PROTEIN){
                *type = KALIGN_TYPE_PROTEIN_PFASUM43;
                return OK;
        }

        /* Use sequence length ratio (max/min) to select matrix.
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
                *type = KALIGN_TYPE_PROTEIN_PFASUM43;
        }else{
                *type = KALIGN_TYPE_PROTEIN_PFASUM60;
        }
        if(!msa->quiet){
                LOG_MSG("Auto matrix: len_ratio=%.2f -> %s",
                        len_ratio,
                        *type == KALIGN_TYPE_PROTEIN_PFASUM60 ? "PFASUM60" : "PFASUM43");
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
        RUN(resolve_pfasum_auto(msa, &type));

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

        RUN(resolve_pfasum_auto(msa, &type));

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

        RUN(resolve_pfasum_auto(msa, &type));

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
        RUN(resolve_pfasum_auto(msa, &type));

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

int kalign_post_realign_probmsa(struct msa *msa, int n_threads, int type,
                                double hmm_delta, double hmm_epsilon,
                                double hmm_tau,
                                double mea_threshold, double mea_gpo, double mea_gpe,
                                double prior_scale)
{
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;
        struct probmsa_params pm;
        float** dm = NULL;
        int si;

        pm.delta         = (hmm_delta > 0.0)     ? hmm_delta     : 0.07;
        pm.epsilon       = (hmm_epsilon > 0.0)    ? hmm_epsilon   : 0.95;
        pm.tau           = (hmm_tau > 0.0)        ? hmm_tau       : 0.2;
        pm.mea_threshold = (mea_threshold >= 0.0) ? mea_threshold : 0.10;
        pm.mea_gpo       = (mea_gpo >= 0.0)       ? mea_gpo       : 0.05;
        pm.mea_gpe       = (mea_gpe >= 0.0)       ? mea_gpe       : 0.0;
        pm.conc_scale    = 20.0;
        pm.prior_scale   = (prior_scale > 0.0)   ? prior_scale   : 0.003;

        ASSERT(msa != NULL, "No MSA");

        /* Detect biotype if not set */
        if(msa->biotype == ALN_BIOTYPE_UNDEF){
                RUN(detect_alphabet(msa));
        }

#ifdef HAVE_OPENMP
        omp_set_num_threads(n_threads);
#endif

        /* The input MSA should already be aligned (e.g. a true/reference alignment).
           Finalize if not already done so we can compute pairwise distances. */
        if(msa->aligned != ALN_STATUS_FINAL){
                /* If aligned but not finalized, finalize now */
                if(msa->aligned == ALN_STATUS_ALIGNED){
                        RUN(finalise_alignment(msa));
                }else{
                        /* Try to detect alignment status from gap characters */
                        int has_gaps = 0;
                        for(si = 0; si < msa->numseq && !has_gaps; si++){
                                for(int r = 0; msa->sequences[si]->seq[r] != '\0'; r++){
                                        if(msa->sequences[si]->seq[r] == '-'){
                                                has_gaps = 1;
                                                break;
                                        }
                                }
                        }
                        if(has_gaps){
                                msa->aligned = ALN_STATUS_FINAL;
                        }else{
                                ERROR_MSG("Input MSA for guided ProbMSA must be already aligned");
                        }
                }
        }

        /* Compute NxN pairwise identity distances from the alignment */
        RUN(compute_aln_pairwise_dist(msa, &dm));

        /* Strip gap characters from seq->seq and fix seq->len.
           Assign ranks to preserve input order for msa_sort_rank(). */
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
                seq->rank = si;
                for(r = 0; r <= w; r++){
                        seq->gaps[r] = 0;
                }
        }
        msa->aligned = ALN_STATUS_UNALIGNED;

        /* Do NOT sort — distance matrix indices must match sequence order.
           (Unlike kalign_run_probmsa which sorts for the k-means tree builder.) */

        /* Convert directly to full alphabet (no reduced alphabet needed
           since we're not using k-means for tree building). */
        if(msa->biotype == ALN_BIOTYPE_DNA){
                msa->L = ALPHA_defDNA;
                RUN(convert_msa_to_internal(msa, ALPHA_defDNA));
        }else if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                msa->L = ALPHA_ambigiousPROTEIN;
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }else{
                ERROR_MSG("Unable to determine what alphabet to use.");
        }

        /* Reset profile tracking */
        RUN(set_sip_nsip(msa));

        /* Build UPGMA tree from alignment-derived distances */
        RUN(alloc_tasks(&tasks, msa->numseq));
        RUN(build_tree_from_pairwise(msa, &tasks, dm));
        free_aln_dm(dm, msa->numseq);
        dm = NULL;

        /* Get substitution matrix */
        RUN(aln_param_init(&ap, msa->biotype, n_threads, type,
                           -1.0f, -1.0f, -1.0f));

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Aligning (ProbMSA guided by reference alignment)");
        }
        START_TIMER(t1);

        /* Run ProbMSA with the guide-derived tree */
        RUN(create_msa_tree_probmsa(msa, tasks, msa->biotype,
                                     ap->subm, (int)msa->L, &pm));
        msa->aligned = ALN_STATUS_ALIGNED;

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
        if(dm) free_aln_dm(dm, msa->numseq);
        aln_param_free(ap);
        if(tasks) free_tasks(tasks);
        return FAIL;
}

int kalign_run_probmsa(struct msa *msa, int n_threads, int type,
                       double hmm_delta, double hmm_epsilon,
                       double hmm_tau,
                       double mea_threshold, double mea_gpo, double mea_gpe,
                       double prior_scale,
                       int use_5state,
                       double delta_s, double epsilon_s,
                       double delta_l, double epsilon_l)
{
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;
        struct probmsa_params pm;

        /* Set up probmsa params with defaults for negative sentinels */
        pm.delta         = (hmm_delta > 0.0)     ? hmm_delta     : 0.07;
        pm.epsilon       = (hmm_epsilon > 0.0)    ? hmm_epsilon   : 0.95;
        pm.tau           = (hmm_tau > 0.0)        ? hmm_tau       : 0.2;
        pm.mea_threshold = (mea_threshold >= 0.0) ? mea_threshold : 0.10;
        pm.mea_gpo       = (mea_gpo >= 0.0)       ? mea_gpo       : 0.05;
        pm.mea_gpe       = (mea_gpe >= 0.0)       ? mea_gpe       : 0.0;
        pm.conc_scale    = 20.0;
        pm.prior_scale   = (prior_scale > 0.0)   ? prior_scale   : 0.003;

        /* 5-state parameters */
        pm.use_5state = use_5state;
        pm.delta_s    = (delta_s > 0.0)   ? delta_s   : 0.10;
        pm.epsilon_s  = (epsilon_s > 0.0) ? epsilon_s : 0.50;
        pm.delta_l    = (delta_l > 0.0)   ? delta_l   : 0.01;
        pm.epsilon_l  = (epsilon_l > 0.0) ? epsilon_l : 0.95;

        RUN(kalign_essential_input_check(msa, 0));

        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }
        RUN(msa_sort_len_name(msa));

        /* Convert to internal representation */
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

        /* Convert to full alphabet for alignment */
        if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }

        RUN(resolve_pfasum_auto(msa, &type));

        /* Get substitution matrix for joint probability conversion */
        RUN(aln_param_init(&ap, msa->biotype, n_threads, type,
                           -1.0f, -1.0f, -1.0f));

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Aligning (ProbMSA %s pair-HMM)",
                        use_5state ? "5-state" : "3-state");
        }
        START_TIMER(t1);

        /* Run ProbMSA progressive alignment */
        RUN(create_msa_tree_probmsa(msa, tasks, msa->biotype,
                                     ap->subm, (int)msa->L, &pm));
        msa->aligned = ALN_STATUS_ALIGNED;

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
