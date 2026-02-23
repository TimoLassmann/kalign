#include <math.h>
#include <string.h>

#include "tldevel.h"
#include "msa_struct.h"
#include "task.h"
#include "weave_alignment.h"

#include "probmsa_submat.h"
#include "probmsa_dirichlet.h"
#include "probmsa_profile.h"
#include "probmsa_pairhmm.h"

#define PROBMSA_IMPORT
#include "probmsa.h"

/* Mapping tables: diri-20 (ACDEFGHIKLMNPQRSTVWY) ↔ kalign-23 (ARNDCQEGHILKMFPSTWYVBZX) */
static const int diri_to_kalign[DIRI_ALPHA_SIZE] = {
         0, /* A(0) → A(0)  */
         4, /* C(1) → C(4)  */
         3, /* D(2) → D(3)  */
         6, /* E(3) → E(6)  */
        13, /* F(4) → F(13) */
         7, /* G(5) → G(7)  */
         8, /* H(6) → H(8)  */
         9, /* I(7) → I(9)  */
        11, /* K(8) → K(11) */
        10, /* L(9) → L(10) */
        12, /* M(10) → M(12)*/
         2, /* N(11) → N(2) */
        14, /* P(12) → P(14)*/
         5, /* Q(13) → Q(5) */
         1, /* R(14) → R(1) */
        15, /* S(15) → S(15)*/
        16, /* T(16) → T(16)*/
        19, /* V(17) → V(19)*/
        17, /* W(18) → W(17)*/
        18, /* Y(19) → Y(18)*/
};

/* Default parameters */
static const struct probmsa_params default_params = {
        .delta          = 0.07,
        .epsilon        = 0.95,
        .tau            = 0.2,    /* semi-global: cheaper leading/trailing gaps for extensions */
        .mea_threshold  = 0.10,
        .mea_gpo        = 0.05,
        .mea_gpe        = 0.0,
        .conc_scale     = 20.0,
        .prior_scale    = 0.003,
        .use_5state     = 0,
        .delta_s        = 0.10,   /* short gap-open (frequent) */
        .epsilon_s      = 0.50,   /* short gap-extend → mean length 2 */
        .delta_l        = 0.01,   /* long gap-open (rare) */
        .epsilon_l      = 0.95,   /* long gap-extend → mean length 20 */
};

/* Compute Dirichlet-derived joint matrix in kalign-23 order.
 * Uses probmsa_diri_joint_matrix (diri-20 order) and remaps. */
static int compute_diri_joint(const struct probmsa_diri_mix* diri,
                               double joint[][PROBMSA_ALPHA_MAX],
                               double* background,
                               int alpha_size)
{
        double dj[DIRI_ALPHA_SIZE][DIRI_ALPHA_SIZE];
        double dbg[DIRI_ALPHA_SIZE];
        int s, t;

        RUN(probmsa_diri_joint_matrix(diri, dj, dbg));

        /* Zero kalign-size arrays */
        for(s = 0; s < PROBMSA_ALPHA_MAX; s++){
                background[s] = 0.0;
                for(t = 0; t < PROBMSA_ALPHA_MAX; t++){
                        joint[s][t] = 0.0;
                }
        }

        /* Remap diri-20 → kalign-23 */
        for(s = 0; s < DIRI_ALPHA_SIZE; s++){
                int ks = diri_to_kalign[s];
                background[ks] = dbg[s];
                for(t = 0; t < DIRI_ALPHA_SIZE; t++){
                        int kt = diri_to_kalign[t];
                        joint[ks][kt] = dj[s][t];
                }
        }

        /* B(20) = avg(N,D) = avg(kalign 2, 3) */
        if(alpha_size > 20){
                background[20] = 0.5 * (background[2] + background[3]);
                for(t = 0; t < alpha_size; t++){
                        joint[20][t] = 0.5 * (joint[2][t] + joint[3][t]);
                        joint[t][20] = 0.5 * (joint[t][2] + joint[t][3]);
                }
        }
        /* Z(21) = avg(Q,E) = avg(kalign 5, 6) */
        if(alpha_size > 21){
                background[21] = 0.5 * (background[5] + background[6]);
                for(t = 0; t < alpha_size; t++){
                        joint[21][t] = 0.5 * (joint[5][t] + joint[6][t]);
                        joint[t][21] = 0.5 * (joint[t][5] + joint[t][6]);
                }
        }
        /* X(22) = uniform */
        if(alpha_size > 22){
                double u = 1.0 / (double)DIRI_ALPHA_SIZE;
                background[22] = u;
                for(t = 0; t < alpha_size; t++){
                        joint[22][t] = u * background[t];
                        joint[t][22] = background[t] * u;
                }
        }

        return OK;
ERROR:
        return FAIL;
}

int create_msa_tree_probmsa(struct msa* msa, struct aln_tasks* t,
                              int biotype, float** subm, int alpha_size,
                              const struct probmsa_params* params)
{
        (void)biotype;  /* alpha_size already encodes the relevant info */
        struct probmsa_profile** profiles = NULL;
        struct probmsa_dp* dp = NULL;
        struct probmsa_diri_mix* diri = NULL;
        int i, j, g, c;
        int n_profiles;
        int max_len;
        int use_dirichlet;

        if(!params) params = &default_params;

        ASSERT(msa != NULL, "NULL msa");
        ASSERT(t != NULL, "NULL tasks");

        RUN(sort_tasks(t, TASK_ORDER_TREE));

        n_profiles = msa->num_profiles;

        /* Initialize Dirichlet prior for protein (alpha_size > 5).
         * For DNA/RNA, use legacy joint-matrix emission. */
        use_dirichlet = (alpha_size > 5) ? 1 : 0;
        if(use_dirichlet){
                RUN(probmsa_diri_init(&diri, params->prior_scale));
        }

        /* Allocate profile array */
        MMALLOC(profiles, sizeof(struct probmsa_profile*) * n_profiles);
        for(i = 0; i < n_profiles; i++){
                profiles[i] = NULL;
        }

        /* Create leaf profiles from sequences */
        for(i = 0; i < msa->numseq; i++){
                RUN(probmsa_profile_from_seq(msa->sequences[i]->s,
                                              msa->sequences[i]->len,
                                              alpha_size, diri, &profiles[i]));
        }

        /* Find max sequence length for DP allocation */
        max_len = 0;
        for(i = 0; i < msa->numseq; i++){
                if(msa->sequences[i]->len > max_len){
                        max_len = msa->sequences[i]->len;
                }
        }
        max_len = max_len * 2 + 100;

        /* Allocate DP workspace */
        if(params->use_5state){
                RUN(probmsa_dp_alloc_5state(&dp, max_len, max_len, alpha_size));
                RUN(probmsa_dp_set_transitions_5state(dp,
                        params->delta_s, params->epsilon_s,
                        params->delta_l, params->epsilon_l));
        }else{
                RUN(probmsa_dp_alloc(&dp, max_len, max_len, alpha_size));
                RUN(probmsa_dp_set_transitions(dp, params->delta, params->epsilon, params->tau));
        }

        /* Set MEA parameters */
        dp->mea_threshold = params->mea_threshold;
        dp->mea_gpo = params->mea_gpo;
        dp->mea_gpe = params->mea_gpe;

        if(use_dirichlet){
                /* Dirichlet mode: Dirichlet profiles for uncertainty,
                 * BLOSUM-derived joint matrix for substitution rates,
                 * position-specific gaps from concentration. */
                dp->use_dotproduct = 0;
                dp->conc_scale = params->use_5state ? 0.0 : params->conc_scale;

                RUN(probmsa_submat_from_kalign(subm, alpha_size,
                                                dp->joint, dp->background));
        }else{
                /* Legacy mode: joint matrix emission, global gaps */
                dp->use_dotproduct = 0;
                dp->conc_scale = 0.0;

                RUN(probmsa_submat_from_kalign(subm, alpha_size,
                                                dp->joint, dp->background));
        }

        /* Walk task tree sequentially (bottom-up) */
        for(c = 0; c < t->n_tasks; c++){
                int a_idx = t->list[c]->a;
                int b_idx = t->list[c]->b;
                int c_idx = t->list[c]->c;
                struct probmsa_profile* prof_a;
                struct probmsa_profile* prof_b;
                struct probmsa_profile* prof_c = NULL;
                int* path = NULL;
                double f_score, b_score;
                int len_a, len_b;

                prof_a = profiles[a_idx];
                prof_b = profiles[b_idx];

                ASSERT(prof_a != NULL, "NULL profile for node %d", a_idx);
                ASSERT(prof_b != NULL, "NULL profile for node %d", b_idx);

                len_a = prof_a->length;
                len_b = prof_b->length;

                /* Reallocate DP if needed */
                if(len_a > dp->alloc_a || len_b > dp->alloc_b){
                        int new_a = (len_a > dp->alloc_a) ? len_a + 100 : dp->alloc_a;
                        int new_b = (len_b > dp->alloc_b) ? len_b + 100 : dp->alloc_b;
                        double mea_t = dp->mea_threshold;
                        double mea_o = dp->mea_gpo;
                        double mea_e = dp->mea_gpe;
                        int dotprod = dp->use_dotproduct;
                        double cs = dp->conc_scale;
                        int is5 = dp->use_5state;
                        probmsa_dp_free(dp);
                        dp = NULL;
                        if(is5){
                                RUN(probmsa_dp_alloc_5state(&dp, new_a, new_b, alpha_size));
                                RUN(probmsa_dp_set_transitions_5state(dp,
                                        params->delta_s, params->epsilon_s,
                                        params->delta_l, params->epsilon_l));
                        }else{
                                RUN(probmsa_dp_alloc(&dp, new_a, new_b, alpha_size));
                                RUN(probmsa_dp_set_transitions(dp, params->delta, params->epsilon, params->tau));
                        }
                        dp->mea_threshold = mea_t;
                        dp->mea_gpo = mea_o;
                        dp->mea_gpe = mea_e;
                        dp->use_dotproduct = dotprod;
                        dp->conc_scale = cs;
                        if(use_dirichlet){
                                RUN(compute_diri_joint(diri, dp->joint,
                                                        dp->background, alpha_size));
                        }else{
                                RUN(probmsa_submat_from_kalign(subm, alpha_size,
                                                                dp->joint, dp->background));
                        }
                }

                /* Forward/Backward/Posterior: dispatch by model */
                if(params->use_5state){
                        RUN(probmsa_forward_5state(dp, prof_a, prof_b, &f_score));
                        RUN(probmsa_backward_5state(dp, prof_a, prof_b, &b_score));
                }else{
                        RUN(probmsa_forward(dp, prof_a, prof_b, &f_score));
                        RUN(probmsa_backward(dp, prof_a, prof_b, &b_score));
                }

                /* Consistency check */
                if(fabs(f_score - b_score) > 0.01 * (fabs(f_score) + 1.0)){
                        LOG_MSG("ProbMSA: forward/backward mismatch: f=%f b=%f (diff=%e)",
                                f_score, b_score, fabs(f_score - b_score));
                }

                /* Posterior */
                if(params->use_5state){
                        RUN(probmsa_posterior_5state(dp, prof_a, prof_b, len_a, len_b, f_score));
                }else{
                        RUN(probmsa_posterior(dp, prof_a, prof_b, len_a, len_b, f_score));
                }

                /* MEA traceback -> kalign path */
                RUN(probmsa_mea_traceback(dp, len_a, len_b, &path));

                /* Compute marginal gap posteriors from match posterior matrix */
                {
                        double* gpa = NULL;
                        double* gpb = NULL;
                        int ii, jj;

                        MMALLOC(gpa, sizeof(double) * len_a);
                        MMALLOC(gpb, sizeof(double) * len_b);

                        for(ii = 0; ii < len_a; ii++){
                                double msum = 0.0;
                                for(jj = 1; jj <= len_b; jj++){
                                        msum += dp->Y[ii + 1][jj];
                                }
                                gpa[ii] = 1.0 - msum;
                                if(gpa[ii] < 0.0) gpa[ii] = 0.0;
                        }
                        for(jj = 0; jj < len_b; jj++){
                                double msum = 0.0;
                                for(ii = 1; ii <= len_a; ii++){
                                        msum += dp->Y[ii][jj + 1];
                                }
                                gpb[jj] = 1.0 - msum;
                                if(gpb[jj] < 0.0) gpb[jj] = 0.0;
                        }

                        /* Merge profiles with Dirichlet regularization */
                        RUN(probmsa_profile_merge(prof_a, prof_b, path,
                                                   dp->Y,  /* match posteriors */
                                                   gpa, gpb,
                                                   alpha_size, diri, &prof_c));

                        MFREE(gpa);
                        MFREE(gpb);
                }

                /* Update kalign gap tracking */
                RUN(make_seq(msa, a_idx, b_idx, path));

                /* Update kalign profile metadata */
                msa->plen[c_idx] = path[0];

                msa->nsip[c_idx] = msa->nsip[a_idx] + msa->nsip[b_idx];
                MREALLOC(msa->sip[c_idx],
                         sizeof(int) * (msa->nsip[a_idx] + msa->nsip[b_idx]));

                g = 0;
                for(j = msa->nsip[a_idx]; j--;){
                        msa->sip[c_idx][g] = msa->sip[a_idx][j];
                        g++;
                }
                for(j = msa->nsip[b_idx]; j--;){
                        msa->sip[c_idx][g] = msa->sip[b_idx][j];
                        g++;
                }

                /* Store confidence */
                t->list[c]->confidence = 0.5f;

                /* Free old profiles, store new */
                probmsa_profile_free(profiles[a_idx]);
                profiles[a_idx] = NULL;
                probmsa_profile_free(profiles[b_idx]);
                profiles[b_idx] = NULL;
                profiles[c_idx] = prof_c;

                MFREE(path);
        }

        /* Cleanup */
        for(i = 0; i < n_profiles; i++){
                probmsa_profile_free(profiles[i]);
        }
        MFREE(profiles);
        probmsa_dp_free(dp);
        probmsa_diri_free(diri);

        return OK;
ERROR:
        if(profiles){
                for(i = 0; i < n_profiles; i++){
                        probmsa_profile_free(profiles[i]);
                }
                MFREE(profiles);
        }
        probmsa_dp_free(dp);
        probmsa_diri_free(diri);
        return FAIL;
}
