#include <math.h>

#include "tldevel.h"
#include "probmsa_logmath.h"

#define PROBMSA_EMISSION_IMPORT
#include "probmsa_emission.h"

double probmsa_emit_match(const struct probmsa_profile* a, int col_a,
                           const struct probmsa_profile* b, int col_b,
                           const double joint[][PROBMSA_ALPHA_MAX],
                           int alpha_size)
{
        const struct probmsa_column* ca = &a->cols[col_a];
        const struct probmsa_column* cb = &b->cols[col_b];

        double sum = 0.0;
        int s;

        if(joint == NULL){
                /* Log-odds dot-product emission (Dirichlet mode):
                 * S = sum_s e_A[s] * e_B[s] / bg[s]
                 * Dividing by background upweights rare amino acid
                 * agreements (HHalign-style profile-profile score).
                 * Returns log(S); gap emission returns 0 in this mode. */

                /* Need background — stored in the profile's alpha_size range.
                 * We access it via the static rr_freq table. */
                extern const double rr_freq_protein[];
                const double* bg;
                double bg_uniform[PROBMSA_ALPHA_MAX];

                if(alpha_size > 5){
                        bg = rr_freq_protein;
                }else{
                        for(s = 0; s < alpha_size; s++){
                                bg_uniform[s] = 1.0 / (double)alpha_size;
                        }
                        bg = bg_uniform;
                }

                for(s = 0; s < alpha_size; s++){
                        if(ca->prob[s] <= 0.0) continue;
                        if(cb->prob[s] <= 0.0) continue;
                        if(bg[s] <= 0.0) continue;
                        sum += ca->prob[s] * cb->prob[s] / bg[s];
                }
        }else{
                /* Joint matrix emission (legacy/DNA mode):
                 * P(col_A, col_B | Match) = sum_s sum_t p_A[s] * p_B[t] * joint[s][t] */
                int t;
                for(s = 0; s < alpha_size; s++){
                        if(ca->prob[s] <= 0.0) continue;
                        for(t = 0; t < alpha_size; t++){
                                if(cb->prob[t] <= 0.0) continue;
                                if(joint[s][t] <= 0.0) continue;
                                sum += ca->prob[s] * cb->prob[t] * joint[s][t];
                        }
                }
        }

        if(sum <= 0.0){
                return -INFINITY;
        }
        return log(sum);
}

double probmsa_emit_gap(const struct probmsa_profile* p, int col,
                         const double* background, int alpha_size)
{
        /* In log-odds mode (background == NULL), gap emission is 0.0
         * (neutral in log-odds space). All gap penalty comes from
         * transition probabilities. */
        if(background == NULL){
                return 0.0;
        }

        const struct probmsa_column* c = &p->cols[col];

        double sum = 0.0;
        int s;

        for(s = 0; s < alpha_size; s++){
                if(c->prob[s] <= 0.0) continue;
                if(background[s] <= 0.0) continue;
                sum += c->prob[s] * background[s];
        }

        if(sum <= 0.0){
                return -INFINITY;
        }
        return log(sum);
}
