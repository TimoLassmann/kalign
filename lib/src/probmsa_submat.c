#include <math.h>
#include <float.h>

#include "tldevel.h"

#define PROBMSA_SUBMAT_IMPORT
#include "probmsa_submat.h"

/* Robinson-Robinson amino acid background frequencies
   in kalign's ARNDCQEGHILKMFPSTWYVBZX order */
const double rr_freq_protein[23] = {
        0.07805,  /* A */
        0.05129,  /* R */
        0.04487,  /* N */
        0.05364,  /* D */
        0.01925,  /* C */
        0.03926,  /* Q */
        0.06295,  /* E */
        0.07377,  /* G */
        0.02199,  /* H */
        0.05966,  /* I */
        0.09660,  /* L */
        0.05841,  /* K */
        0.02243,  /* M */
        0.03856,  /* F */
        0.04730,  /* P */
        0.06880,  /* S */
        0.05560,  /* T */
        0.01330,  /* W */
        0.02935,  /* Y */
        0.06726,  /* V */
        0.04926,  /* B = avg(N,D) */
        0.05111,  /* Z = avg(Q,E) */
        0.05000   /* X = uniform */
};

int probmsa_submat_from_kalign(float** subm, int alpha_size,
                                double joint[][PROBMSA_ALPHA_MAX],
                                double* background)
{
        double bg[PROBMSA_ALPHA_MAX];
        int i, j;
        double lambda;
        double lo, hi, mid;
        int iter;

        ASSERT(alpha_size >= 1 && alpha_size <= PROBMSA_ALPHA_MAX,
               "Invalid alphabet size: %d", alpha_size);

        /* Set background frequencies */
        if(alpha_size <= 5){
                /* DNA/RNA: uniform over 4 bases, small for N */
                double base = 0.2499;
                double n_freq = 0.0004;
                for(i = 0; i < alpha_size && i < 4; i++){
                        bg[i] = base;
                }
                if(alpha_size == 5){
                        bg[4] = n_freq;
                }
                /* Normalize */
                double sum = 0.0;
                for(i = 0; i < alpha_size; i++) sum += bg[i];
                for(i = 0; i < alpha_size; i++) bg[i] /= sum;
        }else{
                /* Protein: Robinson-Robinson frequencies */
                for(i = 0; i < alpha_size; i++){
                        bg[i] = rr_freq_protein[i];
                }
                /* Normalize */
                double sum = 0.0;
                for(i = 0; i < alpha_size; i++) sum += bg[i];
                for(i = 0; i < alpha_size; i++) bg[i] /= sum;
        }

        /* Find lambda by bisection such that
           sum_{s,t} bg[s] * bg[t] * exp(lambda * S[s][t]) = 1.0 */
        lo = 0.001;
        hi = 5.0;

        /* Check that lo and hi bracket the root */
        double sum_lo = 0.0, sum_hi = 0.0;
        for(i = 0; i < alpha_size; i++){
                for(j = 0; j < alpha_size; j++){
                        sum_lo += bg[i] * bg[j] * exp(lo * (double)subm[i][j]);
                        sum_hi += bg[i] * bg[j] * exp(hi * (double)subm[i][j]);
                }
        }

        /* If the range doesn't bracket 1.0, adjust */
        if(sum_lo > 1.0){
                lo = 0.0001;
        }
        if(sum_hi < 1.0){
                hi = 10.0;
        }

        /* Bisection: find lambda such that sum = 1.0 */
        lambda = 0.3;  /* reasonable starting point */
        for(iter = 0; iter < 200; iter++){
                mid = (lo + hi) / 2.0;
                double sum = 0.0;
                for(i = 0; i < alpha_size; i++){
                        for(j = 0; j < alpha_size; j++){
                                sum += bg[i] * bg[j] * exp(mid * (double)subm[i][j]);
                        }
                }
                if(sum > 1.0){
                        hi = mid;
                }else{
                        lo = mid;
                }
                if(fabs(sum - 1.0) < 1e-12){
                        break;
                }
                lambda = mid;
        }

        /* Compute joint probability matrix */
        double total = 0.0;
        for(i = 0; i < alpha_size; i++){
                for(j = 0; j < alpha_size; j++){
                        joint[i][j] = bg[i] * bg[j] * exp(lambda * (double)subm[i][j]);
                        total += joint[i][j];
                }
        }

        /* Final normalization to ensure exact sum = 1.0 */
        if(total > 0.0){
                for(i = 0; i < alpha_size; i++){
                        for(j = 0; j < alpha_size; j++){
                                joint[i][j] /= total;
                        }
                }
        }

        /* Zero out entries beyond alpha_size */
        for(i = 0; i < PROBMSA_ALPHA_MAX; i++){
                for(j = 0; j < PROBMSA_ALPHA_MAX; j++){
                        if(i >= alpha_size || j >= alpha_size){
                                joint[i][j] = 0.0;
                        }
                }
        }

        /* Store background frequencies */
        for(i = 0; i < PROBMSA_ALPHA_MAX; i++){
                if(i < alpha_size){
                        background[i] = bg[i];
                }else{
                        background[i] = 0.0;
                }
        }

        return OK;
ERROR:
        return FAIL;
}
