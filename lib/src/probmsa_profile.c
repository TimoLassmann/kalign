#include <math.h>
#include <string.h>

#include "tldevel.h"

#define PROBMSA_PROFILE_IMPORT
#include "probmsa_profile.h"

/* Alphabet remapping tables.
 * Kalign alignment alphabet (ALPHA_ambigiousPROTEIN): ARNDCQEGHILKMFPSTWYVBZX
 * Dirichlet alphabet (Sjölander):                     ACDEFGHIKLMNPQRSTVWY
 * These orderings are DIFFERENT and must be mapped for Dirichlet operations. */

/* kalign-23 index → diri-20 index (-1 = no direct mapping, e.g. B/Z/X) */
static const int kalign_to_diri[23] = {
         0, /* A(0) → A(0)  */
        14, /* R(1) → R(14) */
        11, /* N(2) → N(11) */
         2, /* D(3) → D(2)  */
         1, /* C(4) → C(1)  */
        13, /* Q(5) → Q(13) */
         3, /* E(6) → E(3)  */
         5, /* G(7) → G(5)  */
         6, /* H(8) → H(6)  */
         7, /* I(9) → I(7)  */
         9, /* L(10) → L(9) */
         8, /* K(11) → K(8) */
        10, /* M(12) → M(10)*/
         4, /* F(13) → F(4) */
        12, /* P(14) → P(12)*/
        15, /* S(15) → S(15)*/
        16, /* T(16) → T(16)*/
        18, /* W(17) → W(18)*/
        19, /* Y(18) → Y(19)*/
        17, /* V(19) → V(17)*/
        -1, /* B(20) */
        -1, /* Z(21) */
        -1, /* X(22) */
};

/* diri-20 index → kalign-23 index */
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

int probmsa_profile_from_seq(const uint8_t* seq, int len, int alpha_size,
                              const struct probmsa_diri_mix* diri,
                              struct probmsa_profile** out)
{
        struct probmsa_profile* p = NULL;
        int i, j;
        (void)diri;  /* Dirichlet regularization only used at merge time */

        ASSERT(seq != NULL, "NULL sequence");
        ASSERT(len > 0, "Zero-length sequence");

        MMALLOC(p, sizeof(struct probmsa_profile));
        p->cols = NULL;
        MMALLOC(p->cols, sizeof(struct probmsa_column) * len);

        for(i = 0; i < len; i++){
                int code = (int)seq[i];

                for(j = 0; j < PROBMSA_ALPHA_MAX; j++){
                        p->cols[i].prob[j] = 0.0;
                }

                if(code >= 0 && code < alpha_size){
                        /* Point estimate (DNA or ambiguous codes) */
                        p->cols[i].prob[code] = 1.0;
                        p->cols[i].concentration = 1.0;
                }else{
                        /* Unknown residue: uniform */
                        for(j = 0; j < alpha_size; j++){
                                p->cols[i].prob[j] = 1.0 / (double)alpha_size;
                        }
                        p->cols[i].concentration = 0.0;
                }
                p->cols[i].quality = 40.0;
                p->cols[i].depth = 1;
                p->cols[i].gap_frac = 0.0;
        }

        p->length = len;
        p->n_sequences = 1;
        p->alpha_size = alpha_size;

        *out = p;
        return OK;
ERROR:
        probmsa_profile_free(p);
        return FAIL;
}

int probmsa_profile_merge(struct probmsa_profile* a,
                           struct probmsa_profile* b,
                           const int* path,
                           double** posteriors,
                           const double* gap_post_a,
                           const double* gap_post_b,
                           int alpha_size,
                           const struct probmsa_diri_mix* diri,
                           struct probmsa_profile** out)
{
        struct probmsa_profile* p = NULL;
        int aln_len;
        int i, k, c;
        int pos_a, pos_b;

        aln_len = path[0];
        ASSERT(aln_len > 0, "Empty path");

        MMALLOC(p, sizeof(struct probmsa_profile));
        p->cols = NULL;
        MMALLOC(p->cols, sizeof(struct probmsa_column) * aln_len);

        pos_a = 0;  /* 0-indexed into profile a */
        pos_b = 0;  /* 0-indexed into profile b */

        c = 1;
        for(i = 0; i < aln_len; i++){
                int op = path[c];
                struct probmsa_column* col = &p->cols[i];

                if(op == 0){
                        /* Match: merge columns from a and b */
                        int da = a->cols[pos_a].depth;
                        int db = b->cols[pos_b].depth;

                        if(diri != NULL && (da + db) >= 2){
                                /* Dirichlet mixture regularization:
                                 * Combine effective counts from both profiles,
                                 * then compute posterior mean.
                                 * Profiles are in kalign-23 order; must remap
                                 * to diri-20 for the Dirichlet computation. */
                                double counts[DIRI_ALPHA_SIZE];
                                double diri_prob[DIRI_ALPHA_SIZE];
                                int ds = da + db;
                                for(k = 0; k < DIRI_ALPHA_SIZE; k++){
                                        counts[k] = 0.0;
                                }
                                /* Remap kalign-23 → diri-20 */
                                for(k = 0; k < 20 && k < alpha_size; k++){
                                        int dk = kalign_to_diri[k];
                                        if(dk >= 0){
                                                counts[dk] = a->cols[pos_a].prob[k] * da
                                                           + b->cols[pos_b].prob[k] * db;
                                        }
                                }
                                RUN(probmsa_diri_posterior_mean(diri, counts, (double)ds,
                                                                diri_prob, NULL));
                                /* Remap diri-20 → kalign-23 */
                                for(k = 0; k < PROBMSA_ALPHA_MAX; k++){
                                        col->prob[k] = 0.0;
                                }
                                for(k = 0; k < DIRI_ALPHA_SIZE; k++){
                                        col->prob[diri_to_kalign[k]] = diri_prob[k];
                                }
                        }else{
                                /* No Dirichlet: weighted average by depth*quality */
                                double wa = (double)da * a->cols[pos_a].quality;
                                double wb = (double)db * b->cols[pos_b].quality;
                                double wt = wa + wb;

                                if(wt > 0.0){
                                        for(k = 0; k < alpha_size; k++){
                                                col->prob[k] = (wa * a->cols[pos_a].prob[k] +
                                                                 wb * b->cols[pos_b].prob[k]) / wt;
                                        }
                                }else{
                                        for(k = 0; k < alpha_size; k++){
                                                col->prob[k] = 0.5 * (a->cols[pos_a].prob[k] +
                                                                       b->cols[pos_b].prob[k]);
                                        }
                                }
                                for(k = alpha_size; k < PROBMSA_ALPHA_MAX; k++){
                                        col->prob[k] = 0.0;
                                }
                        }

                        /* Quality: posterior * average quality */
                        double avg_qual = 0.5 * (a->cols[pos_a].quality + b->cols[pos_b].quality);
                        double post = 1.0;
                        /* posteriors are 1-indexed: [pos_a+1][pos_b+1] */
                        if(posteriors != NULL){
                                post = posteriors[pos_a + 1][pos_b + 1];
                        }
                        col->quality = post * avg_qual;
                        col->depth = a->cols[pos_a].depth + b->cols[pos_b].depth;
                        col->gap_frac = (a->cols[pos_a].gap_frac * a->n_sequences +
                                          b->cols[pos_b].gap_frac * b->n_sequences) /
                                         (double)(a->n_sequences + b->n_sequences);
                        col->concentration = a->cols[pos_a].concentration
                                           + b->cols[pos_b].concentration;

                        pos_a++;
                        pos_b++;
                }else if(op == 1){
                        /* Gap in A: copy from b, A has a gap */
                        for(k = 0; k < alpha_size; k++){
                                col->prob[k] = b->cols[pos_b].prob[k];
                        }
                        for(k = alpha_size; k < PROBMSA_ALPHA_MAX; k++){
                                col->prob[k] = 0.0;
                        }

                        /* Quality: use actual gap posterior, capped at 0.5 */
                        {
                                double pg = (gap_post_b != NULL) ? gap_post_b[pos_b] : 0.3;
                                if(pg > 0.5) pg = 0.5;
                                col->quality = b->cols[pos_b].quality * pg;
                        }
                        col->depth = b->cols[pos_b].depth + a->n_sequences;
                        col->gap_frac = ((double)a->n_sequences +
                                          b->cols[pos_b].gap_frac * b->n_sequences) /
                                         (double)(a->n_sequences + b->n_sequences);
                        col->concentration = b->cols[pos_b].concentration;

                        pos_b++;
                }else if(op == 2){
                        /* Gap in B: copy from a, B has a gap */
                        for(k = 0; k < alpha_size; k++){
                                col->prob[k] = a->cols[pos_a].prob[k];
                        }
                        for(k = alpha_size; k < PROBMSA_ALPHA_MAX; k++){
                                col->prob[k] = 0.0;
                        }

                        {
                                double pg = (gap_post_a != NULL) ? gap_post_a[pos_a] : 0.3;
                                if(pg > 0.5) pg = 0.5;
                                col->quality = a->cols[pos_a].quality * pg;
                        }
                        col->depth = a->cols[pos_a].depth + b->n_sequences;
                        col->gap_frac = (a->cols[pos_a].gap_frac * a->n_sequences +
                                          (double)b->n_sequences) /
                                         (double)(a->n_sequences + b->n_sequences);
                        col->concentration = a->cols[pos_a].concentration;

                        pos_a++;
                }
                c++;
        }

        p->length = aln_len;
        p->n_sequences = a->n_sequences + b->n_sequences;
        p->alpha_size = alpha_size;

        *out = p;
        return OK;
ERROR:
        probmsa_profile_free(p);
        return FAIL;
}

void probmsa_profile_free(struct probmsa_profile* p)
{
        if(p){
                if(p->cols){
                        MFREE(p->cols);
                }
                MFREE(p);
        }
}
