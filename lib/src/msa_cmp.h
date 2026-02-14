#ifndef MSA_CMP_H
#define MSA_CMP_H

#include <stdint.h>

#ifdef MSA_CMP_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;

struct poar_score {
        double recall;       /* |common| / |ref_scored_pairs| — bali_score SP equivalent */
        double precision;    /* |common| / |test_pairs| — fraction of test POARs correct */
        double f1;           /* 2*recall*precision / (recall+precision) */
        double tc;           /* fraction of scored columns that are completely correct */
        int64_t ref_pairs;   /* total scored reference POARs */
        int64_t test_pairs;  /* total test POARs */
        int64_t common;      /* POARs in both ref and test */
};

EXTERN int kalign_msa_compare(struct msa *r, struct msa *t, float *score);

EXTERN int kalign_msa_compare_detailed(struct msa *r, struct msa *t,
                                       float max_gap_frac,
                                       struct poar_score *out);

EXTERN int kalign_msa_compare_with_mask(struct msa *r, struct msa *t,
                                        int *scored_cols, int n_cols,
                                        struct poar_score *out);

#undef MSA_CMP_IMPORT
#undef EXTERN


#endif
