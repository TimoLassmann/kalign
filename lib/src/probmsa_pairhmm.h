#ifndef PROBMSA_PAIRHMM_H
#define PROBMSA_PAIRHMM_H

#include "probmsa_profile.h"
#include "probmsa_submat.h"

#ifdef PROBMSA_PAIRHMM_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct probmsa_dp {
        double** M;    double** X;    double** Y;     /* forward matrices (3-state) or M/Xs/Ys (5-state) */
        double** bM;   double** bX;   double** bY;    /* backward matrices (3-state) or bM/bXs/bYs (5-state) */
        double** Xl;   double** Yl;                    /* 5-state: long gap forward (NULL for 3-state) */
        double** bXl;  double** bYl;                   /* 5-state: long gap backward (NULL for 3-state) */
        double MM, MX, MY, XM, XX, YM, YY;            /* 3-state log transition probs (global fallback) */
        double SX, SY;                                  /* log Start→X, Start→Y */
        /* 5-state transition probabilities (log space) */
        double t5_MM;                           /* M→M = log(1 - 2*δs - 2*δl) */
        double t5_MXs, t5_MYs;                 /* M→Xs, M→Ys = log(δs) */
        double t5_MXl, t5_MYl;                 /* M→Xl, M→Yl = log(δl) */
        double t5_XsM, t5_XsXs;               /* Xs→M = log(1-εs), Xs→Xs = log(εs) */
        double t5_YsM, t5_YsYs;               /* Ys→M, Ys→Ys */
        double t5_XlM, t5_XlXl;               /* Xl→M = log(1-εl), Xl→Xl = log(εl) */
        double t5_YlM, t5_YlYl;               /* Yl→M, Yl→Yl */
        double t5_SXs, t5_SYs, t5_SXl, t5_SYl; /* Start transitions (5-state) */
        int semiglobal;                                 /* 1: free trailing gaps */
        double joint[PROBMSA_ALPHA_MAX][PROBMSA_ALPHA_MAX];
        double background[PROBMSA_ALPHA_MAX];
        double mea_threshold;   /* posterior threshold for MEA (default 0) */
        double mea_gpo;         /* MEA gap-open cost (default 0) */
        double mea_gpe;         /* MEA gap-extend cost (default 0) */
        int alpha_size;
        int alloc_a, alloc_b;                           /* allocated dimensions */
        int use_dotproduct;     /* 1: dot-product emission (Dirichlet mode) */
        double delta_base;      /* base gap-open probability for position-specific modulation */
        double epsilon;         /* gap-extend probability (stored for reuse) */
        double conc_scale;      /* concentration scale: delta_eff = delta_base * cs / (conc + cs) */
        int use_5state;         /* 1: use 5-state model */
};

/* Allocate DP workspace. */
EXTERN int probmsa_dp_alloc(struct probmsa_dp** dp, int max_a, int max_b,
                             int alpha_size);

/* Set transition probabilities from gap-open, gap-extend, and terminal gap-open.
   tau <= 0: global alignment (terminal gaps same as internal).
   tau > 0: semi-global (leading gaps use tau, trailing gaps free). */
EXTERN int probmsa_dp_set_transitions(struct probmsa_dp* dp,
                                       double delta, double epsilon,
                                       double tau);

/* Free DP workspace. */
EXTERN void probmsa_dp_free(struct probmsa_dp* dp);

/* Forward algorithm. Returns total log-score. */
EXTERN int probmsa_forward(struct probmsa_dp* dp,
                            const struct probmsa_profile* a,
                            const struct probmsa_profile* b,
                            double* score);

/* Backward algorithm. Returns total log-score. */
EXTERN int probmsa_backward(struct probmsa_dp* dp,
                              const struct probmsa_profile* a,
                              const struct probmsa_profile* b,
                              double* score);

/* Compute posterior probabilities: overwrite forward M with log posteriors.
   Requires profiles to recompute emissions (avoids double-counting). */
EXTERN int probmsa_posterior(struct probmsa_dp* dp,
                              const struct probmsa_profile* a,
                              const struct probmsa_profile* b,
                              int len_a, int len_b, double total);

/* MEA traceback -> kalign path format.
   path[0] = alignment length, path[1..len] = 0/1/2, path[len+1] = 3.
   Match posteriors are available in dp->Y[i][j] as probabilities. */
EXTERN int probmsa_mea_traceback(struct probmsa_dp* dp, int len_a, int len_b,
                                   int** path);

/* --- 5-state pair-HMM (M, Xs, Ys, Xl, Yl) --- */

/* Allocate DP workspace for 5-state model (10 matrices). */
EXTERN int probmsa_dp_alloc_5state(struct probmsa_dp** dp, int max_a, int max_b,
                                    int alpha_size);

/* Set 5-state transitions from short/long gap-open and gap-extend probs. */
EXTERN int probmsa_dp_set_transitions_5state(struct probmsa_dp* dp,
                                              double delta_s, double epsilon_s,
                                              double delta_l, double epsilon_l);

/* 5-state forward algorithm. */
EXTERN int probmsa_forward_5state(struct probmsa_dp* dp,
                                    const struct probmsa_profile* a,
                                    const struct probmsa_profile* b,
                                    double* score);

/* 5-state backward algorithm. */
EXTERN int probmsa_backward_5state(struct probmsa_dp* dp,
                                     const struct probmsa_profile* a,
                                     const struct probmsa_profile* b,
                                     double* score);

/* 5-state posterior: overwrite forward M with log P(Match|data). */
EXTERN int probmsa_posterior_5state(struct probmsa_dp* dp,
                                      const struct probmsa_profile* a,
                                      const struct probmsa_profile* b,
                                      int len_a, int len_b, double total);

#undef PROBMSA_PAIRHMM_IMPORT
#undef EXTERN

#endif
