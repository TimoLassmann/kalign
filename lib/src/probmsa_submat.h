#ifndef PROBMSA_SUBMAT_H
#define PROBMSA_SUBMAT_H

#ifdef PROBMSA_SUBMAT_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

#define PROBMSA_ALPHA_MAX 23

/* Convert kalign's integer log-odds substitution matrix (float** subm)
   to a joint probability matrix suitable for the pair-HMM.
   alpha_size: 5 for DNA, 23 for protein.
   joint[s][t] = P(s,t | Match), sums to ~1.0.
   background[s] = marginal frequency of residue s. */
EXTERN int probmsa_submat_from_kalign(float** subm, int alpha_size,
                                       double joint[][PROBMSA_ALPHA_MAX],
                                       double* background);

#undef PROBMSA_SUBMAT_IMPORT
#undef EXTERN

#endif
