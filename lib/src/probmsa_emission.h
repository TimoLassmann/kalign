#ifndef PROBMSA_EMISSION_H
#define PROBMSA_EMISSION_H

#include "probmsa_profile.h"
#include "probmsa_submat.h"

#ifdef PROBMSA_EMISSION_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

/* Match emission: log P(col_A, col_B | Match)
   = log(sum_s sum_t p_A[s] * p_B[t] * joint[s][t])
   col_a and col_b are 0-indexed into profiles. */
EXTERN double probmsa_emit_match(const struct probmsa_profile* a, int col_a,
                                  const struct probmsa_profile* b, int col_b,
                                  const double joint[][PROBMSA_ALPHA_MAX],
                                  int alpha_size);

/* Gap emission: log P(col | Gap)
   = log(sum_s p[s] * background[s])
   col is 0-indexed into profile. */
EXTERN double probmsa_emit_gap(const struct probmsa_profile* p, int col,
                                const double* background, int alpha_size);

#undef PROBMSA_EMISSION_IMPORT
#undef EXTERN

#endif
