#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <stdint.h>

#ifdef ENSEMBLE_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;

EXTERN int kalign_ensemble(struct msa* msa, int n_threads, int type,
                           int n_runs, float gpo, float gpe, float tgpe,
                           uint64_t seed);

#undef ENSEMBLE_IMPORT
#undef EXTERN

#endif
