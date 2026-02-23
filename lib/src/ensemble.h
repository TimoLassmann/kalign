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
                           uint64_t seed, int min_support,
                           const char* save_poar_path,
                           int refine, float dist_scale, float vsm_amax,
                           int realign, float use_seq_weights,
                           int consistency_anchors, float consistency_weight);

EXTERN int kalign_consensus_from_poar(struct msa* msa,
                                      const char* poar_path,
                                      int min_support);

#undef ENSEMBLE_IMPORT
#undef EXTERN

#endif
