#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <stdint.h>
#include <kalign/kalign_config.h>

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

EXTERN int kalign_ensemble_custom(struct msa* msa, int n_threads, int type,
                                  int n_runs,
                                  const float* run_gpo,
                                  const float* run_gpe,
                                  const float* run_tgpe,
                                  const int* run_types,
                                  const float* run_noise,
                                  uint64_t seed, int min_support,
                                  int refine, float vsm_amax,
                                  int realign, float use_seq_weights,
                                  int consistency_anchors, float consistency_weight);

EXTERN int kalign_ensemble_from_configs(struct msa* msa,
                                        const struct kalign_run_config* runs,
                                        int n_runs,
                                        const struct kalign_ensemble_config* ens,
                                        int n_threads);

EXTERN int kalign_generate_ensemble_runs(const struct kalign_run_config* base,
                                         int n_runs, uint64_t seed,
                                         struct kalign_run_config* out);

EXTERN int kalign_consensus_from_poar(struct msa* msa,
                                      const char* poar_path,
                                      int min_support);

#undef ENSEMBLE_IMPORT
#undef EXTERN

#endif
