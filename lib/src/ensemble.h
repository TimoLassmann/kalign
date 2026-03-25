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
