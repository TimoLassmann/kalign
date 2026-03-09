#ifndef ALN_WRAP_H
#define ALN_WRAP_H

#include <stdint.h>
#include <kalign/kalign_config.h>

#ifdef ALN_WRAP_IMPORT
   #define EXTERN
#else
   #ifndef EXTERN
      #ifdef __cplusplus
         #define EXTERN extern "C"
      #else
         #define EXTERN extern
      #endif
   #endif
#endif

struct msa;

EXTERN int kalign_run(struct msa *msa, int n_threads, int type, float gpo, float gpe, float tgpe, int refine, int adaptive_budget);
EXTERN int kalign_run_seeded(struct msa *msa, int n_threads, int type,
                             float gpo, float gpe, float tgpe,
                             int refine, int adaptive_budget,
                             uint64_t tree_seed, float tree_noise,
                             float dist_scale, float vsm_amax,
                             float use_seq_weights,
                             int consistency_anchors, float consistency_weight);
EXTERN int kalign_run_realign(struct msa *msa, int n_threads, int type,
                              float gpo, float gpe, float tgpe,
                              int refine, int adaptive_budget,
                              float dist_scale, float vsm_amax,
                              int realign_iterations,
                              float use_seq_weights,
                              int consistency_anchors, float consistency_weight);

EXTERN int kalign_run_dist_scale(struct msa *msa, int n_threads, int type,
                                  float gpo, float gpe, float tgpe,
                                  int refine, int adaptive_budget,
                                  float dist_scale, float vsm_amax,
                                  float use_seq_weights);

EXTERN int kalign_post_realign(struct msa *msa, int n_threads, int type,
                               float gpo, float gpe, float tgpe,
                               int refine, int adaptive_budget,
                               float dist_scale, float vsm_amax,
                               int realign_iterations,
                               float use_seq_weights);

EXTERN struct kalign_run_config kalign_run_config_defaults(void);
EXTERN struct kalign_ensemble_config kalign_ensemble_config_defaults(void);

EXTERN int kalign_align_full(struct msa* msa,
                             const struct kalign_run_config* runs,
                             int n_runs,
                             const struct kalign_ensemble_config* ens,
                             int n_threads);

#undef ALN_WRAP_IMPORT
#undef EXTERN


#endif
