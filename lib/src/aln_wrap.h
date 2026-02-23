#ifndef ALN_WRAP_H
#define ALN_WRAP_H

#include <stdint.h>

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

EXTERN int kalign_post_realign(struct msa *msa, int n_threads, int type,
                               float gpo, float gpe, float tgpe,
                               int refine, int adaptive_budget,
                               float dist_scale, float vsm_amax,
                               int realign_iterations,
                               float use_seq_weights);

EXTERN int kalign_run_probmsa(struct msa *msa, int n_threads, int type,
                              double hmm_delta, double hmm_epsilon,
                              double hmm_tau,
                              double mea_threshold, double mea_gpo, double mea_gpe,
                              double prior_scale,
                              int use_5state,
                              double delta_s, double epsilon_s,
                              double delta_l, double epsilon_l);

EXTERN int kalign_post_realign_probmsa(struct msa *msa, int n_threads, int type,
                                       double hmm_delta, double hmm_epsilon,
                                       double hmm_tau,
                                       double mea_threshold, double mea_gpo, double mea_gpe,
                                       double prior_scale);

#undef ALN_WRAP_IMPORT
#undef EXTERN


#endif
