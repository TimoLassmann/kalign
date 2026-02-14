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
                             uint64_t tree_seed, float tree_noise);
EXTERN int kalign_run_realign(struct msa *msa, int n_threads, int type,
                              float gpo, float gpe, float tgpe,
                              int refine, int adaptive_budget,
                              float dist_scale, float vsm_amax,
                              int realign_iterations);

EXTERN int kalign_post_realign(struct msa *msa, int n_threads, int type,
                               float gpo, float gpe, float tgpe,
                               int refine, int adaptive_budget,
                               float dist_scale, float vsm_amax,
                               int realign_iterations);

#undef ALN_WRAP_IMPORT
#undef EXTERN


#endif
