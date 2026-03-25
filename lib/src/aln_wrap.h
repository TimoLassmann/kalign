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

/* Internal single-run entry point.  Caller must set msa->pool before
   calling (if USE_THREADPOOL is enabled).  Does NOT create or destroy
   the threadpool.  Used by kalign_align_full and ensemble.c. */
EXTERN int kalign_single_run(struct msa *msa,
                             const struct kalign_run_config *cfg,
                             int n_threads);

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
