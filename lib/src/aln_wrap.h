#ifndef ALN_WRAP_H
#define ALN_WRAP_H

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

EXTERN int kalign_run(struct msa *msa, int n_threads, int type, float gpo, float gpe, float tgpe);

#undef ALN_WRAP_IMPORT
#undef EXTERN


#endif
