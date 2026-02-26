#ifndef ALPHABET_H
#define ALPHABET_H

#include <inttypes.h>

#ifdef ALPHABET_IMPORT
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

#define ALPHA_defPROTEIN 21
#define ALPHA_ambigiousPROTEIN 23
#define ALPHA_redPROTEIN 13
#define ALPHA_redPROTEIN2 8
#define ALPHA_defDNA 5

#define ALPHA_UNKNOWN 255
#define ALPHA_UNDEFINED -1

struct alphabet{
        int8_t to_internal[128];
        int8_t to_external[32];
        int type;
        int L;
};

EXTERN struct alphabet* create_alphabet(int type);
EXTERN int switch_alphabet(struct alphabet* a, int type);


#undef ALPHABET_IMPORT
#undef EXTERN

#endif
