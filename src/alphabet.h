#ifndef ALPHABET_H
#define ALPHABET_H

#include "global.h"

#define defPROTEIN 1
#define redPROTEIN 2
#define defDNA 3


struct alphabet{
        int_fast8_t to_internal[128];
        int_fast8_t to_external[32];
        int type;
        int L;
};

struct alphabet* create_alphabet(int type);




#endif
