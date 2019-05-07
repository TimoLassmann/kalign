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

extern struct alphabet* create_alphabet(int type);
extern int switch_alphabet(struct alphabet* a, int type);



#endif
