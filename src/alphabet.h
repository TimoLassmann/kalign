#ifndef ALPHABET_H
#define ALPHABET_H

#include "global.h"

#define defPROTEIN 21
#define redPROTEIN 13
#define defDNA 5


struct alphabet{
        int8_t to_internal[128];
        int8_t to_external[32];
        int type;
        int L;
};

extern struct alphabet* create_alphabet(int type);
extern int switch_alphabet(struct alphabet* a, int type);



#endif
