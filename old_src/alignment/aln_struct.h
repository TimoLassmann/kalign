#ifndef ALN_STRUCT_H
#define ALN_STRUCT_H

#include "stdint.h"

#define ALN_MODE_SCORE_ONLY 2
#define ALN_MODE_FULL 1


struct states{
        float a;
        float ga;
        float gb;
        // float x;
};

struct aln_mem{
        const float* prof1;
        const float* prof2;
        const uint8_t* seq1;
        const uint8_t* seq2;
        struct aln_param* ap;
        struct states* f;
        struct states* b;
        int* path;
        int* tmp_path;
        int alloc_path_len;
        float score;

        int starta;
        int starta_2;
        int startb;
        int enda;
        int enda_2;
        int endb;
        int size;
        int len_a;
        int len_b;

        int sip;
        int mode;
};

#endif
