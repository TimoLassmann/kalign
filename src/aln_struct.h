#ifndef ALN_STRUCT_H
#define ALN_STRUCT_H


struct states{
        float a;
        float ga;
        float gb;
        // float x;
};

struct aln_mem{
        struct states* f;
        struct states* b;
        int starta;
        int startb;
        int enda;
        int endb;
        int size;
        int len_a;
        int len_b;
};

#endif
