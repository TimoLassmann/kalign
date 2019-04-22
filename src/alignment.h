#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "global.h"

#include <float.h>
struct states{
        float a;
        float ga;
        float gb;
        float x;
};

struct hirsch_mem{
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

struct dp_matrix{
        struct states* s;
        void* tb_mem;
        char** tb;
        int x;
        int y;
};


/* Memory allocation for forward and backward slices  */
struct hirsch_mem* hirsch_mem_alloc(struct hirsch_mem* hm,int x);
int hirsch_mem_realloc(struct hirsch_mem* hm,int x);
void hirsch_mem_free(struct hirsch_mem* hm);

extern int align(struct alignment* aln, struct aln_param* ap);

#endif
