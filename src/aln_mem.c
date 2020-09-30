#include "tldevel.h"

#include "aln_struct.h"

#define ALN_MEM_IMPORT
#include "aln_mem.h"


int alloc_aln_mem(struct aln_mem** mem, int x)
{

        struct aln_mem* m = NULL;

        ASSERT(x>=1, "Given size %d is too small",x);
        // a=((typeof(a))(((int)(((void *)malloc(c+15))+15))&-16)).
        MMALLOC(m,sizeof(struct aln_mem));

        m->seq1 = NULL;
        m->seq2 = NULL;
        m->prof1 = NULL;
        m->prof2 = NULL;
        m->sip = 0;
        m->mode = ALN_MODE_FULL;

        m->ap = NULL;

        m->starta = 0;
        m->startb = 0;
        m->enda = 0;
        m->endb = 0;
        m->size = x+1;
        m->len_a = 0;
        m->len_b = 0;
        m->f = NULL;
        m->b = NULL;
        MMALLOC(m->f,sizeof(struct states)* (x+1));
        MMALLOC(m->b,sizeof(struct states)* (x+1));
        *mem = m;
        return OK;
ERROR:
        free_aln_mem(m);
        return FAIL;
}

int resize_aln_mem(struct aln_mem* m,int x)
{
        if((x+1) > m->size){
                m->size = x+1;

                MREALLOC(m->f,sizeof(struct states)* (x+1));
                MREALLOC(m->b,sizeof(struct states)* (x+1));
        }
        return OK;
ERROR:
        free_aln_mem(m);
        return FAIL;
}

void free_aln_mem(struct aln_mem* m)
{
        if(m){
                MFREE(m->f);
                MFREE(m->b);
                MFREE(m);
        }
}
