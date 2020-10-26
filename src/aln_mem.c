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

        m->score = 0.0F;
        m->ap = NULL;

        m->starta = 0;
        m->startb = 0;
        m->enda = 0;
        m->endb = 0;
        m->size = x;
        m->len_a = 0;
        m->len_b = 0;
        m->f = NULL;
        m->b = NULL;
        m->path = NULL;
        m->tmp_path = NULL;
        m->alloc_path_len = x;
        MMALLOC(m->f,sizeof(struct states)* m->size);
        MMALLOC(m->b,sizeof(struct states)* m->size);
        MMALLOC(m->path, sizeof(int) * m->alloc_path_len);
        MMALLOC(m->tmp_path, sizeof(int) * m->alloc_path_len);
        *mem = m;
        return OK;
ERROR:
        free_aln_mem(m);
        return FAIL;
}

int resize_aln_mem(struct aln_mem* m)
{

        int g;

        /* For dynamic programming I only need a slice of the dyn. prog. matrix.
           To be precise, a slide of the shorter sequence */
        g = MACRO_MAX(m->len_a, m->len_b) + 2;

        if(g > m->size){
                while(m->size < g){
                        m->size = m->size + m->size / 2;
                }
                MREALLOC(m->f,sizeof(struct states)* m->size);
                MREALLOC(m->b,sizeof(struct states)* m->size);

        }

        /* For the alignment path I need at most:  */
        g = m->len_a+ m->len_b + 2;
        /* memory */
        if(g > m->alloc_path_len){
                while(m->alloc_path_len < g){
                        m->alloc_path_len = m->alloc_path_len + m->alloc_path_len / 2;
                }
                MREALLOC(m->path, sizeof(int) * m->alloc_path_len);
                MREALLOC(m->tmp_path, sizeof(int) * m->alloc_path_len);
        }
        return OK;
ERROR:
        free_aln_mem(m);
        return FAIL;
}

void free_aln_mem(struct aln_mem* m)
{
        if(m){
                MFREE(m->tmp_path);
                MFREE(m->path);
                MFREE(m->f);
                MFREE(m->b);
                MFREE(m);
        }
}
