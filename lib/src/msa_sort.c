#include "tldevel.h"
#include "tlrng.h"

#include "msa_struct.h"

#include <string.h>

#define MSA_SORT_IMPORT
#include "msa_sort.h"

static int sort_by_len_name(const void *a, const void *b);
static int sort_by_rank(const void *a, const void *b);

int msa_sort_len_name(struct msa *m)
{
        ASSERT(m != NULL, "No alignment");

        qsort(m->sequences, m->numseq, sizeof(struct msa_seq*),sort_by_len_name);

        return OK;
ERROR:
        return FAIL;
}

int msa_sort_rank(struct msa *m)
{
        ASSERT(m != NULL, "No alignment");
        /* for(int i = 0; i < m->numseq;i++){ */
        /*         fprintf(stdout,"%d %s\n", */
        /*                 m->sequences[i]->rank, */
        /*                 m->sequences[i]->name */
        /*                 ); */
        /* } */
        qsort(m->sequences, m->numseq, sizeof(struct msa_seq*),sort_by_rank);
        /* for(int i = 0; i < m->numseq;i++){ */
        /*         fprintf(stdout,"%d %s - sorted\n", */
        /*                 m->sequences[i]->rank, */
        /*                 m->sequences[i]->name */
        /*                 ); */
        /* } */
        return OK;
ERROR:
        return FAIL;
}

int msa_shuffle_seq(struct msa *m, struct rng_state* rng)
{
        int r;
        int i,j;
        struct msa_seq* tmp;
        int n = m->numseq;
        for (i = 0; i < n - 1; i++) {
                r = tl_random_int(rng,n);
                j = i +  r % (n-i);
                tmp =  m->sequences[j];
                m->sequences[j] =  m->sequences[i];
                m->sequences[i] = tmp;
        }
        return OK;
}

int sort_by_len_name(const void *a, const void *b)
{
        struct msa_seq* const *one = a;
        struct msa_seq* const *two = b;

        if((*one)->len > (*two)->len){
                return -1;
        }else if((*one)->len == (*two)->len){
                int c = strncmp((*one)->name, (*two)->name, MSA_NAME_LEN);
                if(c < 0){
                        return -1;
                }else{
                        return 1;
                }

        }else{
                return 1;
        }
}

int sort_by_rank(const void *a, const void *b)
{
        struct msa_seq* const *one = a;
        struct msa_seq* const *two = b;

        if((*one)->rank > (*two)->rank){
                return 1;
        }else{
                return -1;
        }
}
