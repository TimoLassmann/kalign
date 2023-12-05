#include "tldevel.h"
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

        qsort(m->sequences, m->numseq, sizeof(struct msa_seq*),sort_by_rank);

        return OK;
ERROR:
        return FAIL;
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
