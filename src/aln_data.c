
#include "aln_data.h"


struct kalign_sequence* kalign_seq_alloc(void);
void kalign_seq_free(struct kalign_sequence* ks);

struct kalign_alignmment* kalign_aln_alloc(void)
{
        struct kalign_alignmment* aln = NULL;
        int i;
        MMALLOC(aln, sizeof(aln));
        aln->alloc_numseq = 256;
        aln->numseq = 0;
        aln->s_arr = NULL;

        MMALLOC(aln->s_arr, sizeof(struct kalign_sequence*) * aln->alloc_numseq);
        for(i = 0; i < aln->alloc_numseq;i++){
                aln->s_arr[i] = NULL;
                RUNP(aln->s_arr[i] = kalign_seq_alloc());
                aln->s_arr[i]->id = i;
        }
        return aln;
ERROR:
        return NULL;
}


int kalign_alignment_resize(struct kalign_alignmment* aln)
{
        int i;
        int c;
        ASSERT(aln!= NULL, "No aln.");
        c = aln->alloc_numseq;
        aln->alloc_numseq = aln->alloc_numseq << 1;

        MREALLOC(aln->s_arr, sizeof(struct kalign_sequence*) * aln->alloc_numseq);
        for(i = c; i < aln->alloc_numseq;i++){
                aln->s_arr[i] = NULL;
                RUNP(aln->s_arr[i] = kalign_seq_alloc());
                aln->s_arr[i]->id = i;
        }

        return OK;
ERROR:
        return FAIL;
}


void kalign_alignmment_free(struct kalign_alignmment* aln)
{
        int i;
        if(aln){
                if(aln->s_arr){
                        for(i = 0; i < aln->alloc_numseq;i++){
                                kalign_seq_free(aln->s_arr[i]);
                        }
                        MFREE(aln->s_arr);
                }
                MFREE(aln);
        }
}

struct kalign_sequence* kalign_seq_alloc(void)
{
        struct kalign_sequence* ks = NULL;

        MMALLOC(ks, sizeof(struct kalign_sequence));
        ks->alloc_seq_len = 1024;
        ks->id = 0;
        ks->len = 0;
        ks->s =NULL;
        ks->seq = NULL;
        ks->name = NULL;
        MMALLOC(ks->s, sizeof(int_fast8_t) * ks->alloc_seq_len);
        MMALLOC(ks->seq, sizeof(char) * ks->alloc_seq_len);
        MMALLOC(ks->name, sizeof(char) * BUFFER_LEN);
        return ks;
ERROR:
        kalign_seq_free(ks);
        return NULL;
}

int kalign_seq_resize(struct kalign_sequence* ks)
{
        ASSERT(ks != NULL,"No seq.");
        ks->alloc_seq_len = ks->alloc_seq_len << 1;
        MREALLOC(ks->s, sizeof(int_fast8_t) * ks->alloc_seq_len);
        MREALLOC(ks->seq, sizeof(char) * ks->alloc_seq_len);
        MREALLOC(ks->name, sizeof(char) * BUFFER_LEN);
        return OK;
ERROR:
        kalign_seq_free(ks);
        return FAIL;
}

void kalign_seq_free(struct kalign_sequence* ks)
{
        if(ks){
                if(ks->s){
                        MFREE(ks->s);
                }
                if(ks->seq){
                        MFREE(ks->seq);
                }
                if(ks->name){
                        MFREE(ks->name);
                }
        }
}
