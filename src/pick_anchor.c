
#include "pick_anchor.h"

struct sort_struct{
        int len;
        int id;
};

int sort_by_len(const void *a, const void *b);

int* select_seqs(struct alignment* aln, int num_anchor);

int* pick_anchor(struct alignment* aln, int* n)
{


        int* anchors = NULL;
        int num_anchor = 0;


        ASSERT(aln != NULL, "No alignment.");

        num_anchor = MACRO_MAX(MACRO_MIN(8, aln->numseq), (int) pow(log2((double) aln->numseq), 2.0));
        RUNP(anchors = select_seqs(aln, num_anchor));
        *n = num_anchor;
        return anchors;
ERROR:
        return NULL;
}

int* select_seqs(struct alignment* aln, int num_anchor)
{
        struct sort_struct** seq_sort = NULL;
        int* anchors = NULL;
        int i,stride;

        MMALLOC(seq_sort, sizeof(struct sort_struct*) * aln->numseq);
        for(i = 0; i < aln->numseq;i++){
                seq_sort[i] = NULL;
                MMALLOC(seq_sort[i], sizeof(struct sort_struct));
                seq_sort[i]->id = i;
                seq_sort[i]->len = aln->sl[i];
        }

        qsort(seq_sort, aln->numseq, sizeof(struct sort_struct*),sort_by_len);
        //for(i = 0; i < aln->numseq;i++){
        //fprintf(stdout,"%d\t%d\n", seq_sort[i]->id,seq_sort[i]->len);
        //}


        //fprintf(stdout,"%d\t seeds\n", num_anchor);

        MMALLOC(anchors, sizeof(int) * num_anchor);
        stride = aln->numseq / num_anchor;
//        fprintf(stdout,"%d\tstride\n", stride);
        //c = 0;
        for(i = 0; i < num_anchor;i++){
                anchors[i] = seq_sort[i*stride]->id;
        }
        ASSERT(i == num_anchor,"Cound not select all anchors\tnum_anchor:%d\t numseq:%d",num_anchor,aln->numseq);

        for(i = 0; i < aln->numseq;i++){
                MFREE(seq_sort[i]);
        }
        MFREE(seq_sort);
        return anchors;
ERROR:
        return NULL;
}

int sort_by_len(const void *a, const void *b)
{
        struct sort_struct* const *one = a;
        struct sort_struct* const *two = b;

        if((*one)->len > (*two)->len){
                return -1;
        }else{
                return 1;
        }
}
