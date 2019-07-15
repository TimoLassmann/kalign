/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/


#include "pick_anchor.h"

struct sort_struct{
        int len;
        int id;
};

int sort_by_len(const void *a, const void *b);

int* select_seqs(struct msa* msa, int num_anchor);

int* pick_anchor(struct msa* msa, int* n)
{


        int* anchors = NULL;
        int num_anchor = 0;




        ASSERT(msa != NULL, "No alignment.");

        num_anchor = MACRO_MAX(MACRO_MIN(32, msa->numseq), (int) pow(log2((double) msa->numseq), 2.0));
        RUNP(anchors = select_seqs(msa, num_anchor));
        *n = num_anchor;
        return anchors;
ERROR:
        return NULL;
}

int* select_seqs(struct msa* msa, int num_anchor)
{
        struct sort_struct** seq_sort = NULL;
        int* anchors = NULL;
        int i,stride;

        MMALLOC(seq_sort, sizeof(struct sort_struct*) * msa->numseq);
        for(i = 0; i < msa->numseq;i++){
                seq_sort[i] = NULL;
                MMALLOC(seq_sort[i], sizeof(struct sort_struct));
                seq_sort[i]->id = i;
                seq_sort[i]->len = msa->sequences[i]->len;//  aln->sl[i];
        }

        qsort(seq_sort, msa->numseq, sizeof(struct sort_struct*),sort_by_len);
        //for(i = 0; i < aln->numseq;i++){
        //fprintf(stdout,"%d\t%d\n", seq_sort[i]->id,seq_sort[i]->len);
        //}


        //fprintf(stdout,"%d\t seeds\n", num_anchor);

        MMALLOC(anchors, sizeof(int) * num_anchor);
        stride = msa->numseq / num_anchor;
//        fprintf(stdout,"%d\tstride\n", stride);
        //c = 0;
        for(i = 0; i < num_anchor;i++){
                anchors[i] = seq_sort[i*stride]->id;
        }
        ASSERT(i == num_anchor,"Cound not select all anchors\tnum_anchor:%d\t numseq:%d",num_anchor,msa->numseq);

        for(i = 0; i < msa->numseq;i++){
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
