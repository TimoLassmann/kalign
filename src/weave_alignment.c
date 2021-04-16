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

#include "weave_alignment.h"


//struct alignment* make_seq(struct alignment* aln,int a,int b,int* path);

//int update_gaps(int old_len,int*gis,int new_len,int *newgaps);
int update_gaps(int old_len,int*gis,int *newgaps);

int clean_aln(struct msa* msa)
{
        int i,j;
        int* p = NULL;
        for (i = 0; i < msa->numseq;i++){
                p = msa->sequences[i]->gaps;//  aln->gaps[i];
                for (j = 0; j <= msa->sequences[i]->len;j++){
                        p[j] = 0;
                }
        }
        for(i =0;i < msa->numseq;i++){
                msa->nsip[i] = 1;
                msa->sip[i][0] = i;
        }

        for (i = msa->numseq;i < msa->num_profiles ;i++){
                if(msa->sip[i]){
                        MFREE(msa->sip[i]);

                        msa->sip[i] = NULL;
                }
                msa->nsip[i] =0;
        }




        return OK;

}

int make_seq(struct msa* msa,int a,int b,int* path)
{
        int* gap_a = NULL;
        int* gap_b = NULL;

        int c;
        int i;
        int posa = 0;
        int posb = 0;

        MMALLOC(gap_a,(path[0]+1)*sizeof(int));
        MMALLOC(gap_b,(path[0]+1)*sizeof(int));

        for (i = path[0]+1;i--;){
                gap_a[i] = 0;
                gap_b[i] = 0;
        }
        c = 1;
        while(path[c] != 3){

                if (!path[c]){
                        posa++;
                        posb++;
                }else
                if (path[c] & 1){
                        gap_a[posa] += 1;
                        posb++;
                }else
                if (path[c] & 2){
                        gap_b[posb] += 1;
                        posa++;
                }
                c++;
        }

        for (i = msa->nsip[a];i--;){

                RUN(update_gaps(msa->sequences[msa->sip[a][i]]->len, msa->sequences[msa->sip[a][i]]->gaps,gap_a));
                //RUN(update_gaps(aln->sl[aln->sip[a][i]],aln->gaps[aln->sip[a][i]],path[0],gap_a));
        }
        for (i = msa->nsip[b];i--;){
                RUN(update_gaps(msa->sequences[msa->sip[b][i]]->len,msa->sequences[msa->sip[b][i]]->gaps,gap_b));
                //RUN(update_gaps(aln->sl[aln->sip[b][i]],aln->gaps[aln->sip[b][i]],path[0],gap_b));
        }
        MFREE(gap_a);
        MFREE(gap_b);
        return OK;
ERROR:
        if(gap_a){
                MFREE(gap_a);
        }
        if(gap_b){
                MFREE(gap_b);
        }
        return FAIL;
}

int update_gaps(int old_len,int*gis,int *newgaps)
{
        int i,j;
        int add = 0;
        int rel_pos = 0;
        for (i = 0; i <= old_len;i++){
                add = 0;
                for (j = rel_pos;j <= rel_pos + gis[i];j++){
                        if (newgaps[j] != 0){
                                add += newgaps[j];
                        }
                }
                rel_pos += gis[i]+1;
                gis[i] += add;
        }
        return OK;
}
