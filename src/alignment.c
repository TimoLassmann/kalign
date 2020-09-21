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

#include <xmmintrin.h>
#include "alignment.h"


#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)

struct states{
        float a;
        float ga;
        float gb;
        // float x;
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
struct hirsch_mem* hirsch_mem_alloc(int x);
int hirsch_mem_realloc(struct hirsch_mem* hm,int x);
void hirsch_mem_free(struct hirsch_mem* hm);

/* setting up fast data structures for alignment */
float* make_profile(struct aln_param* ap,const uint8_t* seq,const int len);
int set_gap_penalties(float* prof,int len,int nsip);


/* Main dyn programming functions */

//int hirsch_ss_dyn_score(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path,float* score);

int hirsch_ss_dyn_score(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, float* score);
int hirsch_align_two_ss_vector_score(const struct aln_param* ap,struct hirsch_mem* hm,int old_cor[],float* score);

//int hirsch_align_two_ss_vector_score(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm,int old_cor[],float* score);
/* Align 2 sequences  */
int hirsch_ss_dyn(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path);
int hirsch_align_two_ss_vector(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int foward_hirsch_ss_dyn(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm);
int backward_hirsch_ss_dyn(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm);

/* Align one sequence to a profile */
int hirsch_ps_dyn(const struct aln_param* ap, const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip);
int hirsch_align_two_ps_vector(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip);
int foward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip);
int backward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip);

/* Align 2 profiles  */
int hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
int hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
int backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);

/* auxiliary dyn prog functions  */
int* mirror_hirsch_path(int* hirsch_path,int len_a,int len_b);
int* add_gap_info_to_hirsch_path(int* hirsch_path,int len_a,int len_b);


//int update(const float* profa, const float* profb,float* newp,const int* path);
int update(const float* profa, const float* profb,float* newp, struct aln_param*ap, int* path,int sipa,int sipb);
int calc_seq_id(const int* path,int* a, int*b,float* dist);




int** hirschberg_alignment(struct msa* msa, struct aln_param* ap)
{
        struct hirsch_mem* hm = NULL;
        int i,j,g,a,b,c;
        int len_a;
        int len_b;
        float** profile = NULL;

        int** map = NULL;
        int* tree = NULL;
        int numseq;

        g = msa->num_profiles;
        numseq = msa->numseq;

        tree = ap->tree;
        MMALLOC(profile,sizeof(float*)*g);
        MMALLOC(map,sizeof(int*)*g);
        for ( i = 0;i< g;i++){
                profile[i] = NULL;
                map[i] = NULL;
        }

        RUNP(hm = hirsch_mem_alloc(2048));




        for (i = 0; i < (numseq-1);i++){
                a = tree[i*3];
                b = tree[i*3+1];
                c = tree[i*3+2];
                //fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)numseq * 100);
                //fprintf(stdout,"Aligning:%d %d->%d	done:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
                if(a < numseq){
                        len_a = msa->sequences[a]->len;//  aln->sl[a];
                }else{
                        len_a = msa->plen[a];
                }
                if(b < numseq){

                        len_b = msa->sequences[b]->len;// aln->sl[b];
                }else{
                        len_b = msa->plen[b];
                }


                g = (len_a > len_b)? len_a:len_b;
                MMALLOC(map[c],sizeof(int) * (g+2));

                RUN(hirsch_mem_realloc(hm, g));

                for (j = 0; j < (g+2);j++){
                        map[c][j] = -1;
                }

                if (a < numseq){
                        RUNP(profile[a] = make_profile(ap,msa->sequences[a]->s,len_a));
                }else{
                        RUN(set_gap_penalties(profile[a],len_a,msa->nsip[b]));
                        //smooth_gaps(profile[a],len_a,window,strength);

                        //increase_gaps(profile[a],len_a,window,strength);
                }
                if (b < numseq){
                        RUNP(profile[b] = make_profile(ap,msa->sequences[b]->s,len_b));
                }else{
                        RUN(set_gap_penalties(profile[b],len_b,msa->nsip[a]));
                        //smooth_gaps(profile[b],len_b,window,strength);
                        //increase_gaps(profile[b],len_b,window,strength);
                }

                hm->starta = 0;
                hm->startb = 0;
                hm->enda = len_a;
                hm->endb = len_b;
                hm->len_a = len_a;
                hm->len_b = len_b;

                hm->f[0].a = 0.0;
                hm->f[0].ga =  -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = 0.0;
                hm->b[0].ga =  -FLT_MAX;
                hm->b[0].gb =  -FLT_MAX;
                //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
                if(a < numseq){
                        if(b < numseq){
                                hirsch_ss_dyn(ap,msa->sequences[a]->s, msa->sequences[b]->s,hm,map[c]);
                        }else{
                                hm->enda = len_b;
                                hm->endb = len_a;
                                hm->len_a = len_b;
                                hm->len_b = len_a;
                                hirsch_ps_dyn(ap,profile[b], msa->sequences[a]->s,hm,map[c],msa->nsip[b]);
                                RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                        }
                }else{
                        if(b < numseq){
                                hirsch_ps_dyn(ap,profile[a],msa->sequences[b]->s ,hm,map[c],msa->nsip[a]);
                        }else{
                                if(len_a < len_b){
                                        hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
                                }else{
                                        hm->enda = len_b;
                                        hm->endb = len_a;
                                        hm->len_a = len_b;
                                        hm->len_b = len_a;
                                        hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
                                        RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                                }
                        }
                }

                map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

                if(i != numseq-2){
                        //MREALLOC(profile_ptr, sizeof(float)*64*(map[c][0]+2));
                        MMALLOC(profile[c],sizeof(float)*64*(map[c][0]+2));
                        //update(profile[a],profile[b],profile[c],map[c]);
                        update(profile[a],profile[b],profile[c],ap,map[c],msa->nsip[a],msa->nsip[b]);
                }

                msa->plen[c] = map[c][0];

                msa->nsip[c] = msa->nsip[a] + msa->nsip[b];
                MMALLOC(msa->sip[c],sizeof(int)*(msa->nsip[a] + msa->nsip[b]));
                g =0;
                for (j = msa->nsip[a];j--;){
                        msa->sip[c][g] = msa->sip[a][j];
                        g++;
                }
                for (j = msa->nsip[b];j--;){
                        msa->sip[c][g] = msa->sip[b][j];
                        g++;
                }


                MFREE(profile[a]);
                MFREE(profile[b]);

        }
        MFREE(profile);

        hirsch_mem_free(hm);
        return map;
ERROR:
        return NULL;
}


int hirsch_ss_dyn(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path)
{
        int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
        float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
        int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

        if(hm->starta  >= hm->enda){
                return OK;//hirsch_path;
        }
        if(hm->startb  >= hm->endb){
                return OK;///hirsch_path;
        }


        hm->enda = mid;

        //fprintf(stderr,"Forward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
        foward_hirsch_ss_dyn(ap,seq1,seq2,hm);

        hm->starta = mid;
        hm->enda = old_cor[1];
        //fprintf(stderr,"Backward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
        backward_hirsch_ss_dyn(ap,seq1,seq2,hm);


        hirsch_align_two_ss_vector(ap,seq1,seq2,hm,hirsch_path,input_states,old_cor);

        return  OK;
}

int hirsch_ss_dyn_score(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, float* score)
{
        int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
        //float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
        int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

        if(hm->starta  >= hm->enda){
                return OK;//hirsch_path;
        }
        if(hm->startb  >= hm->endb){
                return OK;///hirsch_path;
        }


        hm->enda = mid;

        //fprintf(stderr,"Forward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
        foward_hirsch_ss_dyn(ap,seq1,seq2,hm);

        hm->starta = mid;
        hm->enda = old_cor[1];
        //fprintf(stderr,"Backward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
        backward_hirsch_ss_dyn(ap,seq1,seq2,hm);

        hirsch_align_two_ss_vector_score(ap,hm,old_cor,score);


        return  OK;
}


int hirsch_align_two_ss_vector_score(const struct aln_param* ap,struct hirsch_mem* hm,int old_cor[],float* score)
{
        struct states* f = hm->f;
        struct states* b = hm->b;
        int i;

        float gpo,gpe,tgpe;

        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;


        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;

        //int max = -FLT_MAX;
        float max = -FLT_MAX;
        //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
        float sub = 0.0;

        //i = hm->startb;
        i = old_cor[2];
        //c = -1;
        //for(i = hm->startb; i < hm->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){

                sub = fabsf(middle -i);

                sub /= 1000;
                //	fprintf(stderr,"%d-%d	%f\n",hm->startb,hm->endb,sub);
                if(f[i].a+b[i].a-sub > max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);

                        //              c = i;
                }
                if(f[i].a+b[i].ga-gpo-sub > max){
                        max = f[i].a+b[i].ga-gpo-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);

                        //c = i;
                }
                if(f[i].a+b[i].gb -gpo-sub > max){
                        max = f[i].a+b[i].gb - gpo-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);

                        //c = i;
                }
                if(f[i].ga+b[i].a - gpo-sub > max){
                        max = f[i].ga+b[i].a - gpo-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);

                        //c = i;
                }


                if(hm->startb == 0){
                        if(f[i].gb+b[i].gb - tgpe-sub > max){
                                max = f[i].gb+b[i].gb -tgpe-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);

                                //      c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb - gpe -sub> max){
                                max = f[i].gb+b[i].gb - gpe-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);

                                //c = i;
                        }
                }
                if(f[i].gb+b[i].a - gpo-sub > max){
                        max = f[i].gb+b[i].a - gpo-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);

                        //c = i;
                }
        }
        //i = hm->endb;
        i = old_cor[3];
        sub = fabsf(middle -i);
        sub /= 1000;

        if(f[i].a+b[i].gb-gpo-sub > max){
                max = f[i].a+b[i].gb - gpo-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);

                //c = i;
        }
        if(hm->endb == hm->len_b){
                if(f[i].gb+b[i].gb -tgpe-sub > max){
                        max = f[i].gb+b[i].gb - tgpe-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);

                        //c = i;
                }
        }else{
                if(f[i].gb+b[i].gb - gpe-sub > max){
                        max = f[i].gb+b[i].gb - gpe-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);

                        //c = i;
                }
        }
        *score  = max;
        return OK;
}


int hirsch_align_two_ss_vector(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
{
        struct states* f = hm->f;
        struct states* b = hm->b;
        int i,c;
        int transition = -1;

        float gpo,gpe,tgpe;

        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;


        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;

        //int max = -FLT_MAX;
        float max = -FLT_MAX;
        //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
        float sub = 0.0;

        //i = hm->startb;
        i = old_cor[2];
        c = -1;
        //for(i = hm->startb; i < hm->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){

                sub = fabsf(middle -i);
                sub /= 1000;
                //	fprintf(stderr,"%d-%d	%f\n",hm->startb,hm->endb,sub);
                if(f[i].a+b[i].a-sub > max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga-gpo-sub > max){
                        max = f[i].a+b[i].ga-gpo-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb -gpo-sub > max){
                        max = f[i].a+b[i].gb - gpo-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a - gpo-sub > max){
                        max = f[i].ga+b[i].a - gpo-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(hm->startb == 0){
                        if(f[i].gb+b[i].gb - tgpe-sub > max){
                                max = f[i].gb+b[i].gb -tgpe-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb - gpe -sub> max){
                                max = f[i].gb+b[i].gb - gpe-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a - gpo-sub > max){
                        max = f[i].gb+b[i].a - gpo-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        //i = hm->endb;
        i = old_cor[3];
        sub = fabsf(middle -i);
        sub /= 1000;

        if(f[i].a+b[i].gb-gpo-sub > max){
                max = f[i].a+b[i].gb - gpo-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(hm->endb == hm->len_b){
                if(f[i].gb+b[i].gb -tgpe-sub > max){
                        max = f[i].gb+b[i].gb - tgpe-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb - gpe-sub > max){
                        max = f[i].gb+b[i].gb - gpe-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }


        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
        //LOG_MSG("MAX: %f",max);
        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                hirsch_path[old_cor[4]] = c;
                hirsch_path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;
                //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
                break;
        case 2:// a -> ga = 2

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;


                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);


                //backward:
                hm->starta = old_cor[4];
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = 0.0;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
                break;
        case 3:// a -> gb = 3

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
                break;
        case 5://ga -> a = 5
                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = 0.0;
                hm->b[0].gb = -FLT_MAX;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4];

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
                break;
        case 6://gb->gb = 6;

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
                break;
        case 7://gb->a = 7;

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
                break;
        default:
                break;

        }

        return OK;//hirsch_path;
//ERROR:
        //      return NULL;
}



int foward_hirsch_ss_dyn(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm)
{
        struct states* s = hm->f;
        float *subp = 0;
        const int starta = hm->starta;
        const int enda = hm->enda;
        const int startb =hm->startb;
        const int endb = hm->endb;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register float xa = 0;
        register float xga = 0;
        register int i = 0;
        register int j = 0;

        float gpo,gpe,tgpe;
        float** subm = NULL;
        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;
        subm = ap->subm;


        s[startb].a = s[0].a;
        s[startb].ga = s[0].ga;
        s[startb].gb = s[0].gb;
        if(startb){
                for (j = startb+1; j < endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga - gpe,s[j-1].a-gpo);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for (j = startb+1; j < endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpe;
                        s[j].gb = -FLT_MAX;
                }
        }
        s[endb].a = -FLT_MAX;
        s[endb].ga = -FLT_MAX;
        s[endb].gb = -FLT_MAX;

        seq2--;
        for (i = starta;i < enda;i++){
                subp = subm[seq1[i]];

                pa = s[startb].a;
                pga = s[startb].ga;
                pgb = s[startb].gb;
                s[startb].a = -FLT_MAX;
                s[startb].ga = -FLT_MAX;

                xa = s[startb].a;
                xga = s[startb].ga;

                if(startb){
                        s[startb].gb = MAX(pgb - gpe,pa - gpo);
                }else{
                        s[startb].gb = MAX(pgb,pa) - tgpe;
                }
                for (j = startb+1; j < endb;j++){
                        ca = s[j].a;
                        pa = MAX3(pa,pga-gpo,pgb-gpo);
                        pa += subp[seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;
                        //s[j].ga = MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
                        s[j].ga = MAX(xga-gpe,xa-gpo);

                        pgb = s[j].gb;
                        s[j].gb = MAX(pgb-gpe ,ca-gpo);

                        pa = ca;

                        xa = s[j].a;
                        xga = s[j].ga;

                }
                ca = s[j].a;
                pa = MAX3(pa,pga-gpo,pgb-gpo);
                pa += subp[seq2[j]];

                s[j].a = pa;

                s[j].ga = -FLT_MAX;//MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
                if (endb != hm->len_b){
                        s[j].gb = MAX(s[j].gb-gpe ,ca-gpo);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)-tgpe;
                }

        }
        return OK;
}

int backward_hirsch_ss_dyn(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm)
{

        struct states* s = hm->b;
        float** subm = NULL;
        float *subp = NULL;
        const int starta = hm->starta;
        const int enda = hm->enda;
        const int startb =hm->startb;
        const int endb = hm->endb;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;
        float gpo,gpe,tgpe;

        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;
        subm = ap->subm;

        s[endb].a = s[0].a ;
        s[endb].ga = s[0].ga;
        s[endb].gb = s[0].gb;


        //init of first row;

        //j = endb-startb;
        if(endb != hm->len_b){
                for(j = endb-1;j > startb;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for(j = endb-1;j > startb;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpe;
                        s[j].gb = -FLT_MAX;
                }
        }


        s[startb].a = -FLT_MAX;
        s[startb].ga = -FLT_MAX;
        s[startb].gb = -FLT_MAX;

        i = enda-starta;
        seq1+= starta;
        while(i--){
                subp = subm[seq1[i]];
                pa = s[endb].a;
                pga = s[endb].ga;
                pgb = s[endb].gb;
                s[endb].a = -FLT_MAX;
                s[endb].ga = -FLT_MAX;

                xa = s[endb].a;
                xga = s[endb].ga;

                if(endb != hm->len_b){
                        s[endb].gb = MAX(pgb-gpe,pa-gpo);
                }else{
                        s[endb].gb = MAX(pgb,pa)-tgpe;
                }

                for(j = endb-1;j > startb;j--){

                        ca = s[j].a;

                        pa = MAX3(pa,pga - gpo,pgb-gpo);

                        pa += subp[seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

                        s[j].ga = MAX(xga-gpe,xa-gpo);

                        pgb = s[j].gb;
                        s[j].gb = MAX(pgb-gpe,ca-gpo);

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                ca = s[j].a;

                pa = MAX3(pa,pga - gpo,pgb-gpo);

                pa += subp[seq2[j]];

                s[j].a = pa;

                s[j].ga = -FLT_MAX;//MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

                if(startb){
                        s[j].gb = MAX(s[j].gb-gpe,ca-gpo);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)-tgpe;
                }
        }
        return OK;
}


int hirsch_ps_dyn(const struct aln_param* ap, const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip)
{
        int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
        float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
        int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};


        if(hm->starta  >= hm->enda){
                return OK;
        }
        if(hm->startb  >= hm->endb){
                return OK;
        }

        hm->enda = mid;
        foward_hirsch_ps_dyn(ap,prof1,seq2,hm,sip);

        /*int i;
          fprintf(stderr,"FOWARD\n");
          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
          }*/

        hm->starta = mid;
        hm->enda = old_cor[1];
        backward_hirsch_ps_dyn(ap,prof1,seq2,hm,sip);

        /*fprintf(stderr,"BaCKWARD\n");
          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
          }*/

        hirsch_align_two_ps_vector(ap,prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);
        return OK;
}



int hirsch_align_two_ps_vector(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)
{
        struct states* f = hm->f;
        struct states* b = hm->b;
        int i,c;
        int transition = -1;

        const float open = ap->gpo * sip;


        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;

        //int max = -FLT_MAX;
        float max = -FLT_MAX;
        //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
        float sub = 0.0;


        prof1+= ((old_cor[4]+1)<<6);

        //i = hm->startb;
        i = old_cor[2];
        c = -1;
        //for(i = hm->startb; i < hm->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){
                sub = fabsf(middle -i);
                sub /= 1000;
                if(f[i].a+b[i].a-sub> max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga-open-sub > max){
                        max = f[i].a+b[i].ga-open-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb+prof1[27]-sub > max){
                        max = f[i].a+b[i].gb+prof1[27]-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a-open-sub > max){
                        max = f[i].ga+b[i].a-open-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(hm->startb == 0){
                        if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[29]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[28]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a+prof1[-37]-sub > max){
                        max = f[i].gb+b[i].a+prof1[-37]-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        //i = hm->endb;
        i = old_cor[3];

        sub = fabsf(middle -i);
        sub /= 1000;
        if(f[i].a+b[i].gb+prof1[27]-sub > max){
                max = f[i].a+b[i].gb+prof1[27]-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(hm->endb == hm->len_b){
                if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[29]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[28]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }



        prof1-= ((old_cor[4]+1)<<6);

        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                hirsch_path[old_cor[4]] = c;
                hirsch_path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;
                //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
                break;
        case 2:// a -> ga = 2

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;


                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4];
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = 0.0;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
                break;
        case 3:// a -> gb = 3

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
                break;
        case 5://ga -> a = 5

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = 0.0;
                hm->b[0].gb = -FLT_MAX;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4];

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
                break;
        case 6://gb->gb = 6;

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);


                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
                break;
        case 7://gb->a = 7;

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
                break;
        default:
                break;

        }

        return OK;
}

int foward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip)
{
        struct states* s = hm->f;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;

        const float open = ap->gpo * sip;
        const float ext = ap->gpe *sip;
        const float text = ap->tgpe * sip;



        prof1 += (hm->starta)<< 6;
        s[hm->startb].a = s[0].a;
        s[hm->startb].ga = s[0].ga;
        s[hm->startb].gb = s[0].gb;
        if(hm->startb){
                for (j = hm->startb+1; j < hm->endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for (j = hm->startb+1; j < hm->endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;
                        s[j].gb = -FLT_MAX;
                }

        }


        s[hm->endb].a = -FLT_MAX;
        s[hm->endb].ga = -FLT_MAX;
        s[hm->endb].gb = -FLT_MAX;
        seq2--;

        for (i = hm->starta;i < hm->enda;i++){
                prof1 += 64;

                pa = s[hm->startb].a;
                pga = s[hm->startb].ga;
                pgb = s[hm->startb].gb;
                s[hm->startb].a = -FLT_MAX;
                s[hm->startb].ga = -FLT_MAX;

                xa = s[hm->startb].a;
                xga = s[hm->startb].ga;


                if(hm->startb){
                        s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
                }else{
                        s[hm->startb].gb = MAX(pgb,pa)+prof1[29];
                }
                for (j = hm->startb+1; j < hm->endb;j++){
                        ca = s[j].a;

                        pa = MAX3(pa,pga -open,pgb + prof1[-37]);

                        pa += prof1[32 + seq2[j]];


                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
                        s[j].ga = MAX(xga-ext,xa-open);


                        pgb = s[j].gb;

                        s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;

                }
                ca = s[j].a;

                pa = MAX3(pa,pga -open,pgb + prof1[-37]);

                pa += prof1[32 + seq2[j]];


                s[j].a = pa;

                s[j].ga = -FLT_MAX;//MAX(s[j-1].ga-ext,s[j-1].a-open);

                if (hm->endb != hm->len_b){
                        s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
                }

        }
        prof1 -= hm->enda << 6;
        return OK;
}

int backward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip)
{
        struct states* s = hm->b;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;

        const float open = ap->gpo * sip;
        const float ext = ap->gpe *sip;
        const float text = ap->tgpe * sip;


        prof1 += (hm->enda+1) << 6;

        s[hm->endb].a = s[0].a;
        s[hm->endb].ga = s[0].ga;
        s[hm->endb].gb = s[0].gb;

        if(hm->endb != hm->len_b){
                for(j = hm->endb-1;j > hm->startb;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for(j = hm->endb-1;j > hm->startb;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;
                        s[j].gb = -FLT_MAX;
                }
        }

        s[hm->startb].a = -FLT_MAX;
        s[hm->startb].ga = -FLT_MAX;
        s[hm->startb].gb = -FLT_MAX;

        i = hm->enda-hm->starta;
        while(i--){
                prof1 -= 64;
                pa = s[hm->endb].a;
                pga = s[hm->endb].ga;
                pgb = s[hm->endb].gb;
                s[hm->endb].a = -FLT_MAX;
                s[hm->endb].ga = -FLT_MAX;

                xa = s[hm->endb].a;
                xga = s[hm->endb].ga;


                if(hm->endb != hm->len_b){
                        s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
                }else{
                        s[hm->endb].gb = MAX(pgb,pa) +prof1[29];
                }

                for(j = hm->endb-1;j > hm->startb;j--){
                        ca = s[j].a;

                        pa = MAX3(pa,pga - open,pgb +prof1[91]);
                        pa += prof1[32 + seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
                        s[j].ga = MAX(xga-ext,xa-open);

                        pgb = s[j].gb;

                        s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;


                }
                ca = s[j].a;

                pa = MAX3(pa,pga - open,pgb +prof1[91]);
                pa += prof1[32 + seq2[j]];

                s[j].a = pa;


                s[j].ga = -FLT_MAX;//MAX(s[j+1].ga-ext,s[j+1].a-open);
                if(hm->startb){
                        s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)+prof1[29];
                }

        }
        return OK;
}

int hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)
{
        int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
        float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
        int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};


        //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);


        if(hm->starta  >= hm->enda){
                return OK;
        }
        if(hm->startb  >= hm->endb){
                return OK;
        }

        hm->enda = mid;
        foward_hirsch_pp_dyn(prof1,prof2,hm);
        /*int i;
          fprintf(stderr,"FOWARD\n");
          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
          }*/

        hm->starta = mid;
        hm->enda = old_cor[1];
        backward_hirsch_pp_dyn(prof1,prof2,hm);
        /*fprintf(stderr,"BaCKWARD\n");

          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
          }*/

        hirsch_align_two_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
        return OK;
}



int hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
{
        struct states* f = hm->f;
        struct states* b = hm->b;
        int i,c;
        int transition = -1;


        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;

        //int max = -FLT_MAX;
        float max = -FLT_MAX;
        //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
        float sub = 0.0;


        prof1+= ((old_cor[4]+1) << 6);
        //prof2 += 64 * (hm->startb);
        //i = hm->startb;
        prof2 += old_cor[2] << 6;
        i = old_cor[2];
        c = -1;
        //for(i = hm->startb; i < hm->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){
                sub = fabsf(middle -i);
                sub /= 1000;
                prof2 += 64;
                //fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
                if(f[i].a+b[i].a-sub > max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga+prof2[27]-sub > max){
                        max = f[i].a+b[i].ga+prof2[27]-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb+prof1[27] -sub> max){
                        max = f[i].a+b[i].gb+prof1[27]-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a+prof2[-37]-sub > max){
                        max = f[i].ga+b[i].a+prof2[-37]-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(hm->startb == 0){
                        if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[29]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[28]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a+prof1[-37]-sub > max){
                        max = f[i].gb+b[i].a+prof1[-37]-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        //i = hm->endb;
        i = old_cor[3];
        sub = fabsf(middle -i);
        sub /= 1000;
        if(f[i].a+b[i].gb+prof1[27]-sub > max){
                max = f[i].a+b[i].gb+prof1[27]-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(hm->endb == hm->len_b){
                if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[29]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[28]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }



        prof1-= (old_cor[4]+1)<<6;
        //prof2 -= hm->endb << 6;
        prof2 -= old_cor[3] << 6;

        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
        //if(transition == -1){
        //	exit(0);
        //}

        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                hirsch_path[old_cor[4]] = c;
                hirsch_path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;
                //fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 2:// a -> ga = 2

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;


                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4];
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = 0.0;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 3:// a -> gb = 3

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = -FLT_MAX;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 5://ga -> a = 5
                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = 0.0;
                hm->b[0].gb = -FLT_MAX;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4];

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 6://gb->gb = 6;

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLT_MAX;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 7://gb->a = 7;

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLT_MAX;
                hm->b[0].ga = -FLT_MAX;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        default:
                break;

        }

        return OK;
}

int foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
        unsigned int freq[23];

        struct states* s = hm->f;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;
        register int c = 0;

        register int f = 0;
        prof1 += (hm->starta) << 6;
        prof2 +=  (hm->startb) << 6;
        s[hm->startb].a = s[0].a;
        s[hm->startb].ga = s[0].ga;
        s[hm->startb].gb = s[0].gb;
        if(hm->startb){
                for (j = hm->startb+1; j < hm->endb;j++){
                        prof2+=64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
                        s[j].gb = -FLT_MAX;
                }
                prof2+=64;
        }else{
                for (j = hm->startb+1; j < hm->endb;j++){
                        prof2+=64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];
                        s[j].gb = -FLT_MAX;
                }
                prof2+=64;
        }

        prof2 -= (hm->endb-hm->startb) << 6;

        s[hm->endb].a = -FLT_MAX;
        s[hm->endb].ga = -FLT_MAX;
        s[hm->endb].gb = -FLT_MAX;


        for (i = hm->starta;i < hm->enda;i++){
                prof1 += 64;
                //c = 1;
                f = 0;
                for (j = 0;j < 23; j++){
                        if(prof1[j]){
                                freq[f] = j;
                                f++;
                        }
                }
                f--;
                //freq[0] = c;

                pa = s[hm->startb].a;
                pga = s[hm->startb].ga;
                pgb = s[hm->startb].gb;
                s[hm->startb].a = -FLT_MAX;
                s[hm->startb].ga = -FLT_MAX;

                xa = s[hm->startb].a;
                xga = s[hm->startb].ga;


                if(hm->startb){
                        s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
                }else{
                        s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];
                }
                for (j = hm->startb+1; j < hm->endb;j++){
                        prof2 += 64;
                        ca = s[j].a;

                        pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);

                        prof2 += 32;
                        for (c = f;c >= 0;c--){
                                //for (c = 0;c < f;c++){
                                //for (c = 1;c < freq[0];c++){
                                pa += prof1[freq[c]]*prof2[freq[c]];
                        }
                        prof2 -= 32;

                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
                        s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);

                        pgb = s[j].gb;

                        s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);

                        pa = ca;


                        xa = s[j].a;
                        xga = s[j].ga;
                }
                prof2 += 64;
                ca = s[j].a;

                pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);

                prof2 += 32;
                for (c = f;c >= 0;c--){
                        //for (c = 0;c < f;c++){
                        pa += prof1[freq[c]]*prof2[freq[c]];
                }
                prof2 -= 32;

                s[j].a = pa;

                s[j].ga = -FLT_MAX;

                if (hm->endb != hm->len_b){
                        s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
                }
                prof2 -= (hm->endb-hm->startb) << 6;

        }
        prof1 -=  (hm->enda) << 6;
        return OK;
}

int backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
        unsigned int freq[23];
        struct states* s = hm->b;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;
        register int c = 0;

        register int f = 0;
        prof1 += (hm->enda+1) << 6;
        prof2 += (hm->endb+1) << 6;
        s[hm->endb].a = s[0].a;
        s[hm->endb].ga = s[0].ga;
        s[hm->endb].gb = s[0].gb;
        if(hm->endb != hm->len_b){
                for(j = hm->endb-1;j > hm->startb;j--){
                        prof2 -= 64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);
                        s[j].gb = -FLT_MAX;
                }
                prof2 -= 64;
        }else{
                for(j = hm->endb-1;j > hm->startb;j--){
                        prof2 -= 64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];
                        s[j].gb = -FLT_MAX;
                }
                prof2 -= 64;
        }

        s[hm->startb].a = -FLT_MAX;
        s[hm->startb].ga = -FLT_MAX;
        s[hm->startb].gb = -FLT_MAX;

        i = hm->enda-hm->starta;
        while(i--){
                prof1 -= 64;

                //c = 1;
                f = 0;
                for (j = 0;j < 23; j++){
                        if(prof1[j]){
                                freq[f] = j;
                                f++;
                        }
                }
                f--;
                //freq[0] = c;

                pa = s[hm->endb].a;
                pga = s[hm->endb].ga;
                pgb = s[hm->endb].gb;
                s[hm->endb].a = -FLT_MAX;
                s[hm->endb].ga = -FLT_MAX;

                xa = s[hm->endb].a;
                xga = s[hm->endb].ga;

                if(hm->endb != hm->len_b){
                        s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);
                }else{
                        s[hm->endb].gb = MAX(pgb,pa)+prof1[29];
                }

                prof2 += (hm->endb-hm->startb) << 6;
                for(j = hm->endb-1;j > hm->startb;j--){
                        prof2 -= 64;
                        ca = s[j].a;

                        pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);

                        prof2 += 32;
                        for (c = f;c >= 0;c--){
                                //for (c = 0;c < f;c++){
                                //for (c = 1;c < freq[0];c++){
                                pa += prof1[freq[c]]*prof2[freq[c]];
                        }
                        prof2 -= 32;

                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
                        s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);

                        pgb = s[j].gb;

                        s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                prof2 -= 64;
                ca = s[j].a;

                pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
                prof2 += 32;
                for (c = f;c >= 0;c--){
                        //for (c = 0;c < f;c++){
                        pa += prof1[freq[c]]*prof2[freq[c]];
                }
                prof2 -= 32;
                s[j].a = pa;

                //pga = s[j].ga;
                s[j].ga = -FLT_MAX;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

                //pgb = s[j].gb;
                if(hm->startb){
                        s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)+prof1[29];
                }

                //pa = ca;
        }
        return OK;
}




int* mirror_hirsch_path(int* hirsch_path,int len_a,int len_b)
{
        int* np =NULL;

        int i;

        MMALLOC(np,sizeof(int)*(len_a+2));
        for(i =0; i < len_a+2;i++){
                np[i] = -1;
        }

        for(i = 1; i <= len_b;i++){
                if(hirsch_path[i] != -1){
                        np[hirsch_path[i]] = i;
                }
        }

        MFREE(hirsch_path);
        return np;
ERROR:
        return NULL;
}

int* add_gap_info_to_hirsch_path(int* hirsch_path,int len_a,int len_b)
{
        int i,j;
        int a = 0;
        int b = 0;

        int* np = NULL;
        MMALLOC(np,sizeof(int)*(len_a+len_b+2));
        for(i =0; i < len_a+len_b+2;i++){
                np[i] = 0;
        }

        j = 1;
        b = -1;
        if(hirsch_path[1] == -1){
                np[j] = 2;
                j++;
        }else{
                if(hirsch_path[1] != 1){
                        for ( a = 0;a < hirsch_path[1] -1;a++){
                                np[j] = 1;
                                j++;
                        }
                        np[j] = 0;
                        j++;
                }else{
                        np[j] = 0;
                        j++;
                }
        }
        b = hirsch_path[1];

        /*for ( i= 0;i <= len_a;i++){
          fprintf(stderr,"%d,",hirsch_path[i]);
          }
          fprintf(stderr,"\n");*/

        for(i = 2; i <= len_a;i++){

                if(hirsch_path[i] == -1){
                        np[j] = 2;
                        j++;
                }else{
                        if(hirsch_path[i]-1 != b && b != -1){
                                for ( a = 0;a < hirsch_path[i] - b-1;a++){
                                        np[j] = 1;
                                        j++;
                                }
                                np[j] = 0;
                                j++;
                        }else{
                                np[j] = 0;
                                j++;
                        }
                }
                b = hirsch_path[i];
        }





        if(hirsch_path[len_a] < len_b && hirsch_path[len_a] != -1){
                //	fprintf(stderr,"WARNING:%d	%d\n",hirsch_path[len_a],len_b);
                for ( a = 0;a < len_b - hirsch_path[len_a];a++){
                        np[j] = 1;
                        j++;
                }
        }
        np[0] = j-1;
        np[j] = 3;

        MREALLOC(np,sizeof(int)* (np[0]+2));
        //for ( i= 0;i <= np[0];i++){
        //	fprintf(stderr,"%d,",np[i]);
        //}
        //fprintf(stderr,"\n");

        MFREE(hirsch_path);

        //add gap info..
        i = 2;
        while(np[i] != 3){
                if ((np[i-1] &3) && !(np[i] & 3)){
                        if(np[i-1] & 8){
                                np[i-1] += 8;
                        }else{
                                np[i-1] |= 16;
                        }
                }else if (!(np[i-1] & 3) &&(np[i] &3)){
                        np[i] |= 4;
                }else if ((np[i-1] & 1) && (np[i] & 1)){
                        np[i] |= 8;
                }else if ((np[i-1] & 2) && (np[i] & 2)){
                        np[i] |= 8;
                }
                i++;
        }
        //add terminal gap...
        i = 1;
        while(np[i] != 0){
                np[i] |= 32;
                i++;
        }
        j = i;
        i = np[0];
        while(np[i] != 0){
                np[i] |= 32;
                i--;
        }
        //for ( i= 0;i <= np[0];i++){
        //	fprintf(stderr,"%d,",np[i]);
        //}
        //fprintf(stderr,"\n");
        return np;
ERROR:
        return NULL;
}


float* make_profile(struct aln_param* ap,const uint8_t* seq,const int len)
{
        int i,j,c;
        float** subm = NULL;
        float* prof = NULL;
        float gpo,gpe,tgpe;

        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;
        subm = ap->subm;

        MMALLOC(prof,sizeof(float)*(len+2)*64);
        prof +=  (64 *(len+1));

        for (i = 0;i < 64;i++){
                prof[i] = 0;
        }
        prof[23+32] = -gpo;
        prof[24+32] = -gpe;
        prof[25+32] = -tgpe;

        i = len;
        while(i--){
                prof -= 64;

                for (j = 0;j < 64;j++){
                        prof[j] = 0;
                }
                c = seq[i];

                prof[c] += 1;

                prof += 32;

                for(j = 21;j--;){
                        prof[j] = subm[c][j];
                }
                prof[23] = -gpo;
                prof[24] = -gpe;
                prof[25] = -tgpe;

                prof -= 32;
        }
        prof -= 64;
        for (i = 0;i < 64;i++){
                prof[i] = 0;
        }
        prof[23+32] = -gpo;
        prof[24+32] = -gpe;
        prof[25+32] = -tgpe;
        return prof;
ERROR:
        return NULL;
}


int set_gap_penalties(float* prof,int len,int nsip)
{
        int i;
        //int j;

        prof +=  (64 *(len+1));

        prof[27] = prof[55]*nsip;//gap open or close  23
        prof[28] = prof[56]*nsip;//gap extention 24

        prof[29] = prof[57]*nsip;//gap open or close 25
        i = len+1;
        while(i--){
                prof -= 64;
                prof[27] = prof[55]*nsip;//gap open or close
                prof[28] = prof[56]*nsip;//gap extention

                prof[29] = prof[57]*nsip;//gap open or close
        }
        return OK;
}


int update(const float* profa, const float* profb,float* newp, struct aln_param*ap, int* path,int sipa,int sipb)
{
        int i,j,c;
        for (i = 64; i--;){
                newp[i] = profa[i] + profb[i];
        }

        profa += 64;
        profb += 64;
        newp += 64;

        c = 1;

        while(path[c] != 3){
                //Idea: limit the 'virtual' number of residues of one type to x.
                // i.e. only allow a maximum of 10 alanines to be registered in each column
                // the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
                // the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase
                // with the number of sequences. -> see Durbin pp 140

                if (!path[c]){
                        //fprintf(stderr,"Align	%d\n",c);
                        for (i = 64; i--;){
                                newp[i] = profa[i] + profb[i];
                        }
                        profa += 64;
                        profb += 64;
                }

                if (path[c] & 1){
                        //fprintf(stderr,"Gap_A:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
                        for (i = 64; i--;){
                                newp[i] = profb[i];
                        }
                        profb += 64;
                        if(!(path[c] & 20)){
                                if(path[c] & 32){
                                        newp[25] += sipa;//1;
                                        i = ap->tgpe*sipa;
                                }else{
                                        newp[24] += sipa;//1;
                                        i = ap->gpe*sipa;
                                }

                                for (j = 32; j < 55;j++){
                                        newp[j] -=i;
                                }
                        }else{
                                if (path[c] & 16){
                                        //			fprintf(stderr,"close_open");
                                        if(path[c] & 32){
                                                newp[25] += sipa;//1;
                                                i = ap->tgpe*sipa;
                                                newp[23] += sipa;//1;
                                                i += ap->gpo*sipa;
                                        }else{
                                                newp[23] += sipa;//1;
                                                i = ap->gpo*sipa;
                                        }

                                        for (j = 32; j < 55;j++){
                                                newp[j] -=i;
                                        }
                                }
                                if (path[c] & 4){
                                        //			fprintf(stderr,"Gap_open");
                                        if(path[c] & 32){
                                                newp[25] += sipa;//1;
                                                i = ap->tgpe*sipa;
                                                newp[23] += sipa;//1;
                                                i += ap->gpo*sipa;
                                        }else{
                                                newp[23] += sipa;//1;
                                                i = ap->gpo*sipa;
                                        }
                                        for (j = 32; j < 55;j++){
                                                newp[j] -=i;
                                        }
                                }
                        }
                }
                if (path[c] & 2){
                        //fprintf(stderr,"Gap_B:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
                        for (i = 64; i--;){
                                newp[i] = profa[i];
                        }
                        profa+=64;
                        if(!(path[c] & 20)){
                                if(path[c] & 32){
                                        newp[25] += sipb;//1;
                                        i = ap->tgpe*sipb;
                                }else{
                                        newp[24] += sipb;//1;
                                        i = ap->gpe*sipb;
                                }
                                for (j = 32; j < 55;j++){
                                        newp[j] -=i;
                                }
                        }else{
                                if (path[c] & 16){
                                        //			fprintf(stderr,"close_open");
                                        if(path[c] & 32){
                                                newp[25] += sipb;//1;
                                                i =  ap->tgpe*sipb;
                                                newp[23] += sipb;//1;
                                                i +=  ap->gpo*sipb;
                                        }else{
                                                newp[23] += sipb;//1;
                                                i =  ap->gpo*sipb;
                                        }
                                        for (j = 32; j < 55;j++){
                                                newp[j] -=i;
                                        }
                                }
                                if (path[c] & 4){
                                        //			fprintf(stderr,"Gap_open");
                                        if(path[c] & 32){
                                                newp[25] += sipb;//1;
                                                i = ap->tgpe*sipb;
                                                newp[23] += sipb;//1;
                                                i += ap->gpo*sipb;
                                        }else{
                                                newp[23] += sipb;//1;
                                                i = ap->gpo*sipb;
                                        }

                                        for (j = 32; j < 55;j++){
                                                newp[j] -=i;
                                        }
                                }
                        }
                }
                newp += 64;
                c++;
        }
        for (i = 64; i--;){
                newp[i] =  profa[i] + profb[i];
        }
        newp -= (path[0]+1) *64;
        return OK;
}

/*int update(const float* profa, const float* profb,float* newp,const int* path)
{
        int i,c;
        //float gpo,gpe,tgpe;


        //gpo = ap->gpo;
        //gpe = ap->gpe;
        //tgpe = ap->tgpe;


        for (i = 64; i--;){
                newp[i] = profa[i] + profb[i];
        }

        profa += 64;
        profb += 64;
        newp += 64;

        c = 1;

        while(path[c] != 3){
                //Idea: limit the 'virtual' number of residues of one type to x.
                // i.e. only allow a maximum of 10 alanines to be registered in each column
                // the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
                // the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase
                // with the number of sequences. -> see Durbin pp 140

                if (!path[c]){
                        //fprintf(stderr,"Align	%d\n",c);
                        for (i = 64; i--;){
                                newp[i] = profa[i] + profb[i];
                        }


                        profa += 64;
                        profb += 64;
                }

                if (path[c] & 1){
                        //fprintf(stderr,"Gap_A:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
                        for (i = 64; i--;){
                                newp[i] = profb[i];
                        }
                        profb += 64;


                }
                if (path[c] & 2){
                        //fprintf(stderr,"Gap_B:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
                        for (i = 64; i--;){
                                newp[i] = profa[i];
                        }
                        profa+=64;
                }
                newp += 64;
                c++;
        }
        for (i = 64; i--;){
                newp[i] =  profa[i] + profb[i];
        }
        newp -= (path[0]+1) *64;
        return OK;
        }*/

int calc_seq_id(const int* path,int* a, int*b,float* dist)
{
        int i,j,c;
        float id = 0.0f;
        float len = 0.0f;
        //float gpo,gpe,tgpe;


        //gpo = ap->gpo;
        //gpe = ap->gpe;
        //tgpe = ap->tgpe;
        i = 0;
        j = 0;
        c = 1;

        while(path[c] != 3){
                //Idea: limit the 'virtual' number of residues of one type to x.
                // i.e. only allow a maximum of 10 alanines to be registered in each column
                // the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
                // the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase
                // with the number of sequences. -> see Durbin pp 140
                //


                //fprintf(stdout,"%d\t%d",a[i],b[j]);
                if (!path[c]){
                        if(a[i] == b[j]){
                                id++;
                                //              fprintf(stdout,"\t%f",id);
                        }

                        i++;
                        j++;
                }

                if (path[c] & 1){
                        //fprintf(stderr,"Gap_A:%d\n",c);
                        j++;
                }
                if (path[c] & 2){
                        i++;
                //fprintf(stderr,"Gap_B:%d\n",c);
                }
                len += 1.0f;
                //fprintf(stdout,"\n");
                c++;
        }
        *dist = id / len;
        return OK;
}





struct hirsch_mem* hirsch_mem_alloc(int x)
{

        struct hirsch_mem* hm = NULL;
        // a=((typeof(a))(((int)(((void *)malloc(c+15))+15))&-16)).
        MMALLOC(hm,sizeof(struct hirsch_mem));
        hm->starta = 0;
        hm->startb = 0;
        hm->enda = 0;
        hm->endb = 0;
        hm->size = x+1;
        hm->len_a = 0;
        hm->len_b = 0;
        hm->f = NULL;
        hm->b = NULL;
        MMALLOC(hm->f,sizeof(struct states)* (x+1));
        MMALLOC(hm->b,sizeof(struct states)* (x+1));
        return hm;
ERROR:
        return NULL;
}

int hirsch_mem_realloc(struct hirsch_mem* hm,int x)
{
        if((x+1) > hm->size){
                hm->size = x+1;

                MREALLOC(hm->f,sizeof(struct states)* (x+1));
                MREALLOC(hm->b,sizeof(struct states)* (x+1));
        }
        return OK;
ERROR:
        return FAIL;
}

void hirsch_mem_free(struct hirsch_mem* hm)
{
        if(hm){
                MFREE(hm->f);
                MFREE(hm->b);
                MFREE(hm);
        }
}
