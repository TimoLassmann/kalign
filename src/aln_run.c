#include "tldevel.h"
#include "tlrng.h"

#include "msa.h"
#include "alignment_parameters.h"

#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"

#define ALN_RUN_IMPORT
#include "aln_run.h"

static int score_aln(struct aln_mem* m,struct aln_param* ap,float** profile, struct msa* msa, int a,int b,int numseq,float* score);
static int SampleWithoutReplacement(struct rng_state* rng, int N, int n,int* samples);
static int int_cmp(const void *a, const void *b);

int** create_chaos_msa(struct msa* msa, struct aln_param* ap)
{
        struct aln_mem* m = NULL;
        int i,j,g,f,a,b,c,l;
        int best_a, best_b;
        int len_a;
        int len_b;
        float** profile = NULL;
        int* samples = NULL;
        int** map = NULL;
        int* tree = NULL;

        int* active = NULL;
        float max_score;
        float score;
        int numseq;

        ap->mode = ALN_MODE_FULL;

        g = msa->num_profiles;
        numseq = msa->numseq;
        MMALLOC(profile,sizeof(float*)*g);
        MMALLOC(map,sizeof(int*)*g);
        for ( i = 0;i< g;i++){
                profile[i] = NULL;
                map[i] = NULL;
        }

        RUN(alloc_aln_mem(&m, 2048));
        MMALLOC(samples,sizeof(int) * 6);
        MMALLOC(active, sizeof(int) * numseq);
        for(i = 0; i < numseq;i++){
                active[i] = i;
        }
        qsort(active, numseq, sizeof(int), int_cmp);

        for(i = 0; i < numseq-1;i++){
                /* pick one sequence / profile  */
                max_score = -FLT_MAX;
                l = MACRO_MIN(6, numseq-1);
                SampleWithoutReplacement(ap->rng, numseq-i, l, samples);

                for(g = 0;g < l-1;g++){
                        a = samples[g];
                        for(f = g + 1; f < l;f++){
                                b = samples[f];
                                score_aln(m, ap, profile, msa, active[a], active[b], numseq, &score);
                                //LOG_MSG("TEsting %d %d : %f", a,b, ap->score);
                                if(ap->score > max_score){
                                        best_a = a;
                                        best_b = b;
                                        max_score = ap->score;
                                }
                        }
                }
                /* //LOG_MSG("L:%d", l); */
                /* for(g = 0; g < l;g++){ */
                /*         a = tl_random_int(ap->rng, numseq-i); */
                /*         b = tl_random_int(ap->rng, numseq-i); */
                /*         while(b == a){ */
                /*                 b = tl_random_int(ap->rng, numseq-i); */
                /*         } */
                /*         score_aln(m, ap, profile, msa, active[a], active[b], numseq, &score); */
                /*         //LOG_MSG("TEsting %d %d : %f", a,b, ap->score); */
                /*         if(ap->score > max_score){ */
                /*                 best_a = a; */
                /*                 best_b = b; */
                /*                 max_score = ap->score; */
                /*         } */
                /* } */

                a = best_a;
                b = best_b;
                //LOG_MSG("samples: %d %d", active[a],active[b]);

                ap->tree[i*3] = active[a];
                ap->tree[i*3+1] = active[b];
                ap->tree[i*3+2] = numseq+i;

                active[a] = numseq+i;
                active[b] = -1;
                qsort(active, numseq-i, sizeof(int), int_cmp);

                a = ap->tree[i*3];
                b = ap->tree[i*3+1];
                c = ap->tree[i*3+2];

                //score_aln(m, ap, profile, msa, a, b, numseq, &score);
                //fprintf(stdout,"Aligning:%d %d->%d	done:%f score:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100,score);
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

                RUN(resize_aln_mem(m, g));


                for (j = 0; j < (g+2);j++){
                        map[c][j] = -1;
                }

                if (a < numseq){
                        RUN(make_profile_n(ap, msa->sequences[a]->s,len_a,&profile[a]));
                        //RUNP(profile[a] = make_profile(ap,msa->sequences[a]->s,len_a));
                }else{
                        RUN(set_gap_penalties_n(profile[a],len_a,msa->nsip[b]));
                        //RUN(set_gap_penalties(profile[a],len_a,msa->nsip[b]));
                        //smooth_gaps(profile[a],len_a,window,strength);

                        //increase_gaps(profile[a],len_a,window,strength);
                }
                if (b < numseq){
                        RUN(make_profile_n(ap, msa->sequences[b]->s,len_b,&profile[b]));
                        //RUNP(profile[b] = make_profile(ap,msa->sequences[b]->s,len_b));
                }else{
                        RUN(set_gap_penalties_n(profile[b],len_b,msa->nsip[a]));
                        //smooth_gaps(profile[b],len_b,window,strength);
                        //increase_gaps(profile[b],len_b,window,strength);
                }

                init_alnmem(m, len_a, len_b);
                //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
                if(a < numseq){
                        if(b < numseq){
                                ap->seq1 = msa->sequences[a]->s;
                                ap->seq2 = msa->sequences[b]->s;
                                ap->prof1 = NULL;
                                ap->prof2 = NULL;

                                /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", ap->score); */
                                ap->mode = ALN_MODE_FULL;
                                aln_runner(m, ap, map[c]);
                                //hirsch_ss_dyn(ap,msa->sequences[a]->s, msa->sequences[b]->s,hm,map[c]);
                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                ap->seq1 = NULL;
                                ap->seq2 = msa->sequences[a]->s;
                                ap->prof1 = profile[b];
                                ap->prof2 = NULL;
                                ap->sip = msa->nsip[b];

                                /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", ap->score); */
                                ap->mode = ALN_MODE_FULL;
                                aln_runner(m, ap, map[c]);

                                //hirsch_ps_dyn(ap,profile[b], msa->sequences[a]->s,hm,map[c],msa->nsip[b]);
                                RUN(mirror_path_n(&map[c],len_a,len_b));
                                //RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                        }
                }else{
                        if(b < numseq){
                                ap->seq1 = NULL;
                                ap->seq2 = msa->sequences[b]->s;
                                ap->prof1 = profile[a];
                                ap->prof2 = NULL;
                                ap->sip = msa->nsip[a];
                                /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", ap->score); */
                                ap->mode = ALN_MODE_FULL;
                                aln_runner(m, ap, map[c]);

                                //hirsch_ps_dyn(ap,profile[a],msa->sequences[b]->s ,hm,map[c],msa->nsip[a]);
                        }else{
                                if(len_a < len_b){
                                        ap->seq1 = NULL;
                                        ap->seq2 = NULL;
                                        ap->prof1 = profile[a];
                                        ap->prof2 = profile[b];
                                        /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                        /* aln_runner(m, ap, map[c]); */
                                        /* LOG_MSG("SCORE: %f", ap->score); */
                                        ap->mode = ALN_MODE_FULL;
                                        aln_runner(m, ap, map[c]);
                                        //hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
                                }else{
                                        m->enda = len_b;
                                        m->endb = len_a;
                                        m->len_a = len_b;
                                        m->len_b = len_a;

                                        ap->seq1 = NULL;
                                        ap->seq2 = NULL;
                                        ap->prof1 = profile[b];
                                        ap->prof2 = profile[a];
                                        /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                        /* aln_runner(m, ap, map[c]); */
                                        /* LOG_MSG("SCORE: %f", ap->score); */
                                        ap->mode = ALN_MODE_FULL;
                                        aln_runner(m, ap, map[c]);
                                        //hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
                                        RUN(mirror_path_n(&map[c],len_a,len_b));
                                        //RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                                }
                        }
                }

                RUN(add_gap_info_to_path_n(&map[c], len_a, len_b));
                //map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

                if(i != numseq-2){
                        //MREALLOC(profile_ptr, sizeof(float)*64*(map[c][0]+2));
                        MMALLOC(profile[c],sizeof(float)*64*(map[c][0]+2));
                        //update(profile[a],profile[b],profile[c],map[c]);
                        update_n(profile[a],profile[b],profile[c],ap,map[c],msa->nsip[a],msa->nsip[b]);
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



        return map;
ERROR:
        if(m){
                free_aln_mem(m);
        }
        return NULL;
}


int** create_msa(struct msa* msa, struct aln_param* ap)
{
        struct aln_mem* m = NULL;
        int i,j,g,a,b,c;
        int len_a;
        int len_b;
        float** profile = NULL;

        int** map = NULL;
        int* tree = NULL;
        int numseq;

        ap->mode = ALN_MODE_FULL;

        g = msa->num_profiles;
        numseq = msa->numseq;

        tree = ap->tree;
        MMALLOC(profile,sizeof(float*)*g);
        MMALLOC(map,sizeof(int*)*g);
        for ( i = 0;i< g;i++){
                profile[i] = NULL;
                map[i] = NULL;
        }

        RUN(alloc_aln_mem(&m, 2048));

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

                RUN(resize_aln_mem(m, g));


                for (j = 0; j < (g+2);j++){
                        map[c][j] = -1;
                }

                if (a < numseq){
                        RUN(make_profile_n(ap, msa->sequences[a]->s,len_a,&profile[a]));
                        //RUNP(profile[a] = make_profile(ap,msa->sequences[a]->s,len_a));
                }else{
                        RUN(set_gap_penalties_n(profile[a],len_a,msa->nsip[b]));
                        //RUN(set_gap_penalties(profile[a],len_a,msa->nsip[b]));
                        //smooth_gaps(profile[a],len_a,window,strength);
                        //increase_gaps(profile[a],len_a,window,strength);
                }
                if (b < numseq){
                        RUN(make_profile_n(ap, msa->sequences[b]->s,len_b,&profile[b]));
                        //RUNP(profile[b] = make_profile(ap,msa->sequences[b]->s,len_b));
                }else{
                        RUN(set_gap_penalties_n(profile[b],len_b,msa->nsip[a]));
                        //smooth_gaps(profile[b],len_b,window,strength);
                        //increase_gaps(profile[b],len_b,window,strength);
                }

                init_alnmem(m, len_a, len_b);
                //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
                if(a < numseq){
                        if(b < numseq){
                                ap->seq1 = msa->sequences[a]->s;
                                ap->seq2 = msa->sequences[b]->s;
                                ap->prof1 = NULL;
                                ap->prof2 = NULL;

                                /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", ap->score); */
                                ap->mode = ALN_MODE_FULL;
                                aln_runner(m, ap, map[c]);
                                //hirsch_ss_dyn(ap,msa->sequences[a]->s, msa->sequences[b]->s,hm,map[c]);
                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                ap->seq1 = NULL;
                                ap->seq2 = msa->sequences[a]->s;
                                ap->prof1 = profile[b];
                                ap->prof2 = NULL;
                                ap->sip = msa->nsip[b];

                                /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", ap->score); */
                                ap->mode = ALN_MODE_FULL;
                                aln_runner(m, ap, map[c]);

                                //hirsch_ps_dyn(ap,profile[b], msa->sequences[a]->s,hm,map[c],msa->nsip[b]);
                                RUN(mirror_path_n(&map[c],len_a,len_b));
                                //RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                        }
                }else{
                        if(b < numseq){
                                ap->seq1 = NULL;
                                ap->seq2 = msa->sequences[b]->s;
                                ap->prof1 = profile[a];
                                ap->prof2 = NULL;
                                ap->sip = msa->nsip[a];
                                /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", ap->score); */
                                ap->mode = ALN_MODE_FULL;
                                aln_runner(m, ap, map[c]);

                                //hirsch_ps_dyn(ap,profile[a],msa->sequences[b]->s ,hm,map[c],msa->nsip[a]);
                        }else{
                                if(len_a < len_b){
                                        ap->seq1 = NULL;
                                        ap->seq2 = NULL;
                                        ap->prof1 = profile[a];
                                        ap->prof2 = profile[b];
                                        /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                        /* aln_runner(m, ap, map[c]); */
                                        /* LOG_MSG("SCORE: %f", ap->score); */
                                        ap->mode = ALN_MODE_FULL;
                                        aln_runner(m, ap, map[c]);
                                        //hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
                                }else{
                                        m->enda = len_b;
                                        m->endb = len_a;
                                        m->len_a = len_b;
                                        m->len_b = len_a;

                                        ap->seq1 = NULL;
                                        ap->seq2 = NULL;
                                        ap->prof1 = profile[b];
                                        ap->prof2 = profile[a];
                                        /* ap->mode = ALN_MODE_SCORE_ONLY; */
                                        /* aln_runner(m, ap, map[c]); */
                                        /* LOG_MSG("SCORE: %f", ap->score); */
                                        ap->mode = ALN_MODE_FULL;
                                        aln_runner(m, ap, map[c]);
                                        //hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
                                        RUN(mirror_path_n(&map[c],len_a,len_b));
                                        //RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                                }
                        }
                }

                RUN(add_gap_info_to_path_n(&map[c], len_a, len_b));
                //map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

                if(i != numseq-2){
                        //MREALLOC(profile_ptr, sizeof(float)*64*(map[c][0]+2));
                        MMALLOC(profile[c],sizeof(float)*64*(map[c][0]+2));
                        //update(profile[a],profile[b],profile[c],map[c]);
                        update_n(profile[a],profile[b],profile[c],ap,map[c],msa->nsip[a],msa->nsip[b]);
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

        free_aln_mem(m);
        return map;
ERROR:
        if(m){
                free_aln_mem(m);
        }
        return NULL;
}




int int_cmp(const void *a, const void *b)
{
    const int *ia = (const int *)a; // casting pointer types
    const int *ib = (const int *)b;
    return *ib  - *ia;

}

int SampleWithoutReplacement(struct rng_state* rng, int N, int n,int* samples)
{

        int t = 0; // total input records dealt with
        int m = 0; // number of items selected so far
        double u;

        while (m < n)
        {
                u = tl_random_double(rng);
                //u = GetUniform(); // call a uniform(0,1) random number generator

                if ( (N - t)*u >= n - m ){
                        t++;
                }else{
                        samples[m] = t;
                        t++;
                        m++;
                }
        }
        return OK;
}


int score_aln(struct aln_mem* m,struct aln_param* ap,float** profile, struct msa* msa, int a,int b,int numseq,float* score)
{
        int i,j,g;
        int len_a;
        int len_b;
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
        ap->mode = ALN_MODE_SCORE_ONLY;

        g = (len_a > len_b)? len_a:len_b;

        RUN(resize_aln_mem(m, g));


        if (a < numseq){
                RUN(make_profile_n(ap, msa->sequences[a]->s,len_a,&profile[a]));
                //RUNP(profile[a] = make_profile(ap,msa->sequences[a]->s,len_a));
        }else{
                RUN(set_gap_penalties_n(profile[a],len_a,msa->nsip[b]));
                //RUN(set_gap_penalties(profile[a],len_a,msa->nsip[b]));
                //smooth_gaps(profile[a],len_a,window,strength);

                //increase_gaps(profile[a],len_a,window,strength);
        }
        if (b < numseq){
                RUN(make_profile_n(ap, msa->sequences[b]->s,len_b,&profile[b]));
                //RUNP(profile[b] = make_profile(ap,msa->sequences[b]->s,len_b));
        }else{
                RUN(set_gap_penalties_n(profile[b],len_b,msa->nsip[a]));
                //smooth_gaps(profile[b],len_b,window,strength);
                //increase_gaps(profile[b],len_b,window,strength);
        }

        init_alnmem(m, len_a, len_b);
        //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
        if(a < numseq){
                if(b < numseq){
                        ap->seq1 = msa->sequences[a]->s;
                        ap->seq2 = msa->sequences[b]->s;
                        ap->prof1 = NULL;
                        ap->prof2 = NULL;

                        aln_runner(m, ap, NULL);
                }else{
                        m->enda = len_b;
                        m->endb = len_a;
                        m->len_a = len_b;
                        m->len_b = len_a;

                        ap->seq1 = NULL;
                        ap->seq2 = msa->sequences[a]->s;
                        ap->prof1 = profile[b];
                        ap->prof2 = NULL;
                        ap->sip = msa->nsip[b];

                        aln_runner(m, ap, NULL);
                        ap->score = ap->score / msa->nsip[b];

                }
        }else{
                if(b < numseq){
                        ap->seq1 = NULL;
                        ap->seq2 = msa->sequences[b]->s;
                        ap->prof1 = profile[a];
                        ap->prof2 = NULL;
                        ap->sip = msa->nsip[a];
                        aln_runner(m, ap, NULL);
                        ap->score = ap->score / msa->nsip[a];
                }else{
                        if(len_a < len_b){
                                ap->seq1 = NULL;
                                ap->seq2 = NULL;
                                ap->prof1 = profile[a];
                                ap->prof2 = profile[b];
                                aln_runner(m, ap, NULL);

                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                ap->seq1 = NULL;
                                ap->seq2 = NULL;
                                ap->prof1 = profile[b];
                                ap->prof2 = profile[a];
                                aln_runner(m, ap, NULL);

                        }
                        ap->score = ap->score / (msa->nsip[a] * msa->nsip[b]);
                }
        }
        *score = ap->score;
        return OK;
ERROR:
        return FAIL;
}
