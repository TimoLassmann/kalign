#include "tldevel.h"
#include "tlrng.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "msa.h"
#include "alignment_parameters.h"

#include "aln_task.h"
#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"

#define ALN_RUN_IMPORT
#include "aln_run.h"

static int score_aln(struct aln_mem* m,float** profile, struct msa* msa, int a,int b,int numseq,float* score);
static int SampleWithoutReplacement(struct rng_state* rng, int N, int n,int* samples);
static int int_cmp(const void *a, const void *b);
//static int do_align(struct msa* msa, struct aln_param* ap,struct aln_mem* m, int a,int b, int c);

static int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);

int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id)
{
        int a,b,c;
        int len_a;
        int len_b;
        int j,g;
        int numseq;

        a = t->list[task_id]->a;
        b = t->list[task_id]->b;
        c = t->list[task_id]->c;

        numseq = msa->numseq;

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


        MMALLOC(t->map[c],sizeof(int) * (g+2));

        RUN(resize_aln_mem(m, g));

        /* I should not need to do this */
        for (j = 0; j < (g+2);j++){
                t->map[c][j] = -1;
        }

        if (a < numseq){
                RUN(make_profile_n(m->ap, msa->sequences[a]->s,len_a,&t->profile[a]));
        }else{
                RUN(set_gap_penalties_n(t->profile[a],len_a,msa->nsip[b]));
        }

        if (b < numseq){
                RUN(make_profile_n(m->ap, msa->sequences[b]->s,len_b,&t->profile[b]));
        }else{
                RUN(set_gap_penalties_n(t->profile[b],len_b,msa->nsip[a]));
        }

        init_alnmem(m, len_a, len_b);
        //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
        if(a < numseq){
                if(b < numseq){

                        m->seq1 = msa->sequences[a]->s;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = NULL;
                        m->prof2 = NULL;
                        /* ap->mode = ALN_MODE_SCORE_ONLY; */
                        /* aln_runner(m, ap, map[c]); */
                        /* LOG_MSG("SCORE: %f", ap->score); */
                        m->mode = ALN_MODE_FULL;
#ifdef HAVE_OPENMP
                        /* omp_set_num_threads(4); */
#pragma omp parallel
                        // Only the first thread will spawn other threads
#pragma omp single nowait
                        {
#endif
                                aln_runner(m, t->map[c]);
                                //hirsch_ss_dyn(ap,msa->sequences[a]->s, msa->sequences[b]->s,hm,map[c]);
#ifdef HAVE_OPENMP
                        }
#endif
                }else{
                        m->enda = len_b;
                        m->endb = len_a;
                        m->len_a = len_b;
                        m->len_b = len_a;

                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[a]->s;
                        m->prof1 = t->profile[b];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[b];

                        /* ap->mode = ALN_MODE_SCORE_ONLY; */
                        /* aln_runner(m, ap, map[c]); */
                        /* LOG_MSG("SCORE: %f", ap->score); */
                        m->mode = ALN_MODE_FULL;
#ifdef HAVE_OPENMP
                        /* omp_set_num_threads(4); */
#pragma omp parallel
                        // Only the first thread will spawn other threads
#pragma omp single nowait
                        {
#endif
                                aln_runner(m,t->map[c]);
#ifdef HAVE_OPENMP
                        }
#endif
                        //hirsch_ps_dyn(ap,profile[b], msa->sequences[a]->s,hm,map[c],msa->nsip[b]);
                        RUN(mirror_path_n(&t->map[c],len_a,len_b));
                        //RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                }
        }else{
                if(b < numseq){
                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = t->profile[a];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[a];
                        /* m->mode = ALN_MODE_SCORE_ONLY; */
                        /* aln_runner(m, ap, map[c]); */
                        /* LOG_MSG("SCORE: %f", m->score); */
                        m->mode = ALN_MODE_FULL;
#ifdef HAVE_OPENMP
                        /* omp_set_num_threads(4); */
#pragma omp parallel
                        // Only the first thread will spawn other threads
#pragma omp single nowait
                        {
#endif

                                aln_runner(m,t->map[c]);
#ifdef HAVE_OPENMP
                        }
#endif

                        //hirsch_ps_dyn(ap,profile[a],msa->sequences[b]->s ,hm,map[c],msa->nsip[a]);
                }else{
                        if(len_a < len_b){
                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[a];
                                m->prof2 = t->profile[b];
                                /* m->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", m->score); */
                                m->mode = ALN_MODE_FULL;
                                aln_runner(m, t->map[c]);
                                //hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[b];
                                m->prof2 = t->profile[a];
                                /* m->mode = ALN_MODE_SCORE_ONLY; */
                                /* aln_runner(m, ap, map[c]); */
                                /* LOG_MSG("SCORE: %f", m->score); */
                                m->mode = ALN_MODE_FULL;
#ifdef HAVE_OPENMP
                                /* omp_set_num_threads(4); */
#pragma omp parallel
                                // Only the first thread will spawn other threads
#pragma omp single nowait
                                {
#endif

                                        aln_runner(m, t->map[c]);
#ifdef HAVE_OPENMP
                                }
#endif

                                //hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
                                RUN(mirror_path_n(&t->map[c],len_a,len_b));
                                //RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                        }
                }
        }

        RUN(add_gap_info_to_path_n(&t->map[c], len_a, len_b));
        //map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);


        if(task_id != t->n_tasks-1){
                //if(i != numseq-2){
                //MREALLOC(profile_ptr, sizeof(float)*64*(map[c][0]+2));
                MMALLOC(t->profile[c],sizeof(float)*64*(t->map[c][0]+2));
                //update(profile[a],profile[b],profile[c],map[c]);
                update_n(t->profile[a],t->profile[b],t->profile[c],m->ap,t->map[c],msa->nsip[a],msa->nsip[b]);
        }

        msa->plen[c] = t->map[c][0];

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
        MFREE(t->profile[a]);
        MFREE(t->profile[b]);
        return OK;
ERROR:
        return FAIL;
}

#ifdef HAVE_OPENMP
int create_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{

        int i,j,g,s;


        struct aln_mem** m = NULL;
        int n_threads =  omp_get_max_threads();

        MMALLOC(m, sizeof(struct aln_mem*) * n_threads);
        for(i = 0; i < n_threads;i++){
                m[i] = NULL;
                RUN(alloc_aln_mem(&m[i], 2048));
                m[i]->ap = ap;
                m[i]->mode = ALN_MODE_FULL;
        }
        //LOG_MSG(" Allocated %d threads", n_threads);
        s = 0;

        g = t->list[0]->p;
        for(i = 0; i < t->n_tasks;i++){
                if(t->list[i]->p != g){
#pragma omp parallel for shared(msa,t,m,s,i) private(j)
                        for(j = s; j < i;j++){
                                int tid = omp_get_thread_num();
                                do_align(msa,t,m[tid],j);
                        }
                        //fprintf(stdout,"\n");

                        g =t->list[i]->p;
                        s = i;
                }
        }
        for(j = s; j < i;j++){
                //fprintf(stdout,"%3d %3d -> %3d (p: %d)\n", t->list[j]->a, t->list[j]->b, t->list[j]->c, t->list[j]->p);
                do_align(msa,t,m[0],j);
        }

        for(i = 0; i < n_threads;i++){
                free_aln_mem(m[i]);
        }
        MFREE(m);
        return OK;
ERROR:
        if(m){
                for(i = 0; i < n_threads;i++){
                        free_aln_mem(m[i]);
                }
                MFREE(m);

        }
        return FAIL;
}
#endif

int create_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{

        int i,j,g,s;
        struct aln_mem* m = NULL;
        RUN(alloc_aln_mem(&m, 2048));


        m->ap = ap;
        m->mode = ALN_MODE_FULL;

        s = 0;

        g = t->list[0]->p;
        for(i = 0; i < t->n_tasks;i++){
                if(t->list[i]->p != g){
                        for(j = s; j < i;j++){
                                fprintf(stdout,"%3d %3d -> %3d (p: %d)\n", t->list[j]->a, t->list[j]->b, t->list[j]->c, t->list[j]->p);

                                do_align(msa,t,m,j);
                        }
                        fprintf(stdout,"\n");
                        g =t->list[i]->p;
                        s = i;
                }
        }
        for(j = s; j < i;j++){
                fprintf(stdout,"%3d %3d -> %3d (p: %d)\n", t->list[j]->a, t->list[j]->b, t->list[j]->c, t->list[j]->p);
                do_align(msa,t,m,j);
        }

        free_aln_mem(m);

        return OK;
ERROR:
        if(m){
                free_aln_mem(m);
        }
        return FAIL;
}


int create_chaos_msa(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{
        struct aln_mem* m = NULL;
        int i,g,f,a,b,l;
        int best_a, best_b;

        int* samples = NULL;


        int* active = NULL;
        float max_score;
        float score;
        int numseq;

        m->mode = ALN_MODE_FULL;

        g = msa->num_profiles;
        numseq = msa->numseq;

        RUN(alloc_aln_mem(&m, 2048));
        m->ap = ap;
        MMALLOC(samples,sizeof(int) * m->ap->chaos);
        MMALLOC(active, sizeof(int) * numseq);
        for(i = 0; i < numseq;i++){
                active[i] = i;
        }
        qsort(active, numseq, sizeof(int), int_cmp);

        for(i = 0; i < numseq-1;i++){
                /* pick one sequence / profile  */
                max_score = -FLT_MAX;
                l = MACRO_MIN(ap->chaos, numseq-i);
                SampleWithoutReplacement(ap->rng, numseq-i, l, samples);

                for(g = 0;g < l-1;g++){
                        a = samples[g];
                        for(f = g + 1; f < l;f++){
                                b = samples[f];
                                score_aln(m, t->profile, msa, active[a], active[b], numseq, &score);
                                //LOG_MSG("TEsting %d %d : %f", a,b, ap->score);
                                if(m->score > max_score){
                                        best_a = a;
                                        best_b = b;
                                        max_score = m->score;
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
                //exit(0);
                a = best_a;
                b = best_b;
                //LOG_MSG("samples: %d %d", active[a],active[b]);
                t->list[i]->a = active[a];
                t->list[i]->b = active[b];
                t->list[i]->c = numseq+i;

                /* ap->tree[i*3] = active[a]; */
                /* ap->tree[i*3+1] = active[b]; */
                /* ap->tree[i*3+2] = numseq+i; */

                active[a] = numseq+i;
                active[b] = -1;
                qsort(active, numseq-i, sizeof(int), int_cmp);

                /* a = ap->tree[i*3]; */
                /* b = ap->tree[i*3+1]; */
                /* c = ap->tree[i*3+2]; */
                do_align(msa,t,m,i);
                //score_aln(m, ap, profile, msa, a, b, numseq, &score);
                //fprintf(stdout,"Aligning:%d %d->%d	done:%f score:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100,score);

        }

        MFREE(active);
        MFREE(samples);
        free_aln_mem(m);
        return OK;
ERROR:
        if(m){
                free_aln_mem(m);
        }
        return FAIL;
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


int score_aln(struct aln_mem* m,float** profile, struct msa* msa, int a,int b,int numseq,float* score)
{
        int g;
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
        m->mode = ALN_MODE_SCORE_ONLY;

        g = (len_a > len_b)? len_a:len_b;

        RUN(resize_aln_mem(m, g));


        if (a > numseq){
                RUN(set_gap_penalties_n(profile[a],len_a,msa->nsip[b]));
        }
        if (b > numseq){
                RUN(set_gap_penalties_n(profile[b],len_b,msa->nsip[a]));
        }

        init_alnmem(m, len_a, len_b);
        //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
        if(a < numseq){
                if(b < numseq){
                        m->seq1 = msa->sequences[a]->s;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = NULL;
                        m->prof2 = NULL;

                        aln_runner(m,  NULL);
                }else{
                        m->enda = len_b;
                        m->endb = len_a;
                        m->len_a = len_b;
                        m->len_b = len_a;

                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[a]->s;
                        m->prof1 = profile[b];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[b];

                        aln_runner(m,  NULL);
                        m->score = m->score / msa->nsip[b];

                }
        }else{
                if(b < numseq){
                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = profile[a];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[a];
                        aln_runner(m, NULL);
                        m->score = m->score / msa->nsip[a];
                }else{
                        if(len_a < len_b){
                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = profile[a];
                                m->prof2 = profile[b];
                                aln_runner(m, NULL);

                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = profile[b];
                                m->prof2 = profile[a];
                                aln_runner(m, NULL);

                        }
                        m->score = m->score / (msa->nsip[a] * msa->nsip[b]);
                }
        }
        *score = m->score;
        return OK;
ERROR:
        return FAIL;
}
