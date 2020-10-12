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

static int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);
static int do_score(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);

static int SampleWithoutReplacement(struct rng_state* rng, int N, int n,int* samples);
static int int_cmp(const void *a, const void *b);

/* #ifdef HAVE_OPENMP */
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

        RUN(sort_tasks(t, TASK_ORDER_PRIORITY));
/* #pragma omp parallel shared(msa,t,m,s) private(i) */
/* #pragma omp single */

/*         { */
/*                 for(i = 0; i < msa->numseq;i++){ */
/* #pragma omp task firstprivate(i) */
/*                         { */
/* #pragma omp task depend(out: t->profile[i],msa->plen[i]) */
/*                                 t->profile[i] = NULL; */
/*                                 msa->plen[i] = 0; */
/*                         } */
/*                 } */

/*                 for(i = 0; i < t->n_tasks;i++){ */

/* #pragma omp task depend(in:msa->plen[t->list[i]->a],msa->plen[t->list[i]->b],t->profile[t->list[i]->a],t->profile[t->list[i]->b]) depend (out:t->profile[t->list[i]->c], msa->plen[t->list[i]->c]) */
/*                         { */
/*                                 int tid = omp_get_thread_num(); */
/*                                 LOG_MSG("Aligning %d %d -> %d thread %d", t->list[i]->a,t->list[i]->b,t->list[i]->c,tid); */
/*                                 do_align(msa,t,m[tid],i); */
/*                                 LOG_MSG("SONE %d %d -> %d thread %d %d", t->list[i]->a,t->list[i]->b,t->list[i]->c,tid,msa->plen[t->list[i]->c]); */
/*                         } */
/*                 } */

/*         } */

        for(i = 0; i < t->n_tasks;i++){
                if(t->list[i]->p != g){
#pragma omp parallel for shared(msa,t,m,s,i) private(j)
                        for(j = s; j < i;j++){
                                int tid = omp_get_thread_num();
                                do_align(msa,t,m[tid],j);
                        }
                        //fprintf(stdout,"\n");
                        g = t->list[i]->p;
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

int create_chaos_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{
        struct aln_mem** m = NULL;
        struct aln_tasks* t_chaos = NULL;
        int i,j,g,f,a,b,l;
        int best;

        int* samples = NULL;


        int* active = NULL;
        float max_score;

        int numseq;
        int n_threads =  omp_get_max_threads();

        RUN(alloc_tasks(&t_chaos, (ap->chaos * (ap->chaos-1)) / 2));

        MFREE(t_chaos->profile);
        t_chaos->profile = t->profile;


        MMALLOC(m, sizeof(struct aln_mem*) * n_threads);
        for(i = 0; i < n_threads;i++){
                m[i] = NULL;
                RUN(alloc_aln_mem(&m[i], 2048));
                m[i]->ap = ap;
                m[i]->mode = ALN_MODE_FULL;
        }


        g = msa->num_profiles;
        numseq = msa->numseq;

        MMALLOC(samples,sizeof(int) * m[0]->ap->chaos);
        MMALLOC(active, sizeof(int) * numseq);
        for(i = 0; i < numseq;i++){
                active[i] = i;
        }
        qsort(active, numseq, sizeof(int), int_cmp);
        t->n_tasks = 0;
        for(i = 0; i < numseq-1;i++){
                /* pick one sequence / profile  */
                max_score = -FLT_MAX;
                l = MACRO_MIN(ap->chaos, numseq-i);
                SampleWithoutReplacement(ap->rng, numseq-i, l, samples);

                t_chaos->n_tasks = 0;
                /* prepare tasks */
                for(g = 0;g < l-1;g++){
                        for(f = g + 1; f < l;f++){
                                t_chaos->list[t_chaos->n_tasks]->a = active[samples[g]];
                                t_chaos->list[t_chaos->n_tasks]->b = active[samples[f]];
                                /* HACK! -> used below to update the active array */
                                t_chaos->list[t_chaos->n_tasks]->c = samples[g];
                                t_chaos->list[t_chaos->n_tasks]->p = samples[f];
                                t_chaos->list[t_chaos->n_tasks]->score = 0.0F;
                                t_chaos->n_tasks++;
                        }
                }

                /* Run chaos tasks in parallel  */
#pragma omp parallel for shared(msa,t_chaos,m) private(j)
                for(j = 0; j < t_chaos->n_tasks;j++){
                        int tid = omp_get_thread_num();
                        do_score(msa, t_chaos, m[tid], j);
                        //LOG_MSG("%d running %d %d", tid, t_chaos->list[j]->a, t_chaos->list[j]->b);
                }

                max_score = -FLT_MAX;
                best = -1;
                for(j = 0; j < t_chaos->n_tasks;j++){
                        //fprintf(stdout,"%5d\t%5d -> %f\n", t_chaos->list[j]->a,t_chaos->list[j]->b ,t_chaos->list[j]->score);

                        if(t_chaos->list[j]->score > max_score){
                                max_score = t_chaos->list[j]->score;
                                best = j;
                        }
                }


                //LOG_MSG("samples: %d %d", active[a],active[b]);
                t->list[t->n_tasks]->a = t_chaos->list[best]->a;
                t->list[t->n_tasks]->b = t_chaos->list[best]->b;
                t->list[t->n_tasks]->c = numseq+i;

                active[t_chaos->list[best]->c ] = numseq+i;
                active[t_chaos->list[best]->p ] = -1;
                qsort(active, numseq-i, sizeof(int), int_cmp);
                do_align(msa,t,m[0],t->n_tasks);
                t->n_tasks++;
                //score_aln(m, ap, profile, msa, a, b, numseq, &score);
                //fprintf(stdout,"Aligning:%d %d->%d	done:%f score:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100,score);

        }
        for(i = 0; i < n_threads;i++){
                free_aln_mem(m[i]);
        }
        MFREE(m);

        MFREE(active);
        MFREE(samples);
        t_chaos->profile = NULL;
        free_tasks(t_chaos);
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

/* #endif */
static void recursive_aln_openMP(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c);
static void recursive_aln_serial(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c);


int create_msa_serial_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{
        int i,j,c,g,s;
        uint8_t* active = NULL;

        RUN(sort_tasks(t, TASK_ORDER_TREE));
        MMALLOC(active, sizeof(uint8_t)* msa->num_profiles);

        for(i = 0; i < msa->numseq;i++){
                active[i] = 1;
        }
        for(i = msa->numseq; i < msa->num_profiles;i++){
                active[i] = 0;
        }
        /* for(i = 0; i < t->n_tasks;i++){ */
        /*         fprintf(stdout,"%d %3d %3d -> %3d (p: %d)\n",i + msa->numseq, t->list[i]->a, t->list[i]->b, t->list[i]->c, t->list[i]->p); */
        /* } */

#ifdef HAVE_OPENMP
        /* omp_set_nested(1); */
                /* omp_set_num_threads(param->nthreads); */
        omp_set_max_active_levels(4);
        /* omp_set_num_threads(1); */
        /* omp_set_dynamic(1); */
/* #pragma omp parallel num_threads(1) */
/*         { */
/*                 #pragma omp single */
/*                 { */
        /* recursive_aln_serial(msa, t, ap, active, t->n_tasks-1); */
        recursive_aln_openMP(msa, t, ap, active, t->n_tasks-1);
                /* } */

        /* } */
#else
        recursive_aln_serial(msa, t, ap, active, t->n_tasks-1);

#endif

        MFREE(active);

        return OK;
ERROR:
        if(active){
                MFREE(active);
        }
        return FAIL;
}

void recursive_aln_openMP(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c)
{
        struct task* local_t = NULL;

        /* Follow left and right branch until I arrive at sequences / profiles
           ready to align.
        */
        int a;
        int b;
        local_t = t->list[c];

        /* LOG_MSG("Work: %d", local_t->n); */
        /* if(local_t->n < 5){ */
        /*         recursive_aln_serial(msa, t, ap, active, c); */
        /* }else{ */

        a = local_t->a - msa->numseq;
        b = local_t->b - msa->numseq;
#pragma omp parallel num_threads(2)
        {
#pragma omp single nowait
                {

                        if(local_t->a >= msa->numseq){ /* I have an internal node  */
                                if(!active[local_t->a]){
#pragma omp task shared(msa,t,ap,active) firstprivate(a)
                                        {
                                                recursive_aln_openMP(msa, t, ap, active, a);
                                        }
                                }
                                /* I have a lead node - do nothing */
                        }



                        if(local_t->b >= msa->numseq){ /* I have an internal node */
                                if(!active[local_t->b]){

#pragma omp task shared(msa,t,ap,active) firstprivate(b)
                                        {
                                                recursive_aln_openMP(msa, t, ap, active,b);
                                        }

                                }
                                /* I have a lead node - do nothing */
                        }
                }
        }


#pragma omp taskwait

        struct aln_mem* ml = NULL;

        alloc_aln_mem(&ml, 2048);

        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;

        /* if(active[local_t->a] && active[local_t->b]){ */
        /* fprintf(stdout,"THREAD: %d %3d %3d -> %3d (p: %d)\n",tid, t->list[c]->a, t->list[c]->b, t->list[c]->c, t->list[c]->p); */
        do_align(msa,t,ml,c);
        active[local_t->c] = 1;

        free_aln_mem(ml);
        /* } */
}


void recursive_aln_serial(struct msa* msa, struct aln_tasks*t,struct aln_param* ap, uint8_t* active, int c)
{
        struct task* local_t = NULL;
        local_t = t->list[c];
        /* Follow left and right branch until I arrive at sequences / profiles
           ready to align.
        */
        /* LOG_MSG("Work: %d", local_t->n); */
        if(!active[local_t->a]){

                if(local_t->a >= msa->numseq){ /* I have an internal node  */

                        recursive_aln_serial(msa, t, ap, active, local_t->a - msa->numseq);
                }
                /* I have a lead node - do nothing */
        }

        if(!active[local_t->b]){

                if(local_t->b >= msa->numseq){ /* I have an internal node */

                        recursive_aln_serial(msa, t, ap, active, local_t->b - msa->numseq);
                }
                /* I have a lead node - do nothing */
        }


        /* if(active[local_t->a] && active[local_t->b]){ */
        /* fprintf(stdout,"%3d %3d -> %3d (p: %d)\n", t->list[c]->a, t->list[c]->b, t->list[c]->c, t->list[c]->p);
 */
        struct aln_mem* ml = NULL;



        alloc_aln_mem(&ml, 2048);


        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;

        do_align(msa,t,ml,c);
        active[local_t->c] = 1;
        free_aln_mem(ml);
}



int create_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{

        int i,j,g,s;
        struct aln_mem* m = NULL;


        RUN(sort_tasks(t, TASK_ORDER_PRIORITY));

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


int create_chaos_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{
        struct aln_tasks* t_chaos = NULL;
        struct aln_mem* m = NULL;
        int i,g,f,a,b,l;
        int best_a, best_b;

        int* samples = NULL;


        int* active = NULL;
        float max_score;

        int numseq;

        numseq = msa->numseq;


        /* LOG_MSG("Allocating %d",f); */

        RUN(alloc_tasks(&t_chaos,2));

        MFREE(t_chaos->profile);
        t_chaos->profile = t->profile;

        g = msa->num_profiles;


        RUN(alloc_aln_mem(&m, 2048));
        m->ap = ap;
        m->mode = ALN_MODE_FULL;
        MMALLOC(samples,sizeof(int) * m->ap->chaos);
        MMALLOC(active, sizeof(int) * numseq);
        for(i = 0; i < numseq;i++){
                active[i] = i;
        }
        qsort(active, numseq, sizeof(int), int_cmp);
        t->n_tasks = 0;
        for(i = 0; i < numseq-1;i++){
                /* pick one sequence / profile  */
                max_score = -FLT_MAX;
                l = MACRO_MIN(ap->chaos, numseq-i);
                SampleWithoutReplacement(ap->rng, numseq-i, l, samples);

                for(g = 0;g < l-1;g++){
                        a = samples[g];
                        for(f = g + 1; f < l;f++){
                                b = samples[f];
                                t_chaos->list[0]->a = active[a];
                                t_chaos->list[0]->b = active[b];
                                do_score(msa, t_chaos , m, 0);
                                /* score_aln(m, t->profile, msa, active[a], active[b], numseq, &score); */
                                //LOG_MSG("TEsting %d %d : %f", a,b, t_chaos->list[0]->score);
                                if(t_chaos->list[0]->score > max_score){
                                        best_a = a;
                                        best_b = b;
                                        max_score = t_chaos->list[0]->score;
                                }
                        }
                }


                //LOG_MSG("samples: %d %d", active[a],active[b]);
                t->list[t->n_tasks]->a = active[best_a];
                t->list[t->n_tasks]->b = active[best_b];
                t->list[t->n_tasks]->c = numseq+i;
                t->list[t->n_tasks]->p = 0;

                active[best_a] = numseq+i;
                active[best_b] = -1;
                qsort(active, numseq-i, sizeof(int), int_cmp);

                do_align(msa,t,m,t->n_tasks);

                t->n_tasks++;
        }

        MFREE(active);
        MFREE(samples);
        t_chaos->profile = NULL;
        free_tasks(t_chaos);
        free_aln_mem(m);
        return OK;
ERROR:
        if(m){
                free_aln_mem(m);
        }
        if(t_chaos){
                t_chaos->profile = NULL;
                free_tasks(t_chaos);
        }
        return FAIL;
}





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
        /* int tid = omp_get_thread_num(); */
        /* LOG_MSG("Allocing c: %d  tid: %d\n",c,tid); */
        /* LOG_MSG("%d :a", a); */
        /* LOG_MSG("%d :b", b); */
        /* LOG_MSG("%d a", len_a); */
        /* LOG_MSG("%d b", len_b); */
        MMALLOC(t->map[c],sizeof(int) * (g+2));

        RUN(resize_aln_mem(m, g));
        m->mode = ALN_MODE_FULL;
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

/* #ifdef HAVE_OPENMP */
/*                         /\* omp_set_num_threads(4); *\/ */
/* #pragma omp parallel omp_set_num_threads(4) */
/*                         // Only the first thread will spawn other threads */
/* #pragma omp single nowait */
/*                         { */
/* #endif */
                                aln_runner(m, t->map[c]);
                                //hirsch_ss_dyn(ap,msa->sequences[a]->s, msa->sequences[b]->s,hm,map[c]);
/* #ifdef HAVE_OPENMP */
/*                         } */
/* #endif */
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
/* #ifdef HAVE_OPENMP */
/*                         /\* omp_set_num_threads(4); *\/ */
/* #pragma omp parallel */
/*                         // Only the first thread will spawn other threads */
/* #pragma omp single nowait */
/*                         { */
/* #endif */
                                aln_runner(m,t->map[c]);
/* #ifdef HAVE_OPENMP */
/*                         } */
/* #endif */
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
/* #ifdef HAVE_OPENMP */
/*                         /\* omp_set_num_threads(4); *\/ */
/* #pragma omp parallel */
/*                         // Only the first thread will spawn other threads */
/* #pragma omp single nowait */
/*                         { */
/* #endif */

                                aln_runner(m,t->map[c]);
/* #ifdef HAVE_OPENMP */
/*                         } */
/* #endif */

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
/* #ifdef HAVE_OPENMP */
/*                                 /\* omp_set_num_threads(4); *\/ */
/* #pragma omp parallel */
/*                                 // Only the first thread will spawn other threads */
/* #pragma omp single nowait */
/*                                 { */
/* #endif */

                                aln_runner(m, t->map[c]);
/* #ifdef HAVE_OPENMP */
/*                                 } */
/* #endif */

                                //hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
                                RUN(mirror_path_n(&t->map[c],len_a,len_b));
                                //RUNP(map[c] = mirror_hirsch_path(map[c],len_a,len_b));
                        }
                }
        }

        RUN(add_gap_info_to_path_n(&t->map[c], len_a, len_b));
        //map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);



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

        if(task_id != t->n_tasks-1){
                //if(i != numseq-2){
                //MREALLOC(profile_ptr, sizeof(float)*64*(map[c][0]+2));
                MMALLOC(t->profile[c],sizeof(float)*64*(t->map[c][0]+2));
                //update(profile[a],profile[b],profile[c],map[c]);
                update_n(t->profile[a],t->profile[b],t->profile[c],m->ap,t->map[c],msa->nsip[a],msa->nsip[b]);
        }
        MFREE(t->profile[a]);
        MFREE(t->profile[b]);
        return OK;
ERROR:
        return FAIL;
}

int do_score(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id)
/* int score_aln(struct aln_mem* m,float** profile, struct msa* msa, int a,int b,int numseq,float* score) */
{
        int g;
        int len_a;
        int len_b;
        int a,b;
        int numseq;

        numseq = msa->numseq;
        a = t->list[task_id]->a;
        b = t->list[task_id]->b;


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
                RUN(set_gap_penalties_n(t->profile[a],len_a,msa->nsip[b]));
        }
        if (b > numseq){
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

                        aln_runner(m,  NULL);
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

                        aln_runner(m,  NULL);
                        m->score = m->score / (float) msa->nsip[b];

                }
        }else{
                if(b < numseq){
                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = t->profile[a];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[a];
                        aln_runner(m, NULL);
                        m->score = m->score / (float)msa->nsip[a];
                }else{
                        if(len_a < len_b){
                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[a];
                                m->prof2 = t->profile[b];
                                aln_runner(m, NULL);

                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[b];
                                m->prof2 = t->profile[a];
                                aln_runner(m, NULL);

                        }
                        m->score = m->score / (float)(msa->nsip[a] * msa->nsip[b]);
                }
        }
        t->list[task_id]->score = m->score;
        return OK;
ERROR:
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
