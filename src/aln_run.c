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

#include "weave_alignment.h"

#define ALN_RUN_IMPORT
#include "aln_run.h"

static void recursive_aln_openMP(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c);
static void recursive_aln_serial(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c);

static int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);
static int do_align_serial(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);
static int do_score(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);

static int SampleWithoutReplacement(struct rng_state* rng, int N, int n,int* samples);
static int int_cmp(const void *a, const void *b);

int create_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{

        int i,j,g,s;
        struct aln_mem** m = NULL;

        #ifdef HAVE_OPENMP
        int n_threads =  omp_get_max_threads();
        #else
        int n_threads = 1;
        #endif

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

        for(i = 0; i < t->n_tasks;i++){
                if(t->list[i]->p != g){
                        #ifdef HAVE_OPENMP
#pragma omp parallel for shared(msa,t,m,s,i) private(j)
                        #endif
                        for(j = s; j < i;j++){
                                #ifdef HAVE_OPENMP
                                int tid = omp_get_thread_num();
                                #else
                                int tid = 1;
                                #endif
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
        struct rng_state* rng = NULL;
        int i,j,c,g,f,l;
        int best;
        int task_iter;
        int* active = NULL;
        int work;
        float max_score;
        int old_aln_task = 0;
        int numseq;
#ifdef HAVE_OPENMP
        int n_threads =  omp_get_max_threads();
#else
        int n_threads = 1;
#endif

        RUNP(rng = init_rng(0));
        numseq = msa->numseq;
        g = numseq / ap->chaos *  (ap->chaos * (ap->chaos-1)) / 2   + (ap->chaos*2 * (ap->chaos*2-1)) / 2;

        RUN(alloc_tasks(&t_chaos, g));// (ap->chaos * (ap->chaos-1)) / 2));

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

        MMALLOC(active, sizeof(int) * numseq);
        for(i = 0; i < numseq;i++){
                active[i] = i;
        }
        qsort(active, numseq, sizeof(int), int_cmp);
        t->n_tasks = 0;

        work = msa->numseq;

        c = msa->numseq;
        while(work > 1){

                /* create blocks to be aligned.  */

                t_chaos->n_tasks = 0;
                /* pick one sequence / profile  */
                for(l = 0; l+ap->chaos*2 < work; l+= ap->chaos){
                        for(g = 0; g < l+ ap->chaos-1;g++){
                                for(f = g + 1; f < ap->chaos;f++){
                                        t_chaos->list[t_chaos->n_tasks]->a = active[g+l];
                                        t_chaos->list[t_chaos->n_tasks]->b = active[f+l];
                                        t_chaos->list[t_chaos->n_tasks]->c = g+l;
                                        t_chaos->list[t_chaos->n_tasks]->p = f+l;
                                        t_chaos->list[t_chaos->n_tasks]->score = 0.0F;
                                        t_chaos->n_tasks++;
                                }
                        }

                }
                int size = work - l;

                for(g = 0; g < size-1 ;g++){
                        for(f = g + 1; f < size;f++){
                                t_chaos->list[t_chaos->n_tasks]->a = active[g+l];
                                t_chaos->list[t_chaos->n_tasks]->b = active[f+l];
                                t_chaos->list[t_chaos->n_tasks]->c = g+l;
                                t_chaos->list[t_chaos->n_tasks]->p = f+l;
                                t_chaos->list[t_chaos->n_tasks]->score = 0.0F;
                                t_chaos->n_tasks++;

                        }
                }

#ifdef HAVE_OPENMP
#pragma omp parallel for shared(msa,t_chaos,m) private(j) schedule(dynamic)
#endif
                for(j = 0; j < t_chaos->n_tasks;j++){
#ifdef HAVE_OPENMP
                        int tid = omp_get_thread_num();
#else
                        int tid = 0;
#endif
                        do_score(msa, t_chaos, m[tid], j);
                        /* LOG_MSG("%d", tid); */
                }



                task_iter = 0;
                /* for each block find the best pairwise alignment  */

                for(l = 0; l+ap->chaos*2 < work; l+= ap->chaos){

                        max_score = -FLT_MAX;
                        best = -1;
                        for(g = 0; g < l+ ap->chaos-1;g++){
                                for(f = g + 1; f < ap->chaos;f++){
                                        if(t_chaos->list[task_iter]->score > max_score){
                                                max_score = t_chaos->list[task_iter]->score;
                                                best = task_iter;
                                        }
                                        task_iter++;

                                }
                        }

                        t->list[t->n_tasks]->a = t_chaos->list[best]->a;
                        t->list[t->n_tasks]->b = t_chaos->list[best]->b;
                        t->list[t->n_tasks]->c = c;
                        /* LOG_MSG("Processing: %d %d ->%d",t->list[t->n_tasks]->a,t->list[t->n_tasks]->b,t->list[t->n_tasks]->c); */
                        t->n_tasks++;

                        active[t_chaos->list[best]->c ] = c;
                        active[t_chaos->list[best]->p ] = -1;
                        c++;
                        /* fprintf(stdout,"\n"); */
                }
                size = work - l;
                /* for(l = 0; l < numseq-i; l+= ap->chaos){ */

                max_score = -FLT_MAX;
                best = -1;

                for(g = 0; g < size-1 ;g++){
                        for(f = g + 1; f < size;f++){

                                if(t_chaos->list[task_iter]->score > max_score){
                                        max_score = t_chaos->list[task_iter]->score;
                                        best = task_iter;
                                }
                                task_iter++;
                        }
                }

                t->list[t->n_tasks]->a = t_chaos->list[best]->a;
                t->list[t->n_tasks]->b = t_chaos->list[best]->b;
                t->list[t->n_tasks]->c = c;
                /* LOG_MSG("Processing: %d %d ->%d",t->list[t->n_tasks]->a,t->list[t->n_tasks]->b,t->list[t->n_tasks]->c); */
                t->n_tasks++;

                active[t_chaos->list[best]->c ] = c;
                active[t_chaos->list[best]->p ] = -1;
                c++;

#ifdef HAVE_OPENMP
#pragma omp parallel for shared(msa,t,m) private(j)  schedule(dynamic)
#endif
                for(j = old_aln_task;j < t->n_tasks; j++){
#ifdef HAVE_OPENMP
                        int tid = omp_get_thread_num();
#else
                        int tid = 0;
#endif
                        do_align(msa,t,m[tid],j);
                }

                qsort(active, work, sizeof(int), int_cmp);
                work -= t->n_tasks- old_aln_task;
                old_aln_task = t->n_tasks;
                /* for(i = 0; i < msa->numseq;i++){ */

                /*         fprintf(stdout,"%3d ",active[i]); */
                /* } */
                /* fprintf(stdout,"\n"); */
                /* LOG_MSG("Work: %d", work); */
        }
         /*        exit(0); */

/*                 max_score = -FLT_MAX; */
/*                 l = MACRO_MIN(ap->chaos, numseq-i); */
/*                 SampleWithoutReplacement(rng, numseq-i, l, samples); */


/*                 t_chaos->n_tasks = 0; */
/*                 /\* prepare tasks *\/ */
/*                 for(g = 0;g < l-1;g++){ */
/*                         for(f = g + 1; f < l;f++){ */
/*                                 t_chaos->list[t_chaos->n_tasks]->a = active[samples[g]]; */
/*                                 t_chaos->list[t_chaos->n_tasks]->b = active[samples[f]]; */
/*                                 /\* HACK! -> used below to update the active array *\/ */
/*                                 t_chaos->list[t_chaos->n_tasks]->c = samples[g]; */
/*                                 t_chaos->list[t_chaos->n_tasks]->p = samples[f]; */
/*                                 t_chaos->list[t_chaos->n_tasks]->score = 0.0F; */
/*                                 t_chaos->n_tasks++; */
/*                         } */
/*                 } */
/*                 LOG_MSG("%d tasks", t_chaos->n_tasks); */
/*                 /\* Run chaos tasks in parallel  *\/ */
/* #ifdef HAVE_OPENMP */
/* #pragma omp parallel for shared(msa,t_chaos,m) private(j) */
/* #endif */
/*                 for(j = 0; j < t_chaos->n_tasks;j++){ */
/* #ifdef HAVE_OPENMP */
/*                         int tid = omp_get_thread_num(); */
/* #else */
/*                         int tid = 0; */
/* #endif */
/*                         do_score(msa, t_chaos, m[tid], j); */
/*                         LOG_MSG("%d running %d %d", tid, t_chaos->list[j]->a, t_chaos->list[j]->b); */
/*                 } */

/*                 exit(0); */
/*                 max_score = -FLT_MAX; */
/*                 best = -1; */
/*                 for(j = 0; j < t_chaos->n_tasks;j++){ */
/*                         //fprintf(stdout,"%5d\t%5d -> %f\n", t_chaos->list[j]->a,t_chaos->list[j]->b ,t_chaos->list[j]->score); */

/*                         if(t_chaos->list[j]->score > max_score){ */
/*                                 max_score = t_chaos->list[j]->score; */
/*                                 best = j; */
/*                         } */
/*                 } */


/*                 //LOG_MSG("samples: %d %d", active[a],active[b]); */
/*                 t->list[t->n_tasks]->a = t_chaos->list[best]->a; */
/*                 t->list[t->n_tasks]->b = t_chaos->list[best]->b; */
/*                 t->list[t->n_tasks]->c = numseq+i; */

/*                 active[t_chaos->list[best]->c ] = numseq+i; */
/*                 active[t_chaos->list[best]->p ] = -1; */


/*                 qsort(active, numseq-i, sizeof(int), int_cmp); */
/*                 do_align(msa,t,m[0],t->n_tasks); */
/*                 t->n_tasks++; */
/*                 //score_aln(m, ap, profile, msa, a, b, numseq, &score); */
/*                 //fprintf(stdout,"Aligning:%d %d->%d	done:%f score:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100,score); */

/*         } */
        for(i = 0; i < n_threads;i++){
                free_aln_mem(m[i]);
        }
        MFREE(m);

        MFREE(active);
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


int create_msa_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t, int n_threads)
{
        int i;
        uint8_t* active = NULL;

        RUN(sort_tasks(t, TASK_ORDER_TREE));
        MMALLOC(active, sizeof(uint8_t)* msa->num_profiles);

        for(i = 0; i < msa->numseq;i++){
                active[i] = 1;
        }
        for(i = msa->numseq; i < msa->num_profiles;i++){
                active[i] = 0;
        }
#ifdef HAVE_OPENMP

        if(n_threads == 1){
                recursive_aln_serial(msa, t, ap, active, t->n_tasks-1);
        }else{
                recursive_aln_openMP(msa, t, ap, active, t->n_tasks-1);
        }
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

        a = local_t->a - msa->numseq;
        b = local_t->b - msa->numseq;
#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(2)
        {
#pragma omp single nowait
                {
#endif
                        if(!active[local_t->a]){
#ifdef HAVE_OPENMP
#pragma omp task shared(msa,t,ap,active) firstprivate(a)
#endif
                                {
                                        recursive_aln_openMP(msa, t, ap, active, a);
                                }
                        }
                        if(!active[local_t->b]){
#ifdef HAVE_OPENMP
#pragma omp task shared(msa,t,ap,active) firstprivate(b)
#endif
                                {
                                        recursive_aln_openMP(msa, t, ap, active,b);
                                }

                        }
#ifdef HAVE_OPENMP
                }
        }
#pragma omp taskwait
#endif

        struct aln_mem* ml = NULL;

        alloc_aln_mem(&ml, 256);

        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;

        /* if(active[local_t->a] && active[local_t->b]){ */
        /* fprintf(stdout,"THREAD: %d %3d %3d -> %3d (p: %d)\n",tid, t->list[c]->a, t->list[c]->b, t->list[c]->c, t->list[c]->p); */
        do_align(msa,t,ml,c);
        active[local_t->b] = 0;

        free_aln_mem(ml);
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

        alloc_aln_mem(&ml, 256);

        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;

        do_align_serial(msa,t,ml,c);
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
        struct rng_state* rng = NULL;
        int i,g,f,a,b,l;
        int best_a, best_b;

        int* samples = NULL;


        int* active = NULL;
        float max_score;

        int numseq;
        RUNP(rng = init_rng(0));
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
                SampleWithoutReplacement(rng, numseq-i, l, samples);

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
        float* tmp = NULL;
        int a,b,c;
        int len_a;
        int len_b;
        int j,g;

        a = t->list[task_id]->a;
        b = t->list[task_id]->b;
        c = t->list[task_id]->c;

        if(msa->nsip[a] == 1){
                m->len_a = msa->sequences[a]->len;//  aln->sl[a];
                RUN(make_profile_n(m->ap, msa->sequences[a]->s,m->len_a,&t->profile[a]));
        }else{
                m->len_a = msa->plen[a];
                RUN(set_gap_penalties_n(t->profile[a],m->len_a,msa->nsip[b]));
        }

        if(msa->nsip[b] == 1){
                m->len_b = msa->sequences[b]->len;// aln->sl[b];
                RUN(make_profile_n(m->ap, msa->sequences[b]->s,m->len_b,&t->profile[b]));
        }else{
                m->len_b = msa->plen[b];
                RUN(set_gap_penalties_n(t->profile[b],m->len_b,msa->nsip[a]));
        }

        RUN(init_alnmem(m));

        m->mode = ALN_MODE_FULL;
        if(msa->nsip[a] == 1){
                if(msa->nsip[b] == 1){
                        m->seq1 = msa->sequences[a]->s;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = NULL;
                        m->prof2 = NULL;
                        aln_runner(m);
                }else{
                        len_b = m->len_b;
                        len_a = m->len_a;

                        m->enda = len_b;
                        m->endb = len_a;
                        m->len_a = len_b;
                        m->len_b = len_a;

                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[a]->s;
                        m->prof1 = t->profile[b];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[b];

                        aln_runner(m);
                        RUN(mirror_path_n(m, len_a,len_b));
                        m->len_a = len_a;
                        m->len_b = len_b;
                }
        }else{
                if(msa->nsip[b] == 1){
                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = t->profile[a];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[a];
                        aln_runner(m);
                }else{
                        if(m->len_a < m->len_b){
                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[a];
                                m->prof2 = t->profile[b];
                                aln_runner(m);
                        }else{
                                len_b = m->len_b;
                                len_a = m->len_a;

                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[b];
                                m->prof2 = t->profile[a];

                                aln_runner(m);

                                RUN(mirror_path_n(m,len_a,len_b));
                                m->len_a = len_a;
                                m->len_b = len_b;
                        }
                }
        }

        RUN(add_gap_info_to_path_n(m)) ;

        MMALLOC(tmp,sizeof(float)*64*(m->path[0]+2));

        if(task_id != t->n_tasks-1){
                update_n(t->profile[a],t->profile[b],tmp,m->ap,m->path,msa->nsip[a],msa->nsip[b]);
        }

        MFREE(t->profile[a]);
        MFREE(t->profile[b]);

        t->profile[c] = tmp;
        RUN(make_seq(msa,a,b,m->path));

        msa->plen[c] = m->path[0];

        msa->nsip[c] = msa->nsip[a] + msa->nsip[b];

        MREALLOC(msa->sip[c],sizeof(int)*(msa->nsip[a] + msa->nsip[b]));

        g = 0;
        for (j = msa->nsip[a];j--;){
                msa->sip[c][g] = msa->sip[a][j];
                g++;
        }
        for (j = msa->nsip[b];j--;){
                msa->sip[c][g] = msa->sip[b][j];
                g++;
        }

        return OK;
ERROR:
        return FAIL;
}

int do_align_serial(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id)
{
        float* tmp = NULL;
        int a,b,c;
        int len_a;
        int len_b;
        int j,g;

        a = t->list[task_id]->a;
        b = t->list[task_id]->b;
        c = t->list[task_id]->c;

        if(msa->nsip[a] == 1){
                m->len_a = msa->sequences[a]->len;//  aln->sl[a];
                RUN(make_profile_n(m->ap, msa->sequences[a]->s,m->len_a,&t->profile[a]));
        }else{
                m->len_a = msa->plen[a];
                RUN(set_gap_penalties_n(t->profile[a],m->len_a,msa->nsip[b]));
        }

        if(msa->nsip[b] == 1){
                m->len_b = msa->sequences[b]->len;// aln->sl[b];
                RUN(make_profile_n(m->ap, msa->sequences[b]->s,m->len_b,&t->profile[b]));
        }else{
                m->len_b = msa->plen[b];
                RUN(set_gap_penalties_n(t->profile[b],m->len_b,msa->nsip[a]));
        }

        RUN(init_alnmem(m));

        m->mode = ALN_MODE_FULL;
        if(msa->nsip[a] == 1){
                if(msa->nsip[b] == 1){
                        m->seq1 = msa->sequences[a]->s;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = NULL;
                        m->prof2 = NULL;
                        aln_runner_serial(m);
                }else{
                        len_b = m->len_b;
                        len_a = m->len_a;

                        m->enda = len_b;
                        m->endb = len_a;
                        m->len_a = len_b;
                        m->len_b = len_a;

                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[a]->s;
                        m->prof1 = t->profile[b];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[b];

                        aln_runner_serial(m);
                        RUN(mirror_path_n(m, len_a,len_b));
                        m->len_a = len_a;
                        m->len_b = len_b;
                }
        }else{
                if(msa->nsip[b] == 1){
                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = t->profile[a];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[a];
                        aln_runner_serial(m);
                }else{
                        if(m->len_a < m->len_b){
                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[a];
                                m->prof2 = t->profile[b];
                                aln_runner_serial(m);
                        }else{
                                len_b = m->len_b;
                                len_a = m->len_a;

                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[b];
                                m->prof2 = t->profile[a];

                                aln_runner_serial(m);

                                RUN(mirror_path_n(m,len_a,len_b));
                                m->len_a = len_a;
                                m->len_b = len_b;
                        }
                }
        }

        RUN(add_gap_info_to_path_n(m)) ;

        MMALLOC(tmp,sizeof(float)*64*(m->path[0]+2));

        /* LOG_MSG("%d TASK ID", task_id); */
        if(task_id != t->n_tasks-1){
                update_n(t->profile[a],t->profile[b],tmp,m->ap,m->path,msa->nsip[a],msa->nsip[b]);
        }

        MFREE(t->profile[a]);
        MFREE(t->profile[b]);

        t->profile[c] = tmp;

        RUN(make_seq(msa,a,b,m->path));

        msa->plen[c] = m->path[0];

        msa->nsip[c] = msa->nsip[a] + msa->nsip[b];

        MREALLOC(msa->sip[c],sizeof(int)*(msa->nsip[a] + msa->nsip[b]));

        g = 0;
        for (j = msa->nsip[a];j--;){
                msa->sip[c][g] = msa->sip[a][j];
                g++;
        }
        for (j = msa->nsip[b];j--;){
                msa->sip[c][g] = msa->sip[b][j];
                g++;
        }
        return OK;
ERROR:
        return FAIL;
}


int do_score(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id)
/* int score_aln(struct aln_mem* m,float** profile, struct msa* msa, int a,int b,int numseq,float* score) */
{
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

        m->len_a = len_a;
        m->len_b = len_b;
         /* = (len_a > len_b)? len_a:len_b; */

        if (a > numseq){
                RUN(set_gap_penalties_n(t->profile[a],len_a,msa->nsip[b]));
        }
        if (b > numseq){
                RUN(set_gap_penalties_n(t->profile[b],len_b,msa->nsip[a]));
        }
        RUN(resize_aln_mem(m));
        init_alnmem(m);
        //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
        if(a < numseq){
                if(b < numseq){
                        m->seq1 = msa->sequences[a]->s;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = NULL;
                        m->prof2 = NULL;

                        aln_runner(m);
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

                        aln_runner(m);
                        m->score = m->score / (float) msa->nsip[b];

                }
        }else{
                if(b < numseq){
                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = t->profile[a];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[a];
                        aln_runner(m);
                        m->score = m->score / (float)msa->nsip[a];
                }else{
                        if(len_a < len_b){
                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[a];
                                m->prof2 = t->profile[b];
                                aln_runner(m);

                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[b];
                                m->prof2 = t->profile[a];
                                aln_runner(m);

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
