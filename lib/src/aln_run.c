#include "tldevel.h"
#include "tlrng.h"

#include "msa_struct.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "task.h"

#include "aln_param.h"

#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"

#include "weave_alignment.h"
/* #include "weave_alignment.h" */

#define ALN_RUN_IMPORT
#include "aln_run.h"

static void recursive_aln(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c);
/* static void recursive_aln_openMP(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c); */
/* static void recursive_aln_serial(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c); */

static int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);
/* static int do_align_serial(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id); */
/* static int do_score(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id); */

/* static int SampleWithoutReplacement(struct rng_state* rng, int N, int n,int* samples); */
/* static int int_cmp(const void *a, const void *b); */

int create_msa_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
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
        /* LOG_MSG("Setting threads to 1 for debugging!"); */
        /* ap->nthreads = 1; */
        msa->run_parallel = 1;
        if(ap->nthreads == 1){
                msa->run_parallel = 0;
        }

#ifdef HAVE_OPENMP
#pragma omp parallel
#pragma omp single nowait
#endif
        recursive_aln(msa, t, ap, active, t->n_tasks-1);

        MFREE(active);
        return OK;
ERROR:
        if(active){
                MFREE(active);
        }
        return FAIL;
}


void recursive_aln(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c)
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

        if(!active[local_t->a] && local_t->a >= msa->numseq){
#ifdef HAVE_OPENMP
#pragma omp task shared(msa,t,ap,active) firstprivate(a)
#endif
                recursive_aln(msa, t, ap, active, a);
        }
        if(!active[local_t->b] && local_t->b >= msa->numseq){
#ifdef HAVE_OPENMP
#pragma omp task shared(msa,t,ap,active) firstprivate(b)
#endif
                recursive_aln(msa, t, ap, active, b);
        }
#ifdef HAVE_OPENMP
#pragma omp taskwait
#endif

        struct aln_mem* ml = NULL;

        alloc_aln_mem(&ml, 256);

        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;
        do_align(msa,t,ml,c);

        active[local_t->a] = 0;
        active[local_t->b] = 0;
        active[local_t->c] = 1;
        /* LOG_MSG("Local: %d %d %d p:%d", local_t->a, local_t->b, local_t->c, local_t->p); */
        free_aln_mem(ml);
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
                        if(m->len_a < m->len_b){
                                m->seq1 = msa->sequences[a]->s;
                                m->seq2 = msa->sequences[b]->s;
                                /* LOG_MSG("%d %d", m->len_a, m->len_b); */
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

                                m->seq1 = msa->sequences[b]->s;
                                m->seq2 = msa->sequences[a]->s;
                                m->prof1 = NULL;
                                m->prof2 = NULL;

                                aln_runner(m);

                                RUN(mirror_path_n(m,len_a,len_b));
                                m->len_a = len_a;
                                m->len_b = len_b;
                        }
                        /* m->seq1 = msa->sequences[a]->s; */
                        /* m->seq2 = msa->sequences[b]->s; */
                        /* m->prof1 = NULL; */
                        /* m->prof2 = NULL; */
                        /* aln_runner(m); */
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
        /* LOG_MSG("Aligned %d and %d (len %d %d) -> path is of length: %d",a,b, m->len_a,m->len_b, 64*(m->path[0]+2)); */

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
