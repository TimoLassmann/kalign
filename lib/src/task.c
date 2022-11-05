#include "tldevel.h"

#define ALN_TASK_IMPORT
#include "task.h"

static int sort_tasks_by_priority(const void *a, const void *b);
static int sort_tasks_by_c(const void *a, const void *b);

#ifdef TASKWRITETEST
#include "tlrng.h"
int main(void)
{
        struct rng_state* rng = NULL;
        struct aln_tasks *t = NULL;
        rng = init_rng(0);

        int n_tasks = 54;

        alloc_tasks(&t, n_tasks);

        for(int i = 0; i < n_tasks;i++){
                t->list[i]->score = 0.0;
                t->list[i]->a = tl_random_int(rng, 1000);
                t->list[i]->b = tl_random_int(rng, 1000);
                t->list[i]->c = tl_random_int(rng, 1000);
                t->list[i]->p = tl_random_int(rng, 1000);
                t->list[i]->n = tl_random_int(rng, 1000);
                t->n_tasks++;
        }


        RUN(write_tasks(t, "task_write_test.txt"));


        free_tasks(t);
        t = 0;

        RUN(read_tasks(&t, "task_write_test.txt" ));

        for(int i = 0; i < t->n_tasks;i++){
                struct task* a = t->list[i];
                fprintf(stdout,"%d %d %d %d %d\n",a->a,a->b,a->c,a->p,a->n);
        }


        free_rng(rng);
        return EXIT_SUCCESS;
ERROR:
        if(t){
                free_tasks(t);
        }
        if(rng){
                free_rng(rng);
        }
        return EXIT_FAILURE;
}

#endif

int write_tasks(struct aln_tasks *t, char *filename)
{

        FILE* f_ptr = NULL;

        RUNP(f_ptr = fopen(filename, "w"));

        fprintf(f_ptr,"%d\n", t->n_tasks);

        for(int i = 0; i < t->n_tasks;i++){
                struct task* a = t->list[i];
                fprintf(f_ptr,"%d,%d,%d,%d,%d\n",a->a,a->b,a->c,a->p,a->n);
        }
        fclose(f_ptr);
        return OK;

ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}

int read_tasks(struct aln_tasks **tasks, char *filename)
{
        struct aln_tasks *t = NULL;

        FILE* f_ptr = NULL;

        RUNP(f_ptr = fopen(filename, "r"));
        int n_tasks = 0;



        fscanf( f_ptr, "%d", &n_tasks);

        RUN(alloc_tasks(&t, n_tasks));
        for(int i = 0; i < n_tasks;i++){
                struct task* a = t->list[i];
                fscanf(f_ptr,"%d,%d,%d,%d,%d\n",&a->a,&a->b,&a->c,&a->p,&a->n);
                t->n_tasks++;
        }
        fclose(f_ptr);

        *tasks = t;
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}


int sort_tasks(struct aln_tasks* t , int order)
{
        ASSERT(t != NULL, "No tasks");
        ASSERT(t->n_tasks != 0, "No tasks");

        switch (order) {
        case TASK_ORDER_PRIORITY: {
                qsort(t->list, t->n_tasks, sizeof(struct task*), sort_tasks_by_priority);
                break;
        }
        case TASK_ORDER_TREE: {
                qsort(t->list, t->n_tasks, sizeof(struct task*), sort_tasks_by_c);
                break;
        }
        default:
                ERROR_MSG("Task ordering %d not recognised.");
                break;
        }

        return OK;
ERROR:
        return FAIL;
}


int sort_tasks_by_priority(const void *a, const void *b)
{
        struct task* const *one = a;
        struct task* const *two = b;

        if((*one)->p >= (*two)->p){
                return 1;
        }else{
                return -1;
        }
}

int sort_tasks_by_c(const void *a, const void *b)
{
        struct task* const *one = a;
        struct task* const *two = b;

        if((*one)->c >= (*two)->c){
                return 1;
        }else{
                return -1;
        }
}

int alloc_tasks(struct aln_tasks** tasks,int numseq)
{
        struct aln_tasks* t = NULL;
        int np;
        int i;

        MMALLOC(t, sizeof(struct aln_tasks));

        t->n_tasks = 0;
        t->n_alloc_tasks = numseq;
        t->list = NULL;
        t->profile = NULL;
        np = (numseq << 1) - 1;

        MMALLOC(t->profile,sizeof(float*)*np);
        for(i = 0; i < np;i++){
                t->profile[i] = NULL;
        }

        MMALLOC(t->list, sizeof(struct task*) * t->n_alloc_tasks);
        for(i = 0; i < t->n_alloc_tasks;i++){
                t->list[i] = NULL;
                MMALLOC(t->list[i], sizeof(struct task));
        }

        *tasks = t;
        return OK;
ERROR:
        free_tasks(t);
        return FAIL;
}

void free_tasks(struct aln_tasks* t)
{
        if(t){
                for(int i = 0; i < t->n_alloc_tasks;i++){
                        MFREE(t->list[i]);
                }
                if(t->profile){
                        int np = t->n_alloc_tasks;
                        np =  (np << 1) - 1;

                        for(int i = 0; i < np;i++){
                                if(t->profile[i]){
                                        MFREE(t->profile[i]);
                                }
                        }
                        MFREE(t->profile);
                }
                MFREE(t->list);
                MFREE(t);
        }
}
