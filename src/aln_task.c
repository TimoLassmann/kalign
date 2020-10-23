#include "tldevel.h"

#define ALN_TASK_IMPORT
#include "aln_task.h"

static int sort_tasks_by_priority(const void *a, const void *b);
static int sort_tasks_by_c(const void *a, const void *b);



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
        /* t->map = NULL; */

        np =  (numseq << 1) - 1;
        MMALLOC(t->profile,sizeof(float*)*np);
        /* MMALLOC(t->map,sizeof(int*)*np); */

        for(i = 0; i < np;i++){
                t->profile[i] = NULL;
                /* t->map[i] = NULL; */
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
                int i;
                int np;
                for(i = 0; i < t->n_alloc_tasks;i++){
                        MFREE(t->list[i]);
                }
                np = t->n_alloc_tasks;
                np =  (np << 1) - 1;
                /* for(i = 0; i < np;i++){ */
                /*         if(t->map[i]){ */
                /*                 MFREE(t->map[i]); */
                /*         } */
                /* } */


                MFREE(t->list);
                /* MFREE(t->map); */
                if(t->profile){
                        MFREE(t->profile);
                }
                MFREE(t);
        }

}
