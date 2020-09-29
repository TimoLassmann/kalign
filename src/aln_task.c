#include "tldevel.h"

#define ALN_TASK_IMPORT
#include "aln_task.h"


int alloc_tasks(struct aln_tasks** tasks,int numseq)
{
        struct aln_tasks* t = NULL;
        int np;
        int i;

        MMALLOC(t, sizeof(struct aln_tasks));

        t->n_tasks = 0;
        t->n_alloc_tasks = numseq-1;
        t->list = NULL;
        t->profile = NULL;
        t->map = NULL;

        np =  (numseq << 1) - 1;
        MMALLOC(t->profile,sizeof(float*)*np);
        MMALLOC(t->map,sizeof(int*)*np);
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
                for(i = 0; i < t->n_alloc_tasks;i++){

                        MFREE(t->list[i]);
                }
                MFREE(t->list);
                MFREE(t->map);
                MFREE(t->profile);
                MFREE(t);
        }

}
