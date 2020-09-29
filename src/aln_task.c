#include "tldevel.h"

#define ALN_TASK_IMPORT
#include "aln_task.h"


int alloc_task_list(struct aln_task_list** tasks,int numseq)
{
        struct aln_task_list* t = NULL;

        t->n_tasks = 0;
        t->n_alloc_tasks = numseq-1;
        t->list = NULL;
        MMALLOC(t->list, sizeof(struct task*) * t->n_alloc_tasks);
        for(i = 0; i < t->n_alloc_tasks;i++){
                t->list[i] = NULL;
                MMALLOC(t->list[i], sizeof(struct task));
        }

        *tasks = t;

        return OK;
ERROR:
        return FAIL;
}

void free_task_list(struct aln_task_list* tasks)
{
        if(tasks){

                MFREE(tasks);
        }

}
