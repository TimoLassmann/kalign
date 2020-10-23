#ifndef ALN_TASK_H
#define ALN_TASK_H

struct task{
        float score;            /* score of output alignment */
        int a;                  /* input 1 */
        int b;                  /* input 2 */
        int c;                  /* output  */
        int p;                  /* priority */
        int n;                  /* amount of work */
};

struct aln_tasks{
        struct task** list;     /* list of pairwise alignments and their priority */
        float** profile;        /* buffer to hold output profiles */
        /* int** map;              /\* traceback paths *\/ */
        int n_tasks;
        int n_alloc_tasks;
};

#ifdef ALN_TASK_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define TASK_ORDER_PRIORITY 1
#define TASK_ORDER_TREE 2

EXTERN int sort_tasks(struct aln_tasks* t , int order);
EXTERN int alloc_tasks(struct aln_tasks** tasks,int numseq);
EXTERN void free_tasks(struct aln_tasks* tasks);

#undef ALN_TASK_IMPORT
#undef EXTERN
#endif
