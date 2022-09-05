#ifndef ALN_TASK_H
#define ALN_TASK_H

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif


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


#define TASK_ORDER_PRIORITY 1
#define TASK_ORDER_TREE 2

kalign_extern int sort_tasks(struct aln_tasks* t , int order);
kalign_extern int alloc_tasks(struct aln_tasks** tasks,int numseq);
kalign_extern void free_tasks(struct aln_tasks* tasks);


#endif
