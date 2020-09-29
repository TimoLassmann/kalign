#ifndef ALN_TASK_H
#define ALN_TASK_H

struct task{
        int a;                  /* input 1 */
        int b;                  /* input 2 */
        int c;                  /* output  */
        int p;                  /* priority */
};

struct aln_task_list{
        struct task** list;     /* list of pairwise alignments and their priority */
        float** profile;        /* buffer to hold output profiles */
        int** map;              /* traceback paths */
        int n_tasks;
        int n_alloc_tasks;
};

#ifdef ALN_TASK_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int alloc_task_list(struct aln_task_list** tasks,int numseq);
EXTERN void free_task_list(stuct aln_task_list* tasks);

#undef ALN_TASK_IMPORT
#undef EXTERN
#endif
