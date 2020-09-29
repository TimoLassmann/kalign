#ifndef ALN_RUN_H
#define ALN_RUN_H

#ifdef ALN_RUN_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN int** create_msa(struct msa* msa, struct aln_param* ap,struct aln_task_list* t);
/* EXTERN int** create_msa(struct msa* msa, struct aln_param* ap); */
EXTERN int** create_chaos_msa(struct msa* msa, struct aln_param* ap);
#undef ALN_RUN_IMPORT
#undef EXTERN
#endif
