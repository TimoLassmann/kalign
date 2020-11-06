#ifndef ALN_RUN_H
#define ALN_RUN_H

#ifdef ALN_RUN_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct aln_tasks;

/* EXTERN int create_msa(struct msa* msa, struct aln_param* ap,struct aln_tasks* t); */

/* EXTERN int create_msa_serial_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t); */
EXTERN int create_msa_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t, int n_threads);
EXTERN int create_chaos_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);
EXTERN int create_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);

EXTERN int create_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);
EXTERN int create_chaos_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);

#undef ALN_RUN_IMPORT
#undef EXTERN
#endif
