#ifndef ALN_RUN_H
#define ALN_RUN_H


#ifdef ALN_RUN_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif


struct aln_param;
struct aln_tasks;
struct msa;
/* EXTERN int create_msa(struct msa* msa, struct aln_param* ap,struct aln_tasks* t); */

EXTERN int create_msa_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);
EXTERN int create_chaos_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);
EXTERN int create_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);

/* EXTERN int create_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t); */
EXTERN int create_chaos_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);

#undef ALN_RUN_IMPORT
#undef EXTERN


#endif
