#ifndef ALN_RUN_H
#define ALN_RUN_H

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif

struct aln_param;
struct aln_tasks;
struct msa;
/* EXTERN int create_msa(struct msa* msa, struct aln_param* ap,struct aln_tasks* t); */

/* kalign_extern int create_msa_serial_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t); */
kalign_extern int create_msa_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t, int n_threads);
kalign_extern int create_chaos_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);
kalign_extern int create_msa_serial(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);

kalign_extern int create_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);
kalign_extern int create_chaos_msa_openMP(struct msa* msa, struct aln_param* ap,struct aln_tasks* t);

#endif
