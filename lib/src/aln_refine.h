#ifndef ALN_REFINE_H
#define ALN_REFINE_H

#ifdef ALN_REFINE_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

#define KALIGN_REFINE_NONE 0
#define KALIGN_REFINE_ALL 1
#define KALIGN_REFINE_CONFIDENT 2
#define KALIGN_REFINE_INLINE 3

struct msa;
struct aln_param;
struct aln_tasks;

struct aln_mem;
EXTERN int dispatch_alignment(struct msa* msa, struct aln_mem* ml, float* prof_a, float* prof_b, int a, int b, int len_a, int len_b);
EXTERN int refine_alignment(struct msa* msa, struct aln_param* ap, struct aln_tasks* t, int refine_mode);

#undef ALN_REFINE_IMPORT
#undef EXTERN

#endif
