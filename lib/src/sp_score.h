#ifndef SP_SCORE_H
#define SP_SCORE_H

#ifdef SP_SCORE_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;
struct aln_param;

EXTERN int compute_sp_score(struct msa* msa, struct aln_param* ap,
                            int* path, int* sip_a, int nsip_a,
                            int* sip_b, int nsip_b, float* score);

#undef SP_SCORE_IMPORT
#undef EXTERN

#endif
