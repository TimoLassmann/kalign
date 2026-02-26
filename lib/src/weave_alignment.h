#ifndef WEAVE_ALIGNMENT_H
#define WEAVE_ALIGNMENT_H

#ifdef WEAVE_ALIGNMENT_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;

EXTERN int make_seq(struct msa* msa,int a,int b,int* path);
EXTERN int clean_aln(struct msa* msa);

#undef WEAVE_ALIGNMENT_IMPORT
#undef EXTERN

#endif
