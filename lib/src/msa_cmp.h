#ifndef MSA_CMP_H
#define MSA_CMP_H

#ifdef MSA_CMP_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

EXTERN int kalign_msa_compare(struct msa* r, struct msa* t, int cmp_type, float* score);

#undef MSA_CMP_IMPORT
#undef EXTERN


#endif
