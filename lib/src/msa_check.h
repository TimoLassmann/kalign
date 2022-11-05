#ifndef MSA_CHECK_H
#define MSA_CHECK_H

#ifdef MSA_CHECK_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;

EXTERN int kalign_essential_input_check(struct msa *msa, int exit_on_error);
EXTERN int kalign_check_msa(struct msa* msa, int exit_on_error);
EXTERN int kalign_sort_msa(struct msa *msa);
#undef MSA_CHECK_IMPORT
#undef EXTERN


#endif
