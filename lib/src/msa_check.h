#ifndef MSA_CHECK_H
#define MSA_CHECK_H

#ifdef MSA_CHECK_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct msa;
int check_msa(struct msa* msa);


#undef MSA_CHECK_IMPORT
#undef EXTERN


#endif
