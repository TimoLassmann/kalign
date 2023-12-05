#ifndef MSA_SORT_H
#define MSA_SORT_H

struct msa;

#ifdef MSA_SORT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int msa_sort_len_name(struct msa *m);
EXTERN int msa_sort_rank(struct msa *m);

#undef MSA_SORT_IMPORT
#undef EXTERN


#endif
