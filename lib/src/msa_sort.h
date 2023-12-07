#ifndef MSA_SORT_H
#define MSA_SORT_H

struct msa;
struct rng_state;

#ifdef MSA_SORT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int msa_sort_len_name(struct msa *m);
EXTERN int msa_sort_rank(struct msa *m);
EXTERN int msa_shuffle_seq(struct msa *m, struct rng_state* rng);
#undef MSA_SORT_IMPORT
#undef EXTERN


#endif
