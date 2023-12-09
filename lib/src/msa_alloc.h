#ifndef MSA_ALLOC_H
#define MSA_ALLOC_H

#ifdef MSA_ALLOC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct msa;
struct msa_seq;

EXTERN int alloc_msa(struct msa** msa, int numseq);

EXTERN int resize_msa(struct msa* msa);
EXTERN void kalign_free_msa(struct msa* msa);

EXTERN int alloc_msa_seq(struct msa_seq** s);
EXTERN int resize_msa_seq(struct msa_seq* seq);
EXTERN void free_msa_seq(struct msa_seq* seq);


#undef MSA_ALLOC_IMPORT
#undef EXTERN


#endif
