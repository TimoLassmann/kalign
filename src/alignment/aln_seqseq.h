#ifndef ALN_SEQSEQ_H
#define ALN_SEQSEQ_H

#ifdef ALN_SEQSEQ_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int aln_seqseq_foward(struct aln_mem* m);
EXTERN int aln_seqseq_backward(struct aln_mem* m);
EXTERN int aln_seqseq_meetup(struct aln_mem* m,int old_cor[],int* meet,int* t,float* score);
#undef ALN_SEQSEQ_IMPORT
#undef EXTERN
#endif
