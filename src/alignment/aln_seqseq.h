#ifndef ALN_SEQSEQ_H
#define ALN_SEQSEQ_H

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif

kalign_extern int aln_seqseq_foward(struct aln_mem* m);
kalign_extern int aln_seqseq_backward(struct aln_mem* m);
kalign_extern int aln_seqseq_meetup(struct aln_mem* m,int old_cor[],int* meet,int* t,float* score);

#endif
