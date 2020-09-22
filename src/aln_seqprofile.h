#ifndef ALN_SEQPROF_H
#define ALN_SEQPROF_H

#ifdef ALN_SEQPROF_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int aln_seqprofile_foward(struct aln_mem* m,const struct aln_param* ap);
EXTERN int aln_seqprofile_backward(struct aln_mem* m,const struct aln_param* ap);
EXTERN int aln_seqprofile_meetup(struct aln_mem* m,struct aln_param* ap,int old_cor[],int* meet,int* t,float* score);
#undef ALN_SEQPROF_IMPORT
#undef EXTERN
#endif
