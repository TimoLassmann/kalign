#ifndef ALN_SEQPROF_H
#define ALN_SEQPROF_H

#ifdef ALN_SEQPROFILE_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif


struct aln_mem;

EXTERN int aln_seqprofile_foward(struct aln_mem* m);
EXTERN int aln_seqprofile_backward(struct aln_mem* m);
EXTERN int aln_seqprofile_meetup(struct aln_mem* m,int old_cor[],int* meet,int* t,float* score);



#undef ALN_SEQPROFILE_IMPORT
#undef EXTERN

#endif
