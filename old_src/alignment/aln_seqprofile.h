#ifndef ALN_SEQPROF_H
#define ALN_SEQPROF_H

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif

struct aln_mem;

kalign_extern int aln_seqprofile_foward(struct aln_mem* m);
kalign_extern int aln_seqprofile_backward(struct aln_mem* m);
kalign_extern int aln_seqprofile_meetup(struct aln_mem* m,int old_cor[],int* meet,int* t,float* score);


#endif
