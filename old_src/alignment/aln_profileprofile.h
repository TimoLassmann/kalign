#ifndef ALN_PROFILEPROFILE_H
#define ALN_PROFILEPROFILE_H

#ifdef ALN_PROFILEPROFILE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int aln_profileprofile_foward(struct aln_mem* m);
EXTERN int aln_profileprofile_backward(struct aln_mem* m);
EXTERN int aln_profileprofile_meetup(struct aln_mem* m,int old_cor[], int* meet,int* t,float* score);


#undef ALN_PROFILEPROFILE_IMPORT
#undef EXTERN
#endif
