#ifndef ALN_SETUP_H
#define ALN_SETUP_H

#ifdef ALN_SETUP_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int init_alnmem(struct aln_mem* m);
EXTERN int make_profile_n(struct aln_param* ap,const uint8_t* seq,const int len, float** p);
EXTERN int set_gap_penalties_n(float* prof,int len,int nsip);
EXTERN int add_gap_info_to_path_n(struct aln_mem* m);
EXTERN int update_n(const float* profa, const float* profb,float* newp, struct aln_param*ap, int* path,int sipa,int sipb);


EXTERN int mirror_path_n(struct aln_mem* m,int len_a,int len_b);

#undef ALN_SETUP_IMPORT
#undef EXTERN
#endif
