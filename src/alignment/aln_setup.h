#ifndef ALN_SETUP_H
#define ALN_SETUP_H

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif

kalign_extern int init_alnmem(struct aln_mem* m);
kalign_extern int make_profile_n(struct aln_param* ap,const uint8_t* seq,const int len, float** p);
kalign_extern int set_gap_penalties_n(float* prof,int len,int nsip);
kalign_extern int add_gap_info_to_path_n(struct aln_mem* m);
kalign_extern int update_n(const float* profa, const float* profb,float* newp, struct aln_param*ap, int* path,int sipa,int sipb);


kalign_extern int mirror_path_n(struct aln_mem *m, int len_a, int len_b);

#endif
