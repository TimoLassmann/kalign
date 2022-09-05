#ifndef ALN_PARAM_H
#define ALN_PARAM_H

#ifdef ALN_PARAM_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct aln_param{
        int nthreads;
        /* actual parameters  */
        float** subm;
        float gpo;
        float gpe;
        float tgpe;
        float score;
};

EXTERN int aln_param_init(struct aln_param **aln_param,int biotype , int type,int n_threads, float gpo, float gpe, float tgpe);

EXTERN void aln_param_free(struct aln_param* ap);

#undef ALN_PARAM_IMPORT
#undef EXTERN


#endif
