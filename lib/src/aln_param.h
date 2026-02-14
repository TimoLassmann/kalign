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

/* #define KALIGN_DNA 0 */
/* #define KALIGN_DNA_INTERNAL 1 */
/* #define KALIGN_RNA 2 */
/* #define KALIGN_PROTEIN 3 */

struct aln_param{
        int nthreads;
        /* actual parameters  */
        float** subm;
        float gpo;
        float gpe;
        float tgpe;
        float score;
        float dist_scale;       /* distance-dependent gap scaling: 0=off, >0 scales gap penalties down for divergent pairs */
        float vsm_amax;         /* variable scoring matrix: 0=off, >0 subtracts a(d)=max(0,amax-d) from subm scores */
        float subm_offset;      /* computed per alignment step: amount to subtract from substitution scores */
        int adaptive_budget;    /* 0=off, 1=scale trial count by uncertainty */
};

EXTERN int aln_param_init(struct aln_param **aln_param,int biotype , int n_threads, int type, float gpo, float gpe, float tgpe);

EXTERN void aln_param_free(struct aln_param* ap);

#undef ALN_PARAM_IMPORT
#undef EXTERN


#endif
