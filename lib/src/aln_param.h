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

#ifdef USE_THREADPOOL
typedef struct threadpool threadpool_t;
#endif

struct aln_param{
        int nthreads;
#ifdef USE_THREADPOOL
        threadpool_t *pool;   /* shared pool for all parallel work */
#endif
        /* actual parameters  */
        float** subm;
        float gpo;
        float gpe;
        float tgpe;
        float dist_scale;       /* distance-dependent gap scaling: 0=off, >0 scales gap penalties down for divergent pairs */
        float vsm_amax;         /* variable scoring matrix: 0=off, >0 subtracts a(d)=max(0,amax-d) from subm scores */
        float subm_offset;      /* computed per alignment step: amount to subtract from substitution scores */
        int adaptive_budget;    /* 0=off, 1=scale trial count by uncertainty */
        float use_seq_weights;    /* 0=off, >0=pseudocount for profile rebalancing */
        int consistency_anchors;  /* 0=off, >0=number of anchor sequences K for consistency */
        float consistency_weight; /* bonus scale for consistency (default: 2.0) */
};

EXTERN int aln_param_init(struct aln_param **aln_param,int biotype , int n_threads, int type, float gpo, float gpe, float tgpe);

EXTERN void aln_param_free(struct aln_param* ap);

#undef ALN_PARAM_IMPORT
#undef EXTERN


#endif
