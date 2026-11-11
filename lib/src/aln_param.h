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

typedef struct threadpool threadpool_t;

struct aln_param{
        int nthreads;
        threadpool_t *pool;   /* shared pool for all parallel work */
        /* actual parameters  */
        float** subm;
        float gpo;
        float gpe;
        float tgpe;
        float dist_scale;       /* distance-dependent gap scaling: 0=off, >0 scales gap penalties down for divergent pairs */
        float vsm_amax;         /* variable scoring matrix: 0=off, >0 subtracts a(d)=max(0,amax-d) from subm scores */
        float subm_offset;      /* computed per alignment step: amount to subtract from substitution scores */
        int vsm_mode;           /* KALIGN_VSM_LEGACY (0) or KALIGN_VSM_MERGE (1), see aln_vsm.h */
        float vsm_kappa;        /* MERGE: fraction (0..1) of the offset bound used for identical groups */
        float vsm_dzero;        /* MERGE: corrected distance (subst/site) at which the offset reaches 0 */
        int adaptive_budget;    /* 0=off, 1=scale trial count by uncertainty */
        float use_seq_weights;    /* 0=off, >0=pseudocount for profile rebalancing */
        int consistency_anchors;  /* 0=off, >0=number of anchor sequences K for consistency */
        float consistency_weight; /* bonus scale for consistency (default: 2.0) */
        int gap_model;            /* 0 = Kalign gap rule, 1 = member gap model (aln_gapmodel.c) */
};

EXTERN int aln_param_init(struct aln_param **aln_param,int biotype , int n_threads, int type, float gpo, float gpe, float tgpe);

EXTERN void aln_param_free(struct aln_param* ap);

#undef ALN_PARAM_IMPORT
#undef EXTERN


#endif
