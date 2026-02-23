#ifndef PROBMSA_H
#define PROBMSA_H

#ifdef PROBMSA_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;
struct aln_tasks;

struct probmsa_params {
        double delta;           /* HMM gap-open probability (default 0.01) */
        double epsilon;         /* HMM gap-extend probability (default 0.2) */
        double tau;             /* terminal gap-open probability; <=0 = use delta (global) */
        double mea_threshold;   /* posterior threshold for MEA (default 0.0) */
        double mea_gpo;         /* MEA gap-open cost (default 0.0) */
        double mea_gpe;         /* MEA gap-extend cost (default 0.0) */
        double conc_scale;      /* concentration scale for position-specific gap penalties
                                   delta_eff = delta * cs/(conc + cs). 0 = global delta. */
        double prior_scale;     /* Dirichlet prior strength (default 0.1).
                                   Multiplies Sjölander alpha pseudocounts.
                                   1.0 = full prior, 0.1 = gentle regularization. */
        /* 5-state pair-HMM parameters */
        int use_5state;         /* 0 = 3-state (default), 1 = 5-state */
        double delta_s;         /* short gap-open probability (5-state, default 0.10) */
        double epsilon_s;       /* short gap-extend probability (5-state, default 0.50) */
        double delta_l;         /* long gap-open probability (5-state, default 0.01) */
        double epsilon_l;       /* long gap-extend probability (5-state, default 0.95) */
};

/* Progressive pair-HMM alignment using forward-backward MEA.
   Replaces create_msa_tree() for probabilistic alignment.
   biotype: ALN_BIOTYPE_DNA or ALN_BIOTYPE_PROTEIN.
   subm: 23x23 substitution matrix from aln_param.
   alpha_size: actual alphabet size (5 for DNA, 23 for protein).
   params: if NULL, uses defaults. */
EXTERN int create_msa_tree_probmsa(struct msa* msa, struct aln_tasks* t,
                                    int biotype, float** subm, int alpha_size,
                                    const struct probmsa_params* params);

#undef PROBMSA_IMPORT
#undef EXTERN

#endif
