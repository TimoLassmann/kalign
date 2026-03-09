#ifndef KALIGN_CONFIG_H
#define KALIGN_CONFIG_H

#include <stdint.h>
#include <stddef.h>

/* Per-run alignment configuration.
   Each field controls one aspect of a single alignment run.
   Use kalign_run_config_defaults() to get a config with all sentinels/defaults,
   then override only the fields you care about. */
struct kalign_run_config {
        int type;                    /* KALIGN_TYPE_* constant (UNDEFINED = auto-detect) */
        float gpo;                   /* gap open penalty (-1.0 = use matrix default) */
        float gpe;                   /* gap extend penalty (-1.0 = use matrix default) */
        float tgpe;                  /* terminal gap extend penalty (-1.0 = use matrix default) */
        float vsm_amax;              /* variable scoring matrix amplitude (-1.0 = biotype default) */
        float dist_scale;            /* distance-dependent gap scaling (0.0 = off) */
        float use_seq_weights;       /* profile rebalancing pseudocount (-1.0 = biotype default) */
        int consistency_anchors;     /* consistency transform: K anchor sequences (0 = off) */
        float consistency_weight;    /* consistency transform: bonus scale (default: 2.0) */
        int refine;                  /* KALIGN_REFINE_* constant (default: NONE) */
        int adaptive_budget;         /* scale refinement trials by uncertainty (0 = off) */
        int realign;                 /* iterative tree-rebuild iterations (0 = off) */
        uint64_t tree_seed;          /* random seed for guide tree perturbation (0 = deterministic) */
        float tree_noise;            /* guide tree perturbation sigma (0.0 = none) */
};

/* Ensemble orchestration configuration.
   Only used when n_runs > 1. */
struct kalign_ensemble_config {
        uint64_t seed;               /* base RNG seed for diversity generation */
        int min_support;             /* POAR consensus threshold (0 = auto) */
        const char* save_poar;       /* path to save POAR table (NULL = don't save) */
};

#endif
