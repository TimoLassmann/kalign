#ifndef KALIGN_CONFIG_H
#define KALIGN_CONFIG_H

#include <stdint.h>
#include <stddef.h>

/* Per-run alignment configuration.
   Each field controls one aspect of a single alignment run.
   All values are concrete — no sentinel values.
   Use kalign_run_config_defaults() for sensible PFASUM43 protein defaults. */
struct kalign_run_config {
        int matrix;                  /* KALIGN_MATRIX_* constant (AUTO = auto-detect) */
        float gpo;                   /* gap open penalty                              */
        float gpe;                   /* gap extend penalty                            */
        float tgpe;                  /* terminal gap extend penalty                   */
        float vsm_amax;              /* variable scoring matrix amplitude (0 = off)   */
        float seq_weights;           /* profile rebalancing pseudo-count (0 = off)    */
        float dist_scale;            /* distance-dependent gap scaling (0 = off)      */
        int refine;                  /* KALIGN_REFINE_* constant (default: NONE)      */
        int adaptive_budget;         /* scale refinement trials by uncertainty (0=off)*/
        int realign;                 /* iterative tree-rebuild iterations (0 = off)   */
        uint64_t tree_seed;          /* random seed for guide tree perturbation       */
        float tree_noise;            /* guide tree perturbation sigma (0.0 = none)    */
        int consistency_anchors;     /* anchor consistency rounds (0 = off)           */
        float consistency_weight;    /* anchor consistency bonus weight (default: 2.0)*/
};

/* Ensemble orchestration configuration.
   Only used when n_runs > 1. */
struct kalign_ensemble_config {
        int min_support;             /* POAR consensus threshold (0 = auto)           */
};

#endif
