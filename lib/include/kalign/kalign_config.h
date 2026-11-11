#ifndef KALIGN_CONFIG_H
#define KALIGN_CONFIG_H

#include <stdint.h>
#include <stddef.h>

/* Guide-tree distance (kalign_run_config.guide_dist) */
#define KALIGN_GUIDE_DIST_BPM     0  /* alignment-free BPM edit distance (shipped)       */
#define KALIGN_GUIDE_DIST_PAIRALN 1  /* Kimura/JC distance from pairwise global alignments */

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
        int vsm_mode;                /* 0 = legacy VSM (vsm_amax), 1 = per-merge VSM  */
        float vsm_kappa;             /* per-merge VSM: fraction of offset bound (0-1) */
        float vsm_dzero;             /* per-merge VSM: distance where offset hits 0   */
        int guide_dist;              /* KALIGN_GUIDE_DIST_* (default: BPM)            */
        int gap_model;               /* 0 = Kalign gap rule (default), 1 = member gap model */
        /* Guide-tree options; -1 = default (the shipped tree).  See struct tree_opts (msa_struct.h). */
        int tree_method;
        int tree_dist;
        int tree_alpha;
        int tree_maxlen;
        int tree_leaf;
        int tree_kmeans_dist;
        int realign_tree;
        float tree_sl_noise;         /* hybrid tree: perturbation sigma of the single-
                                        linkage subtrees, keyed by tree_seed (0 = none) */
};

/* Consensus fill-in modes (kalign_ensemble_config.fill_in) */
#define KALIGN_FILL_OFF     0  /* leave unsupported residues in own columns     */
#define KALIGN_FILL_REALIGN 1  /* realign each inter-anchor segment afresh     */
#define KALIGN_FILL_PROJECT 2  /* place them as in the POAR-best ensemble run  */
#define KALIGN_FILL_ATTACH  3  /* let a still-unaligned residue join a column
                                  on a pair with >= fill_support votes         */

/* Ensemble orchestration configuration.
   Only used when n_runs > 1. */
struct kalign_ensemble_config {
        int min_support;             /* POAR consensus threshold (0 = auto)           */
        int consistency_merge;       /* 0 = POAR (default), 1 = MSA consistency merge */
        float consistency_merge_weight; /* MSA consistency bonus weight (default: 2.0)*/
        int fill_in;                 /* fill residues between consensus anchor
                                        columns: KALIGN_FILL_OFF (default),
                                        KALIGN_FILL_REALIGN, KALIGN_FILL_PROJECT      */
        float fill_anchor_frac;      /* fill-in anchor: column needs at least
                                        max(2, ceil(frac * numseq)) residues
                                        (default 0.0 = any merged column)            */
        int fill_support;            /* KALIGN_FILL_ATTACH vote threshold (default 1) */
};

#endif
