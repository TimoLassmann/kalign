#ifndef MSA_STRUCT_H
#define MSA_STRUCT_H

#include <stdint.h>

typedef struct threadpool threadpool_t;

#ifdef MSA_STRUCT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define MSA_NAME_LEN 256

#define ALN_STATUS_UNALIGNED 1   /* no gaps sequences may or may not have equal lengths  */
#define ALN_STATUS_ALIGNED 2   /* sequences have equal lengths and may or may not contain gaps*/
#define ALN_STATUS_FINAL 3   /* sequences have equal lengths and may or may not contain gaps*/
#define ALN_STATUS_UNKNOWN 3     /* sequences have un-equal length and contain gaps  */

#define ALN_BIOTYPE_PROTEIN 0
#define ALN_BIOTYPE_DNA 1
#define ALN_BIOTYPE_UNDEF 2

struct msa_seq{
        char* name;
        char* seq;
        uint8_t* s;
        int* gaps;
        float* confidence;   /* per-position confidence [0..1], NULL if not computed */
        int rank;
        int len;
        int alloc_len;
};

/* Guide-tree options of one run; set by sl_opts_resolve() from the run config (defaults: the shipped tree). */
struct tree_opts{
        int mode;          /* SL_MODE_*: 0 k-means (shipped), 1 single linkage (Prim), 2 SLINK, 3 hybrid */
        int dist;          /* 0 raw Myers (shipped), 1 normalised Myers, 2 LCS indel ratio, 3 1 - LCS/shorter */
        int full;          /* 0 reduced tree alphabet (shipped), 1 full alphabet */
        int maxlen;        /* residues per sequence for the tree distance (0 = whole sequences, LCS only) */
        int leaf;          /* groups below this size get an exact subtree (shipped: 50) */
        int kmeans_dist;   /* 0 shipped anchor distances, 1 k-means on the tree distance above */
        int realign_tree;  /* 0 realign pass keeps the shipped leaf rule, 1 it follows leaf / hybrid */
        int in_realign;    /* set by kalign_single_run around the realign iterations */
        float sl_noise;    /* hybrid: perturbation sigma of the single-linkage subtrees (0 = none) */
        uint64_t sl_seed;  /* seed of that perturbation (the run's tree_seed) */
};

struct msa{
        struct msa_seq** sequences;
        float* seq_distances;   /* per-sequence mean distance (normalized), set during tree building */
        float* col_confidence;  /* per-column confidence [0..1], NULL if not computed */
        float* seq_weights;     /* per-sequence tree-based weight, NULL when not computed */
        int** sip;
        int* nsip;
        int* plen;
        uint8_t run_parallel;
        threadpool_t *pool;
        int numseq;
        int num_profiles;
        int alloc_numseq;
        int aligned;
        int alnlen;
        int defer_finalise; /* 1: a single run leaves the alignment as gap counts; kalign_write_msa writes FASTA rows from them */
        int letter_freq[128];
        uint8_t L;
        uint8_t biotype;
        int quiet;
        void* consistency_table;    /* struct consistency_table*, NULL when disabled */
        void* poar_consistency;     /* struct poar_consistency_ctx*, NULL when disabled */
        struct tree_opts tree;
};

#undef MSA_STRUCT_IMPORT
#undef EXTERN

#endif
