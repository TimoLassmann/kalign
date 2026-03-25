#ifndef BISECTINGKMEANS_H
#define BISECTINGKMEANS_H

#include <stdint.h>

#ifdef BISECTINGKMEANS_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct aln_tasks;
struct msa;

/* Callback for computing pairwise distances within a leaf cluster.
   Given a subset of sequence indices (samples[0..n-1]), return an n×n
   symmetric distance matrix.  The caller (UPGMA) frees it via gfree().
   Returns NULL on failure. */
typedef float** (*pair_dist_fn)(struct msa* msa, int* samples, int n);

EXTERN int build_tree_kmeans(struct msa* msa, struct aln_tasks** tasks);
EXTERN int build_tree_kmeans_noisy(struct msa* msa, struct aln_tasks** tasks,
                                   uint64_t seed, float noise_sigma);
EXTERN int build_tree_kmeans_from_dm(struct msa* msa, struct aln_tasks** tasks,
                                     float** dm, int num_anchors,
                                     pair_dist_fn leaf_dist);
EXTERN int build_tree_from_pairwise(struct msa* msa, struct aln_tasks** tasks,
                                    float** dm);

#undef BISECTINGKMEANS_IMPORT
#undef EXTERN

#endif
