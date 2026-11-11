#ifndef SL_TREE_H
#define SL_TREE_H
/* Guide-tree options of a run (single-linkage and hybrid trees, tree distances); see sl_tree.c. */
#include <stdint.h>

struct msa;

/* binary guide-tree node (leaves: id = sequence index; internal: id = -1 until labelled) */
struct node{
        struct node* left;
        struct node* right;
        int id;
};

#define SL_MODE_KMEANS 0
#define SL_MODE_PRIM   1
#define SL_MODE_SLINK  2
#define SL_MODE_HYBRID 3

struct tree_opts;
struct kalign_run_config;
int sl_opts_resolve(struct tree_opts* o, const struct kalign_run_config* cfg);
int sl_tree_mode(const struct msa* msa);
int sl_kmeans_dist_override(const struct msa* msa);
int sl_kmeans_leaf(const struct msa* msa);
int sl_kmeans_leaf_is_sl(const struct msa* msa);
/* single-linkage tree of all sequences / of a subset (leaf ids = idx[k]); pool: threadpool_t* or NULL */
int sl_build(struct msa* msa, int mode, float noise, uint64_t seed, void* pool, struct node** ret);
int sl_build_subset(struct msa* msa, const int* idx, int N, int mode, float noise, uint64_t seed, void* pool,
                    struct node** ret);
float** sl_anchor_matrix(struct msa* msa, const int* anchors, int K);
float** sl_pair_matrix(struct msa* msa, int* samples, int n);
#endif
