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
/* int build_tree_kmeans(struct msa* msa, struct aln_param* ap,struct  aln_tasks** task_list); */

EXTERN int build_tree_kmeans(struct msa* msa, struct aln_tasks** tasks);
EXTERN int build_tree_kmeans_noisy(struct msa* msa, struct aln_tasks** tasks,
                                   uint64_t seed, float noise_sigma);
EXTERN int build_tree_from_pairwise(struct msa* msa, struct aln_tasks** tasks,
                                    float** dm);

#undef BISECTINGKMEANS_IMPORT
#undef EXTERN

#endif
