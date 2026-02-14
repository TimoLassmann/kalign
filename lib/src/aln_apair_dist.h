#ifndef ALN_APAIR_DIST_H
#define ALN_APAIR_DIST_H

#ifdef ALN_APAIR_DIST_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;

/* Compute NxN pairwise identity distances from an aligned MSA.
   Distance = 1.0 - (matches / aligned_positions) for each pair.
   The MSA must be finalized (sequences contain gap characters).
   Caller must free the returned matrix with free_dm(). */
EXTERN int compute_aln_pairwise_dist(struct msa* msa, float*** dm_ptr);

/* Free an NxN distance matrix allocated by compute_aln_pairwise_dist. */
EXTERN void free_aln_dm(float** dm, int n);

#undef ALN_APAIR_DIST_IMPORT
#undef EXTERN

#endif
