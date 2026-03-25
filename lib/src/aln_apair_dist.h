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

#define REALIGN_NUM_ANCHORS 32

struct msa;

/* Compute NxN pairwise identity distances from an aligned MSA.
   Distance = 1.0 - (matches / aligned_positions) for each pair.
   The MSA must be finalized (sequences contain gap characters).
   Caller must free the returned matrix with free_aln_dm(). */
EXTERN int compute_aln_pairwise_dist(struct msa* msa, float*** dm_ptr);

/* Free an NxN distance matrix allocated by compute_aln_pairwise_dist. */
EXTERN void free_aln_dm(float** dm, int n);

/* Identity distance between two aligned sequences (0=identical, 1=no matches).
   Only counts columns where both have a residue (no gap). */
EXTERN float pairwise_identity_dist(const char* a, const char* b, int alnlen);

/* Pairwise identity distance callback for bisecting k-means leaf clusters.
   Returns an n×n distance matrix (allocated with galloc, freed by gfree).
   Sequences must be finalized (seq->seq has gap characters, msa->alnlen set). */
EXTERN float** aln_identity_pair_dist(struct msa* msa, int* samples, int n);

/* Select K diverse anchors from a finalized alignment using farthest-first
   traversal with identity distance.  Returns anchor indices and actual K. */
EXTERN int pick_anchor_from_alignment(struct msa* msa, int K,
                                       int** anchors_out, int* K_out);

/* Compute N×K identity distances from each sequence to K anchors.
   Output format matches d_estimation(pair=0): rows AVX2-padded.
   The MSA must be finalized.  Caller frees with free_aln_anchor_dm(). */
EXTERN int compute_aln_anchor_dist(struct msa* msa, int* anchors, int K,
                                    float*** dm_out);

#undef ALN_APAIR_DIST_IMPORT
#undef EXTERN

#endif
