#ifndef GUIDE_DIST_H
#define GUIDE_DIST_H

#include <stdint.h>

#ifdef GUIDE_DIST_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

/* Guide-tree distance from independent pairwise global alignments.

   Each pair of sequences is aligned with Kalign's seq-seq aligner using the
   run's scoring matrix and gap penalties (no VSM offset, no consistency
   bonus).  The distance is the evolutionary correction of the mismatch
   fraction p over aligned residue pairs:
     protein : Kimura (1983)   d = -ln(1 - p - 0.2 p^2)
     DNA/RNA : Jukes-Cantor    d = -3/4 ln(1 - 4/3 p)
   capped at GUIDE_DIST_MAX when the argument of the log falls below 0.01
   (or no residue pair is aligned).  No length term is added.

   The context holds full-alphabet encodings of the sequences, so it can be
   used while msa->sequences[]->s is still in the reduced tree-building
   alphabet. */

struct msa;
struct aln_param;

struct guide_dist_ctx {
        struct aln_param* ap;   /* borrowed, not freed */
        uint8_t** s;            /* full-alphabet encoding per sequence */
        int* len;
        int numseq;
        int biotype;
        /* Reuse of the anchor matrix (filled by guide_dist_anchor_matrix):
           unperturbed distances d(r, anchors[k]) and, per sequence, its
           anchor column (-1 if not an anchor).  A later request for the same
           ordered alignment (same seq1/seq2 after the shorter-first rule)
           returns the stored value instead of realigning. */
        float* anchor_d;        /* numseq x num_anchors, row-major; NULL if unset */
        int* anchor_col;        /* numseq */
        int num_anchors;
};

#define GUIDE_DIST_MAX 4.60517019f   /* -ln(0.01) */

EXTERN int guide_dist_ctx_create(struct msa* msa, struct aln_param* ap,
                                 struct guide_dist_ctx** out);
EXTERN void guide_dist_ctx_free(struct guide_dist_ctx* c);

/* Distance between sequences i and j (indices into msa->sequences). */
EXTERN int guide_dist_pair(const struct guide_dist_ctx* c, int i, int j, float* d);

/* n x n matrix for samples[0..n-1], allocated with galloc (free with gfree).
   Matches pair_dist_fn; ctx is a struct guide_dist_ctx*. */
EXTERN float** guide_dist_pair_matrix(struct msa* msa, int* samples, int n, void* ctx);

/* numseq x num_anchors matrix, rows padded to a multiple of 8 and allocated
   exactly like d_estimation(..., pair = 0) (kalign_aligned_malloc; free the rows with kalign_aligned_free). */
EXTERN float** guide_dist_anchor_matrix(struct msa* msa, struct guide_dist_ctx* c,
                                        const int* anchors, int num_anchors);

#undef GUIDE_DIST_IMPORT
#undef EXTERN

#endif
