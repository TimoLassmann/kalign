#ifndef CONSENSUS_FILL_H
#define CONSENSUS_FILL_H

#ifdef CONSENSUS_FILL_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;
struct kalign_run_config;

/* Fill in a POAR consensus alignment by anchored realignment.

   Columns of the consensus holding at least
   max(2, ceil(anchor_frac * numseq)) residues are anchors and are kept
   exactly as they are.  Every maximal run of non-anchor columns (a
   segment) is dissolved: each sequence's residues in the segment (a
   contiguous stretch of that sequence) are realigned with a plain
   progressive Kalign run using the substitution scores and gap
   penalties of cfg, and the realigned block is spliced back between
   the flanking anchors.

   mode KALIGN_FILL_REALIGN realigns each segment from scratch;
   mode KALIGN_FILL_PROJECT instead places the segment residues where
   the guide alignment (a full progressive run over the same sequences,
   same row order) puts them.

   cons must be a final (gapped, ALN_STATUS_FINAL) alignment whose
   seq->seq strings hold the original residue letters. */
EXTERN int consensus_fill_in(struct msa* cons,
                             const struct kalign_run_config* cfg,
                             const struct msa* guide,
                             int mode,
                             float anchor_frac,
                             int n_threads);

#undef CONSENSUS_FILL_IMPORT
#undef EXTERN

#endif
