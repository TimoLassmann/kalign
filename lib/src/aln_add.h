#ifndef ALN_ADD_H
#define ALN_ADD_H

#ifdef ALN_ADD_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;

/* Add new unaligned sequences to an existing finalized alignment.
 *
 * Each new sequence is independently aligned against a consensus profile
 * built from the existing alignment. The existing sequences are NOT
 * modified — their gaps are preserved exactly.
 *
 * After this call, existing->sequences contains the original sequences
 * (unchanged) plus the new sequences (with gaps inserted). existing->numseq
 * is updated accordingly.
 *
 * existing: Must be a finalized alignment (ALN_STATUS_FINAL).
 * new_seqs: Unaligned sequences to add. Modified in place (seq->seq gets gaps).
 * n_threads: Number of threads for parallel alignment of new sequences.
 *
 * Returns OK on success, FAIL on error. */
EXTERN int kalign_add_sequences(struct msa* existing,
                                 struct msa* new_seqs,
                                 int n_threads);

#undef ALN_ADD_IMPORT
#undef EXTERN

#endif
