#ifndef MSA_OP_H
#define MSA_OP_H

#ifdef MSA_OP_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct msa;
struct msa_seq;

EXTERN int msa_seq_cpy(struct msa_seq *d, struct msa_seq *src);
EXTERN int msa_cpy(struct msa** dest, struct msa* src);
EXTERN int merge_msa(struct msa** dest, struct msa* src);

EXTERN int dealign_msa(struct msa *msa);

EXTERN int detect_alphabet(struct msa *msa);
EXTERN int detect_aligned(struct msa *msa);
EXTERN int set_sip_nsip(struct msa *msa);

EXTERN int reformat_settings_msa(struct msa *msa, int rename, int unalign);

EXTERN int convert_msa_to_internal(struct msa* msa, int type);

/* convert alinged msa sequences to character array */
EXTERN int kalign_msa_to_arr(struct msa *msa, char ***aligned, int *out_aln_len);
/* Used to convert sequences read by non-kalign code into the msa struct.. */
EXTERN int kalign_arr_to_msa(char **input_sequences, int *len, int numseq, struct msa **multiple_aln);

EXTERN int finalise_alignment(struct msa* msa);
EXTERN int make_linear_sequence(struct msa_seq *seq, char *linear_seq, int *out_len);

/* Confidence masking styles */
#define KALIGN_MASK_LOWERCASE 0
#define KALIGN_MASK_REMOVE    1

/* Mask low-confidence alignment columns.
   style: KALIGN_MASK_LOWERCASE (residues → lowercase) or KALIGN_MASK_REMOVE (→ gaps).
   No-op if threshold <= 0 or col_confidence is NULL (non-ensemble modes). */
EXTERN int kalign_mask_by_confidence(struct msa* msa, float threshold, int style);

/* Write per-column confidence scores to a text file (one value per line). */
EXTERN int kalign_write_confidence(struct msa* msa, const char* path);

#undef MSA_OP_IMPORT
#undef EXTERN


#endif
