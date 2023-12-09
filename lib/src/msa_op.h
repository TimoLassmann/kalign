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
EXTERN int make_linear_sequence(struct msa_seq *seq, char *linear_seq);

#undef MSA_OP_IMPORT
#undef EXTERN


#endif
