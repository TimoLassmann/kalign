#ifndef KALIGN_H
#define KALIGN_H

#ifdef KALIGN_IMPORT
   #define EXTERN
#else
   #ifndef EXTERN
      #ifdef __cplusplus
         #define EXTERN extern "C"
      #else
         #define EXTERN extern
      #endif
   #endif
#endif

#define KALIGN_TYPE_DNA 0
#define KALIGN_TYPE_DNA_INTERNAL 1
#define KALIGN_TYPE_RNA 2
#define KALIGN_TYPE_PROTEIN 3
#define KALIGN_TYPE_PROTEIN_DIVERGENT 4
#define KALIGN_TYPE_UNDEFINED 5

struct msa;
/* input output routines  */

EXTERN int kalign_read_input(char* infile, struct msa** msa,int quiet);

EXTERN int kalign_write_msa(struct msa *msa, char *outfile, char *format);

/* EXTERN int kalign_msa_to_arr(struct msa *msa, char ***aligned, int *out_aln_len); */
/* Used to convert sequences read by non-kalign code into the msa struct.. */
/* EXTERN int kalign_arr_to_msa(char **input_sequences, int *len, int numseq, struct msa **multiple_aln); */


EXTERN int kalign(char **seq, int *len, int numseq, int n_threads, int type,
                  float gpo, float gpe, float tgpe, char ***aligned,
                  int *out_aln_len);

EXTERN int kalign_run(struct msa *msa, int n_threads, int type, float gpo, float gpe, float tgpe);
/* Memory */
EXTERN void kalign_free_msa(struct msa* msa);

/* Auxillary...  */
EXTERN int reformat_settings_msa(struct msa *msa, int rename, int unalign);

EXTERN int kalign_check_msa(struct msa* msa, int exit_on_error);

EXTERN int kalign_msa_compare(struct msa *r, struct msa *t, float *score);
#undef KALIGN_IMPORT
#undef EXTERN

#endif
