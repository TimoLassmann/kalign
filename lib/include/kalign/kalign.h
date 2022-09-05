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


struct msa;
/* input output routines  */

EXTERN int kalign_read_input(char* infile, struct msa** msa,int verbose);

EXTERN int kalign_write_msa(struct msa *msa, char *outfile, char *format);

EXTERN int kalign_msa_to_arr(struct msa *msa, char ***aligned, int *out_aln_len);
/* Used to convert sequences read by non-kalign code into the msa struct.. */
EXTERN int kalign_arr_to_msa(char **input_sequences, int *len, int numseq, struct msa **multiple_aln);

EXTERN int kalign_run(struct msa *msa, int n_threads, int type, float gpo, float gpe, float tgpe);

/* Memory */
EXTERN void kalign_free_msa(struct msa* msa);

/* Auxillary...  */
EXTERN int reformat_settings_msa(struct msa *msa, int rename, int unalign);
EXTERN int kalign_check_msa(struct msa* msa);

#undef KALIGN_IMPORT
#undef EXTERN


#endif
