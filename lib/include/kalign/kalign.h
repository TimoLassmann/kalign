#ifndef KALIGN_H
#define KALIGN_H

#include <stdint.h>

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
#define KALIGN_TYPE_PROTEIN_PFASUM43 5
#define KALIGN_TYPE_PROTEIN_PFASUM60 6
#define KALIGN_TYPE_PROTEIN_PFASUM_AUTO 7
#define KALIGN_TYPE_UNDEFINED 8

#define KALIGN_REFINE_NONE 0
#define KALIGN_REFINE_ALL 1
#define KALIGN_REFINE_CONFIDENT 2
#define KALIGN_REFINE_INLINE 3

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

EXTERN int kalign_run(struct msa *msa, int n_threads, int type, float gpo, float gpe, float tgpe, int refine, int adaptive_budget);

EXTERN int kalign_run_seeded(struct msa *msa, int n_threads, int type,
                             float gpo, float gpe, float tgpe,
                             int refine, int adaptive_budget,
                             uint64_t tree_seed, float tree_noise,
                             float dist_scale, float vsm_amax,
                             float use_seq_weights,
                             int consistency_anchors, float consistency_weight);

EXTERN int kalign_run_dist_scale(struct msa *msa, int n_threads, int type,
                                  float gpo, float gpe, float tgpe,
                                  int refine, int adaptive_budget,
                                  float dist_scale, float vsm_amax,
                                  float use_seq_weights);

EXTERN int kalign_run_realign(struct msa *msa, int n_threads, int type,
                              float gpo, float gpe, float tgpe,
                              int refine, int adaptive_budget,
                              float dist_scale, float vsm_amax,
                              int realign_iterations,
                              float use_seq_weights,
                              int consistency_anchors, float consistency_weight);

EXTERN int kalign_post_realign(struct msa *msa, int n_threads, int type,
                               float gpo, float gpe, float tgpe,
                               int refine, int adaptive_budget,
                               float dist_scale, float vsm_amax,
                               int realign_iterations,
                               float use_seq_weights);

EXTERN int kalign_ensemble(struct msa* msa, int n_threads, int type,
                           int n_runs, float gpo, float gpe, float tgpe,
                           uint64_t seed, int min_support,
                           const char* save_poar_path,
                           int refine, float dist_scale, float vsm_amax,
                           int realign, float use_seq_weights,
                           int consistency_anchors, float consistency_weight);

EXTERN int kalign_consensus_from_poar(struct msa* msa,
                                      const char* poar_path,
                                      int min_support);

/* Memory */
EXTERN void kalign_free_msa(struct msa* msa);

/* Auxillary...  */
EXTERN int reformat_settings_msa(struct msa *msa, int rename, int unalign);

EXTERN int kalign_check_msa(struct msa* msa, int exit_on_error);

EXTERN int kalign_msa_compare(struct msa *r, struct msa *t, float *score);

struct poar_score;
EXTERN int kalign_msa_compare_detailed(struct msa *r, struct msa *t,
                                       float max_gap_frac,
                                       struct poar_score *out);

EXTERN int kalign_msa_compare_with_mask(struct msa *r, struct msa *t,
                                        int *scored_cols, int n_cols,
                                        struct poar_score *out);
#undef KALIGN_IMPORT
#undef EXTERN

#endif
