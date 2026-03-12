#ifndef KALIGN_H
#define KALIGN_H

#include <stdint.h>
#include <kalign/kalign_config.h>

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

/* Substitution matrix constants.
   Every value maps to exactly one scoring table. No duplicates.
   Value 4 is reserved for legacy GONNET (KALIGN_TYPE_PROTEIN_DIVERGENT). */
#define KALIGN_MATRIX_AUTO         0  /* auto-select for biotype          */
#define KALIGN_MATRIX_PFASUM43     1  /* 1/3 bit, divergent protein       */
#define KALIGN_MATRIX_PFASUM60     2  /* 1/3 bit, moderate protein        */
#define KALIGN_MATRIX_CORBLOSUM66  3  /* 1/3 bit, close protein           */
#define KALIGN_MATRIX_DNA          5  /* DNA match/mismatch (+5/-4)       */
#define KALIGN_MATRIX_DNA_INTERNAL 6  /* DNA internal (tgpe=8)            */
#define KALIGN_MATRIX_RNA          7  /* RNA RIBOSUM-like (~160-383)      */

/* Backward compatibility — old KALIGN_TYPE_* map to KALIGN_MATRIX_*.
   KALIGN_TYPE_PROTEIN_DIVERGENT stays at 4 (GONNET, dead code). */
#define KALIGN_TYPE_DNA              KALIGN_MATRIX_DNA
#define KALIGN_TYPE_DNA_INTERNAL     KALIGN_MATRIX_DNA_INTERNAL
#define KALIGN_TYPE_RNA              KALIGN_MATRIX_RNA
#define KALIGN_TYPE_PROTEIN          KALIGN_MATRIX_PFASUM43
#define KALIGN_TYPE_PROTEIN_DIVERGENT 4  /* GONNET — dead code, do not use */
#define KALIGN_TYPE_PROTEIN_PFASUM43 KALIGN_MATRIX_PFASUM43
#define KALIGN_TYPE_PROTEIN_PFASUM60 KALIGN_MATRIX_PFASUM60
#define KALIGN_TYPE_PROTEIN_PFASUM_AUTO KALIGN_MATRIX_AUTO
#define KALIGN_TYPE_UNDEFINED        KALIGN_MATRIX_AUTO
#define KALIGN_TYPE_PROTEIN_CORBLOSUM66 KALIGN_MATRIX_CORBLOSUM66

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

EXTERN int kalign_ensemble_custom(struct msa* msa, int n_threads, int type,
                                  int n_runs,
                                  const float* run_gpo,
                                  const float* run_gpe,
                                  const float* run_tgpe,
                                  const int* run_types,
                                  const float* run_noise,
                                  uint64_t seed, int min_support,
                                  int refine, float vsm_amax,
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

/* Unified alignment entry point — all callers should use this. */
EXTERN struct kalign_run_config kalign_run_config_defaults(void);
EXTERN struct kalign_ensemble_config kalign_ensemble_config_defaults(void);

EXTERN int kalign_align_full(struct msa* msa,
                             const struct kalign_run_config* runs,
                             int n_runs,
                             const struct kalign_ensemble_config* ens,
                             int n_threads);

/* Expand a base config into n_runs configs using the built-in diversity table.
   Caller must pass resolved (non-sentinel) gap penalties in base.
   Caller allocates out[n_runs]. */
EXTERN int kalign_generate_ensemble_runs(const struct kalign_run_config* base,
                                         int n_runs, uint64_t seed,
                                         struct kalign_run_config* out);

/* Get a built-in mode preset.
 *
 * Presets were derived from NSGA-III multi-objective optimization
 * (objectives: F1, TC, wall_time) with 5-fold cross-validation on
 * BAliBASE v4 (protein) and BRAliBASE (RNA).
 *
 * mode:    "fast", "default", or "accurate" (case-insensitive).
 *          NULL is treated as "default".
 * biotype: ALN_BIOTYPE_PROTEIN, ALN_BIOTYPE_DNA, or ALN_BIOTYPE_RNA.
 *          Determines which preset grid slot to use.
 * runs:    caller-allocated array of at least KALIGN_MAX_PRESET_RUNS configs.
 * n_runs:  filled with the number of runs in the preset.
 * ens:     filled with ensemble config (only meaningful when *n_runs > 1).
 *
 * Returns 0 on success, -1 if mode is unknown. */
#define KALIGN_MAX_PRESET_RUNS 8

EXTERN int kalign_get_mode_preset(const char *mode,
                                   int biotype,
                                   struct kalign_run_config *runs,
                                   int *n_runs,
                                   struct kalign_ensemble_config *ens);

#undef KALIGN_IMPORT
#undef EXTERN

#endif
