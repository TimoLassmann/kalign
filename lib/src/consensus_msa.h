#ifndef CONSENSUS_MSA_H
#define CONSENSUS_MSA_H

#ifdef CONSENSUS_MSA_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct poar_table;
struct pos_matrix;
struct msa;

typedef struct threadpool threadpool_t;

EXTERN int build_consensus(struct poar_table* table,
                           int* seq_lengths, int numseq,
                           int min_support,
                           int attach_support,
                           struct msa* out_msa,
                           threadpool_t* pool);

EXTERN int score_alignment_poar(struct poar_table* table,
                                struct pos_matrix* pm,
                                int numseq,
                                int n_alignments,
                                double* out_score,
                                threadpool_t* pool);

/* 1 if run k of the (implicit) table is the alignment pm: every residue in the same column. */
EXTERN int poar_run_is_alignment(const struct poar_table* table, const struct pos_matrix* pm, int k);

/* scores[k] = what score_alignment_poar returns for run k of the table, for every run, from one pass over
   the sequence pairs (the same doubles).  Implicit tables only. */
EXTERN int score_runs_poar(struct poar_table* table, int numseq, int n_alignments,
                           double* scores, threadpool_t* pool);

EXTERN int compute_residue_confidence(struct poar_table* table,
                                      struct msa* aligned_msa,
                                      threadpool_t* pool);

#undef CONSENSUS_MSA_IMPORT
#undef EXTERN

#endif
