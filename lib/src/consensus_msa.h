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

EXTERN int build_consensus(struct poar_table* table,
                           int* seq_lengths, int numseq,
                           int min_support,
                           struct msa* out_msa);

EXTERN int score_alignment_poar(struct poar_table* table,
                                struct pos_matrix* pm,
                                int numseq,
                                int n_alignments,
                                double* out_score);

EXTERN int compute_residue_confidence(struct poar_table* table,
                                      struct msa* aligned_msa);

#undef CONSENSUS_MSA_IMPORT
#undef EXTERN

#endif
