#ifndef ANCHOR_CONSISTENCY_H
#define ANCHOR_CONSISTENCY_H

#ifdef ANCHOR_CONSISTENCY_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;
struct aln_param;

struct consistency_table {
        int** pos_maps;     /* [seq_idx * K + anchor_slot][pos] = anchor_pos or -1 */
        int* map_lengths;   /* [seq_idx * K + anchor_slot] = length of seq_idx */
        int* anchor_ids;    /* [0..K-1] = sequence indices of anchors */
        int n_anchors;      /* K */
        int numseq;         /* N */
        float weight;       /* scaling factor for bonus */
};

EXTERN int anchor_consistency_build(struct msa* msa, struct aln_param* ap,
                                    int n_anchors, float weight,
                                    struct consistency_table** ct_out);

EXTERN int anchor_consistency_get_bonus(struct consistency_table* ct,
                                        int seq_a, int len_a,
                                        int seq_b, int len_b,
                                        float** bonus_out);

EXTERN int anchor_consistency_get_bonus_profile(struct consistency_table* ct,
                                                struct msa* msa,
                                                int node_a, int len_a,
                                                int node_b, int len_b,
                                                float** bonus_out);

EXTERN void anchor_consistency_free(struct consistency_table* ct);

#undef ANCHOR_CONSISTENCY_IMPORT
#undef EXTERN

#endif
