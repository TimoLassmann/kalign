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

struct sparse_bonus {
        int*   cols;    /* cols[i * K + k] = column index, or -1 if unused */
        float* vals;    /* vals[i * K + k] = bonus value */
        int    n_rows;  /* = len_a (number of DP rows) */
        int    K;       /* max entries per row (= n_anchors) */
};

static inline float sparse_bonus_lookup(const struct sparse_bonus* sb, int i, int j)
{
        float bonus = 0.0f;
        const int base = i * sb->K;
        int k;
        for(k = 0; k < sb->K; k++){
                if(sb->cols[base + k] < 0) break;
                if(sb->cols[base + k] == j)
                        bonus += sb->vals[base + k];
        }
        return bonus;
}

EXTERN void sparse_bonus_free(struct sparse_bonus* sb);

EXTERN int anchor_consistency_build(struct msa* msa, struct aln_param* ap,
                                    int n_anchors, float weight,
                                    struct consistency_table** ct_out);

EXTERN int anchor_consistency_get_bonus(struct consistency_table* ct,
                                        int seq_a, int len_a,
                                        int seq_b, int len_b,
                                        struct sparse_bonus** bonus_out);

EXTERN int anchor_consistency_get_bonus_profile(struct consistency_table* ct,
                                                struct msa* msa,
                                                int node_a, int len_a,
                                                int node_b, int len_b,
                                                struct sparse_bonus** bonus_out);

EXTERN void anchor_consistency_free(struct consistency_table* ct);

#undef ANCHOR_CONSISTENCY_IMPORT
#undef EXTERN

#endif
