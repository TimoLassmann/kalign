#include "tldevel.h"

#include "msa_struct.h"
#include "anchor_consistency.h"  /* for sparse_bonus, sparse_bonus_free */
#include "poar.h"

#define MSA_CONSISTENCY_IMPORT
#include "msa_consistency.h"

/* Flat index for pair (i,j) where i < j.
   Duplicated from poar.c (static there). */
static inline int poar_pair_idx(int i, int j, int numseq)
{
        return i * numseq - (i * (i + 1)) / 2 + (j - i - 1);
}

/* Build ungapped_pos -> DP_position map for a member sequence.
 * Walks gaps[] to find the DP column for each ungapped residue.
 * Caller must free the returned array. */
static int build_ungapped_to_dp(struct msa* msa, int seq_idx,
                                 int** map_out, int* ungapped_len_out)
{
        int* gaps = msa->sequences[seq_idx]->gaps;
        int seq_len = msa->sequences[seq_idx]->len;
        int* map = NULL;
        int dp, p, g;

        if(seq_len == 0){
                *map_out = NULL;
                *ungapped_len_out = 0;
                return OK;
        }

        MMALLOC(map, sizeof(int) * seq_len);

        dp = 0;
        for(p = 0; p <= seq_len; p++){
                for(g = 0; g < gaps[p]; g++){
                        dp++;
                }
                if(p < seq_len){
                        map[p] = dp;
                        dp++;
                }
        }

        *map_out = map;
        *ungapped_len_out = seq_len;
        return OK;
ERROR:
        if(map) MFREE(map);
        *map_out = NULL;
        *ungapped_len_out = 0;
        return FAIL;
}

#define POAR_BONUS_K 16   /* max distinct target positions per DP row */

int poar_consistency_get_bonus(struct poar_consistency_ctx* ctx,
                                struct msa* msa,
                                int node_a, int len_a,
                                int node_b, int len_b,
                                struct sparse_bonus** bonus_out)
{
        struct sparse_bonus* sb = NULL;
        int* map_a = NULL;
        int* map_b = NULL;
        int K = POAR_BONUS_K;
        int numseq = ctx->poar->numseq;
        int n_members_a = msa->nsip[node_a];
        int n_members_b = msa->nsip[node_b];
        int* members_a = msa->sip[node_a];
        int* members_b = msa->sip[node_b];
        float denom = (float)n_members_a * (float)n_members_b * (float)ctx->n_runs;
        float pair_weight = ctx->weight / denom;
        int ma, mb, e, i;

        /* Allocate sparse bonus */
        MMALLOC(sb, sizeof(struct sparse_bonus));
        sb->cols = NULL;
        sb->vals = NULL;
        sb->n_rows = len_a;
        sb->K = K;

        MMALLOC(sb->cols, sizeof(int) * len_a * K);
        MMALLOC(sb->vals, sizeof(float) * len_a * K);
        for(i = 0; i < len_a * K; i++){
                sb->cols[i] = -1;
                sb->vals[i] = 0.0f;
        }

        /* For each pair of member sequences, look up POAR entries */
        for(ma = 0; ma < n_members_a; ma++){
                int sa = members_a[ma];
                int sa_len = 0;

                RUN(build_ungapped_to_dp(msa, sa, &map_a, &sa_len));

                for(mb = 0; mb < n_members_b; mb++){
                        int sb_seq = members_b[mb];
                        int sb_len = 0;
                        int lo, hi;
                        int pidx;
                        struct poar_pair* pp;
                        int swapped;

                        if(sa == sb_seq) continue;

                        RUN(build_ungapped_to_dp(msa, sb_seq, &map_b, &sb_len));

                        /* POAR pairs are stored with lo < hi */
                        if(sa < sb_seq){
                                lo = sa; hi = sb_seq;
                                swapped = 0;
                        }else{
                                lo = sb_seq; hi = sa;
                                swapped = 1;
                        }

                        pidx = poar_pair_idx(lo, hi, numseq);
                        pp = ctx->poar->pairs[pidx];

                        if(pp != NULL){
                                for(e = 0; e < pp->n_entries; e++){
                                        uint32_t key = pp->entries[e].key;
                                        uint32_t support = pp->entries[e].support;
                                        int pos_lo = (int)(key >> 20);
                                        int pos_hi = (int)(key & 0xFFFFF);
                                        int pos_sa, pos_sb;
                                        int dp_row, dp_col;
                                        int count, base, slot, s;
                                        float val;

                                        if(swapped){
                                                pos_sa = pos_hi;
                                                pos_sb = pos_lo;
                                        }else{
                                                pos_sa = pos_lo;
                                                pos_sb = pos_hi;
                                        }

                                        if(pos_sa >= sa_len || pos_sb >= sb_len) continue;
                                        if(map_a == NULL || map_b == NULL) continue;

                                        dp_row = map_a[pos_sa];
                                        dp_col = map_b[pos_sb];

                                        if(dp_row >= len_a || dp_col >= len_b) continue;

                                        count = __builtin_popcount(support);
                                        val = pair_weight * (float)count;

                                        /* Insert into sparse bonus */
                                        base = dp_row * K;
                                        slot = -1;
                                        for(s = 0; s < K; s++){
                                                if(sb->cols[base + s] == dp_col){
                                                        slot = s; break;
                                                }
                                                if(sb->cols[base + s] < 0){
                                                        slot = s; break;
                                                }
                                        }
                                        if(slot >= 0){
                                                sb->vals[base + slot] += val;
                                                sb->cols[base + slot] = dp_col;
                                        }
                                }
                        }

                        if(map_b){ MFREE(map_b); map_b = NULL; }
                }

                if(map_a){ MFREE(map_a); map_a = NULL; }
        }

        *bonus_out = sb;
        return OK;
ERROR:
        sparse_bonus_free(sb);
        if(map_a) MFREE(map_a);
        if(map_b) MFREE(map_b);
        *bonus_out = NULL;
        return FAIL;
}
