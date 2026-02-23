#include <string.h>
#include <float.h>

#include "tldevel.h"

#include "msa_struct.h"
#include "aln_param.h"
#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"

#define ANCHOR_CONSISTENCY_IMPORT
#include "anchor_consistency.h"

/* Run a pairwise seq-seq alignment and extract position map.
   Returns map[0..len_i-1] where map[p] = position in seq_j that
   aligns to position p of seq_i, or -1 if gapped. */
static int pairwise_align_map(struct aln_param* ap,
                               const uint8_t* s_i, int len_i,
                               const uint8_t* s_j, int len_j,
                               int** map_out)
{
        struct aln_mem* m = NULL;
        int* posmap = NULL;
        int c, pos_a, pos_b;
        int swapped = 0;

        MMALLOC(posmap, sizeof(int) * len_i);
        for(c = 0; c < len_i; c++){
                posmap[c] = -1;
        }

        RUN(alloc_aln_mem(&m, 256));
        m->ap = ap;
        m->mode = ALN_MODE_FULL;
        m->run_parallel = 0;
        m->flip_threshold = 0.0F;
        m->flip_trial = 0;
        m->flip_stride = 1;
        m->flip_counter = 0;
        m->flip_mask = 0;
        m->margin_sum = 0.0F;
        m->margin_count = 0;

        /* Kalign requires seq1 = shorter, seq2 = longer */
        if(len_i <= len_j){
                m->len_a = len_i;
                m->len_b = len_j;
                m->enda = len_i;
                m->endb = len_j;
                m->seq1 = s_i;
                m->seq2 = s_j;
                swapped = 0;
        }else{
                m->len_a = len_j;
                m->len_b = len_i;
                m->enda = len_j;
                m->endb = len_i;
                m->seq1 = s_j;
                m->seq2 = s_i;
                swapped = 1;
        }
        m->prof1 = NULL;
        m->prof2 = NULL;

        m->f[0].a = 0.0F;
        m->f[0].ga = -FLT_MAX;
        m->f[0].gb = -FLT_MAX;
        m->b[0].a = 0.0F;
        m->b[0].ga = -FLT_MAX;
        m->b[0].gb = -FLT_MAX;

        RUN(init_alnmem(m));
        aln_runner(m);

        if(swapped){
                RUN(mirror_path_n(m, len_i, len_j));
                m->len_a = len_i;
                m->len_b = len_j;
        }

        RUN(add_gap_info_to_path_n(m));

        /* Extract position map from alignment path.
           After add_gap_info_to_path_n:
           path[0] = alignment length
           path[c]: 0=match (both advance), 1=gap in a (b advances), 2=gap in b (a advances), 3=end */
        pos_a = 0;
        pos_b = 0;
        c = 1;
        while(m->path[c] != 3){
                if(m->path[c] == 0){
                        /* match: a[pos_a] aligns to b[pos_b] */
                        if(pos_a < len_i){
                                posmap[pos_a] = pos_b;
                        }
                        pos_a++;
                        pos_b++;
                }else if(m->path[c] & 1){
                        /* gap in a: only b advances */
                        pos_b++;
                }else if(m->path[c] & 2){
                        /* gap in b: only a advances (no partner) */
                        if(pos_a < len_i){
                                posmap[pos_a] = -1;
                        }
                        pos_a++;
                }
                c++;
        }

        free_aln_mem(m);
        *map_out = posmap;
        return OK;
ERROR:
        if(posmap) MFREE(posmap);
        if(m) free_aln_mem(m);
        return FAIL;
}

/* Select K diverse anchor sequences using farthest-first traversal
   on msa->seq_distances[] (per-sequence mean distance). */
static int select_anchors(struct msa* msa, int K, int* anchor_ids)
{
        int i, k;
        float* min_dist = NULL;  /* min distance to any selected anchor */
        int N = msa->numseq;

        if(K > N) K = N;

        MMALLOC(min_dist, sizeof(float) * N);

        /* Pick first anchor: sequence with median seq_distances[] value.
           Use a simple approach: pick the one closest to the mean. */
        {
                float sum = 0.0f;
                float best_diff = FLT_MAX;
                int best_idx = 0;
                for(i = 0; i < N; i++){
                        sum += msa->seq_distances[i];
                }
                float mean = sum / (float)N;
                for(i = 0; i < N; i++){
                        float diff = msa->seq_distances[i] - mean;
                        if(diff < 0) diff = -diff;
                        if(diff < best_diff){
                                best_diff = diff;
                                best_idx = i;
                        }
                }
                anchor_ids[0] = best_idx;
        }

        /* Initialize min_dist: distance from each seq to first anchor.
           We use |seq_distances[i] - seq_distances[anchor]| as a proxy
           since we only have per-sequence mean distances, not pairwise. */
        for(i = 0; i < N; i++){
                float d = msa->seq_distances[i] - msa->seq_distances[anchor_ids[0]];
                if(d < 0) d = -d;
                min_dist[i] = d;
        }

        /* Farthest-first: pick remaining K-1 anchors */
        for(k = 1; k < K; k++){
                float best_min = -1.0f;
                int best_idx = 0;
                for(i = 0; i < N; i++){
                        /* Skip already-selected anchors */
                        int skip = 0;
                        for(int j = 0; j < k; j++){
                                if(anchor_ids[j] == i){ skip = 1; break; }
                        }
                        if(skip) continue;

                        if(min_dist[i] > best_min){
                                best_min = min_dist[i];
                                best_idx = i;
                        }
                }
                anchor_ids[k] = best_idx;

                /* Update min_dist with new anchor */
                for(i = 0; i < N; i++){
                        float d = msa->seq_distances[i] - msa->seq_distances[best_idx];
                        if(d < 0) d = -d;
                        if(d < min_dist[i]){
                                min_dist[i] = d;
                        }
                }
        }

        MFREE(min_dist);
        return OK;
ERROR:
        if(min_dist) MFREE(min_dist);
        return FAIL;
}

int anchor_consistency_build(struct msa* msa, struct aln_param* ap,
                             int n_anchors, float weight,
                             struct consistency_table** ct_out)
{
        struct consistency_table* ct = NULL;
        int N = msa->numseq;
        int K = n_anchors;
        int i, k;

        if(K <= 0 || N < 3){
                *ct_out = NULL;
                return OK;
        }
        if(K > N) K = N;

        /* Check that seq_distances are available */
        if(msa->seq_distances == NULL){
                *ct_out = NULL;
                return OK;
        }

        MMALLOC(ct, sizeof(struct consistency_table));
        ct->pos_maps = NULL;
        ct->map_lengths = NULL;
        ct->anchor_ids = NULL;
        ct->n_anchors = K;
        ct->numseq = N;
        ct->weight = weight;

        MMALLOC(ct->anchor_ids, sizeof(int) * K);
        MMALLOC(ct->pos_maps, sizeof(int*) * N * K);
        MMALLOC(ct->map_lengths, sizeof(int) * N * K);

        for(i = 0; i < N * K; i++){
                ct->pos_maps[i] = NULL;
                ct->map_lengths[i] = 0;
        }

        /* Select diverse anchors */
        RUN(select_anchors(msa, K, ct->anchor_ids));

        if(!msa->quiet){
                LOG_MSG("Anchor consistency: K=%d, weight=%.1f", K, weight);
        }

        /* Build position maps: for each sequence i, for each anchor k */
        for(i = 0; i < N; i++){
                int len_i = msa->sequences[i]->len;
                for(k = 0; k < K; k++){
                        int ak = ct->anchor_ids[k];
                        ct->map_lengths[i * K + k] = len_i;

                        if(i == ak){
                                /* Identity map for anchor itself */
                                int p;
                                MMALLOC(ct->pos_maps[i * K + k], sizeof(int) * len_i);
                                for(p = 0; p < len_i; p++){
                                        ct->pos_maps[i * K + k][p] = p;
                                }
                        }else{
                                /* Pairwise alignment: seq_i vs anchor_k */
                                RUN(pairwise_align_map(ap,
                                                       msa->sequences[i]->s, len_i,
                                                       msa->sequences[ak]->s, msa->sequences[ak]->len,
                                                       &ct->pos_maps[i * K + k]));
                        }
                }
        }

        *ct_out = ct;
        return OK;
ERROR:
        anchor_consistency_free(ct);
        *ct_out = NULL;
        return FAIL;
}

int anchor_consistency_get_bonus(struct consistency_table* ct,
                                 int seq_a, int len_a,
                                 int seq_b, int len_b,
                                 float** bonus_out)
{
        float* bonus = NULL;
        int* inv_b = NULL;
        int K = ct->n_anchors;
        int k, i;
        float per_anchor_weight = ct->weight / (float)K;

        MMALLOC(bonus, sizeof(float) * len_a * len_b);
        memset(bonus, 0, sizeof(float) * len_a * len_b);

        /* For each anchor, build inverse map for seq_b (anchor_pos → seq_b_pos),
           then scan seq_a's map to accumulate bonus. O(K * (La + Lb + max_anchor_len)) */
        for(k = 0; k < K; k++){
                int* map_a = ct->pos_maps[seq_a * K + k];
                int* map_b = ct->pos_maps[seq_b * K + k];
                int anchor_len = 0;
                int j;

                if(map_a == NULL || map_b == NULL) continue;

                /* Determine anchor length for inverse map */
                /* Find max mapped position to size the inverse map */
                anchor_len = 0;
                for(i = 0; i < len_a; i++){
                        if(map_a[i] >= anchor_len) anchor_len = map_a[i] + 1;
                }
                for(j = 0; j < len_b; j++){
                        if(map_b[j] >= anchor_len) anchor_len = map_b[j] + 1;
                }
                if(anchor_len == 0) continue;

                /* Build inverse map: anchor_pos → seq_b_pos */
                MMALLOC(inv_b, sizeof(int) * anchor_len);
                for(j = 0; j < anchor_len; j++){
                        inv_b[j] = -1;
                }
                for(j = 0; j < len_b; j++){
                        if(map_b[j] >= 0 && map_b[j] < anchor_len){
                                inv_b[map_b[j]] = j;
                        }
                }

                /* Accumulate bonus */
                for(i = 0; i < len_a; i++){
                        int ak_pos = map_a[i];
                        if(ak_pos >= 0 && ak_pos < anchor_len){
                                int bj = inv_b[ak_pos];
                                if(bj >= 0){
                                        bonus[i * len_b + bj] += per_anchor_weight;
                                }
                        }
                }

                MFREE(inv_b);
                inv_b = NULL;
        }

        *bonus_out = bonus;
        return OK;
ERROR:
        if(bonus) MFREE(bonus);
        if(inv_b) MFREE(inv_b);
        *bonus_out = NULL;
        return FAIL;
}

/* Compute consensus anchor positions for a profile node.
   For a leaf (nsip==1), positions come directly from pos_maps.
   For a profile (nsip>1), positions are derived by majority vote
   across member sequences, using gaps[] to map profile columns
   to ungapped positions. */
static int get_node_anchor_positions(struct consistency_table* ct,
                                      struct msa* msa,
                                      int node, int dp_len, int anchor_k,
                                      int* positions, float* confidence)
{
        int K = ct->n_anchors;
        int i, c, p;

        if(msa->nsip[node] == 1){
                /* Leaf: direct lookup from pos_maps */
                int* map = ct->pos_maps[node * K + anchor_k];
                int seq_len = msa->sequences[node]->len;
                if(map == NULL){
                        for(i = 0; i < dp_len; i++){
                                positions[i] = -1;
                                confidence[i] = 0.0f;
                        }
                        return OK;
                }
                for(i = 0; i < dp_len && i < seq_len; i++){
                        positions[i] = map[i];
                        confidence[i] = (map[i] >= 0) ? 1.0f : 0.0f;
                }
                for(; i < dp_len; i++){
                        positions[i] = -1;
                        confidence[i] = 0.0f;
                }
        }else{
                /* Profile: consensus from member sequences */
                int n_members = msa->nsip[node];
                int* members = msa->sip[node];
                int* col_to_ungapped = NULL;
                int* best_pos = NULL;
                int* agree_count = NULL;
                int* total_count = NULL;
                int mi, col;

                MMALLOC(col_to_ungapped, sizeof(int) * (dp_len + 1));
                MMALLOC(best_pos, sizeof(int) * dp_len);
                MMALLOC(agree_count, sizeof(int) * dp_len);
                MMALLOC(total_count, sizeof(int) * dp_len);

                for(c = 0; c < dp_len; c++){
                        best_pos[c] = -1;
                        agree_count[c] = 0;
                        total_count[c] = 0;
                }

                for(mi = 0; mi < n_members; mi++){
                        int si = members[mi];
                        int* map;
                        int seq_len;
                        int* gaps;
                        int apos;

                        if(si >= ct->numseq) continue;
                        map = ct->pos_maps[si * K + anchor_k];
                        if(map == NULL) continue;

                        seq_len = msa->sequences[si]->len;
                        gaps = msa->sequences[si]->gaps;

                        /* Build col_to_ungapped from gaps[] */
                        col = 0;
                        for(p = 0; p <= seq_len && col < dp_len; p++){
                                int g;
                                for(g = 0; g < gaps[p] && col < dp_len; g++){
                                        col_to_ungapped[col] = -1;
                                        col++;
                                }
                                if(p < seq_len && col < dp_len){
                                        col_to_ungapped[col] = p;
                                        col++;
                                }
                        }
                        while(col < dp_len){
                                col_to_ungapped[col] = -1;
                                col++;
                        }

                        /* Accumulate votes for each profile column */
                        for(c = 0; c < dp_len; c++){
                                int ugp = col_to_ungapped[c];
                                if(ugp < 0 || ugp >= seq_len) continue;
                                apos = map[ugp];
                                if(apos < 0) continue;

                                total_count[c]++;
                                if(best_pos[c] < 0){
                                        best_pos[c] = apos;
                                        agree_count[c] = 1;
                                }else if(apos == best_pos[c]){
                                        agree_count[c]++;
                                }
                        }
                }

                for(c = 0; c < dp_len; c++){
                        if(total_count[c] > 0 && agree_count[c] > 0){
                                positions[c] = best_pos[c];
                                confidence[c] = (float)agree_count[c] / (float)total_count[c];
                        }else{
                                positions[c] = -1;
                                confidence[c] = 0.0f;
                        }
                }

                MFREE(col_to_ungapped);
                MFREE(best_pos);
                MFREE(agree_count);
                MFREE(total_count);
        }
        return OK;
ERROR:
        return FAIL;
}

int anchor_consistency_get_bonus_profile(struct consistency_table* ct,
                                          struct msa* msa,
                                          int node_a, int len_a,
                                          int node_b, int len_b,
                                          float** bonus_out)
{
        float* bonus = NULL;
        int* apos_a = NULL;
        float* conf_a = NULL;
        int* apos_b = NULL;
        float* conf_b = NULL;
        int* inv_b = NULL;
        float* inv_conf_b = NULL;
        int K = ct->n_anchors;
        int k, i, j;
        float per_anchor_weight = ct->weight / (float)K;

        MMALLOC(bonus, sizeof(float) * len_a * len_b);
        memset(bonus, 0, sizeof(float) * len_a * len_b);

        MMALLOC(apos_a, sizeof(int) * len_a);
        MMALLOC(conf_a, sizeof(float) * len_a);
        MMALLOC(apos_b, sizeof(int) * len_b);
        MMALLOC(conf_b, sizeof(float) * len_b);

        for(k = 0; k < K; k++){
                int anchor_len;

                RUN(get_node_anchor_positions(ct, msa, node_a, len_a, k,
                                               apos_a, conf_a));
                RUN(get_node_anchor_positions(ct, msa, node_b, len_b, k,
                                               apos_b, conf_b));

                /* Find max anchor position to size the inverse map */
                anchor_len = 0;
                for(i = 0; i < len_a; i++){
                        if(apos_a[i] >= anchor_len) anchor_len = apos_a[i] + 1;
                }
                for(j = 0; j < len_b; j++){
                        if(apos_b[j] >= anchor_len) anchor_len = apos_b[j] + 1;
                }
                if(anchor_len == 0) continue;

                /* Build inverse map: anchor_pos -> node_b DP col */
                MMALLOC(inv_b, sizeof(int) * anchor_len);
                MMALLOC(inv_conf_b, sizeof(float) * anchor_len);
                for(j = 0; j < anchor_len; j++){
                        inv_b[j] = -1;
                        inv_conf_b[j] = 0.0f;
                }
                for(j = 0; j < len_b; j++){
                        if(apos_b[j] >= 0 && apos_b[j] < anchor_len){
                                inv_b[apos_b[j]] = j;
                                inv_conf_b[apos_b[j]] = conf_b[j];
                        }
                }

                /* Accumulate bonus */
                for(i = 0; i < len_a; i++){
                        int ak_pos = apos_a[i];
                        if(ak_pos >= 0 && ak_pos < anchor_len){
                                int bj = inv_b[ak_pos];
                                if(bj >= 0){
                                        bonus[i * len_b + bj] +=
                                                per_anchor_weight * conf_a[i] * inv_conf_b[ak_pos];
                                }
                        }
                }

                MFREE(inv_b);
                inv_b = NULL;
                MFREE(inv_conf_b);
                inv_conf_b = NULL;
        }

        MFREE(apos_a);
        MFREE(conf_a);
        MFREE(apos_b);
        MFREE(conf_b);

        *bonus_out = bonus;
        return OK;
ERROR:
        if(bonus) MFREE(bonus);
        if(apos_a) MFREE(apos_a);
        if(conf_a) MFREE(conf_a);
        if(apos_b) MFREE(apos_b);
        if(conf_b) MFREE(conf_b);
        if(inv_b) MFREE(inv_b);
        if(inv_conf_b) MFREE(inv_conf_b);
        *bonus_out = NULL;
        return FAIL;
}

void anchor_consistency_free(struct consistency_table* ct)
{
        if(ct){
                if(ct->pos_maps){
                        int total = ct->numseq * ct->n_anchors;
                        for(int i = 0; i < total; i++){
                                if(ct->pos_maps[i]){
                                        MFREE(ct->pos_maps[i]);
                                }
                        }
                        MFREE(ct->pos_maps);
                }
                if(ct->map_lengths) MFREE(ct->map_lengths);
                if(ct->anchor_ids) MFREE(ct->anchor_ids);
                MFREE(ct);
        }
}
