#include "tldevel.h"
#include "msa_struct.h"
#include "bpm.h"

#include "threadpool.h"

struct anchor_dist_ctx {
        struct msa_seq** seqs;
        float* min_dist;
        uint8_t* anchor_s;
        int anchor_len;
        float** dist;   /* optional N x K cache of the BPM distances, or NULL */
        int k;
};

static void anchor_init_fn(int start, int end, void *arg)
{
        struct anchor_dist_ctx *c = (struct anchor_dist_ctx *)arg;
        for (int i = start; i < end; i++) {
                uint8_t *si = c->seqs[i]->s;
                int li = c->seqs[i]->len;
                if (li > c->anchor_len)
                        c->min_dist[i] = (float)BPM(si, c->anchor_s, li, c->anchor_len);
                else
                        c->min_dist[i] = (float)BPM(c->anchor_s, si, c->anchor_len, li);
                if (c->dist) c->dist[i][c->k] = c->min_dist[i];
        }
}

static void anchor_update_fn(int start, int end, void *arg)
{
        struct anchor_dist_ctx *c = (struct anchor_dist_ctx *)arg;
        for (int i = start; i < end; i++) {
                if (c->min_dist[i] < 0.0f) continue;
                uint8_t *si = c->seqs[i]->s;
                int li = c->seqs[i]->len;
                float d;
                if (li > c->anchor_len)
                        d = (float)BPM(si, c->anchor_s, li, c->anchor_len);
                else
                        d = (float)BPM(c->anchor_s, si, c->anchor_len, li);
                if (c->dist) c->dist[i][c->k] = d;
                if (d < c->min_dist[i])
                        c->min_dist[i] = d;
        }
}

#define PICK_ANCHOR_IMPORT
#include "pick_anchor.h"

#define DEFAULT_NUM_ANCHORS 32

/* Farthest-first anchor selection.
 *
 * Picks K diverse anchor sequences by greedily selecting the sequence
 * most distant (by BPM edit distance) from all already-selected anchors.
 * This ensures the anchors span the full diversity of the dataset.
 *
 * Cost: K * N BPM calls — microseconds for typical inputs.
 */
/* dist (optional): numseq x K matrix, filled with -1 by the caller; on
   return dist[i][k] holds the BPM distance between sequence i and anchor k
   wherever it was computed here, called exactly as calc_distance() calls it,
   so d_estimation can reuse it instead of recomputing it. */
static int* pick_anchor_farthest_first(struct msa* msa, int K, int* n_out, float** dist)
{
        int numseq = msa->numseq;
        int* anchors = NULL;
        float* min_dist = NULL;   /* min distance from each seq to any anchor */
        int i, k;

        if(K > numseq) K = numseq;
        if(K < 1) K = 1;

        MMALLOC(anchors, sizeof(int) * K);
        MMALLOC(min_dist, sizeof(float) * numseq);

        /* Pick first anchor: median-length sequence.
           Sort would be overkill — just find the sequence closest
           to the mean length. */
        {
                float mean_len = 0.0f;
                float best_diff = 1e30f;
                int best_idx = 0;
                for(i = 0; i < numseq; i++){
                        mean_len += (float)msa->sequences[i]->len;
                }
                mean_len /= (float)numseq;
                for(i = 0; i < numseq; i++){
                        float diff = (float)msa->sequences[i]->len - mean_len;
                        if(diff < 0) diff = -diff;
                        if(diff < best_diff){
                                best_diff = diff;
                                best_idx = i;
                        }
                }
                anchors[0] = best_idx;
        }

        /* Initialize min_dist: BPM distance from each seq to first anchor */
        {
                uint8_t* anchor_s = msa->sequences[anchors[0]]->s;
                int anchor_len = msa->sequences[anchors[0]]->len;
                if(msa->pool && numseq >= KALIGN_DIST_MIN_SEQS){
                        struct anchor_dist_ctx ctx = { msa->sequences, min_dist, anchor_s, anchor_len, dist, 0 };
                        tp_parallel_for_chunked(msa->pool, 0, numseq, KALIGN_PFOR_MIN_CHUNK, anchor_init_fn, &ctx);
                }else{
                for(i = 0; i < numseq; i++){
                        uint8_t* si = msa->sequences[i]->s;
                        int li = msa->sequences[i]->len;
                        if(li > anchor_len){
                                min_dist[i] = (float)BPM(si, anchor_s, li, anchor_len);
                        }else{
                                min_dist[i] = (float)BPM(anchor_s, si, anchor_len, li);
                        }
                        if(dist) dist[i][0] = min_dist[i];
                }
                }
                min_dist[anchors[0]] = -1.0f;  /* mark as selected */
        }

        /* Farthest-first: pick remaining K-1 anchors */
        for(k = 1; k < K; k++){
                /* Find the sequence with largest min_dist */
                float best_min = -1.0f;
                int best_idx = 0;
                for(i = 0; i < numseq; i++){
                        if(min_dist[i] > best_min){
                                best_min = min_dist[i];
                                best_idx = i;
                        }
                }
                anchors[k] = best_idx;
                min_dist[best_idx] = -1.0f;  /* mark as selected */

                /* Update min_dist with new anchor */
                uint8_t* anchor_s = msa->sequences[best_idx]->s;
                int anchor_len = msa->sequences[best_idx]->len;
                if(msa->pool && numseq >= KALIGN_DIST_MIN_SEQS){
                        struct anchor_dist_ctx ctx = { msa->sequences, min_dist, anchor_s, anchor_len, dist, k };
                        tp_parallel_for_chunked(msa->pool, 0, numseq, KALIGN_PFOR_MIN_CHUNK, anchor_update_fn, &ctx);
                }else{
                for(i = 0; i < numseq; i++){
                        if(min_dist[i] < 0.0f) continue;  /* already selected */
                        uint8_t* si = msa->sequences[i]->s;
                        int li = msa->sequences[i]->len;
                        float d;
                        if(li > anchor_len){
                                d = (float)BPM(si, anchor_s, li, anchor_len);
                        }else{
                                d = (float)BPM(anchor_s, si, anchor_len, li);
                        }
                        if(dist) dist[i][k] = d;
                        if(d < min_dist[i]){
                                min_dist[i] = d;
                        }
                }
                }
        }

        MFREE(min_dist);
        *n_out = K;
        return anchors;
ERROR:
        if(anchors) MFREE(anchors);
        if(min_dist) MFREE(min_dist);
        return NULL;
}

int* pick_anchor(struct msa* msa, int* n)
{ KPERF_REGION(KP_TREE_PICK_ANCHOR);
        ASSERT(msa != NULL, "No alignment.");
        return pick_anchor_farthest_first(msa, DEFAULT_NUM_ANCHORS, n, NULL);
ERROR:
        return NULL;
}

int* pick_anchor_dist(struct msa* msa, int* n, float*** dist_out)
{ KPERF_REGION(KP_TREE_PICK_ANCHOR);
        float** dist = NULL;
        int* anchors = NULL;
        int i, k;
        int K = DEFAULT_NUM_ANCHORS;
        ASSERT(msa != NULL, "No alignment.");
        if(K > msa->numseq) K = msa->numseq;
        if(K < 1) K = 1;
        MMALLOC(dist, sizeof(float*) * msa->numseq);
        for(i = 0; i < msa->numseq; i++){
                dist[i] = NULL;
        }
        for(i = 0; i < msa->numseq; i++){
                MMALLOC(dist[i], sizeof(float) * K);
                for(k = 0; k < K; k++){
                        dist[i][k] = -1.0F;
                }
        }
        RUNP(anchors = pick_anchor_farthest_first(msa, K, n, dist));
        *dist_out = dist;
        return anchors;
ERROR:
        if(dist){
                for(i = 0; i < msa->numseq; i++){
                        if(dist[i]) MFREE(dist[i]);
                }
                MFREE(dist);
        }
        return NULL;
}

void pick_anchor_dist_free(float** dist, int numseq)
{
        int i;
        if(!dist) return;
        for(i = 0; i < numseq; i++){
                MFREE(dist[i]);
        }
        MFREE(dist);
}

int* pick_anchor_n(struct msa* msa, int requested, int* n)
{ KPERF_REGION(KP_TREE_PICK_ANCHOR);
        ASSERT(msa != NULL, "No alignment.");
        return pick_anchor_farthest_first(msa, requested, n, NULL);
ERROR:
        return NULL;
}
