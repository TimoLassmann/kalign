#include "tldevel.h"
#include "msa_struct.h"

#ifdef HAVE_AVX2
#include <xmmintrin.h>
#include <mm_malloc.h>
#endif

#ifdef USE_THREADPOOL
#include "threadpool.h"
#endif

#define ALN_APAIR_DIST_IMPORT
#include "aln_apair_dist.h"

float pairwise_identity_dist(const char* a, const char* b, int alnlen);

#define REALIGN_NUM_ANCHORS 32

#ifdef USE_THREADPOOL
struct apair_ctx {
        struct msa_seq** seqs;
        float** dm;
        int n;
        int alnlen;
};

static void apair_row_fn(int row_start, int row_end, void *arg)
{
        struct apair_ctx *c = (struct apair_ctx *)arg;
        for (int i = row_start; i < row_end; i++) {
                const char *seq_i = c->seqs[i]->seq;
                for (int j = i + 1; j < c->n; j++) {
                        float d = pairwise_identity_dist(seq_i, c->seqs[j]->seq, c->alnlen);
                        c->dm[i][j] = d;
                        c->dm[j][i] = d;
                }
        }
}
#endif

int compute_aln_pairwise_dist(struct msa* msa, float*** dm_ptr)
{
        float** dm = NULL;
        int n;
        int i, j;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(msa->aligned == ALN_STATUS_FINAL, "MSA must be finalized");

        n = msa->numseq;

        MMALLOC(dm, sizeof(float*) * n);
        for(i = 0; i < n; i++){
                dm[i] = NULL;
        }
        for(i = 0; i < n; i++){
                MMALLOC(dm[i], sizeof(float) * n);
                dm[i][i] = 0.0f;
        }

#ifdef USE_THREADPOOL
        if(msa->pool && n >= KALIGN_DIST_MIN_SEQS){
                struct apair_ctx ctx = { msa->sequences, dm, n, msa->alnlen };
                tp_parallel_for_chunked(msa->pool, 0, n - 1, KALIGN_PFOR_MIN_CHUNK, apair_row_fn, &ctx);
        }else{
#endif
        for(i = 0; i < n - 1; i++){
                const char* seq_i = msa->sequences[i]->seq;
                for(j = i + 1; j < n; j++){
                        float d = pairwise_identity_dist(seq_i, msa->sequences[j]->seq,
                                                         msa->alnlen);
                        dm[i][j] = d;
                        dm[j][i] = d;
                }
        }
#ifdef USE_THREADPOOL
        }
#endif

        *dm_ptr = dm;
        return OK;
ERROR:
        if(dm){
                for(i = 0; i < n; i++){
                        if(dm[i]) MFREE(dm[i]);
                }
                MFREE(dm);
        }
        return FAIL;
}

void free_aln_dm(float** dm, int n)
{
        int i;
        if(dm == NULL) return;
        for(i = 0; i < n; i++){
                if(dm[i]) MFREE(dm[i]);
        }
        MFREE(dm);
}

/* Distance = 1.0 - identity.
   Only counts columns where both sequences have a residue (no gap). */
float pairwise_identity_dist(const char* a, const char* b, int alnlen) {
        int matches = 0;
        int aligned = 0;
        int i;

        for(i = 0; i < alnlen; i++){
                if(a[i] != '-' && b[i] != '-'){
                        aligned++;
                        if(a[i] == b[i]){
                                matches++;
                        }
                }
        }

        if(aligned == 0){
                return 1.0f;
        }
        return 1.0f - (float)matches / (float)aligned;
}

/* ======================================================================== */
/* Pairwise identity distance callback for bisecting k-means leaf clusters. */
/* ======================================================================== */

float** aln_identity_pair_dist(struct msa* msa, int* samples, int n)
{
        float** dm = NULL;
        int i, j;

        RUN(galloc(&dm, n, n));

        for(i = 0; i < n; i++){
                dm[i][i] = 0.0f;
                const char* seq_i = msa->sequences[samples[i]]->seq;
                for(j = i + 1; j < n; j++){
                        const char* seq_j = msa->sequences[samples[j]]->seq;
                        float d = pairwise_identity_dist(seq_i, seq_j, msa->alnlen);
                        /* Add small length-preference bonus matching d_estimation */
                        int avg_len = (msa->sequences[samples[i]]->len +
                                       msa->sequences[samples[j]]->len) / 2;
                        float add = (float)(avg_len < 10000 ? avg_len : 10000) / 10000.0f;
                        d += add;
                        dm[i][j] = d;
                        dm[j][i] = d;
                }
        }

        return dm;
ERROR:
        return NULL;
}

/* ======================================================================== */
/* Anchor selection and N×K distances from aligned sequences.               */
/* Used by the realign loop to replace the O(N²) pairwise + O(N³) UPGMA    */
/* with O(N×K) distances + parallel bisecting k-means.                      */
/* ======================================================================== */

#ifdef USE_THREADPOOL
struct anchor_aln_ctx {
        struct msa_seq** seqs;
        float* min_dist;
        int anchor_idx;
        int alnlen;
};

static void anchor_aln_init_fn(int start, int end, void* arg)
{
        struct anchor_aln_ctx* c = (struct anchor_aln_ctx*)arg;
        const char* anchor_seq = c->seqs[c->anchor_idx]->seq;
        for(int i = start; i < end; i++){
                c->min_dist[i] = pairwise_identity_dist(
                        c->seqs[i]->seq, anchor_seq, c->alnlen);
        }
}

static void anchor_aln_update_fn(int start, int end, void* arg)
{
        struct anchor_aln_ctx* c = (struct anchor_aln_ctx*)arg;
        const char* anchor_seq = c->seqs[c->anchor_idx]->seq;
        for(int i = start; i < end; i++){
                if(c->min_dist[i] < 0.0f) continue;
                float d = pairwise_identity_dist(
                        c->seqs[i]->seq, anchor_seq, c->alnlen);
                if(d < c->min_dist[i]){
                        c->min_dist[i] = d;
                }
        }
}

struct anchor_dm_ctx {
        struct msa_seq** seqs;
        int* anchors;
        int K;
        int alnlen;
        float** dm;
};

static void anchor_dm_row_fn(int start, int end, void* arg)
{
        struct anchor_dm_ctx* c = (struct anchor_dm_ctx*)arg;
        for(int i = start; i < end; i++){
                const char* seq_i = c->seqs[i]->seq;
                for(int k = 0; k < c->K; k++){
                        c->dm[i][k] = pairwise_identity_dist(
                                seq_i, c->seqs[c->anchors[k]]->seq, c->alnlen);
                }
        }
}
#endif

int pick_anchor_from_alignment(struct msa* msa, int K,
                                int** anchors_out, int* K_out)
{
        int numseq = msa->numseq;
        int* anchors = NULL;
        float* min_dist = NULL;
        int i, k;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(msa->aligned == ALN_STATUS_FINAL, "MSA must be finalized");

        if(K > numseq) K = numseq;
        if(K < 1) K = 1;

        MMALLOC(anchors, sizeof(int) * K);
        MMALLOC(min_dist, sizeof(float) * numseq);

        /* Pick first anchor: sequence closest to mean length */
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

        /* Initialize min_dist: identity distance from each seq to first anchor */
#ifdef USE_THREADPOOL
        if(msa->pool && numseq >= KALIGN_DIST_MIN_SEQS){
                struct anchor_aln_ctx ctx = { msa->sequences, min_dist,
                                               anchors[0], msa->alnlen };
                tp_parallel_for_chunked(msa->pool, 0, numseq,
                                         KALIGN_PFOR_MIN_CHUNK,
                                         anchor_aln_init_fn, &ctx);
        }else{
#endif
        for(i = 0; i < numseq; i++){
                min_dist[i] = pairwise_identity_dist(
                        msa->sequences[i]->seq,
                        msa->sequences[anchors[0]]->seq,
                        msa->alnlen);
        }
#ifdef USE_THREADPOOL
        }
#endif
        min_dist[anchors[0]] = -1.0f;

        /* Farthest-first: pick remaining K-1 anchors */
        for(k = 1; k < K; k++){
                float best_min = -1.0f;
                int best_idx = 0;
                for(i = 0; i < numseq; i++){
                        if(min_dist[i] > best_min){
                                best_min = min_dist[i];
                                best_idx = i;
                        }
                }
                anchors[k] = best_idx;
                min_dist[best_idx] = -1.0f;

                /* Update min_dist with new anchor */
#ifdef USE_THREADPOOL
                if(msa->pool && numseq >= KALIGN_DIST_MIN_SEQS){
                        struct anchor_aln_ctx ctx = { msa->sequences, min_dist,
                                                       best_idx, msa->alnlen };
                        tp_parallel_for_chunked(msa->pool, 0, numseq,
                                                 KALIGN_PFOR_MIN_CHUNK,
                                                 anchor_aln_update_fn, &ctx);
                }else{
#endif
                for(i = 0; i < numseq; i++){
                        if(min_dist[i] < 0.0f) continue;
                        float d = pairwise_identity_dist(
                                msa->sequences[i]->seq,
                                msa->sequences[best_idx]->seq,
                                msa->alnlen);
                        if(d < min_dist[i]){
                                min_dist[i] = d;
                        }
                }
#ifdef USE_THREADPOOL
                }
#endif
        }

        MFREE(min_dist);
        *anchors_out = anchors;
        *K_out = K;
        return OK;
ERROR:
        if(anchors) MFREE(anchors);
        if(min_dist) MFREE(min_dist);
        return FAIL;
}

int compute_aln_anchor_dist(struct msa* msa, int* anchors, int K,
                             float*** dm_out)
{
        float** dm = NULL;
        int numseq = msa->numseq;
        int i, j;
        int padded_K;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(msa->aligned == ALN_STATUS_FINAL, "MSA must be finalized");
        ASSERT(anchors != NULL, "No anchors");
        ASSERT(K > 0, "K must be > 0");

        /* Pad K to multiple of 8 for AVX2 alignment (matches d_estimation) */
        padded_K = K / 8;
        if(K % 8) padded_K++;
        padded_K <<= 3;

        MMALLOC(dm, sizeof(float*) * numseq);
        for(i = 0; i < numseq; i++){
                dm[i] = NULL;
        }
        for(i = 0; i < numseq; i++){
#ifdef HAVE_AVX2
                dm[i] = _mm_malloc(sizeof(float) * padded_K, 32);
                if(!dm[i]) goto ERROR;
#else
                MMALLOC(dm[i], sizeof(float) * padded_K);
#endif
                for(j = 0; j < padded_K; j++){
                        dm[i][j] = 0.0f;
                }
        }

        /* Fill N×K: identity distance from each seq to each anchor */
#ifdef USE_THREADPOOL
        if(msa->pool && numseq >= KALIGN_DIST_MIN_SEQS){
                struct anchor_dm_ctx ctx = { msa->sequences, anchors, K,
                                              msa->alnlen, dm };
                tp_parallel_for_chunked(msa->pool, 0, numseq,
                                         KALIGN_PFOR_MIN_CHUNK,
                                         anchor_dm_row_fn, &ctx);
        }else{
#endif
        for(i = 0; i < numseq; i++){
                const char* seq_i = msa->sequences[i]->seq;
                for(j = 0; j < K; j++){
                        dm[i][j] = pairwise_identity_dist(
                                seq_i, msa->sequences[anchors[j]]->seq,
                                msa->alnlen);
                }
        }
#ifdef USE_THREADPOOL
        }
#endif

        *dm_out = dm;
        return OK;
ERROR:
        if(dm){
                for(i = 0; i < numseq; i++){
                        if(dm[i]){
#ifdef HAVE_AVX2
                                _mm_free(dm[i]);
#else
                                MFREE(dm[i]);
#endif
                        }
                }
                MFREE(dm);
        }
        return FAIL;
}
