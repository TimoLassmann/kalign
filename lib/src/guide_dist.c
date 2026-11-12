#include <math.h>
#include <float.h>
#include <string.h>

#include "tldevel.h"

#include "kalign_simd.h"

#include "threadpool/threadpool.h"

#include "msa_struct.h"
#include "alphabet.h"
#include "aln_param.h"
#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"

#define GUIDE_DIST_IMPORT
#include "guide_dist.h"

int guide_dist_ctx_create(struct msa* msa, struct aln_param* ap,
                          struct guide_dist_ctx** out)
{ KPERF_PHASE_REGION(KP_GUIDE_CTX, KPH_TREE);
        struct guide_dist_ctx* c = NULL;
        struct alphabet* a = NULL;
        int i, j;

        ASSERT(msa != NULL, "No msa");
        ASSERT(ap != NULL, "No alignment parameters");

        MMALLOC(c, sizeof(struct guide_dist_ctx));
        c->ap = ap;
        c->s = NULL;
        c->len = NULL;
        c->anchor_d = NULL;
        c->anchor_col = NULL;
        c->num_anchors = 0;
        c->numseq = msa->numseq;
        c->biotype = msa->biotype;

        if(msa->biotype == ALN_BIOTYPE_DNA){
                RUNP(a = create_alphabet(ALPHA_defDNA));
        }else{
                RUNP(a = create_alphabet(ALPHA_ambigiousPROTEIN));
        }

        MMALLOC(c->s, sizeof(uint8_t*) * c->numseq);
        MMALLOC(c->len, sizeof(int) * c->numseq);
        for(i = 0; i < c->numseq; i++){
                c->s[i] = NULL;
        }
        for(i = 0; i < c->numseq; i++){
                struct msa_seq* seq = msa->sequences[i];
                int L = seq->len;
                c->len[i] = L;
                MMALLOC(c->s[i], sizeof(uint8_t) * (L + 1));
                for(j = 0; j < L; j++){
                        int8_t v = a->to_internal[(int) seq->seq[j] & 0x7F];
                        c->s[i][j] = (v < 0) ? 0 : (uint8_t) v;
                }
                c->s[i][L] = 0;
        }
        MFREE(a);
        *out = c;
        return OK;
ERROR:
        if(a){
                MFREE(a);
        }
        guide_dist_ctx_free(c);
        return FAIL;
}

void guide_dist_ctx_free(struct guide_dist_ctx* c)
{
        int i;
        if(!c){
                return;
        }
        if(c->s){
                for(i = 0; i < c->numseq; i++){
                        if(c->s[i]){
                                MFREE(c->s[i]);
                        }
                }
                MFREE(c->s);
        }
        if(c->len){
                MFREE(c->len);
        }
        if(c->anchor_d){
                MFREE(c->anchor_d);
        }
        if(c->anchor_col){
                MFREE(c->anchor_col);
        }
        MFREE(c);
}

static float corrected_distance(int aligned, int diff, int biotype)
{
        double p;
        double q;
        if(aligned <= 0){
                return GUIDE_DIST_MAX;
        }
        p = (double) diff / (double) aligned;
        if(biotype == ALN_BIOTYPE_DNA){
                /* Jukes-Cantor */
                q = 1.0 - (4.0 / 3.0) * p;
                if(q <= 0.01){
                        return GUIDE_DIST_MAX;
                }
                return (float)(-0.75 * log(q));
        }
        /* Kimura (1983) protein approximation */
        q = 1.0 - p - 0.2 * p * p;
        if(q <= 0.01){
                return GUIDE_DIST_MAX;
        }
        return (float)(-log(q));
}

/* Align sequences i and j with the caller's alignment memory m (reused
   across calls; init_alnmem resets everything the aligner reads). */
static int pair_dist_mem(const struct guide_dist_ctx* c, struct aln_mem* m,
                         int i, int j, float* d)
{
        const uint8_t* s1;
        const uint8_t* s2;
        int l1, l2;
        int pos_a, pos_b, k;
        int aligned = 0;
        int diff = 0;

        if(i == j){
                *d = 0.0f;
                return OK;
        }
        /* Kalign's seq-seq aligner requires seq1 = shorter sequence */
        if(c->len[i] <= c->len[j]){
                s1 = c->s[i]; l1 = c->len[i];
                s2 = c->s[j]; l2 = c->len[j];
        }else{
                s1 = c->s[j]; l1 = c->len[j];
                s2 = c->s[i]; l2 = c->len[i];
        }
        if(l1 == 0 || l2 == 0){
                *d = GUIDE_DIST_MAX;
                return OK;
        }

        m->ap = c->ap;
        m->mode = ALN_MODE_FULL;
        m->run_parallel = 0;
        m->consistency = NULL;
        m->prof1 = NULL;
        m->prof2 = NULL;
        m->seq1 = s1;
        m->seq2 = s2;
        m->len_a = l1;
        m->len_b = l2;
        RUN(init_alnmem(m));
        aln_runner(m);
        RUN(add_gap_info_to_path_n(m));

        /* path codes after add_gap_info_to_path_n:
           0 = aligned pair, bit 1 = gap in seq1, bit 2 = gap in seq2, 3 = end */
        pos_a = 0;
        pos_b = 0;
        k = 1;
        while(m->path[k] != 3){
                if(m->path[k] == 0){
                        if(pos_a < l1 && pos_b < l2){
                                aligned++;
                                if(s1[pos_a] != s2[pos_b]){
                                        diff++;
                                }
                        }
                        pos_a++;
                        pos_b++;
                }else if(m->path[k] & 1){
                        pos_b++;
                }else if(m->path[k] & 2){
                        pos_a++;
                }
                k++;
        }
        *d = corrected_distance(aligned, diff, c->biotype);
        return OK;
ERROR:
        return FAIL;
}

/* The anchor matrix entry holding the same ordered alignment as (i, j), if
   any.  Entry (r, k) was computed as pair(r, anchors[k]), i.e. seq1 = r when
   len[r] <= len[anchor], else seq1 = anchor.  Pair (i, j) is x = shorter
   (i on ties), y = the other.  (x, anchor y) always matches; (y, anchor x)
   matches only when len[x] < len[y] (on a tie it would put y first). */
static int cached_dist(const struct guide_dist_ctx* c, int i, int j, float* d)
{
        int x, y;
        if(!c->anchor_d || i == j){
                return 0;
        }
        if(c->len[i] <= c->len[j]){
                x = i; y = j;
        }else{
                x = j; y = i;
        }
        if(c->anchor_col[y] >= 0){
                *d = c->anchor_d[(size_t) x * c->num_anchors + c->anchor_col[y]];
                return 1;
        }
        if(c->anchor_col[x] >= 0 && c->len[x] < c->len[y]){
                *d = c->anchor_d[(size_t) y * c->num_anchors + c->anchor_col[x]];
                return 1;
        }
        return 0;
}

int guide_dist_pair(const struct guide_dist_ctx* c, int i, int j, float* d)
{
        struct aln_mem* m = NULL;
        RUN(alloc_aln_mem(&m, 256));
        RUN(pair_dist_mem(c, m, i, j, d));
        free_aln_mem(m);
        return OK;
ERROR:
        if(m){
                free_aln_mem(m);
        }
        return FAIL;
}

struct pair_mat_ctx {
        const struct guide_dist_ctx* c;
        const int* samples;
        float** dm;
        int n;
        int error;
};

/* Pairs (i < j) of the n x n matrix, enumerated row by row as
   p = 0 .. n(n-1)/2 - 1; a chunk [p_start, p_end) of that sequence. */
static void pair_mat_chunk(int p_start, int p_end, void* arg)
{
        struct pair_mat_ctx* a = (struct pair_mat_ctx*) arg;
        struct aln_mem* m = NULL;
        int i = 0;
        int j;
        int p = p_start;
        int row_len = a->n - 1;

        if(p_start >= p_end){
                return;
        }
        while(p >= row_len){
                p -= row_len;
                i++;
                row_len--;
        }
        j = i + 1 + p;
        if(alloc_aln_mem(&m, 256) != OK){
                a->error = 1;
                return;
        }
        for(p = p_start; p < p_end; p++){
                float d = 0.0f;
                int si = a->samples[i];
                int sj = a->samples[j];
                if(!cached_dist(a->c, si, sj, &d)){
                        if(pair_dist_mem(a->c, m, si, sj, &d) != OK){
                                a->error = 1;
                                break;
                        }
                }
                a->dm[i][j] = d;
                a->dm[j][i] = d;
                j++;
                if(j == a->n){
                        i++;
                        j = i + 1;
                }
        }
        free_aln_mem(m);
}

float** guide_dist_pair_matrix(struct msa* msa, int* samples, int n, void* ctx)
{ KPERF_REGION(KP_TREE_DIST_PAIRALN);
        const struct guide_dist_ctx* c = (const struct guide_dist_ctx*) ctx;
        struct pair_mat_ctx pctx;
        float** dm = NULL;
        int n_pairs;
        int i;

        ASSERT(c != NULL, "No guide distance context");
        RUN(galloc(&dm, n, n));
        for(i = 0; i < n; i++){
                dm[i][i] = 0.0f;
        }
        n_pairs = n * (n - 1) / 2;
        pctx.c = c;
        pctx.samples = samples;
        pctx.dm = dm;
        pctx.n = n;
        pctx.error = 0;
        if(msa->pool){
                tp_parallel_for_chunked(msa->pool, 0, n_pairs, KALIGN_PFOR_MIN_CHUNK,
                                        pair_mat_chunk, &pctx);
        }else{
                pair_mat_chunk(0, n_pairs, &pctx);
        }
        ASSERT(pctx.error == 0, "Pairwise guide-tree distance failed");
        return dm;
ERROR:
        if(dm){
                gfree(dm);
        }
        return NULL;
}

struct anchor_row_ctx {
        const struct guide_dist_ctx* c;
        float** dm;
        const int* anchors;
        int num_anchors;
        int error;
};

/* Entry (r, k) with r = anchors[p] and p > k is the same ordered alignment
   as entry (anchors[k], p) when the two lengths differ; it is copied from
   there after all rows are done instead of being computed twice. */
static int anchor_mirrored(const struct anchor_row_ctx* a, int r, int k)
{
        int p = a->c->anchor_col[r];
        int q = a->anchors[k];
        return p > k && a->c->len[r] != a->c->len[q] && a->anchors[p] == r;
}

static void anchor_rows(int row_start, int row_end, void* arg)
{
        struct anchor_row_ctx* a = (struct anchor_row_ctx*) arg;
        struct aln_mem* m = NULL;
        int r, k;
        if(alloc_aln_mem(&m, 256) != OK){
                a->error = 1;
                return;
        }
        for(r = row_start; r < row_end; r++){
                for(k = 0; k < a->num_anchors; k++){
                        float d = 0.0f;
                        if(anchor_mirrored(a, r, k)){
                                continue;
                        }
                        if(pair_dist_mem(a->c, m, r, a->anchors[k], &d) != OK){
                                a->error = 1;
                                free_aln_mem(m);
                                return;
                        }
                        a->dm[r][k] = d;
                }
        }
        free_aln_mem(m);
}

float** guide_dist_anchor_matrix(struct msa* msa, struct guide_dist_ctx* c,
                                 const int* anchors, int num_anchors)
{ KPERF_REGION(KP_TREE_DIST_PAIRALN);
        float** dm = NULL;
        struct anchor_row_ctx actx;
        int numseq;
        int padded;
        int i, j;

        ASSERT(c != NULL, "No guide distance context");
        numseq = msa->numseq;
        padded = num_anchors / 8;
        if(num_anchors % 8){
                padded++;
        }
        padded = padded << 3;

        MMALLOC(dm, sizeof(float*) * numseq);
        for(i = 0; i < numseq; i++){
                dm[i] = NULL;
        }
        for(i = 0; i < numseq; i++){
                dm[i] = kalign_aligned_malloc(sizeof(float) * padded);
                ASSERT(dm[i] != NULL, "aligned allocation failed");
                for(j = 0; j < padded; j++){
                        dm[i][j] = 0.0f;
                }
        }

        /* anchor column of each sequence (first occurrence) */
        if(c->anchor_col == NULL){
                MMALLOC(c->anchor_col, sizeof(int) * numseq);
        }
        for(i = 0; i < numseq; i++){
                c->anchor_col[i] = -1;
        }
        for(j = 0; j < num_anchors; j++){
                if(c->anchor_col[anchors[j]] < 0){
                        c->anchor_col[anchors[j]] = j;
                }
        }
        if(c->anchor_d){
                MFREE(c->anchor_d);
                c->anchor_d = NULL;
        }
        c->num_anchors = num_anchors;

        actx.c = c;
        actx.dm = dm;
        actx.anchors = anchors;
        actx.num_anchors = num_anchors;
        actx.error = 0;

        if(msa->pool && numseq >= KALIGN_DIST_MIN_SEQS){
                tp_parallel_for_chunked(msa->pool, 0, numseq, KALIGN_PFOR_MIN_CHUNK,
                                        anchor_rows, &actx);
        }else{
                anchor_rows(0, numseq, &actx);
        }
        ASSERT(actx.error == 0, "Pairwise guide-tree distance failed");

        /* copy the mirrored anchor-anchor entries */
        for(i = 0; i < numseq; i++){
                for(j = 0; j < num_anchors; j++){
                        if(anchor_mirrored(&actx, i, j)){
                                dm[i][j] = dm[anchors[j]][c->anchor_col[i]];
                        }
                }
        }

        /* keep an unperturbed copy for reuse by the leaf matrices */
        MMALLOC(c->anchor_d, sizeof(float) * (size_t) numseq * num_anchors);
        for(i = 0; i < numseq; i++){
                for(j = 0; j < num_anchors; j++){
                        c->anchor_d[(size_t) i * num_anchors + j] = dm[i][j];
                }
        }
        return dm;
ERROR:
        if(dm){
                for(i = 0; i < numseq; i++){
                        if(dm[i]){
                                kalign_aligned_free(dm[i]);
                        }
                }
                MFREE(dm);
        }
        return NULL;
}
