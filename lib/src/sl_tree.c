/* Guide-tree options of a run (struct tree_opts in msa_struct.h; run config fields tree_* in kalign_config.h).
   The defaults give the shipped tree (bisecting k-means on Kalign's anchor distances); nothing in this file runs
   unless a run selects another option.

   tree_method       0 k-means (shipped) | 1 prim | 2 slink | 3 hybrid
                     prim / slink: exact single linkage on all pairs (identical trees, see below)
                     hybrid: bisecting k-means down to groups below tree_leaf sequences, each group's subtree by
                     single linkage (Prim) instead of UPGMA
   tree_dist         0 raw (shipped) | 1 norm | 2 lcs_indel | 3 lcs_frac
                     raw  = Myers best-fit edit distance d + min(10000, (li+lj)/2)/10000 (Kalign's pair distance)
                     norm = d / pattern length
                     lcs_indel = (n + m - 2 L) / L, L = longest common subsequence (FAMSA's indel_div_lcs)
                     lcs_frac  = 1 - L / min(n, m)
   tree_alpha        0 reduced (shipped: the 13-letter tree alphabet) | 1 full (letters A-Z, upper-cased; other = 26)
   tree_maxlen       residues used per sequence: Myers pattern cap (default and maximum 1024); LCS: both sequences
                     cut to this length (default 1024; 0 = whole sequences)
   tree_kmeans_dist  0 shipped anchor distances | 1 k-means (anchor distances and UPGMA leaf distances) on the
                     distance chosen above.  The alignment's own use of anchor distances (gap scaling, VSM) is
                     unchanged.
   tree_leaf         groups below this size get an exact subtree (default KALIGN_KMEANS_UPGMA_THRESHOLD = 50;
                     hybrid default 2000)
   realign_tree      1 (default) the realign pass follows tree_leaf / hybrid | 0 it keeps the shipped leaf rule
   tree_sl_noise     hybrid only: sigma of the noise on the single-linkage subtrees' distances, keyed by the run's
                     tree_seed (default 0 = none; ignored when tree_seed is 0 and in the realign pass)

   Myers distance: Kalign's bit-parallel best-fit edit distance (as bpm_block: the shorter sequence, capped, fitted
   into the best-matching part of the longer one), made symmetric by taking the longer sequence as text (equal
   lengths: the higher index).  LCS: bit-parallel (Allison-Dix / Hyyro), exact, symmetric; the shorter sequence is
   the bit-vector.  Accurate mode's noisy trees: d * (1 + sigma * z), z standard normal drawn from (seed, i, j).

   Prim and SLINK compute the same single-linkage hierarchy; the binary tree is built canonically from the merge
   links (groups of equal height, components merged left to right by smallest member), so both give the identical
   tree even with tied distances.  Prim aborts a distance as soon as it can no longer beat the current best link of
   that sequence; rows are split over the thread pool. */

#include "tldevel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "msa_struct.h"
#include "sl_tree.h"
#include "kalign/kalign_config.h"

#include "threadpool.h"
#ifdef HAVE_AVX2
#include <immintrin.h>
#include <mm_malloc.h>
#elif defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#define SL_SIGMA 27
#define SL_W 64
#define SL_MAXM 1024

enum { D_RAW = 0, D_NORM, D_LCS_INDEL, D_LCS_FRAC };

/* The bit-parallel kernels below are compiled once per number of 64-bit blocks (1 .. SL_FIXED): with a constant
   block count the compiler unrolls the block loop and keeps the bit-vectors in registers. */
#define SL_FIXED (SL_MAXM / SL_W)
#if defined(__GNUC__) || defined(__clang__)
#define SL_ALWAYS_INLINE static inline __attribute__((always_inline))
#define SL_NOINLINE static __attribute__((noinline))
#else
#define SL_ALWAYS_INLINE static inline
#define SL_NOINLINE static
#endif

struct sl_seq {
        const uint8_t* s;
        uint8_t* own;   /* full-alphabet encoding (owned) or NULL */
        int n;          /* full length */
        /* Myers */
        int m;          /* pattern length (capped) */
        int b;          /* blocks */
        uint64_t* peq;  /* SL_SIGMA x b, padded */
        /* LCS */
        int ln;         /* length used (cut to maxlen) */
        int lb;         /* blocks */
        uint64_t* lpeq; /* SL_SIGMA x lb, match bits only */
        uint32_t have;  /* letters among the first ln (bit c = letter c) */
};

struct sl_ctx {
        struct sl_seq* q;
        int N;
        int dist;
        int full;
        int maxlen;
        float noise;
        uint64_t seed;
};

/* ---------------------------------------------------------------- settings */
/* Run config value if set (>= 0), else the default (cfg NULL: all defaults).  The defaults are the shipped tree;
   leaf defaults to 2000 for the hybrid tree and to KALIGN_KMEANS_UPGMA_THRESHOLD otherwise; the realign pass
   follows leaf / hybrid unless realign_tree = 0. */
int sl_opts_resolve(struct tree_opts* o, const struct kalign_run_config* cfg)
{
        o->mode = (cfg && cfg->tree_method >= 0) ? cfg->tree_method : SL_MODE_KMEANS;
        o->dist = (cfg && cfg->tree_dist >= 0) ? cfg->tree_dist : D_RAW;
        o->full = (cfg && cfg->tree_alpha >= 0) ? cfg->tree_alpha : 0;
        o->maxlen = (cfg && cfg->tree_maxlen >= 0) ? cfg->tree_maxlen : SL_MAXM;
        o->leaf = (cfg && cfg->tree_leaf >= 0) ? cfg->tree_leaf
                : (o->mode == SL_MODE_HYBRID ? 2000 : KALIGN_KMEANS_UPGMA_THRESHOLD);
        o->kmeans_dist = (cfg && cfg->tree_kmeans_dist >= 0) ? cfg->tree_kmeans_dist : 0;
        o->realign_tree = (cfg && cfg->realign_tree >= 0) ? cfg->realign_tree : 1;
        o->in_realign = 0;
        o->sl_seed = cfg ? cfg->tree_seed : 0;
        o->sl_noise = (cfg && o->sl_seed != 0) ? cfg->tree_sl_noise : 0.0f;
        ASSERT(o->sl_noise >= 0.0f, "tree_sl_noise %f out of range", o->sl_noise);
        ASSERT(o->mode >= SL_MODE_KMEANS && o->mode <= SL_MODE_HYBRID, "tree_method %d out of range", o->mode);
        ASSERT(o->dist >= D_RAW && o->dist <= D_LCS_FRAC, "tree_dist %d out of range", o->dist);
        ASSERT(o->full == 0 || o->full == 1, "tree_alpha %d out of range", o->full);
        ASSERT(o->kmeans_dist == 0 || o->kmeans_dist == 1, "tree_kmeans_dist %d out of range", o->kmeans_dist);
        ASSERT(o->realign_tree == 0 || o->realign_tree == 1, "realign_tree %d out of range", o->realign_tree);
        return OK;
ERROR:
        return FAIL;
}

int sl_tree_mode(const struct msa* msa)
{
        return msa->tree.mode;
}

int sl_kmeans_dist_override(const struct msa* msa)
{
        return msa->tree.kmeans_dist;
}

/* the realign pass keeps the shipped leaf rule when realign_tree = 0 */
static int sl_follow(const struct msa* msa)
{
        return !msa->tree.in_realign || msa->tree.realign_tree;
}

int sl_kmeans_leaf(const struct msa* msa)
{
        return sl_follow(msa) ? msa->tree.leaf : KALIGN_KMEANS_UPGMA_THRESHOLD;
}

int sl_kmeans_leaf_is_sl(const struct msa* msa)
{
        return sl_follow(msa) && msa->tree.mode == SL_MODE_HYBRID;
}

/* ---------------------------------------------------------------- per-sequence tables */
static int sl_prep(struct sl_ctx* c, struct msa* msa, const int* idx, int N, float noise, uint64_t seed)
{
        c->N = N;
        c->dist = msa->tree.dist;
        c->full = msa->tree.full;
        c->maxlen = msa->tree.maxlen;
        c->noise = noise;
        c->seed = seed;
        c->q = NULL;
        MMALLOC(c->q, sizeof(struct sl_seq) * (size_t)N);
        for(int k = 0; k < N; k++){
                c->q[k].own = NULL; c->q[k].peq = NULL; c->q[k].lpeq = NULL;
        }
        for(int k = 0; k < N; k++){
                int i = idx ? idx[k] : k;
                struct sl_seq* x = &c->q[k];
                x->s = msa->sequences[i]->s;
                x->n = msa->sequences[i]->len;
                if(c->full){
                        const char* txt = msa->sequences[i]->seq;
                        x->own = malloc((size_t)(x->n > 0 ? x->n : 1));
                        if(!x->own) ERROR_MSG("SL: allocation failed");
                        for(int p = 0; p < x->n; p++){
                                int ch = (unsigned char)txt[p];
                                if(ch >= 'a' && ch <= 'z') ch -= 32;
                                x->own[p] = (ch >= 'A' && ch <= 'Z') ? (uint8_t)(ch - 'A') : 26;
                        }
                        x->s = x->own;
                }
                if(c->dist == D_RAW || c->dist == D_NORM){
                        int cap = (c->maxlen > 0 && c->maxlen < SL_MAXM) ? c->maxlen : SL_MAXM;
                        x->m = x->n > cap ? cap : x->n;
                        x->b = (x->m + SL_W - 1) / SL_W;
                        if(x->b < 1) x->b = 1;
                        x->peq = calloc((size_t)SL_SIGMA * (size_t)x->b, sizeof(uint64_t));
                        if(!x->peq) ERROR_MSG("SL: Peq allocation failed");
                        for(int p = 0; p < x->m; p++){
                                if(x->s[p] < SL_SIGMA){
                                        x->peq[x->s[p] * x->b + p / SL_W] |= (uint64_t)1 << (p % SL_W);
                                }
                        }
                        for(int p = x->m; p < x->b * SL_W; p++){
                                for(int cc = 0; cc < SL_SIGMA; cc++){
                                        x->peq[cc * x->b + p / SL_W] |= (uint64_t)1 << (p % SL_W);
                                }
                        }
                }else{
                        x->ln = (c->maxlen > 0 && x->n > c->maxlen) ? c->maxlen : x->n;
                        x->lb = (x->ln + SL_W - 1) / SL_W;
                        if(x->lb < 1) x->lb = 1;
                        x->lpeq = calloc((size_t)SL_SIGMA * (size_t)x->lb, sizeof(uint64_t));
                        if(!x->lpeq) ERROR_MSG("SL: LCS table allocation failed");
                        x->have = 0;
                        for(int p = 0; p < x->ln; p++){
                                if(x->s[p] < SL_SIGMA){
                                        x->lpeq[x->s[p] * x->lb + p / SL_W] |= (uint64_t)1 << (p % SL_W);
                                        x->have |= (uint32_t)1 << x->s[p];
                                }
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}

static void sl_free(struct sl_ctx* c)
{
        if(c->q){
                for(int i = 0; i < c->N; i++){
                        free(c->q[i].peq);
                        free(c->q[i].lpeq);
                        free(c->q[i].own);
                }
                MFREE(c->q);
        }
}

/* ---------------------------------------------------------------- Myers best-fit (= scalar bpm_block) */
/* Text s (n letters) against the pattern tables peq (m letters in B blocks, padded with positions that match
   every letter); the text is followed by W = 64 B - m pad letters.  Returns the smallest score of the last
   pattern position over the text, or, once that minimum cannot fall below cut any more, the minimum so far
   (>= cut).  Only the last block's score is used, so only that one is kept. */
SL_ALWAYS_INLINE int sl_bpm_rows(const uint8_t* s, const int n, const uint64_t* peq, const int B, const int m,
                                 const int cut)
{
        const int w = SL_W;
        const int W = w * B - m;
        const int L = n + W;
        uint64_t P[SL_FIXED], M[SL_FIXED];
        int32_t score = B * w;
        int k = m;
        for(int b = 0; b < B; b++){
                P[b] = ~(uint64_t)0;
                M[b] = 0;
        }
        for(int i = 0; i < L; i++){
                const uint8_t c = (i >= n) ? 0 : s[i];
                const int known = c < SL_SIGMA;
                const uint64_t* eqc = peq + (size_t)(known ? c : 0) * (size_t)B;
                uint64_t hp = 0, hn = 0;
                for(int b = 0; b < B; b++){
                        const uint64_t Pv = P[b];
                        const uint64_t Mv = M[b];
                        uint64_t Eq = known ? eqc[b] : 0;
                        const uint64_t Xv = Eq | Mv;
                        Eq |= hn;
                        const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
                        uint64_t Ph = Mv | (~(Xh | Pv));
                        uint64_t Mh = Pv & Xh;
                        const uint64_t op = Ph >> (w - 1);
                        const uint64_t on = Mh >> (w - 1);
                        Ph = (Ph << 1) | hp;
                        Mh = (Mh << 1) | hn;
                        P[b] = Mh | (~(Xv | Ph));
                        M[b] = Ph & Xv;
                        hp = op;
                        hn = on;
                }
                score += (int32_t)hp - (int32_t)hn;
                if(score < k){
                        k = score;
                }
                if(k >= cut && score - (L - 1 - i) >= cut){
                        return k;
                }
        }
        return k;
}

#define SL_BPM_B(x) SL_NOINLINE int sl_bpm_b##x(const uint8_t* s, int n, const uint64_t* peq, int m, int cut) \
        { return sl_bpm_rows(s, n, peq, x, m, cut); }
SL_BPM_B(1) SL_BPM_B(2) SL_BPM_B(3) SL_BPM_B(4) SL_BPM_B(5) SL_BPM_B(6) SL_BPM_B(7) SL_BPM_B(8)
SL_BPM_B(9) SL_BPM_B(10) SL_BPM_B(11) SL_BPM_B(12) SL_BPM_B(13) SL_BPM_B(14) SL_BPM_B(15) SL_BPM_B(16)

typedef int (*sl_bpm_fn)(const uint8_t*, int, const uint64_t*, int, int);
static const sl_bpm_fn sl_bpm_fixed[SL_FIXED + 1] = {
        NULL, sl_bpm_b1, sl_bpm_b2, sl_bpm_b3, sl_bpm_b4, sl_bpm_b5, sl_bpm_b6, sl_bpm_b7, sl_bpm_b8,
        sl_bpm_b9, sl_bpm_b10, sl_bpm_b11, sl_bpm_b12, sl_bpm_b13, sl_bpm_b14, sl_bpm_b15, sl_bpm_b16
};

/* the pattern is capped at SL_MAXM letters (sl_prep), so p->b is 1 .. SL_FIXED */
static int sl_bpm(const struct sl_seq* t, const struct sl_seq* p, int cut)
{
        return sl_bpm_fixed[p->b](t->s, t->n, p->peq, p->m, cut);
}

/* Several pairs at once, one pair per 64-bit lane: AVX2 four lanes, NEON two, otherwise scalar only. */
#if defined(HAVE_AVX2)
#define SLV_N 4
typedef __m256i slv;
#define slv_lanes(r, b) _mm256_set_epi64x((long long)(r)[3][b], (long long)(r)[2][b], \
                                          (long long)(r)[1][b], (long long)(r)[0][b])
#define slv_set1(x) _mm256_set1_epi64x((long long)(x))
#define slv_load(p) _mm256_loadu_si256((const __m256i*)(p))
#define slv_store(p, v) _mm256_storeu_si256((__m256i*)(p), v)
#define slv_and(a, b) _mm256_and_si256(a, b)
#define slv_or(a, b) _mm256_or_si256(a, b)
#define slv_xor(a, b) _mm256_xor_si256(a, b)
#define slv_andnot(a, b) _mm256_andnot_si256(a, b)              /* ~a & b */
#define slv_add(a, b) _mm256_add_epi64(a, b)
#define slv_sub(a, b) _mm256_sub_epi64(a, b)
#define slv_srl63(a) _mm256_srli_epi64(a, 63)
#define slv_sll1(a) _mm256_slli_epi64(a, 1)
#define slv_gt(a, b) _mm256_cmpgt_epi64(a, b)                   /* signed a > b, all ones or zero */
#define slv_blend(a, b, mask) _mm256_blendv_epi8(a, b, mask)    /* mask ? b : a */
#define slv_none(mask) _mm256_testz_si256(mask, mask)
#elif defined(__aarch64__) && defined(__ARM_NEON)
#define SLV_N 2
typedef uint64x2_t slv;
#define slv_lanes(r, b) vcombine_u64(vcreate_u64((r)[0][b]), vcreate_u64((r)[1][b]))
#define slv_set1(x) vdupq_n_u64((uint64_t)(int64_t)(x))
#define slv_load(p) vld1q_u64((const uint64_t*)(p))
#define slv_store(p, v) vst1q_u64((uint64_t*)(p), v)
#define slv_and(a, b) vandq_u64(a, b)
#define slv_or(a, b) vorrq_u64(a, b)
#define slv_xor(a, b) veorq_u64(a, b)
#define slv_andnot(a, b) vbicq_u64(b, a)                        /* ~a & b */
#define slv_add(a, b) vaddq_u64(a, b)
#define slv_sub(a, b) vsubq_u64(a, b)
#define slv_srl63(a) vshrq_n_u64(a, 63)
#define slv_sll1(a) vshlq_n_u64(a, 1)
#define slv_gt(a, b) vcgtq_s64(vreinterpretq_s64_u64(a), vreinterpretq_s64_u64(b))
#define slv_blend(a, b, mask) vbslq_u64(mask, b, a)             /* mask ? b : a */
#define slv_none(mask) ((vgetq_lane_u64(mask, 0) | vgetq_lane_u64(mask, 1)) == 0)
#endif

#ifdef SLV_N
/* one pair of the vector kernel: text, pattern tables and the cut of sl_bpm_rows */
struct sl_bpm_lane {
        const uint8_t* s;
        const uint64_t* peq;
        int n;
        int m;
        int cut;
};

static const uint64_t sl_zero_row[SL_FIXED] = { 0 };

/* sl_bpm_rows on SLV_N pairs whose patterns have the same number of blocks B.  Each lane runs the scalar
   recurrence on its own text and pattern; a lane stops updating its minimum after its own L = n + W columns.
   The loop ends early only when every lane is past its end or meets the scalar exit test (its minimum and
   its score minus the columns left are both >= its cut), so every lane returns either its exact value or,
   as the scalar function may, a value >= its cut when the exact value is >= its cut. */
SL_ALWAYS_INLINE void sl_bpmv_rows(const struct sl_bpm_lane* ln, const int B, int* kout)
{
        const int w = SL_W;
        int64_t Ll[SLV_N], cutl[SLV_N], kl[SLV_N];
        int Lmax = 0;
        slv P[SL_FIXED], M[SL_FIXED];
        const slv ones = slv_set1(-1);
        const slv zero = slv_set1(0);
        for(int l = 0; l < SLV_N; l++){
                const int L = ln[l].n + (w * B - ln[l].m);
                Ll[l] = L;
                cutl[l] = ln[l].cut;
                kl[l] = ln[l].m;
                if(L > Lmax) Lmax = L;
        }
        const slv vL = slv_load(Ll);
        const slv vcut = slv_load(cutl);
        slv k = slv_load(kl);
        slv score = slv_set1(B * w);
        for(int b = 0; b < B; b++){
                P[b] = ones;
                M[b] = zero;
        }
        for(int i = 0; i < Lmax; i++){
                const uint64_t* r[SLV_N];
                for(int l = 0; l < SLV_N; l++){
                        const uint8_t c = (i >= ln[l].n) ? 0 : ln[l].s[i];
                        r[l] = (c < SL_SIGMA) ? ln[l].peq + (size_t)c * (size_t)B : sl_zero_row;
                }
                slv hp = zero, hn = zero;
                for(int b = 0; b < B; b++){
                        const slv Pv = P[b];
                        const slv Mv = M[b];
                        slv Eq = slv_lanes(r, b);
                        const slv Xv = slv_or(Eq, Mv);
                        Eq = slv_or(Eq, hn);
                        const slv Xh = slv_or(slv_xor(slv_add(slv_and(Eq, Pv), Pv), Pv), Eq);
                        slv Ph = slv_or(Mv, slv_andnot(slv_or(Xh, Pv), ones));
                        slv Mh = slv_and(Pv, Xh);
                        const slv op = slv_srl63(Ph);
                        const slv on = slv_srl63(Mh);
                        Ph = slv_or(slv_sll1(Ph), hp);
                        Mh = slv_or(slv_sll1(Mh), hn);
                        P[b] = slv_or(Mh, slv_andnot(slv_or(Xv, Ph), ones));
                        M[b] = slv_and(Ph, Xv);
                        hp = op;
                        hn = on;
                }
                score = slv_sub(slv_add(score, hp), hn);
                {
                        const slv active = slv_gt(vL, slv_set1(i));             /* i < L */
                        const slv left = slv_sub(vL, slv_set1(i + 1));          /* L - 1 - i */
                        k = slv_blend(k, score, slv_and(active, slv_gt(k, score)));
                        /* not finished: still inside its text and (k < cut or score - left < cut) */
                        const slv open = slv_and(active, slv_or(slv_gt(vcut, k),
                                                                 slv_gt(vcut, slv_sub(score, left))));
                        if(slv_none(open)){
                                break;
                        }
                }
        }
        slv_store(kl, k);
        for(int l = 0; l < SLV_N; l++){
                kout[l] = (int)kl[l];
        }
}

#define SL_BPMV_B(x) SL_NOINLINE void sl_bpmv_b##x(const struct sl_bpm_lane* ln, int* kout) \
        { sl_bpmv_rows(ln, x, kout); }
SL_BPMV_B(1) SL_BPMV_B(2) SL_BPMV_B(3) SL_BPMV_B(4) SL_BPMV_B(5) SL_BPMV_B(6) SL_BPMV_B(7) SL_BPMV_B(8)
SL_BPMV_B(9) SL_BPMV_B(10) SL_BPMV_B(11) SL_BPMV_B(12) SL_BPMV_B(13) SL_BPMV_B(14) SL_BPMV_B(15) SL_BPMV_B(16)

typedef void (*sl_bpmv_fn)(const struct sl_bpm_lane*, int*);
static const sl_bpmv_fn sl_bpmv_fixed[SL_FIXED + 1] = {
        NULL, sl_bpmv_b1, sl_bpmv_b2, sl_bpmv_b3, sl_bpmv_b4, sl_bpmv_b5, sl_bpmv_b6, sl_bpmv_b7, sl_bpmv_b8,
        sl_bpmv_b9, sl_bpmv_b10, sl_bpmv_b11, sl_bpmv_b12, sl_bpmv_b13, sl_bpmv_b14, sl_bpmv_b15, sl_bpmv_b16
};
#endif

/* ---------------------------------------------------------------- LCS, bit-parallel */
static int lcs_zeros(const uint64_t* V, int B, int m)
{
        int z = 0;
        for(int b = 0; b < B - 1; b++){
                z += 64 - __builtin_popcountll(V[b]);
        }
        int rest = m - 64 * (B - 1);
        uint64_t mask = rest >= 64 ? ~(uint64_t)0 : (((uint64_t)1 << rest) - 1);
        z += __builtin_popcountll(~V[B - 1] & mask);
        return z;
}

/* zeros among the first j bits of V = LCS of the text read so far and the first j letters of p */
static inline int lcs_zeros_prefix(const uint64_t* V, int j)
{
        int z = 0;
        int b = 0;
        for(; j >= 64; j -= 64, b++){
                z += 64 - __builtin_popcountll(V[b]);
        }
        if(j > 0){
                z += __builtin_popcountll(~V[b] & (((uint64_t)1 << j) - 1));
        }
        return z;
}

/* LCS of a text s (n letters) and a bit-vector side (lpeq, m letters in B blocks).
   need > 0: the caller only uses the value if it is >= need, so any value < need may be returned as soon as the
   LCS cannot reach need.  After row i, with r letters of the text left, an alignment that crosses row i at
   column j scores at most L(i, j) + min(r, m - j), which is largest at j = m - r:
       LCS <= L(i, m - r) + r      (zeros of the first m - r bits; checked every SL_LCS_CHECK rows).
   The function is always inlined with a constant B (sl_lcs_b1 .. sl_lcs_b16 below): the compiler then unrolls
   the block loop and keeps the bit-vector in registers, which is 1.5 to 3 times faster than a loop over a
   bit-vector in memory.  V is copied before it is counted so that it never has to live in memory. */
#define SL_LCS_CHECK 64

SL_ALWAYS_INLINE int sl_lcs_rows(const uint8_t* s, const int n, const uint64_t* lpeq, const int B, const int m,
                                 const int need)
{
        uint64_t V[SL_FIXED];
        uint64_t T[SL_FIXED];
        for(int b = 0; b < B; b++) V[b] = ~(uint64_t)0;
        int i = 0;
        while(i < n){
                const int e = (n - i > SL_LCS_CHECK) ? i + SL_LCS_CHECK : n;
                for(; i < e; i++){
                        const uint8_t c = s[i];
                        if(c < SL_SIGMA){
                                const uint64_t* Mk = lpeq + (size_t)c * (size_t)B;
                                uint64_t carry = 0;
                                for(int b = 0; b < B; b++){
                                        const uint64_t old = V[b];
                                        const uint64_t mk = Mk[b];
                                        const uint64_t s1 = old + (old & mk);
                                        const uint64_t c1 = s1 < old;
                                        const uint64_t s2 = s1 + carry;
                                        const uint64_t c2 = s2 < s1;
                                        V[b] = s2 | (old & ~mk);
                                        carry = c1 | c2;
                                }
                        }
                }
                if(need > 0 && i < n){
                        const int r = n - i;
                        if(m - r > 0){
                                for(int b = 0; b < B; b++) T[b] = V[b];
                                if(lcs_zeros_prefix(T, m - r) + r < need){
                                        return 0;
                                }
                        }
                }
        }
        for(int b = 0; b < B; b++) T[b] = V[b];
        return lcs_zeros(T, B, m);
}

#define SL_LCS_B(k) SL_NOINLINE int sl_lcs_b##k(const uint8_t* s, int n, const uint64_t* lpeq, int m, int need) \
        { return sl_lcs_rows(s, n, lpeq, k, m, need); }
SL_LCS_B(1) SL_LCS_B(2) SL_LCS_B(3) SL_LCS_B(4) SL_LCS_B(5) SL_LCS_B(6) SL_LCS_B(7) SL_LCS_B(8)
SL_LCS_B(9) SL_LCS_B(10) SL_LCS_B(11) SL_LCS_B(12) SL_LCS_B(13) SL_LCS_B(14) SL_LCS_B(15) SL_LCS_B(16)

typedef int (*sl_lcs_fn)(const uint8_t*, int, const uint64_t*, int, int);
static const sl_lcs_fn sl_lcs_fixed[SL_FIXED + 1] = {
        NULL, sl_lcs_b1, sl_lcs_b2, sl_lcs_b3, sl_lcs_b4, sl_lcs_b5, sl_lcs_b6, sl_lcs_b7, sl_lcs_b8,
        sl_lcs_b9, sl_lcs_b10, sl_lcs_b11, sl_lcs_b12, sl_lcs_b13, sl_lcs_b14, sl_lcs_b15, sl_lcs_b16
};

/* any number of blocks (bit-vector sides above 1024 letters: tree_maxlen 0 or > 1024) */
static int sl_lcs_any(const uint8_t* s, const int n, const uint64_t* lpeq, const int B, const int m, const int need)
{
        uint64_t* V = malloc(sizeof(uint64_t) * (size_t)B);
        if(!V) return 0;
        for(int b = 0; b < B; b++) V[b] = ~(uint64_t)0;
        for(int i = 0; i < n; i++){
                const uint8_t c = s[i];
                if(c < SL_SIGMA){
                        const uint64_t* Mk = lpeq + (size_t)c * (size_t)B;
                        uint64_t carry = 0;
                        for(int b = 0; b < B; b++){
                                const uint64_t old = V[b];
                                const uint64_t mk = Mk[b];
                                const uint64_t s1 = old + (old & mk);
                                const uint64_t c1 = s1 < old;
                                const uint64_t s2 = s1 + carry;
                                const uint64_t c2 = s2 < s1;
                                V[b] = s2 | (old & ~mk);
                                carry = c1 | c2;
                        }
                }
                if(need > 0 && (i & (SL_LCS_CHECK - 1)) == SL_LCS_CHECK - 1){
                        const int r = n - 1 - i;
                        if(r > 0 && m - r > 0 && lcs_zeros_prefix(V, m - r) + r < need){
                                free(V);
                                return 0;
                        }
                }
        }
        int L = lcs_zeros(V, B, m);
        free(V);
        return L;
}

/* LCS of t (text) and p (bit-vector side); need: see sl_lcs_rows */
static int sl_lcs(const struct sl_seq* t, const struct sl_seq* p, int need)
{
        if(need > p->ln){
                return 0;       /* the LCS is at most the shorter length */
        }
        if(p->lb <= SL_FIXED){
                return sl_lcs_fixed[p->lb](t->s, t->ln, p->lpeq, p->ln, need);
        }
        return sl_lcs_any(t->s, t->ln, p->lpeq, p->lb, p->ln, need);
}

#ifdef HAVE_AVX2
/* Four LCS at once: one text against four bit-vector sides, one side per 64-bit lane.  tab holds the match
   masks of the four sides interleaved (letter c, block b at tab[c * B + b]); a side with fewer than B blocks
   has zero masks in the blocks it lacks, which leaves those blocks at all ones and passes the carry through.
   The carry out of old + u (u a subset of old) is the top bit of u | (old & ~s1); the carry out of adding
   the incoming carry is the top bit of s1 & ~s2.  Same recurrence as sl_lcs_rows, so the same LCS. */
SL_ALWAYS_INLINE void sl_lcs4_rows(const uint8_t* s, const int n, const __m256i* tab, const int B, uint64_t* out)
{
        __m256i V[SL_FIXED];
        const __m256i ones = _mm256_set1_epi64x(-1);
        for(int b = 0; b < B; b++) V[b] = ones;
        for(int i = 0; i < n; i++){
                const uint8_t c = s[i];
                if(c < SL_SIGMA){
                        const __m256i* Mk = tab + (size_t)c * (size_t)B;
                        __m256i carry = _mm256_setzero_si256();
                        for(int b = 0; b < B; b++){
                                const __m256i old = V[b];
                                const __m256i mk = _mm256_load_si256(Mk + b);
                                const __m256i u = _mm256_and_si256(old, mk);
                                const __m256i s1 = _mm256_add_epi64(old, u);
                                const __m256i c1 = _mm256_or_si256(u, _mm256_andnot_si256(s1, old));
                                const __m256i s2 = _mm256_add_epi64(s1, carry);
                                const __m256i c2 = _mm256_andnot_si256(s2, s1);
                                V[b] = _mm256_or_si256(s2, _mm256_andnot_si256(mk, old));
                                carry = _mm256_srli_epi64(_mm256_or_si256(c1, c2), 63);
                        }
                }
        }
        for(int b = 0; b < B; b++) _mm256_storeu_si256((__m256i*)(out + 4 * b), V[b]);
}

#define SL_LCS4_B(k) SL_NOINLINE void sl_lcs4_b##k(const uint8_t* s, int n, const __m256i* tab, uint64_t* out) \
        { sl_lcs4_rows(s, n, tab, k, out); }
SL_LCS4_B(1) SL_LCS4_B(2) SL_LCS4_B(3) SL_LCS4_B(4) SL_LCS4_B(5) SL_LCS4_B(6) SL_LCS4_B(7) SL_LCS4_B(8)
SL_LCS4_B(9) SL_LCS4_B(10) SL_LCS4_B(11) SL_LCS4_B(12) SL_LCS4_B(13) SL_LCS4_B(14) SL_LCS4_B(15) SL_LCS4_B(16)

typedef void (*sl_lcs4_fn)(const uint8_t*, int, const __m256i*, uint64_t*);
static const sl_lcs4_fn sl_lcs4_fixed[SL_FIXED + 1] = {
        NULL, sl_lcs4_b1, sl_lcs4_b2, sl_lcs4_b3, sl_lcs4_b4, sl_lcs4_b5, sl_lcs4_b6, sl_lcs4_b7, sl_lcs4_b8,
        sl_lcs4_b9, sl_lcs4_b10, sl_lcs4_b11, sl_lcs4_b12, sl_lcs4_b13, sl_lcs4_b14, sl_lcs4_b15, sl_lcs4_b16
};

/* LCS of text t with the four sides p[0..3] (each at most SL_FIXED blocks) */
static void sl_lcs4(const struct sl_seq* t, const struct sl_seq* const* p, int* L)
{
        __m256i tab[SL_SIGMA * SL_FIXED];
        uint64_t out[4 * SL_FIXED];
        uint64_t* w = (uint64_t*)tab;
        int B = 1;
        for(int k = 0; k < 4; k++){
                if(p[k]->lb > B) B = p[k]->lb;
        }
        for(int c = 0; c < SL_SIGMA; c++){
                if(!((t->have >> c) & 1)) continue;     /* rows of letters the text lacks are never read */
                for(int k = 0; k < 4; k++){
                        const int lb = p[k]->lb;
                        const uint64_t* src = p[k]->lpeq + (size_t)c * (size_t)lb;
                        uint64_t* dst = w + (size_t)c * (size_t)B * 4 + k;
                        int b = 0;
                        for(; b < lb; b++) dst[4 * b] = src[b];
                        for(; b < B; b++) dst[4 * b] = 0;
                }
        }
        sl_lcs4_fixed[B](t->s, t->ln, tab, out);
        for(int k = 0; k < 4; k++){
                uint64_t V[SL_FIXED];
                const int lb = p[k]->lb;
                for(int b = 0; b < lb; b++) V[b] = out[4 * b + k];
                L[k] = lcs_zeros(V, lb, p[k]->ln);
        }
}
#endif

/* exact LCS of sequence v with each of the cnt (<= SL_MANY) sequences us[]; the LCS does not depend on which
   sequence is the text, so with AVX2 v is the text of all pairs and the others go four to a vector, grouped by
   their number of blocks */
#define SL_MANY 128
static void sl_lcs_many(const struct sl_ctx* c, int v, const int* us, int cnt, int* L)
{
        const struct sl_seq* a = &c->q[v];
        int k = 0;
#ifdef HAVE_AVX2
        int ord[SL_MANY];
        int start[SL_FIXED + 2];
        int no = 0;
        /* counting sort of the positions by block count; sequences above SL_FIXED blocks stay scalar */
        for(int b = 0; b <= SL_FIXED + 1; b++) start[b] = 0;
        for(int x = 0; x < cnt; x++){
                const int lb = c->q[us[x]].lb;
                if(lb <= SL_FIXED) start[lb + 1]++;
        }
        for(int b = 1; b <= SL_FIXED + 1; b++) start[b] += start[b - 1];
        for(int x = 0; x < cnt; x++){
                const int lb = c->q[us[x]].lb;
                if(lb <= SL_FIXED){
                        ord[start[lb]++] = x;
                        no++;
                }else{
                        const struct sl_seq* b = &c->q[us[x]];
                        L[x] = (a->ln >= b->ln) ? sl_lcs(a, b, 0) : sl_lcs(b, a, 0);
                }
        }
        for(; k + 4 <= no; k += 4){
                const struct sl_seq* p[4] = { &c->q[us[ord[k]]], &c->q[us[ord[k + 1]]],
                                              &c->q[us[ord[k + 2]]], &c->q[us[ord[k + 3]]] };
                int l4[4];
                sl_lcs4(a, p, l4);
                L[ord[k]] = l4[0]; L[ord[k + 1]] = l4[1]; L[ord[k + 2]] = l4[2]; L[ord[k + 3]] = l4[3];
        }
        for(; k < no; k++){
                const struct sl_seq* b = &c->q[us[ord[k]]];
                L[ord[k]] = (a->ln >= b->ln) ? sl_lcs(a, b, 0) : sl_lcs(b, a, 0);
        }
#else
        for(; k < cnt; k++){
                const struct sl_seq* b = &c->q[us[k]];
                L[k] = (a->ln >= b->ln) ? sl_lcs(a, b, 0) : sl_lcs(b, a, 0);
        }
#endif
}

/* ---------------------------------------------------------------- pair distance */
static inline uint64_t sl_mix(uint64_t x)
{
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27; x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return x;
}

static double sl_z(uint64_t seed, int i, int j)
{
        uint64_t a = sl_mix(seed ^ sl_mix(((uint64_t)(uint32_t)i << 32) | (uint32_t)j));
        uint64_t b = sl_mix(a ^ 0x9e3779b97f4a7c15ULL);
        double u1 = ((a >> 11) + 0.5) * (1.0 / 9007199254740992.0);
        double u2 = ((b >> 11) + 0.5) * (1.0 / 9007199254740992.0);
        return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

/* distance of local sequences i, j (noise keyed by the original indices gi, gj);
   cut: values >= cut need not be exact (returned value is then >= cut) */
/* noise factor of the pair with original indices gi, gj (1 without noise) */
static inline double sl_pair_f(const struct sl_ctx* c, int gi, int gj)
{
        int lo = gi < gj ? gi : gj, hi = gi < gj ? gj : gi;
        double f = 1.0;
        if(c->noise > 0.0f){
                f = 1.0 + (double)c->noise * sl_z(c->seed, lo, hi);
                if(f < 0.01) f = 0.01;
        }
        return f;
}

/* LCS distance from the LCS L of two sequences with n >= m letters, times the noise factor f */
static inline double sl_lcs_d(const struct sl_ctx* c, double n, double m, int L, double f)
{
        double d;
        if(c->dist == D_LCS_INDEL){
                d = L > 0 ? (n + m - 2.0 * L) / (double)L : (n + m);
        }else{
                d = m > 0 ? 1.0 - (double)L / m : 1.0;
        }
        return d * f;
}

/* Myers distance of a pair: text t, pattern p, the integer cut for sl_bpm, and distance = (k * scale + add) * f */
struct sl_bpm_plan {
        const struct sl_seq* t;
        const struct sl_seq* p;
        double add;
        double scale;
        double f;
        int icut;
};

static inline void sl_bpm_plan_pair(const struct sl_ctx* c, const struct sl_seq* a, const struct sl_seq* b,
                                    int gi, int gj, double cut, struct sl_bpm_plan* pl)
{
        const double f = sl_pair_f(c, gi, gj);
        const double ccut = cut < 1e300 ? cut / f : INFINITY;
        const struct sl_seq *t, *p;
        if(a->n > b->n || (a->n == b->n && gi > gj)){ t = a; p = b; } else { t = b; p = a; }
        double add = 0.0, scale = 1.0;
        if(c->dist == D_RAW){
                int s = (a->n + b->n) / 2;
                add = (double)(s < 10000 ? s : 10000) / 10000.0;
        }else{
                scale = 1.0 / (double)(p->m > 0 ? p->m : 1);
        }
        int icut = INT32_MAX;
        if(ccut < 1e300){
                double x = (ccut - add) / scale;
                icut = x > 1e9 ? INT32_MAX : (int)ceil(x);
                if(icut < 0) icut = 0;
        }
        pl->t = t; pl->p = p; pl->add = add; pl->scale = scale; pl->f = f; pl->icut = icut;
}

static inline double sl_bpm_d(const struct sl_bpm_plan* pl, int k)
{
        return ((double)k * pl->scale + pl->add) * pl->f;
}

static double sl_dist(const struct sl_ctx* c, int i, int j, int gi, int gj, double cut)
{
        const struct sl_seq* a = &c->q[i];
        const struct sl_seq* b = &c->q[j];
        if(c->dist == D_RAW || c->dist == D_NORM){
                struct sl_bpm_plan pl;
                sl_bpm_plan_pair(c, a, b, gi, gj, cut, &pl);
                return sl_bpm_d(&pl, sl_bpm(pl.t, pl.p, pl.icut));
        }
        double f = sl_pair_f(c, gi, gj);
        double ccut = cut < 1e300 ? cut / f : INFINITY;
        /* LCS: the shorter sequence is the bit-vector */
        const struct sl_seq *t, *p;
        if(a->ln > b->ln || (a->ln == b->ln && gi > gj)){ t = a; p = b; } else { t = b; p = a; }
        const double n = t->ln, m = p->ln;
        int need = 0;
        if(ccut < 1e300){
                double x;
                if(c->dist == D_LCS_INDEL){
                        x = (n + m) / (ccut + 2.0);          /* d < cut  <=>  L > (n+m)/(cut+2) */
                }else{
                        x = m * (1.0 - ccut);                 /* 1 - L/m < cut  <=>  L > m (1 - cut) */
                }
                if(x < 0) x = -1;
                need = x > 1e9 ? INT32_MAX : (int)floor(x) + 1;
                if(need < 0) need = 0;
        }
        int L = sl_lcs(t, p, need);
        return sl_lcs_d(c, n, m, L, f);
}

/* LCS distances of v to the cnt (<= SL_MANY) sequences us[] (no cut: every value exact) */
static void sl_dist_many(const struct sl_ctx* c, int v, int gv, const int* us, const int* gid, int cnt, double* d)
{
        int L[SL_MANY];
        const struct sl_seq* a = &c->q[v];
        sl_lcs_many(c, v, us, cnt, L);
        for(int k = 0; k < cnt; k++){
                const struct sl_seq* b = &c->q[us[k]];
                const double n = a->ln > b->ln ? a->ln : b->ln;
                const double m = a->ln > b->ln ? b->ln : a->ln;
                d[k] = sl_lcs_d(c, n, m, L[k], sl_pair_f(c, gv, gid ? gid[us[k]] : us[k]));
        }
}

/* ---------------------------------------------------------------- canonical tree from merge links */
struct sl_link { double h; int a, b; };

static int cmp_link(const void* x, const void* y)
{
        const struct sl_link* p = x;
        const struct sl_link* q = y;
        if(p->h < q->h) return -1;
        if(p->h > q->h) return 1;
        return 0;
}

static int uf_find_sl(int* par, int x)
{
        while(par[x] != x){ par[x] = par[par[x]]; x = par[x]; }
        return x;
}

struct sl_member { int group; int minidx; int root; };

static int cmp_member(const void* x, const void* y)
{
        const struct sl_member* p = x;
        const struct sl_member* q = y;
        if(p->group != q->group) return p->group < q->group ? -1 : 1;
        return p->minidx < q->minidx ? -1 : (p->minidx > q->minidx);
}

static int cmp_int(const void* x, const void* y)
{
        int a = *(const int*)x, b = *(const int*)y;
        return (a > b) - (a < b);
}

static struct node* sl_new_node(struct node* l, struct node* r, int id)
{
        struct node* n = malloc(sizeof(struct node));
        if(!n) return NULL;
        n->left = l; n->right = r; n->id = id;
        return n;
}

/* links over local indices 0..N-1; leaf ids = idx[local] (or local if idx is NULL) */
static struct node* sl_build_tree(struct sl_link* links, int N, const int* idx)
{
        int* par = malloc(sizeof(int) * (size_t)N);
        int* mn = malloc(sizeof(int) * (size_t)N);
        struct node** nd = malloc(sizeof(struct node*) * (size_t)N);
        int* roots = malloc(sizeof(int) * (size_t)(2 * N));
        int* lpar = malloc(sizeof(int) * (size_t)(2 * N));
        struct sl_member* mem = malloc(sizeof(struct sl_member) * (size_t)(2 * N));
        for(int i = 0; i < N; i++){
                par[i] = i; mn[i] = i;
                nd[i] = sl_new_node(NULL, NULL, idx ? idx[i] : i);
        }
        qsort(links, (size_t)(N - 1), sizeof(struct sl_link), cmp_link);
        for(int g = 0; g < N - 1; ){
                int e = g;
                while(e < N - 1 && links[e].h == links[g].h) e++;
                int nr = 0;
                for(int k = g; k < e; k++){
                        roots[nr++] = uf_find_sl(par, links[k].a);
                        roots[nr++] = uf_find_sl(par, links[k].b);
                }
                qsort(roots, (size_t)nr, sizeof(int), cmp_int);
                int u = 0;
                for(int k = 0; k < nr; k++){
                        if(u == 0 || roots[u - 1] != roots[k]) roots[u++] = roots[k];
                }
                nr = u;
                for(int k = 0; k < nr; k++) lpar[k] = k;
                for(int k = g; k < e; k++){
                        int ra = uf_find_sl(par, links[k].a), rb = uf_find_sl(par, links[k].b);
                        int ia = (int)((int*)bsearch(&ra, roots, (size_t)nr, sizeof(int), cmp_int) - roots);
                        int ib = (int)((int*)bsearch(&rb, roots, (size_t)nr, sizeof(int), cmp_int) - roots);
                        int x = uf_find_sl(lpar, ia), y = uf_find_sl(lpar, ib);
                        if(x != y){ if(x < y) lpar[y] = x; else lpar[x] = y; }
                }
                for(int k = 0; k < nr; k++){
                        mem[k].group = uf_find_sl(lpar, k);
                        mem[k].minidx = mn[roots[k]];
                        mem[k].root = roots[k];
                }
                qsort(mem, (size_t)nr, sizeof(struct sl_member), cmp_member);
                for(int k = 0; k < nr; ){
                        int l = k;
                        while(l < nr && mem[l].group == mem[k].group) l++;
                        if(l - k > 1){
                                struct node* acc = nd[mem[k].root];
                                for(int t = k + 1; t < l; t++){
                                        acc = sl_new_node(acc, nd[mem[t].root], -1);
                                }
                                int keep = mem[k].root;
                                for(int t = k + 1; t < l; t++){
                                        par[mem[t].root] = keep;
                                }
                                nd[keep] = acc;
                                mn[keep] = mem[k].minidx;
                        }
                        k = l;
                }
                g = e;
        }
        struct node* root = nd[uf_find_sl(par, 0)];
        free(par); free(mn); free(nd); free(roots); free(lpar); free(mem);
        return root;
}

/* ---------------------------------------------------------------- Prim */
struct prim_ctx {
        const struct sl_ctx* c;
        const int* gid;
        const uint8_t* in;
        double* best;
        int* parent;
        int v;
};

#ifdef SLV_N
/* Myers distances of v to the cnt (<= SL_MANY) sequences us[], each with its own cut best[u]: pairs are grouped
   by the block count of their pattern and go SLV_N to a vector; a group's last pairs fill the free lanes with
   copies of its last pair, or run scalar if only one is left. */
static void prim_bpm_many(struct prim_ctx* x, const int* us, int cnt)
{
        const struct sl_ctx* c = x->c;
        struct sl_bpm_plan pl[SL_MANY];
        int ord[SL_MANY];
        int kk[SL_MANY];
        int start[SL_FIXED + 2];
        for(int b = 0; b <= SL_FIXED + 1; b++) start[b] = 0;
        for(int i = 0; i < cnt; i++){
                sl_bpm_plan_pair(c, &c->q[x->v], &c->q[us[i]], x->gid[x->v], x->gid[us[i]], x->best[us[i]], &pl[i]);
                start[pl[i].p->b + 1]++;
        }
        for(int b = 1; b <= SL_FIXED + 1; b++) start[b] += start[b - 1];
        for(int i = 0; i < cnt; i++){
                ord[start[pl[i].p->b]++] = i;
        }
        for(int g = 0; g < cnt; ){
                const int B = pl[ord[g]].p->b;
                int e = g;
                while(e < cnt && pl[ord[e]].p->b == B) e++;
                for(; g < e; g += SLV_N){
                        const int have = e - g < SLV_N ? e - g : SLV_N;
                        if(have == 1){
                                const struct sl_bpm_plan* q = &pl[ord[g]];
                                kk[ord[g]] = sl_bpm(q->t, q->p, q->icut);
                        }else{
                                struct sl_bpm_lane ln[SLV_N];
                                int k4[SLV_N];
                                for(int l = 0; l < SLV_N; l++){
                                        const struct sl_bpm_plan* q = &pl[ord[g + (l < have ? l : have - 1)]];
                                        ln[l].s = q->t->s;
                                        ln[l].n = q->t->n;
                                        ln[l].peq = q->p->peq;
                                        ln[l].m = q->p->m;
                                        ln[l].cut = q->icut;
                                }
                                sl_bpmv_fixed[B](ln, k4);
                                for(int l = 0; l < have; l++){
                                        kk[ord[g + l]] = k4[l];
                                }
                        }
                }
                g = e;
        }
        for(int i = 0; i < cnt; i++){
                const double d = sl_bpm_d(&pl[i], kk[i]);
                if(d < x->best[us[i]]){
                        x->best[us[i]] = d;
                        x->parent[us[i]] = x->v;
                }
        }
}
#endif

static void prim_chunk(int start, int end, void* arg)
{
        struct prim_ctx* x = arg;
#ifdef SLV_N
        if(x->c->dist == D_RAW || x->c->dist == D_NORM){
                int us[SL_MANY];
                int cnt = 0;
                for(int u = start; u <= end; u++){
                        if(u < end && !x->in[u]) us[cnt++] = u;
                        if(cnt == SL_MANY || (u == end && cnt > 0)){
                                prim_bpm_many(x, us, cnt);
                                cnt = 0;
                        }
                }
                return;
        }
#endif
#ifdef HAVE_AVX2
        /* LCS distances four at a time; all exact, so the comparison below sees the same or a larger value
           than the scalar path, which may stop early only when the distance cannot beat best[u] */
        if(x->c->dist == D_LCS_INDEL || x->c->dist == D_LCS_FRAC){
                int us[SL_MANY];
                double d[SL_MANY];
                int cnt = 0;
                for(int u = start; u <= end; u++){
                        if(u < end && !x->in[u]) us[cnt++] = u;
                        if(cnt == SL_MANY || (u == end && cnt > 0)){
                                sl_dist_many(x->c, x->v, x->gid[x->v], us, x->gid, cnt, d);
                                for(int k = 0; k < cnt; k++){
                                        if(d[k] < x->best[us[k]]){
                                                x->best[us[k]] = d[k];
                                                x->parent[us[k]] = x->v;
                                        }
                                }
                                cnt = 0;
                        }
                }
                return;
        }
#endif
        for(int u = start; u < end; u++){
                if(x->in[u]) continue;
                double d = sl_dist(x->c, x->v, u, x->gid[x->v], x->gid[u], x->best[u]);
                if(d < x->best[u]){
                        x->best[u] = d;
                        x->parent[u] = x->v;
                }
        }
}

static int sl_prim(struct sl_ctx* c, const int* gid, struct sl_link* links, void* pool_v)
{
        int N = c->N;
        uint8_t* in = calloc((size_t)N, 1);
        double* best = malloc(sizeof(double) * (size_t)N);
        int* parent = malloc(sizeof(int) * (size_t)N);
        for(int i = 0; i < N; i++){ best[i] = INFINITY; parent[i] = -1; }
        struct prim_ctx x = { c, gid, in, best, parent, 0 };
        (void)pool_v;
        in[0] = 1;
        for(int step = 0; step < N - 1; step++){
                if(pool_v && N >= 512){
                        tp_parallel_for_chunked((threadpool_t*)pool_v, 0, N, 128, prim_chunk, &x);
                }else
                {
                        prim_chunk(0, N, &x);
                }
                double bm = INFINITY;
                int ba = -1;
                for(int u = 0; u < N; u++){
                        if(!in[u] && (best[u] < bm || ba < 0)){
                                bm = best[u]; ba = u;
                        }
                }
                in[ba] = 1;
                links[step].h = best[ba];
                links[step].a = parent[ba];
                links[step].b = ba;
                x.v = ba;
        }
        free(in); free(best); free(parent);
        return OK;
}

/* ---------------------------------------------------------------- SLINK (Sibson 1973) */
struct slink_row { const struct sl_ctx* c; const int* gid; double* M; int n; };

static void slink_row_fn(int start, int end, void* arg)
{
        struct slink_row* r = arg;
        for(int i = start; i < end; i++){
                r->M[i] = sl_dist(r->c, i, r->n, r->gid[i], r->gid[r->n], INFINITY);
        }
}

static int sl_slink(struct sl_ctx* c, const int* gid, struct sl_link* links, void* pool_v)
{
        int N = c->N;
        int* pi = malloc(sizeof(int) * (size_t)N);
        double* lambda = malloc(sizeof(double) * (size_t)N);
        double* M = malloc(sizeof(double) * (size_t)N);
        (void)pool_v;
        for(int n = 0; n < N; n++){
                pi[n] = n;
                lambda[n] = INFINITY;
                struct slink_row r = { c, gid, M, n };
                if(pool_v && n >= 512){
                        tp_parallel_for_chunked((threadpool_t*)pool_v, 0, n, 256, slink_row_fn, &r);
                }else
                {
                        slink_row_fn(0, n, &r);
                }
                for(int i = 0; i < n; i++){
                        if(lambda[i] >= M[i]){
                                if(lambda[i] < M[pi[i]]) M[pi[i]] = lambda[i];
                                lambda[i] = M[i];
                                pi[i] = n;
                        }else{
                                if(M[i] < M[pi[i]]) M[pi[i]] = M[i];
                        }
                }
                for(int i = 0; i < n; i++){
                        if(lambda[i] >= lambda[pi[i]]) pi[i] = n;
                }
        }
        for(int i = 0; i < N - 1; i++){
                links[i].h = lambda[i];
                links[i].a = i;
                links[i].b = pi[i];
        }
        free(pi); free(lambda); free(M);
        return OK;
}

/* ---------------------------------------------------------------- entry points */
int sl_build_subset(struct msa* msa, const int* idx, int N, int mode, float noise, uint64_t seed, void* pool,
                    struct node** ret)
{
        struct sl_ctx c;
        struct sl_link* links = NULL;
        int* gid = NULL;
        c.q = NULL;
        c.N = 0;
        if(N == 1){
                *ret = sl_new_node(NULL, NULL, idx ? idx[0] : 0);
                return OK;
        }
        ASSERT(N >= 2, "SL tree needs at least 2 sequences");
        RUN(sl_prep(&c, msa, idx, N, noise, seed));
        MMALLOC(gid, sizeof(int) * (size_t)N);
        for(int k = 0; k < N; k++) gid[k] = idx ? idx[k] : k;
        MMALLOC(links, sizeof(struct sl_link) * (size_t)(N - 1));
        if(mode == SL_MODE_SLINK){
                RUN(sl_slink(&c, gid, links, pool));
        }else{
                RUN(sl_prim(&c, gid, links, pool));
        }
        *ret = sl_build_tree(links, N, idx);
        if(!*ret) ERROR_MSG("SL tree construction failed");
        MFREE(links);
        MFREE(gid);
        sl_free(&c);
        return OK;
ERROR:
        if(links) MFREE(links);
        if(gid) MFREE(gid);
        sl_free(&c);
        return FAIL;
}

int sl_build(struct msa* msa, int mode, float noise, uint64_t seed, void* pool, struct node** ret)
{
        return sl_build_subset(msa, NULL, msa->numseq, mode, noise, seed, pool, ret);
}

/* anchor distance matrix for k-means on the tree distance: numseq rows, padded to a multiple of 8 columns,
   allocated like Kalign's anchor matrices (32-byte aligned under AVX2; freed by free_anchor_dm) */
struct anchor_ctx { const struct sl_ctx* c; const int* anchors; int K; float** dm; };

static void anchor_rows(int start, int end, void* arg)
{
        struct anchor_ctx* a = arg;
        const int lcs = a->c->dist == D_LCS_INDEL || a->c->dist == D_LCS_FRAC;
        for(int i = start; i < end; i++){
                if(lcs){
                        int us[SL_MANY], ks[SL_MANY];
                        double d[SL_MANY];
                        int cnt = 0;
                        for(int k = 0; k <= a->K; k++){
                                if(k < a->K){
                                        if(a->anchors[k] == i){
                                                a->dm[i][k] = 0.0f;
                                        }else{
                                                us[cnt] = a->anchors[k];
                                                ks[cnt++] = k;
                                        }
                                }
                                if(cnt == SL_MANY || (k == a->K && cnt > 0)){
                                        sl_dist_many(a->c, i, i, us, NULL, cnt, d);
                                        for(int x = 0; x < cnt; x++) a->dm[i][ks[x]] = (float)d[x];
                                        cnt = 0;
                                }
                        }
                        continue;
                }
                for(int k = 0; k < a->K; k++){
                        int j = a->anchors[k];
                        a->dm[i][k] = (i == j) ? 0.0f : (float)sl_dist(a->c, i, j, i, j, INFINITY);
                }
        }
}

float** sl_anchor_matrix(struct msa* msa, const int* anchors, int K)
{
        struct sl_ctx c;
        float** dm = NULL;
        int numseq = msa->numseq;
        int padded = ((K + 7) / 8) * 8;
        c.q = NULL;
        c.N = 0;
        RUN(sl_prep(&c, msa, NULL, numseq, 0.0f, 0));
        MMALLOC(dm, sizeof(float*) * (size_t)numseq);
        for(int i = 0; i < numseq; i++){
                dm[i] = NULL;
        }
        for(int i = 0; i < numseq; i++){
#ifdef HAVE_AVX2
                dm[i] = _mm_malloc(sizeof(float) * (size_t)padded, 32);
                ASSERT(dm[i] != NULL, "_mm_malloc failed");
#else
                MMALLOC(dm[i], sizeof(float) * (size_t)padded);
#endif
                for(int k = 0; k < padded; k++) dm[i][k] = 0.0f;
        }
        {
                struct anchor_ctx a = { &c, anchors, K, dm };
                if(msa->pool && numseq >= 512){
                        tp_parallel_for_chunked(msa->pool, 0, numseq, 64, anchor_rows, &a);
                }else
                {
                        anchor_rows(0, numseq, &a);
                }
        }
        sl_free(&c);
        return dm;
ERROR:
        sl_free(&c);
        return NULL;
}

/* pair matrix (galloc layout, freed with gfree) of the given samples on the tree distance */
float** sl_pair_matrix(struct msa* msa, int* samples, int n)
{
        struct sl_ctx c;
        float** dm = NULL;
        c.q = NULL;
        c.N = 0;
        RUN(sl_prep(&c, msa, samples, n, 0.0f, 0));
        RUN(galloc(&dm, n, n));
        for(int i = 0; i < n; i++){
                dm[i][i] = 0.0f;
                for(int j = 0; j < i; j++){
                        float d = (float)sl_dist(&c, i, j, samples[i], samples[j], INFINITY);
                        dm[i][j] = d;
                        dm[j][i] = d;
                }
        }
        sl_free(&c);
        return dm;
ERROR:
        sl_free(&c);
        return NULL;
}
