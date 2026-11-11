/* bpm_proto -- PROTOTYPE: exact faster variants of bpm_block (lib/src/bpm.c).
 *
 *   ref     lib's bpm_block (linked)
 *   simple  same arithmetic; the block-count adjustments removed (with
 *           maxd = m the last block's score never reaches maxd + w and the
 *           block count never grows past b_max - 1, so y is constant)
 *   x2      two independent (text, pattern) pairs advanced in one loop, so
 *           the two serial carry chains overlap
 *
 * Reads a FASTA file, encodes residues to 0..12 (any fixed map works for an
 * equality test), checks all three on every ordered pair (longer = text) and
 * times them.
 * usage: bpm_proto in.fasta [reps]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int bpm_block(const uint8_t* t, const uint8_t* p, int n, int m);

#define SIGMA 13
#define W64 64

struct bpm_pat {
        uint64_t peq[SIGMA][16];
        int b_max;
        int m;
        int pad;        /* W: padding characters after the text */
};

static inline void bpm_pat_init(struct bpm_pat* q, const uint8_t* p, int m)
{
        if(m > 1024) m = 1024;
        q->m = m;
        q->b_max = (m == 0) ? 1 : (m + W64 - 1) / W64;
        q->pad = W64 * q->b_max - m;
        memset(q->peq, 0, sizeof(q->peq));
        for(int i = 0; i < m; i++){
                if(p[i] < SIGMA) q->peq[p[i]][i / W64] |= (uint64_t)1 << (i % W64);
        }
        for(int i = m; i < q->b_max * W64; i++){
                for(int c = 0; c < SIGMA; c++) q->peq[c][i / W64] |= (uint64_t)1 << (i % W64);
        }
}

/* one block step, exactly the arithmetic of bpm_block */
#define BPM_STEP(Pv_, Mv_, Eq_, carry_, hout_)                                  \
        do {                                                                    \
                uint64_t Eq = (Eq_);                                            \
                uint64_t Xv = Eq | (Mv_);                                       \
                if((carry_) < 0) Eq |= 1;                                       \
                uint64_t Xh = (((Eq & (Pv_)) + (Pv_)) ^ (Pv_)) | Eq;            \
                uint64_t Ph = (Mv_) | ~(Xh | (Pv_));                            \
                uint64_t Mh = (Pv_) & Xh;                                       \
                int ho_ = 0;                                                     \
                if(Ph & ((uint64_t)1 << 63)) ho_ += 1;                           \
                if(Mh & ((uint64_t)1 << 63)) ho_ -= 1;                           \
                Ph <<= 1;                                                       \
                Mh <<= 1;                                                       \
                if((carry_) < 0) Mh |= 1;                                       \
                else if((carry_) > 0) Ph |= 1;                                  \
                (Pv_) = Mh | ~(Xv | Ph);                                        \
                (Mv_) = Ph & Xv;                                                \
                (hout_) = ho_;                                                   \
        } while(0)

static int bpm_simple(const uint8_t* t, const uint8_t* p, int n, int m)
{
        struct bpm_pat q;
        uint64_t P[16], M[16];
        int32_t score[16];
        bpm_pat_init(&q, p, m);
        const int B = q.b_max;
        int k = q.m;
        for(int b = 0; b < B; b++){
                P[b] = ~(uint64_t)0;
                M[b] = 0;
                score[b] = (b + 1) * W64;
        }
        for(int i = 0; i < n + q.pad; i++){
                const uint8_t c = (i >= n) ? 0 : t[i];
                int carry = 0;
                for(int b = 0; b < B; b++){
                        int ho;
                        BPM_STEP(P[b], M[b], q.peq[c][b], carry, ho);
                        carry = ho;
                        score[b] += carry;
                }
                if(score[B - 1] < k) k = score[B - 1];
        }
        return k;
}


/* Branch-free block step: the horizontal delta of the previous block is
   carried as two bits, hp (= +1) and hn (= -1), never both set (Ph and Mh
   are disjoint).  Same arithmetic as BPM_STEP, no data-dependent branches. */
#define BPM_STEP_BL(Pv_, Mv_, Eq_, hp_, hn_, sc_)                               \
        do {                                                                    \
                uint64_t Eq = (Eq_);                                            \
                uint64_t Xv = Eq | (Mv_);                                       \
                Eq |= (hn_);                                                    \
                uint64_t Xh = (((Eq & (Pv_)) + (Pv_)) ^ (Pv_)) | Eq;            \
                uint64_t Ph = (Mv_) | ~(Xh | (Pv_));                            \
                uint64_t Mh = (Pv_) & Xh;                                       \
                uint64_t op = Ph >> 63, on = Mh >> 63;                          \
                Ph = (Ph << 1) | (hp_);                                         \
                Mh = (Mh << 1) | (hn_);                                         \
                (Pv_) = Mh | ~(Xv | Ph);                                        \
                (Mv_) = Ph & Xv;                                                \
                (hp_) = op;                                                     \
                (hn_) = on;                                                     \
                (sc_) += (int32_t)op - (int32_t)on;                             \
        } while(0)

static int bpm_bl(const uint8_t* t, const uint8_t* p, int n, int m)
{
        struct bpm_pat q;
        uint64_t P[16], M[16];
        int32_t score[16];
        bpm_pat_init(&q, p, m);
        const int B = q.b_max;
        int k = q.m;
        for(int b = 0; b < B; b++){
                P[b] = ~(uint64_t)0;
                M[b] = 0;
                score[b] = (b + 1) * W64;
        }
        for(int i = 0; i < n + q.pad; i++){
                const uint8_t c = (i >= n) ? 0 : t[i];
                uint64_t hp = 0, hn = 0;
                for(int b = 0; b < B; b++){
                        BPM_STEP_BL(P[b], M[b], q.peq[c][b], hp, hn, score[b]);
                }
                k = score[B - 1] < k ? score[B - 1] : k;
        }
        return k;
}


/* Wavefront over (character, block): at step s, block b processes
   character s - b.  Its inputs are its own state after character s-b-1
   (step s-1) and block b-1's carry for character s-b (also step s-1), so the
   blocks of one step are independent and run in parallel on an out-of-order
   core.  Blocks are visited in descending order so that hp/hn[b-1] still hold
   the previous step's carry when block b reads them.  Per block, the
   arithmetic and its order are exactly bpm_block's. */
static int bpm_wf(const uint8_t* t, const uint8_t* p, int n, int m)
{
        struct bpm_pat q;
        uint64_t P[16], M[16], hp[17], hn[17];
        int32_t score[16];
        bpm_pat_init(&q, p, m);
        const int B = q.b_max;
        const int L = n + q.pad;
        int k = q.m;
        for(int b = 0; b < B; b++){
                P[b] = ~(uint64_t)0;
                M[b] = 0;
                score[b] = (b + 1) * W64;
                hp[b] = 0;
                hn[b] = 0;
        }
        /* hp/hn[b] = carry out of block b - 1 (index shifted: [0] is always 0) */
        uint64_t* cp = hp;
        uint64_t* cn = hn;
        cp[0] = 0; cn[0] = 0;
        for(int s = 0; s < L + B - 1; s++){
                int bhi = s < B - 1 ? s : B - 1;
                int blo = s - L + 1 > 0 ? s - L + 1 : 0;
                for(int b = bhi; b >= blo; b--){
                        const int i = s - b;
                        const uint8_t c = (i >= n) ? 0 : t[i];
                        uint64_t xp = b ? cp[b] : 0;
                        uint64_t xn = b ? cn[b] : 0;
                        BPM_STEP_BL(P[b], M[b], q.peq[c][b], xp, xn, score[b]);
                        cp[b + 1] = xp;
                        cn[b + 1] = xn;
                }
                if(bhi == B - 1 && s >= B - 1){
                        k = score[B - 1] < k ? score[B - 1] : k;
                }
        }
        return k;
}


/* ---- SIMD block wavefront: lane b = block b, character s - b ---- */
#if defined(__ARM_NEON)
#include <arm_neon.h>
#define BPM_SIMD 1
typedef struct { uint64x2_t lo, hi; } u64x4;
static inline u64x4 u4_set(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
{
        u64x4 r;
        r.lo = vcombine_u64(vcreate_u64(a), vcreate_u64(b));
        r.hi = vcombine_u64(vcreate_u64(c), vcreate_u64(d));
        return r;
}
static inline u64x4 u4_dup(uint64_t a) { u64x4 r = { vdupq_n_u64(a), vdupq_n_u64(a) }; return r; }
#define U4_BIN(name, op) static inline u64x4 name(u64x4 a, u64x4 b) { u64x4 r = { op(a.lo, b.lo), op(a.hi, b.hi) }; return r; }
U4_BIN(u4_and, vandq_u64)
U4_BIN(u4_or, vorrq_u64)
U4_BIN(u4_xor, veorq_u64)
U4_BIN(u4_add, vaddq_u64)
U4_BIN(u4_sub, vsubq_u64)
static inline u64x4 u4_not(u64x4 a) { return u4_xor(a, u4_dup(~(uint64_t)0)); }
static inline u64x4 u4_shl1(u64x4 a) { u64x4 r = { vshlq_n_u64(a.lo, 1), vshlq_n_u64(a.hi, 1) }; return r; }
static inline u64x4 u4_shr63(u64x4 a) { u64x4 r = { vshrq_n_u64(a.lo, 63), vshrq_n_u64(a.hi, 63) }; return r; }
/* [prev3, cur0, cur1, cur2] */
static inline u64x4 u4_shift_in(u64x4 prev, u64x4 cur)
{
        u64x4 r = { vextq_u64(prev.hi, cur.lo, 1), vextq_u64(cur.lo, cur.hi, 1) };
        return r;
}
/* lanes with mask all-ones take a, others b */
static inline u64x4 u4_sel(u64x4 mask, u64x4 a, u64x4 b)
{
        u64x4 r = { vbslq_u64(mask.lo, a.lo, b.lo), vbslq_u64(mask.hi, a.hi, b.hi) };
        return r;
}
static inline uint64_t u4_lane(u64x4 a, int l)
{
        switch(l){
        case 0: return vgetq_lane_u64(a.lo, 0);
        case 1: return vgetq_lane_u64(a.lo, 1);
        case 2: return vgetq_lane_u64(a.hi, 0);
        default: return vgetq_lane_u64(a.hi, 1);
        }
}
#elif defined(__AVX2__)
#include <immintrin.h>
#define BPM_SIMD 1
typedef __m256i u64x4;
static inline u64x4 u4_set(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
{
        return _mm256_set_epi64x((long long)d, (long long)c, (long long)b, (long long)a);
}
static inline u64x4 u4_dup(uint64_t a) { return _mm256_set1_epi64x((long long)a); }
static inline u64x4 u4_and(u64x4 a, u64x4 b) { return _mm256_and_si256(a, b); }
static inline u64x4 u4_or(u64x4 a, u64x4 b) { return _mm256_or_si256(a, b); }
static inline u64x4 u4_xor(u64x4 a, u64x4 b) { return _mm256_xor_si256(a, b); }
static inline u64x4 u4_add(u64x4 a, u64x4 b) { return _mm256_add_epi64(a, b); }
static inline u64x4 u4_sub(u64x4 a, u64x4 b) { return _mm256_sub_epi64(a, b); }
static inline u64x4 u4_not(u64x4 a) { return _mm256_xor_si256(a, _mm256_set1_epi64x(-1)); }
static inline u64x4 u4_shl1(u64x4 a) { return _mm256_slli_epi64(a, 1); }
static inline u64x4 u4_shr63(u64x4 a) { return _mm256_srli_epi64(a, 63); }
static inline u64x4 u4_shift_in(u64x4 prev, u64x4 cur)
{
        const u64x4 c = _mm256_permute4x64_epi64(cur, _MM_SHUFFLE(2, 1, 0, 3));
        const u64x4 p = _mm256_permute4x64_epi64(prev, _MM_SHUFFLE(2, 1, 0, 3));
        return _mm256_blend_epi32(c, p, 0x03);
}
static inline u64x4 u4_sel(u64x4 mask, u64x4 a, u64x4 b) { return _mm256_blendv_epi8(b, a, mask); }
static inline uint64_t u4_lane(u64x4 a, int l)
{
        uint64_t t[4];
        _mm256_storeu_si256((__m256i*)t, a);
        return t[l];
}
#else
#define BPM_SIMD 0
#endif

#if BPM_SIMD
static int bpm_simd(const uint8_t* t, const uint8_t* p, int n, int m)
{
        struct bpm_pat q;
        bpm_pat_init(&q, p, m);
        const int B = q.b_max;
        const int L = n + q.pad;
        const int V = (B + 3) / 4;
        int k = q.m;
        /* text padded by 16 zero characters on both sides: char(i) = tp[i] */
        uint8_t tbuf[2048 + 64];
        uint8_t* tp;
        uint8_t* heap = NULL;
        if(L + 32 <= (int)sizeof(tbuf)){
                tp = tbuf + 16;
        }else{
                heap = malloc((size_t)L + 32);
                tp = heap + 16;
        }
        memset(tp - 16, 0, 16);
        memcpy(tp, t, (size_t)n);
        memset(tp + n, 0, (size_t)(L - n) + 16);

        u64x4 P[4], M[4], SC[4], HP[4], HN[4];
        for(int v = 0; v < V; v++){
                P[v] = u4_dup(~(uint64_t)0);
                M[v] = u4_dup(0);
                SC[v] = u4_set(64 * (4 * v + 1), 64 * (4 * v + 2), 64 * (4 * v + 3), 64 * (4 * v + 4));
                HP[v] = u4_dup(0);
                HN[v] = u4_dup(0);
        }
        const int lv = (B - 1) / 4, ll = (B - 1) % 4;
        const u64x4 zero = u4_dup(0);
        for(int s = 0; s < L + B - 1; s++){
                u64x4 prevp = zero, prevn = zero;       /* carry into lane 0: none */
                for(int v = 0; v < V; v++){
                        const int b0 = 4 * v;
                        /* carry from lane b-1 at the previous step */
                        const u64x4 hp = u4_shift_in(prevp, HP[v]);
                        const u64x4 hn = u4_shift_in(prevn, HN[v]);
                        prevp = HP[v];
                        prevn = HN[v];
                        u64x4 Eq = u4_set(q.peq[tp[s - b0]][b0], q.peq[tp[s - b0 - 1]][b0 + 1],
                                          q.peq[tp[s - b0 - 2]][b0 + 2], q.peq[tp[s - b0 - 3]][b0 + 3]);
                        const u64x4 Pv = P[v], Mv = M[v];
                        const u64x4 Xv = u4_or(Eq, Mv);
                        Eq = u4_or(Eq, hn);
                        const u64x4 Xh = u4_or(u4_xor(u4_add(u4_and(Eq, Pv), Pv), Pv), Eq);
                        u64x4 Ph = u4_or(Mv, u4_not(u4_or(Xh, Pv)));
                        u64x4 Mh = u4_and(Pv, Xh);
                        const u64x4 op = u4_shr63(Ph);
                        const u64x4 on = u4_shr63(Mh);
                        Ph = u4_or(u4_shl1(Ph), hp);
                        Mh = u4_or(u4_shl1(Mh), hn);
                        u64x4 nP = u4_or(Mh, u4_not(u4_or(Xv, Ph)));
                        u64x4 nM = u4_and(Ph, Xv);
                        u64x4 nS = u4_sub(u4_add(SC[v], op), on);
                        if(s < b0 + 3){
                                /* prologue: lanes b > s have not started; keep their state */
                                const u64x4 act = u4_set(s >= b0 ? ~(uint64_t)0 : 0, s >= b0 + 1 ? ~(uint64_t)0 : 0,
                                                         s >= b0 + 2 ? ~(uint64_t)0 : 0, s >= b0 + 3 ? ~(uint64_t)0 : 0);
                                nP = u4_sel(act, nP, Pv);
                                nM = u4_sel(act, nM, Mv);
                                nS = u4_sel(act, nS, SC[v]);
                        }
                        P[v] = nP;
                        M[v] = nM;
                        SC[v] = nS;
                        HP[v] = op;
                        HN[v] = on;
                }
                if(s >= B - 1){
                        const int sc = (int)(int64_t)u4_lane(SC[lv], ll);
                        k = sc < k ? sc : k;
                }
        }
        free(heap);
        return k;
}
#endif


/* ---- two computations per loop (independent chains overlap) ---- */
#if BPM_SIMD
struct wstate {
        u64x4 P[4], M[4], SC[4], HP[4], HN[4];
        uint64_t peq[SIGMA][16];
        uint8_t tb[2048 + 64];
        const uint8_t* tp;
        int B, L, V, k, lv, ll, end;
};

static void ws_init(struct wstate* w, const uint8_t* t, const uint8_t* p, int n, int m)
{
        struct bpm_pat q;
        bpm_pat_init(&q, p, m);
        memcpy(w->peq, q.peq, sizeof(q.peq));
        w->B = q.b_max;
        w->L = n + q.pad;
        w->V = (w->B + 3) / 4;
        w->k = q.m;
        w->lv = (w->B - 1) / 4;
        w->ll = (w->B - 1) % 4;
        w->end = w->L + w->B - 1;
        uint8_t* tp = w->tb + 16;
        memset(w->tb, 0, 16);
        memcpy(tp, t, (size_t)n);
        memset(tp + n, 0, (size_t)(w->L - n) + 16);
        w->tp = tp;
        for(int v = 0; v < w->V; v++){
                w->P[v] = u4_dup(~(uint64_t)0);
                w->M[v] = u4_dup(0);
                w->SC[v] = u4_set(64 * (4 * v + 1), 64 * (4 * v + 2), 64 * (4 * v + 3), 64 * (4 * v + 4));
                w->HP[v] = u4_dup(0);
                w->HN[v] = u4_dup(0);
        }
}

static inline __attribute__((always_inline)) void ws_step(struct wstate* w, int s)
{
        const uint8_t* tp = w->tp;
        u64x4 prevp = u4_dup(0), prevn = u4_dup(0);
        for(int v = 0; v < w->V; v++){
                const int b0 = 4 * v;
                const u64x4 hp = u4_shift_in(prevp, w->HP[v]);
                const u64x4 hn = u4_shift_in(prevn, w->HN[v]);
                prevp = w->HP[v];
                prevn = w->HN[v];
                u64x4 Eq = u4_set(w->peq[tp[s - b0]][b0], w->peq[tp[s - b0 - 1]][b0 + 1],
                                  w->peq[tp[s - b0 - 2]][b0 + 2], w->peq[tp[s - b0 - 3]][b0 + 3]);
                const u64x4 Pv = w->P[v], Mv = w->M[v];
                const u64x4 Xv = u4_or(Eq, Mv);
                Eq = u4_or(Eq, hn);
                const u64x4 Xh = u4_or(u4_xor(u4_add(u4_and(Eq, Pv), Pv), Pv), Eq);
                u64x4 Ph = u4_or(Mv, u4_not(u4_or(Xh, Pv)));
                u64x4 Mh = u4_and(Pv, Xh);
                const u64x4 op = u4_shr63(Ph);
                const u64x4 on = u4_shr63(Mh);
                Ph = u4_or(u4_shl1(Ph), hp);
                Mh = u4_or(u4_shl1(Mh), hn);
                u64x4 nP = u4_or(Mh, u4_not(u4_or(Xv, Ph)));
                u64x4 nM = u4_and(Ph, Xv);
                u64x4 nS = u4_sub(u4_add(w->SC[v], op), on);
                if(s < b0 + 3){
                        const u64x4 act = u4_set(s >= b0 ? ~(uint64_t)0 : 0, s >= b0 + 1 ? ~(uint64_t)0 : 0,
                                                 s >= b0 + 2 ? ~(uint64_t)0 : 0, s >= b0 + 3 ? ~(uint64_t)0 : 0);
                        nP = u4_sel(act, nP, Pv);
                        nM = u4_sel(act, nM, Mv);
                        nS = u4_sel(act, nS, w->SC[v]);
                }
                w->P[v] = nP;
                w->M[v] = nM;
                w->SC[v] = nS;
                w->HP[v] = op;
                w->HN[v] = on;
        }
        if(s >= w->B - 1){
                const int sc = (int)(int64_t)u4_lane(w->SC[w->lv], w->ll);
                w->k = sc < w->k ? sc : w->k;
        }
}

static void bpm_simd_x2(const uint8_t* t0, const uint8_t* p0, int n0, int m0,
                        const uint8_t* t1, const uint8_t* p1, int n1, int m1, int out[2])
{
        static struct wstate w0, w1;
        ws_init(&w0, t0, p0, n0, m0);
        ws_init(&w1, t1, p1, n1, m1);
        const int e = w0.end > w1.end ? w0.end : w1.end;
        for(int s = 0; s < e; s++){
                if(s < w0.end) ws_step(&w0, s);
                if(s < w1.end) ws_step(&w1, s);
        }
        out[0] = w0.k;
        out[1] = w1.k;
}
#endif

/* two pairs at once; returns the results in out[0], out[1] */
static void bpm_x2(const uint8_t* t0, const uint8_t* p0, int n0, int m0,
                   const uint8_t* t1, const uint8_t* p1, int n1, int m1, int out[2])
{
        struct bpm_pat q0, q1;
        uint64_t P0[16], M0[16], P1[16], M1[16];
        int32_t s0[16], s1[16];
        bpm_pat_init(&q0, p0, m0);
        bpm_pat_init(&q1, p1, m1);
        const int B0 = q0.b_max, B1 = q1.b_max;
        const int L0 = n0 + q0.pad, L1 = n1 + q1.pad;
        int k0 = q0.m, k1 = q1.m;
        for(int b = 0; b < B0; b++){ P0[b] = ~(uint64_t)0; M0[b] = 0; s0[b] = (b + 1) * W64; }
        for(int b = 0; b < B1; b++){ P1[b] = ~(uint64_t)0; M1[b] = 0; s1[b] = (b + 1) * W64; }
        const int Lc = L0 < L1 ? L0 : L1;
        const int Bc = B0 < B1 ? B0 : B1;
        int i;
        for(i = 0; i < Lc; i++){
                const uint8_t c0 = (i >= n0) ? 0 : t0[i];
                const uint8_t c1 = (i >= n1) ? 0 : t1[i];
                int cy0 = 0, cy1 = 0;
                int b;
                for(b = 0; b < Bc; b++){
                        int h0, h1;
                        BPM_STEP(P0[b], M0[b], q0.peq[c0][b], cy0, h0);
                        BPM_STEP(P1[b], M1[b], q1.peq[c1][b], cy1, h1);
                        cy0 = h0; s0[b] += cy0;
                        cy1 = h1; s1[b] += cy1;
                }
                for(int bb = b; bb < B0; bb++){
                        int h0;
                        BPM_STEP(P0[bb], M0[bb], q0.peq[c0][bb], cy0, h0);
                        cy0 = h0; s0[bb] += cy0;
                }
                for(int bb = b; bb < B1; bb++){
                        int h1;
                        BPM_STEP(P1[bb], M1[bb], q1.peq[c1][bb], cy1, h1);
                        cy1 = h1; s1[bb] += cy1;
                }
                if(s0[B0 - 1] < k0) k0 = s0[B0 - 1];
                if(s1[B1 - 1] < k1) k1 = s1[B1 - 1];
        }
        for(int j = i; j < L0; j++){
                const uint8_t c0 = (j >= n0) ? 0 : t0[j];
                int cy0 = 0;
                for(int b = 0; b < B0; b++){
                        int h0;
                        BPM_STEP(P0[b], M0[b], q0.peq[c0][b], cy0, h0);
                        cy0 = h0; s0[b] += cy0;
                }
                if(s0[B0 - 1] < k0) k0 = s0[B0 - 1];
        }
        for(int j = i; j < L1; j++){
                const uint8_t c1 = (j >= n1) ? 0 : t1[j];
                int cy1 = 0;
                for(int b = 0; b < B1; b++){
                        int h1;
                        BPM_STEP(P1[b], M1[b], q1.peq[c1][b], cy1, h1);
                        cy1 = h1; s1[b] += cy1;
                }
                if(s1[B1 - 1] < k1) k1 = s1[B1 - 1];
        }
        out[0] = k0;
        out[1] = k1;
}

static double now(void)
{
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char** argv)
{
        static uint8_t* s[20000];
        static int len[20000];
        int n = 0;
        char line[1 << 16];
        FILE* f = fopen(argv[1], "r");
        int reps = argc > 2 ? atoi(argv[2]) : 3;
        if(!f){ perror(argv[1]); return 1; }
        while(fgets(line, sizeof(line), f)){
                if(line[0] == '>'){
                        s[n] = malloc(1 << 16);
                        len[n] = 0;
                        n++;
                }else if(n){
                        for(char* c = line; *c; c++){
                                if(*c >= 'A' && *c <= 'Z') s[n - 1][len[n - 1]++] = (uint8_t)((*c - 'A') % SIGMA);
                        }
                }
        }
        fclose(f);
        /* pairs (i, j), i < j; longer = text */
        long np = (long)n * (n - 1) / 2;
        int (*pr)[2] = malloc(sizeof(int[2]) * np);
        long z = 0;
        for(int i = 0; i < n; i++) for(int j = i + 1; j < n; j++){ pr[z][0] = i; pr[z][1] = j; z++; }
        int* r_ref = malloc(sizeof(int) * np);
        int* r_s = malloc(sizeof(int) * np);
        int* r_x2 = malloc(sizeof(int) * np);
        int* r_bl = malloc(sizeof(int) * np);
#define TXT(k) (len[pr[k][0]] > len[pr[k][1]] ? pr[k][0] : pr[k][1])
#define PAT(k) (len[pr[k][0]] > len[pr[k][1]] ? pr[k][1] : pr[k][0])
        double best[4] = {1e30, 1e30, 1e30, 1e30};
        for(int rep = 0; rep < reps; rep++){
                double t0 = now();
                for(long k = 0; k < np; k++) r_ref[k] = bpm_block(s[TXT(k)], s[PAT(k)], len[TXT(k)], len[PAT(k)]);
                double t1 = now();
                for(long k = 0; k < np; k++) r_s[k] = bpm_simple(s[TXT(k)], s[PAT(k)], len[TXT(k)], len[PAT(k)]);
                double t2 = now();
                long k;
                for(k = 0; k + 1 < np; k += 2){
                        int o[2];
                        bpm_simd_x2(s[TXT(k)], s[PAT(k)], len[TXT(k)], len[PAT(k)],
                               s[TXT(k + 1)], s[PAT(k + 1)], len[TXT(k + 1)], len[PAT(k + 1)], o);
                        r_x2[k] = o[0];
                        r_x2[k + 1] = o[1];
                }
                if(k < np) r_x2[k] = bpm_simple(s[TXT(k)], s[PAT(k)], len[TXT(k)], len[PAT(k)]);
                double t3 = now();
                for(long kk = 0; kk < np; kk++) r_bl[kk] = bpm_simd(s[TXT(kk)], s[PAT(kk)], len[TXT(kk)], len[PAT(kk)]);
                double t4 = now();
                if(t4 - t3 < best[3]) best[3] = t4 - t3;
                if(t1 - t0 < best[0]) best[0] = t1 - t0;
                if(t2 - t1 < best[1]) best[1] = t2 - t1;
                if(t3 - t2 < best[2]) best[2] = t3 - t2;
        }
        long bad_s = 0, bad_x2 = 0, bad_bl = 0;
        for(long k = 0; k < np; k++) bad_bl += r_bl[k] != r_ref[k];
        for(long k = 0; k < np; k++){
                bad_s += r_s[k] != r_ref[k];
                bad_x2 += r_x2[k] != r_ref[k];
        }
        printf("%ld pairs | mismatches simple %ld x2 %ld bl %ld | us/pair ref %.2f simple %.2f x2 %.2f bl %.2f | speedup simple %.2f x2 %.2f bl %.2f\n",
               np, bad_s, bad_x2, bad_bl, 1e6 * best[0] / np, 1e6 * best[1] / np, 1e6 * best[2] / np, 1e6 * best[3] / np,
               best[0] / best[1], best[0] / best[2], best[0] / best[3]);
        return 0;
}
