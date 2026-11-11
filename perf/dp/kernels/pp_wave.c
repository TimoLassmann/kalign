/* pp_wave -- inter-row wavefront profile-profile forward and backward passes
 * (NEON or AVX2 via lib/src/wave_simd.h), bit-identical to
 * lib/src/aln_profileprofile.c.
 *
 * Same scheme as ss_wave.c (strips of 8 rows, lane q at column t = s - q).
 * Profile-specific parts:
 *  - The reference adds the substitution terms straight into the cell
 *    score, pa += p1[c] * p2[32+c], over the non-zero residues of the row
 *    profile in descending c.  Each lane does the same with vf8_muladd over
 *    the union U of the strip rows' non-zero residues (descending); a
 *    residue a row lacks has weight p1[c] == 0 and adds +0, so the per-lane
 *    sequence of roundings is the reference one.  vf8_muladd rounds like the
 *    scalar statement as the compiler builds it (fused where the compiler
 *    contracts it); the scalar parts here use the reference's own expression.
 *  - Lanes read consecutive columns (t - q), so with the column profile
 *    transposed per residue (and stored reversed for the forward pass) each
 *    term is one contiguous vector load instead of a per-lane gather.  The
 *    same holds for the column gap terms.
 *  - Cell (row r, residue y): D = max3(a, ga + O(y), gb + GO1(r)), with
 *    O(y) = P2[y+1][27], GO1(r) = P1[r+1][27];
 *    ga(t) = max(ga(t-1) + E(t), a(t-1) + O(t)), E = P2[y+1][28] forward,
 *    P2[y+2][28] backward; gb = max(gb + GE1(r), a_up + GO1(r)), GE1 =
 *    P1[r+1][28] forward, P1[r+2][28] backward.
 *
 * Column arrays (transposed scores, gaps) are built once per problem in
 * .prepare.  Kernel: pp_wave8.
 */
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "kdp_kernel.h"
#include "aln_profileprofile.h"
#include "wave_simd.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MAX3(a, b, c) MAX(MAX(a, b), c)

#if WAVE_SIMD

#ifndef PW_MIN_W
#define PW_MIN_W 16
#endif
#define QPAD 32
#ifndef WAVE_PART_MIN
#define WAVE_PART_MIN 3
#endif
#define WL 8

/* Column arrays, indexed by residue y (0..len_b-1), padded by QPAD.
   Forward copies are reversed: index u = (len_b - 1) - y. */
struct pw_state {
        int n;          /* len_b + 2*QPAD */
        float* tf;      /* 23 x n, forward (reversed) P2[y+1][32+c] */
        float* tb;      /* 23 x n, backward */
        float* of;      /* P2[y+1][27] reversed */
        float* ob;
        float* ef;      /* forward E = P2[y+1][28], reversed */
        float* eb;      /* backward E = P2[y+2][28] */
};

static void pw_release(void* state);

static int pw_prepare(struct aln_mem* m, int kind, void** state)
{
        struct pw_state* st;
        const int lb = m->len_b;
        int c, y;
        *state = NULL;
        if(kind != KDP_PP) return 0;
        st = calloc(1, sizeof(*st));
        if(!st) return 1;
        st->n = lb + 2 * QPAD;
        st->tf = calloc((size_t)23 * st->n, sizeof(float));
        st->tb = calloc((size_t)23 * st->n, sizeof(float));
        st->of = calloc(st->n, sizeof(float));
        st->ob = calloc(st->n, sizeof(float));
        st->ef = calloc(st->n, sizeof(float));
        st->eb = calloc(st->n, sizeof(float));
        if(!st->tf || !st->tb || !st->of || !st->ob || !st->ef || !st->eb){
                pw_release(st);
                return 1;
        }
        for(y = 0; y < lb; y++){
                const float* p = m->prof2 + (size_t)(y + 1) * 64;
                const int u = lb - 1 - y;
                for(c = 0; c < 23; c++){
                        st->tf[(size_t)c * st->n + QPAD + u] = p[32 + c];
                        st->tb[(size_t)c * st->n + QPAD + y] = p[32 + c];
                }
                st->of[QPAD + u] = p[27];
                st->ob[QPAD + y] = p[27];
                st->ef[QPAD + u] = p[28];
                st->eb[QPAD + y] = p[64 + 28];
        }
        *state = st;
        return 0;
}

static void pw_release(void* state)
{
        struct pw_state* st = state;
        if(!st) return;
        free(st->tf); free(st->tb); free(st->of); free(st->ob); free(st->ef); free(st->eb);
        free(st);
}

struct pw_buf {
        float* ra;
        float* rga;
        float* rgb;
        float* rd;
        float* bl;
        size_t cap;
};

static struct pw_buf pb = {0};

static void pw_reserve(int W)
{
        size_t n = (size_t)W + 2 * QPAD;
        if(n > pb.cap){
                free(pb.ra); free(pb.rga); free(pb.rgb); free(pb.rd); free(pb.bl);
                pb.ra = malloc(n * sizeof(float));
                pb.rga = malloc(n * sizeof(float));
                pb.rgb = malloc(n * sizeof(float));
                pb.rd = malloc(n * sizeof(float));
                pb.bl = calloc(n * WL, sizeof(float));
                if(!pb.ra || !pb.rga || !pb.rgb || !pb.rd || !pb.bl) abort();
                pb.cap = n;
        }
}

/* Substitution score of row profile p1 (non-zero residues nz[0..f],
   ascending) against column profile p2, added to pa exactly as the
   reference statement does (same expression, so the same contraction). */
static inline float pp_term(float pa, const float* p1, const int* nz, int f, const float* p2)
{
        for(int c = f; c >= 0; c--){
                pa += p1[nz[c]] * p2[32 + nz[c]];
        }
        return pa;
}

static inline __attribute__((always_inline)) int pw_sweep(struct aln_mem* m, const int dir)
{
        const struct pw_state* st = kdp_kernel_state;
        struct states* s = dir > 0 ? m->f : m->b;
        const float* P1 = m->prof1;
        const float* P2 = m->prof2;
        const int startb = m->startb;
        const int endb = m->endb;
        const int lb = m->len_b;
        const int W = endb - startb - 1;
        const struct sparse_bonus* cons = m->consistency;
        const int lcol = dir > 0 ? startb : endb;
        const int rcol = dir > 0 ? endb : startb;
        const int lterm = dir > 0 ? (startb == 0) : (endb == m->len_b);
        const int rterm = dir > 0 ? (endb == m->len_b) : (startb == 0);
        const int y0 = dir > 0 ? startb : endb - 1;     /* y(t) = y0 + dir*t */
        const int yr = dir > 0 ? endb - 1 : startb;
        const int bstride = W + 2 * QPAD;
        /* lane q at step s needs y = y0 + dir*(s - q); in the direction's
           column arrays that is index c0 - s + q (lanes ascending):
           forward (reversed arrays): u = lb-1-y = (lb-1-y0) - s + q;
           backward: y = y0 - s + q. */
        const int c0 = dir > 0 ? (lb - 1 - y0) : y0;
        const float* T = (dir > 0 ? st->tf : st->tb) + QPAD;
        const float* OC = (dir > 0 ? st->of : st->ob) + QPAD;
        const float* EC = (dir > 0 ? st->ef : st->eb) + QPAD;
        const int n = st->n;
        int r, rend, h, t, q, sstep;

        pw_reserve(W);
        float* restrict RA = pb.ra + QPAD;
        float* restrict RGA = pb.rga + QPAD;
        float* restrict RGB = pb.rgb + QPAD;
        float* restrict RD = pb.rd + QPAD;

#define YT(tt) (y0 + dir * (tt))
#define OCOL(tt) (P2[(size_t)(YT(tt) + 1) * 64 + 27])
#define ECOL(tt) (dir > 0 ? P2[(size_t)(YT(tt) + 1) * 64 + 28] : P2[(size_t)(YT(tt) + 2) * 64 + 28])
#define TCOL(tt) (P2[(size_t)(YT(tt) + 1) * 64 + 29])
#define ROWP(rr) (P1 + (size_t)((rr) + 1) * 64)
#define GE1(rr) (dir > 0 ? P1[(size_t)((rr) + 1) * 64 + 28] : P1[(size_t)((rr) + 2) * 64 + 28])

        /* ---- first row ---- */
        s[lcol].a = s[0].a;
        s[lcol].ga = s[0].ga;
        s[lcol].gb = s[0].gb;
        RA[-1] = s[lcol].a;
        RGA[-1] = s[lcol].ga;
        RGB[-1] = s[lcol].gb;
        {
                float pga = RGA[-1], pa = RA[-1];
                for(t = 0; t < W; t++){
                        float g = lterm ? MAX(pga, pa) + TCOL(t) : MAX(pga + ECOL(t), pa + OCOL(t));
                        RA[t] = -FLT_MAX;
                        RGA[t] = g;
                        RGB[t] = -FLT_MAX;
                        pga = g;
                        pa = -FLT_MAX;
                }
        }
        /* RD of the row before the first swept row: its GO1 is that row's */
        {
                const float go = ROWP(dir > 0 ? m->starta - 1 : m->enda_2)[27];
                RD[-1] = MAX3(RA[-1], RGA[-1] + OCOL(-1), RGB[-1] + go);
                for(t = 0; t < W; t++){
                        RD[t] = MAX3(RA[t], RGA[t] + OCOL(t), RGB[t] + go);
                }
        }
        for(t = W; t < W + WL + 1; t++){
                RA[t] = -FLT_MAX; RGA[t] = -FLT_MAX; RGB[t] = -FLT_MAX; RD[t] = -FLT_MAX;
        }
        float ra = -FLT_MAX, rga = -FLT_MAX, rgb = -FLT_MAX;

        const vf8 vneg = vf8_set1(-FLT_MAX);

        if(dir > 0){ r = m->starta; rend = m->enda; h = m->enda - m->starta; }
        else { r = m->enda_2 - 1; rend = m->starta_2 - 1; h = m->enda_2 - m->starta_2; }

        for(int hs = 0; h >= WAVE_PART_MIN; h -= hs, r += dir * hs){
                hs = h < WL ? h : WL;
                float lgb[WL], dleft[WL], rbon[WL], dw[WL];
                float dcap[WL][WL];
                float go1[WL], ge1[WL];
                const float* bp[WL];
                uint32_t umask = 0;
                int ulist[23], nu = 0;

                /* left column chain */
                {
                        float pa0 = RA[-1], pgb0 = RGB[-1];
                        for(q = 0; q < WL; q++){
                                const int rq = r + dir * (q < hs ? q : hs - 1);
                                const float* p1 = ROWP(rq);
                                go1[q] = p1[27];
                                ge1[q] = GE1(rq);
                                lgb[q] = lterm ? MAX(pgb0, pa0) + p1[29] : MAX(pgb0 + ge1[q], pa0 + p1[27]);
                                dleft[q] = MAX3(-FLT_MAX, -FLT_MAX + OCOL(-1), lgb[q] + go1[q]);
                                pa0 = -FLT_MAX;
                                pgb0 = lgb[q];
                                for(int c = 0; c < 23; c++) if(p1[c]) umask |= 1u << c;
                        }
                }
                for(int c = 22; c >= 0; c--) if(umask & (1u << c)) ulist[nu++] = c;

                const vf8 vgo1 = vf8_loadu(go1);
                const vf8 vge1 = vf8_loadu(ge1);
                vf8 vw[23];
                for(int k = 0; k < nu; k++){
                        float tmp[WL];
                        for(q = 0; q < WL; q++) tmp[q] = ROWP(r + dir * (q < hs ? q : hs - 1))[ulist[k]];
                        vw[k] = vf8_loadu(tmp);
                }
                for(q = 0; q < WL; q++){
                        const int rq = r + dir * (q < hs ? q : hs - 1);
                        rbon[q] = 0.0F;
                        bp[q] = pb.bl + (size_t)q * bstride + QPAD - q;
                        if(cons && q < hs){
                                const int bb = rq * cons->K;
                                for(int e = 0; e < cons->K; e++){
                                        int c = cons->cols[bb + e];
                                        if(c < 0) break;
                                        if(c == yr){
                                                rbon[q] += cons->vals[bb + e];
                                        }else{
                                                int tt = (c - y0) * dir;
                                                if(tt >= 0 && tt < W) pb.bl[(size_t)q * bstride + QPAD + tt] += cons->vals[bb + e];
                                        }
                                }
                        }
                }
                const float rdw_prev = RD[W - 1];

                vf8 pA = vf8_set_lane0(vneg, -FLT_MAX);
                vf8 pGA = vf8_set_lane0(vneg, -FLT_MAX);
                vf8 pGB = vf8_set_lane0(vneg, lgb[0]);
                vf8 pD = vf8_set_lane0(vneg, dleft[0]);
                vf8 p2D = vneg;

                const int send = W + hs - 1;
                for(sstep = 0; sstep < send; sstep++){
                        const int ci = c0 - sstep;      /* lane q: ci + q */
                        vf8 a = vf8_shift_in(p2D, RD[sstep - 1]);
                        for(int k = 0; k < nu; k++){
                                a = vf8_muladd(a, vw[k], vf8_loadu(T + (size_t)ulist[k] * n + ci));
                        }
                        const vf8 upA = vf8_shift_in(pA, RA[sstep]);
                        const vf8 upGB = vf8_shift_in(pGB, RGB[sstep]);
                        const vf8 oc = vf8_loadu(OC + ci);
                        const vf8 ec = vf8_loadu(EC + ci);
                        if(cons){
                                a = vf8_add(a, vf8_gather(bp, sstep));
                        }
                        vf8 ga = vf8_max(vf8_add(pGA, ec), vf8_add(pA, oc));
                        vf8 gb = vf8_max(vf8_add(upGB, vge1), vf8_add(upA, vgo1));
                        vf8 d = vf8_max(vf8_max(a, vf8_add(ga, oc)), vf8_add(gb, vgo1));
                        if(sstep < WL - 1){
                                const int ql = sstep + 1;
                                a = vf8_with_lane(a, ql, -FLT_MAX);
                                ga = vf8_with_lane(ga, ql, -FLT_MAX);
                                gb = vf8_with_lane(gb, ql, lgb[ql]);
                                d = vf8_with_lane(d, ql, dleft[ql]);
                        }
                        if(sstep >= hs - 1){
                                const int tw = sstep - (hs - 1);
                                if(hs == WL){
                                        RA[tw] = vf8_last(a);
                                        RGA[tw] = vf8_last(ga);
                                        RGB[tw] = vf8_last(gb);
                                        RD[tw] = vf8_last(d);
                                }else{
                                        RA[tw] = vf8_lane(a, hs - 1);
                                        RGA[tw] = vf8_lane(ga, hs - 1);
                                        RGB[tw] = vf8_lane(gb, hs - 1);
                                        RD[tw] = vf8_lane(d, hs - 1);
                                }
                        }
                        if(sstep >= W - 1){
                                vf8_storeu(dcap[sstep - (W - 1)], d);
                        }
                        p2D = pD;
                        pD = d;
                        pA = a;
                        pGA = ga;
                        pGB = gb;
                }
                for(q = 0; q < hs; q++) dw[q] = dcap[q][q];

                /* right column (scalar) */
                {
                        const float* p2r = P2 + (size_t)(yr + 1) * 64;
                        for(q = 0; q < hs; q++){
                                const float* p1 = ROWP(r + dir * q);
                                int nz[23], f = 0;
                                for(int c = 0; c < 23; c++) if(p1[c]) nz[f++] = c;
                                const float ca = ra;
                                float p = pp_term(q ? dw[q - 1] : rdw_prev, p1, nz, f - 1, p2r);
                                if(cons) p += rbon[q];
                                ra = p;
                                rga = -FLT_MAX;
                                rgb = rterm ? MAX(rgb, ca + p1[27]) + p1[29] : MAX(rgb + ge1[q], ca + p1[27]);
                        }
                }
                RA[-1] = -FLT_MAX;
                RGA[-1] = -FLT_MAX;
                RGB[-1] = lgb[hs - 1];
                RD[-1] = dleft[hs - 1];
                if(cons){
                        for(q = 0; q < hs; q++){
                                const int bb = (r + dir * q) * cons->K;
                                for(int e = 0; e < cons->K; e++){
                                        int c = cons->cols[bb + e];
                                        if(c < 0) break;
                                        int tt = (c - y0) * dir;
                                        if(c != yr && tt >= 0 && tt < W) pb.bl[(size_t)q * bstride + QPAD + tt] = 0.0F;
                                }
                        }
                }
        }

        /* ---- leftover rows: scalar ---- */
        for(; r != rend; r += dir){
                const float* p1 = ROWP(r);
                const float go = p1[27], ge = GE1(r);
                int nz[23], f = 0;
                for(int c = 0; c < 23; c++) if(p1[c]) nz[f++] = c;
                f--;
                float pa0 = RA[-1], pgb0 = RGB[-1];
                float pa = RD[-1];      /* diag for t = 0 */
                float xa = -FLT_MAX, xga = -FLT_MAX;
                RA[-1] = -FLT_MAX;
                RGA[-1] = -FLT_MAX;
                RGB[-1] = lterm ? MAX(pgb0, pa0) + p1[29] : MAX(pgb0 + ge, pa0 + go);
                RD[-1] = MAX3(RA[-1], RGA[-1] + OCOL(-1), RGB[-1] + go);
                for(t = 0; t < W; t++){
                        const float ca = RA[t], cgb = RGB[t];
                        const float dnext = RD[t];
                        const float* p2 = P2 + (size_t)(YT(t) + 1) * 64;
                        float a = pp_term(pa, p1, nz, f, p2);
                        if(cons) a += sparse_bonus_lookup(cons, r, YT(t));
                        RA[t] = a;
                        RGA[t] = MAX(xga + ECOL(t), xa + OCOL(t));
                        RGB[t] = MAX(cgb + ge, ca + go);
                        RD[t] = MAX3(RA[t], RGA[t] + OCOL(t), RGB[t] + go);
                        xa = RA[t];
                        xga = RGA[t];
                        pa = dnext;
                }
                {
                        const float ca = ra;
                        float p = pp_term(pa, p1, nz, f, P2 + (size_t)(yr + 1) * 64);
                        if(cons) p += sparse_bonus_lookup(cons, r, yr);
                        ra = p;
                        rga = -FLT_MAX;
                        rgb = rterm ? MAX(rgb, ca + go) + p1[29] : MAX(rgb + ge, ca + go);
                }
        }

        s[lcol].a = RA[-1];
        s[lcol].ga = RGA[-1];
        s[lcol].gb = RGB[-1];
        for(t = 0; t < W; t++){
                const int c = lcol + dir * (t + 1);
                s[c].a = RA[t];
                s[c].ga = RGA[t];
                s[c].gb = RGB[t];
        }
        s[rcol].a = ra;
        s[rcol].ga = rga;
        s[rcol].gb = rgb;
        return 0;
#undef YT
#undef OCOL
#undef ECOL
#undef TCOL
#undef ROWP
#undef GE1
}

static int pw_fwd(struct aln_mem* m)
{
        if(!kdp_kernel_state || m->endb - m->startb - 1 < PW_MIN_W || m->enda <= m->starta)
                return aln_profileprofile_foward(m);
        return pw_sweep(m, 1);
}

static int pw_bwd(struct aln_mem* m)
{
        if(!kdp_kernel_state || m->endb - m->startb - 1 < PW_MIN_W || m->enda_2 <= m->starta_2)
                return aln_profileprofile_backward(m);
        return pw_sweep(m, -1);
}

static const struct kdp_kernel pp_wave_kernel = {
        .name = "pp_wave8",
        .description = "SIMD (NEON/AVX2) inter-row wavefront profile-profile fwd/bwd, 8 rows (exact)",
        .fwd = { [KDP_PP] = pw_fwd },
        .bwd = { [KDP_PP] = pw_bwd },
        .prepare = pw_prepare,
        .release = pw_release,
};
KDP_REGISTER_KERNEL(pp_wave_kernel)

#endif
