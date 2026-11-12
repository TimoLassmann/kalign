/* ss_wave -- inter-row wavefront seq-seq forward and backward passes
 * (NEON or AVX2 via lib/src/wave_simd.h), bit-identical to
 * lib/src/aln_seqseq.c.
 *
 * Rows are processed in strips of L = 8 rows.  At step s, lane q works on
 * row q of the strip, logical column t = s - q, so every dependency of a
 * cell is either the lane's own previous step (row gap ga: (q, t-1)), lane
 * q-1 one step ago (column gap gb and the cell above: (q-1, t)) or lane q-1
 * two steps ago (diagonal: (q-1, t-1)).  Each cell evaluates exactly the
 * reference expressions, so no fix-up is needed and the floats match.
 *
 * Lane 0 reads the row before the strip from SoA row arrays (RA, RGA, RGB
 * and RD = max3(a, ga - gpo, gb - gpo)); the strip's last row (lane hs-1)
 * writes back.  The left boundary column is injected as t = -1 during the
 * strip prologue; the right edge column (terminal rules) is scalar.  The
 * last strip of a call may hold hs < L rows (down to WAVE_PART_MIN): its
 * lanes hs..L-1 repeat row hs-1 and are discarded (lanes only feed higher
 * lanes).  Fewer rows use a scalar row routine on the same arrays.
 *
 * Substitution scores come from QL[x][y] = subm[x][seq2[y]] - soff, built
 * once per problem in .prepare (23 x len_b floats; in Kalign once per
 * pairwise alignment), so Hirschberg sub-calls need no per-call set-up
 * beyond copying their row in and out.
 *
 * Kernel: ss_wave8.
 */
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "kdp_kernel.h"
#include "aln_seqseq.h"
#include "wave_simd.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MAX3(a, b, c) MAX(MAX(a, b), c)

#if WAVE_SIMD

#ifndef SW_MIN_W
#define SW_MIN_W 16
#endif
#define QPAD 32                 /* padding of QL rows and bonus rows (>= L + 2) */
#ifndef WAVE_PART_MIN
#define WAVE_PART_MIN 3
#endif
#define WL 8                    /* rows per strip = vf8 lanes */

struct sw_state {
        float* ql;      /* 23 rows of (len_b + 2*QPAD) floats */
        int stride;
};

static int sw_prepare(struct aln_mem* m, int kind, void** state)
{
        struct sw_state* st;
        int x, y;
        *state = NULL;
        if(kind != KDP_SS) return 0;
        st = calloc(1, sizeof(*st));
        if(!st) return 1;
        st->stride = m->len_b + 2 * QPAD;
        st->ql = malloc(sizeof(float) * 23 * (size_t)st->stride);
        if(!st->ql){
                free(st);
                return 1;
        }
        for(x = 0; x < 23; x++){
                float* row = st->ql + (size_t)x * st->stride + QPAD;
                for(y = -QPAD; y < m->len_b + QPAD; y++){
                        row[y] = (y >= 0 && y < m->len_b) ? m->ap->subm[x][m->seq2[y]] - m->ap->subm_offset : 0.0F;
                }
        }
        *state = st;
        return 0;
}

static void sw_release(void* state)
{
        struct sw_state* st = state;
        if(!st) return;
        free(st->ql);
        free(st);
}

struct sw_buf {
        float* ra;      /* row arrays, index t + QPAD (t = -1 is the left column) */
        float* rga;
        float* rgb;
        float* rd;
        float* bl;      /* WL bonus rows of (W + 2*QPAD) */
        size_t cap;
};

static struct sw_buf wb = {0};

static void sw_reserve(int W)
{
        size_t n = (size_t)W + 2 * QPAD;
        if(n > wb.cap){
                free(wb.ra); free(wb.rga); free(wb.rgb); free(wb.rd); free(wb.bl);
                wb.ra = malloc(n * sizeof(float));
                wb.rga = malloc(n * sizeof(float));
                wb.rgb = malloc(n * sizeof(float));
                wb.rd = malloc(n * sizeof(float));
                wb.bl = calloc(n * WL, sizeof(float));
                if(!wb.ra || !wb.rga || !wb.rgb || !wb.rd || !wb.bl) abort();
                wb.cap = n;
        }
}

static inline __attribute__((always_inline)) WAVE_TARGET int sw_sweep(struct aln_mem* m, const int dir)
{
        const struct sw_state* st = kdp_kernel_state;
        struct states* s = dir > 0 ? m->f : m->b;
        const uint8_t* seq1 = m->seq1;
        const int startb = m->startb;
        const int endb = m->endb;
        const int W = endb - startb - 1;
        const float gpo = m->ap->gpo;
        const float gpe = m->ap->gpe;
        const float tgpe = m->ap->tgpe;
        const struct sparse_bonus* cons = m->consistency;
        const int lcol = dir > 0 ? startb : endb;
        const int rcol = dir > 0 ? endb : startb;
        const int lterm = dir > 0 ? (startb == 0) : (endb == m->len_b);
        const int rterm = dir > 0 ? (endb == m->len_b) : (startb == 0);
        const int y0 = dir > 0 ? startb : endb - 1;     /* y(t) = y0 + dir*t */
        const int yr = dir > 0 ? endb - 1 : startb;
        const int bstride = W + 2 * QPAD;
        int r, rend, h, t, q, sstep;

        sw_reserve(W);
        float* restrict RA = wb.ra + QPAD;      /* RA[-1] = left column */
        float* restrict RGA = wb.rga + QPAD;
        float* restrict RGB = wb.rgb + QPAD;
        float* restrict RD = wb.rd + QPAD;

        /* ---- first row (identical to the reference init) ---- */
        s[lcol].a = s[0].a;
        s[lcol].ga = s[0].ga;
        s[lcol].gb = s[0].gb;
        RA[-1] = s[lcol].a;
        RGA[-1] = s[lcol].ga;
        RGB[-1] = s[lcol].gb;
        {
                float pga = RGA[-1], pa = RA[-1];
                for(t = 0; t < W; t++){
                        float g = lterm ? MAX(pga, pa) - tgpe : MAX(pga - gpe, pa - gpo);
                        RA[t] = -FLT_MAX;
                        RGA[t] = g;
                        RGB[t] = -FLT_MAX;
                        pga = g;
                        pa = -FLT_MAX;
                }
        }
        for(t = -1; t < W; t++){
                RD[t] = MAX3(RA[t], RGA[t] - gpo, RGB[t] - gpo);
        }
        for(t = W; t < W + WL + 1; t++){
                RA[t] = -FLT_MAX; RGA[t] = -FLT_MAX; RGB[t] = -FLT_MAX; RD[t] = -FLT_MAX;
        }
        float ra = -FLT_MAX, rga = -FLT_MAX, rgb = -FLT_MAX;   /* right column */

        const vf8 vgpo = vf8_set1(gpo);
        const vf8 vgpe = vf8_set1(gpe);
        const vf8 vneg = vf8_set1(-FLT_MAX);

        if(dir > 0){ r = m->starta; rend = m->enda; h = m->enda - m->starta; }
        else { r = m->enda_2 - 1; rend = m->starta_2 - 1; h = m->enda_2 - m->starta_2; }

        const float* qlbase = st->ql + QPAD;

        for(int hs = 0; h >= WAVE_PART_MIN; h -= hs, r += dir * hs){
                hs = h < WL ? h : WL;
                const float* qp[WL];
                const float* bp[WL];
                float lgb[WL], dleft[WL], rbon[WL], dw[WL];
                float dcap[WL][WL];

                /* left column of each strip row (scalar chain) */
                {
                        float pa0 = RA[-1], pgb0 = RGB[-1];
                        for(q = 0; q < WL; q++){
                                lgb[q] = lterm ? MAX(pgb0, pa0) - tgpe : MAX(pgb0 - gpe, pa0 - gpo);
                                dleft[q] = MAX3(-FLT_MAX, -FLT_MAX - gpo, lgb[q] - gpo);
                                pa0 = -FLT_MAX;
                                pgb0 = lgb[q];
                        }
                }
                for(q = 0; q < WL; q++){
                        const int rq = r + dir * (q < hs ? q : hs - 1);
                        /* lane q at step s reads y(s - q) = y0 + dir*(s - q) */
                        qp[q] = qlbase + (size_t)seq1[rq] * st->stride + y0 - dir * q;
                        rbon[q] = 0.0F;
                        bp[q] = wb.bl + (size_t)q * bstride + QPAD - q;
                        if(cons && q < hs){
                                const int bb = rq * cons->K;
                                for(int e = 0; e < cons->K; e++){
                                        int c = cons->cols[bb + e];
                                        if(c < 0) break;
                                        if(c == yr){
                                                rbon[q] += cons->vals[bb + e];
                                        }else{
                                                int tt = (c - y0) * dir;
                                                if(tt >= 0 && tt < W) wb.bl[(size_t)q * bstride + QPAD + tt] += cons->vals[bb + e];
                                        }
                                }
                        }
                }
                const float rdw_prev = RD[W - 1];

                /* step -1: lane 0 holds the left column of strip row 0 */
                vf8 pA = vf8_set_lane0(vneg, -FLT_MAX);
                vf8 pGA = vf8_set_lane0(vneg, -FLT_MAX);
                vf8 pGB = vf8_set_lane0(vneg, lgb[0]);
                vf8 pD = vf8_set_lane0(vneg, dleft[0]);
                vf8 p2D = vneg;

                const int send = W + hs - 1;
                for(sstep = 0; sstep < send; sstep++){
                        const vf8 upA = vf8_shift_in(pA, RA[sstep]);
                        const vf8 upGB = vf8_shift_in(pGB, RGB[sstep]);
                        const vf8 dg = vf8_shift_in(p2D, RD[sstep - 1]);
                        vf8 a = vf8_add(dg, vf8_gather(qp, dir * sstep));
                        if(cons){
                                a = vf8_add(a, vf8_gather(bp, sstep));
                        }
                        vf8 ga = vf8_max(vf8_sub(pGA, vgpe), vf8_sub(pA, vgpo));
                        vf8 gb = vf8_max(vf8_sub(upGB, vgpe), vf8_sub(upA, vgpo));
                        vf8 d = vf8_max(vf8_max(a, vf8_sub(ga, vgpo)), vf8_sub(gb, vgpo));
                        /* prologue: lane sstep+1 enters at t = -1 (left column) */
                        if(sstep < WL - 1){
                                const int ql = sstep + 1;
                                a = vf8_with_lane(a, ql, -FLT_MAX);
                                ga = vf8_with_lane(ga, ql, -FLT_MAX);
                                gb = vf8_with_lane(gb, ql, lgb[ql]);
                                d = vf8_with_lane(d, ql, dleft[ql]);
                        }
                        /* the strip's last row finished column t = sstep - (hs - 1) */
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

                /* right column, rows in order (scalar, as the reference) */
                for(q = 0; q < hs; q++){
                        const float ca = ra;
                        float p = (q ? dw[q - 1] : rdw_prev) + qp[q][dir * (W + q)];
                        if(cons) p += rbon[q];
                        ra = p;
                        rga = -FLT_MAX;
                        rgb = rterm ? MAX(rgb, ca - gpo) - tgpe : MAX(rgb - gpe, ca - gpo);
                }
                /* left column after the strip = strip row hs-1 */
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
                                        if(c != yr && tt >= 0 && tt < W) wb.bl[(size_t)q * bstride + QPAD + tt] = 0.0F;
                                }
                        }
                }
        }

        /* ---- leftover rows: scalar on the row arrays ---- */
        for(; r != rend; r += dir){
                const float* ql = qlbase + (size_t)seq1[r] * st->stride + y0;
                float pa0 = RA[-1], pga0 = RGA[-1], pgb0 = RGB[-1];
                float pa = MAX3(pa0, pga0 - gpo, pgb0 - gpo);   /* diag for t = 0 */
                float xa = -FLT_MAX, xga = -FLT_MAX;
                RA[-1] = -FLT_MAX;
                RGA[-1] = -FLT_MAX;
                RGB[-1] = lterm ? MAX(pgb0, pa0) - tgpe : MAX(pgb0 - gpe, pa0 - gpo);
                RD[-1] = MAX3(RA[-1], RGA[-1] - gpo, RGB[-1] - gpo);
                for(t = 0; t < W; t++){
                        const float ca = RA[t], cga = RGA[t], cgb = RGB[t];
                        const float dnext = MAX3(ca, cga - gpo, cgb - gpo);
                        float a = pa + ql[dir * t];
                        if(cons) a += sparse_bonus_lookup(cons, r, y0 + dir * t);
                        RA[t] = a;
                        RGA[t] = MAX(xga - gpe, xa - gpo);
                        RGB[t] = MAX(cgb - gpe, ca - gpo);
                        RD[t] = MAX3(RA[t], RGA[t] - gpo, RGB[t] - gpo);
                        xa = RA[t];
                        xga = RGA[t];
                        pa = dnext;
                }
                {
                        const float ca = ra;
                        float p = pa + ql[dir * W];
                        if(cons) p += sparse_bonus_lookup(cons, r, yr);
                        ra = p;
                        rga = -FLT_MAX;
                        rgb = rterm ? MAX(rgb, ca - gpo) - tgpe : MAX(rgb - gpe, ca - gpo);
                }
        }

        /* write back */
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
}

static WAVE_TARGET int sw_fwd(struct aln_mem* m)
{
        if(!kalign_simd_on() || !kdp_kernel_state || m->endb - m->startb - 1 < SW_MIN_W || m->enda <= m->starta)
                return aln_seqseq_foward(m);
        return sw_sweep(m, 1);
}

static WAVE_TARGET int sw_bwd(struct aln_mem* m)
{
        if(!kalign_simd_on() || !kdp_kernel_state || m->endb - m->startb - 1 < SW_MIN_W || m->enda_2 <= m->starta_2)
                return aln_seqseq_backward(m);
        return sw_sweep(m, -1);
}

static const struct kdp_kernel ss_wave_kernel = {
        .name = "ss_wave8",
        .description = "SIMD (NEON/AVX2) inter-row wavefront seq-seq fwd/bwd, 8 rows (exact)",
        .fwd = { [KDP_SS] = sw_fwd },
        .bwd = { [KDP_SS] = sw_bwd },
        .prepare = sw_prepare,
        .release = sw_release,
};
KDP_REGISTER_KERNEL(ss_wave_kernel)

#endif
