/* sp_wave -- inter-row wavefront seq-profile forward and backward passes
 * (NEON or AVX2 via lib/src/wave_simd.h), bit-identical to
 * lib/src/aln_seqprofile.c.
 *
 * Same scheme as ss_wave.c.  Rows are the profile (position-specific gaps
 * GO1(r) = P1[r+1][27], GE1 = P1[r+1][28] forward / P1[r+2][28] backward,
 * T1 = P1[r+1][29]), columns the sequence (scalar gaps gpo/gpe/tgpe * sip).
 * Cell: D = max3(a, ga - open, gb + GO1(r)); a = D_diag + P1[r+1][32 + seq2[y]];
 * ga = max(ga - ext, a - open); gb = max(gb + GE1, a_up + GO1).
 * The substitution score is a per-lane gather (row table x column residue).
 * Kernel: sp_wave8.
 */
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "kdp_kernel.h"
#include "aln_seqprofile.h"
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

static inline __attribute__((always_inline)) int sw_sweep(struct aln_mem* m, const int dir)
{
        struct states* s = dir > 0 ? m->f : m->b;
        const uint8_t* seq1 = m->seq1;
        const int startb = m->startb;
        const int endb = m->endb;
        const int W = endb - startb - 1;
        const float* P1 = m->prof1;
        const uint8_t* seq2 = m->seq2;
        const float gpo = m->ap->gpo * m->sip;          /* "open" */
        const float gpe = m->ap->gpe * m->sip;          /* "ext"  */
        const float tgpe = m->ap->tgpe * m->sip;        /* "text" */
        const struct sparse_bonus* cons = m->consistency;
#define ROWP(rr) (P1 + (size_t)((rr) + 1) * 64)
#define GE1(rr) (dir > 0 ? P1[(size_t)((rr) + 1) * 64 + 28] : P1[(size_t)((rr) + 2) * 64 + 28])
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
        {
                const float go = ROWP(dir > 0 ? m->starta - 1 : m->enda_2)[27];
                for(t = -1; t < W; t++){
                        RD[t] = MAX3(RA[t], RGA[t] - gpo, RGB[t] + go);
                }
        }
        for(t = W; t < W + WL + 1; t++){
                RA[t] = -FLT_MAX; RGA[t] = -FLT_MAX; RGB[t] = -FLT_MAX; RD[t] = -FLT_MAX;
        }
        float ra = -FLT_MAX, rga = -FLT_MAX, rgb = -FLT_MAX;   /* right column */

        const vf8 vgpo = vf8_set1(gpo);
        const vf8 vgpe = vf8_set1(gpe);
        (void)seq1;
        const vf8 vneg = vf8_set1(-FLT_MAX);

        if(dir > 0){ r = m->starta; rend = m->enda; h = m->enda - m->starta; }
        else { r = m->enda_2 - 1; rend = m->starta_2 - 1; h = m->enda_2 - m->starta_2; }


        for(int hs = 0; h >= WAVE_PART_MIN; h -= hs, r += dir * hs){
                hs = h < WL ? h : WL;
                const float* p1q[WL];
                const float* bp[WL];
                float lgb[WL], dleft[WL], rbon[WL], dw[WL], go1[WL], ge1[WL];
                float dcap[WL][WL];

                /* left column of each strip row (scalar chain) */
                {
                        float pa0 = RA[-1], pgb0 = RGB[-1];
                        for(q = 0; q < WL; q++){
                                const int rq = r + dir * (q < hs ? q : hs - 1);
                                p1q[q] = ROWP(rq);
                                go1[q] = p1q[q][27];
                                ge1[q] = GE1(rq);
                                lgb[q] = lterm ? MAX(pgb0, pa0) + p1q[q][29] : MAX(pgb0 + ge1[q], pa0 + go1[q]);
                                dleft[q] = MAX3(-FLT_MAX, -FLT_MAX - gpo, lgb[q] + go1[q]);
                                pa0 = -FLT_MAX;
                                pgb0 = lgb[q];
                        }
                }
                const vf8 vgo1 = vf8_loadu(go1);
                const vf8 vge1 = vf8_loadu(ge1);
                for(q = 0; q < WL; q++){
                        const int rq = r + dir * (q < hs ? q : hs - 1);
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
                        float subl[WL];
                        for(int z = 0; z < WL; z++){
                                const int tt = sstep - z;
                                subl[z] = (tt >= 0 && tt < W) ? p1q[z][32 + seq2[y0 + dir * tt]] : 0.0F;
                        }
                        vf8 a = vf8_add(dg, vf8_loadu(subl));
                        if(cons){
                                a = vf8_add(a, vf8_gather(bp, sstep));
                        }
                        vf8 ga = vf8_max(vf8_sub(pGA, vgpe), vf8_sub(pA, vgpo));
                        vf8 gb = vf8_max(vf8_add(upGB, vge1), vf8_add(upA, vgo1));
                        vf8 d = vf8_max(vf8_max(a, vf8_sub(ga, vgpo)), vf8_add(gb, vgo1));
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
                        float p = (q ? dw[q - 1] : rdw_prev) + p1q[q][32 + seq2[yr]];
                        if(cons) p += rbon[q];
                        ra = p;
                        rga = -FLT_MAX;
                        rgb = rterm ? MAX(rgb, ca + go1[q]) + p1q[q][29] : MAX(rgb + ge1[q], ca + go1[q]);
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
                const float* p1 = ROWP(r);
                const float go = p1[27], ge = GE1(r);
                float pa0 = RA[-1], pgb0 = RGB[-1];
                float pa = RD[-1];
                float xa = -FLT_MAX, xga = -FLT_MAX;
                RA[-1] = -FLT_MAX;
                RGA[-1] = -FLT_MAX;
                RGB[-1] = lterm ? MAX(pgb0, pa0) + p1[29] : MAX(pgb0 + ge, pa0 + go);
                RD[-1] = MAX3(RA[-1], RGA[-1] - gpo, RGB[-1] + go);
                for(t = 0; t < W; t++){
                        const float ca = RA[t], cgb = RGB[t];
                        const float dnext = RD[t];
                        float a = pa + p1[32 + seq2[y0 + dir * t]];
                        if(cons) a += sparse_bonus_lookup(cons, r, y0 + dir * t);
                        RA[t] = a;
                        RGA[t] = MAX(xga - gpe, xa - gpo);
                        RGB[t] = MAX(cgb + ge, ca + go);
                        RD[t] = MAX3(RA[t], RGA[t] - gpo, RGB[t] + go);
                        xa = RA[t];
                        xga = RGA[t];
                        pa = dnext;
                }
                {
                        const float ca = ra;
                        float p = pa + p1[32 + seq2[yr]];
                        if(cons) p += sparse_bonus_lookup(cons, r, yr);
                        ra = p;
                        rga = -FLT_MAX;
                        rgb = rterm ? MAX(rgb, ca + go) + p1[29] : MAX(rgb + ge, ca + go);
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
#undef ROWP
#undef GE1
}

static int spw_fwd(struct aln_mem* m)
{
        if(m->endb - m->startb - 1 < SW_MIN_W || m->enda <= m->starta)
                return aln_seqprofile_foward(m);
        return sw_sweep(m, 1);
}

static int spw_bwd(struct aln_mem* m)
{
        if(m->endb - m->startb - 1 < SW_MIN_W || m->enda_2 <= m->starta_2)
                return aln_seqprofile_backward(m);
        return sw_sweep(m, -1);
}

static const struct kdp_kernel sp_wave_kernel = {
        .name = "sp_wave8",
        .description = "SIMD (NEON/AVX2) inter-row wavefront seq-profile fwd/bwd, 8 rows (exact)",
        .fwd = { [KDP_SP] = spw_fwd },
        .bwd = { [KDP_SP] = spw_bwd },
};
KDP_REGISTER_KERNEL(sp_wave_kernel)

#endif
