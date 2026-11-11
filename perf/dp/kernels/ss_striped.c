/* ss_striped -- PROTOTYPE: striped (Farrar-style) NEON seq-seq forward and
 * backward passes, bit-identical to lib/src/aln_seqseq.c.
 *
 * Why it can be exact: the match state of a row only reads the previous row,
 * and the row-gap state (ga, the in-row dependency) never feeds back into the
 * match state of the same row.  So a row is
 *   1. a = diag(max3(a, ga - gpo, gb - gpo)) + (subm[x][y] - soff) [+ bonus]
 *   2. gb = max(gb - gpe, a_prev - gpo)                (vertical, per column)
 *   3. ga = max(ga[t-1] - gpe, a[t-1] - gpo)           (prefix scan in t)
 * Steps 1-3 run in one striped loop; step 3 then gets Farrar's lazy
 * cross-lane fix-up.  max is exact and x - c is monotone, so the fixed-up
 * values are the same floats the scalar recurrence produces.
 *
 * Layout: L = 4*NV lanes, S = ceil(W/L) segments; logical t = lane*S + k is
 * stored at k*L + lane.  Forward t -> column startb+1+t, backward t ->
 * column endb-1-t.  The two edge columns (terminal gap rules) are scalar,
 * exactly as in the reference.
 *
 * Kernels: ss_striped4 / ss_striped8 / ss_striped16 (NV = 1, 2, 4).
 * Calls with fewer than SS_MIN_W interior columns or SS_MIN_H rows use the
 * reference code.
 */
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "kdp_kernel.h"
#include "aln_seqseq.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MAX3(a, b, c) MAX(MAX(a, b), c)

#if defined(__ARM_NEON)
#include <arm_neon.h>

#ifndef SS_MIN_W
#define SS_MIN_W 16
#endif
#ifndef SS_MIN_H
#define SS_MIN_H 1
#endif

struct ss_buf {
        float* A;
        float* GA;
        float* GB;
        float* QP;      /* 23 x [S*L] */
        float* BON;
        float* lin;
        size_t cap;     /* floats per array */
};

static struct ss_buf sb = {0};

#ifdef SS_STATS
#include <stdio.h>
static long long st_rows, st_main_k, st_fix_k, st_calls, st_cells, st_hist_cells[8];
__attribute__((destructor)) static void ss_stats_print(void)
{
        if(!st_calls) return;
        fprintf(stderr, "ss_stats calls %lld rows %lld cells %lld fixup/main k-steps %.3f\n",
                st_calls, st_rows, st_cells, (double)st_fix_k / (double)st_main_k);
        const char* lab[8] = {"h<4", "h<8", "h<16", "h<32", "h<64", "h<128", "h<256", "h>=256"};
        for(int i = 0; i < 8; i++) fprintf(stderr, "  %-7s %5.1f%% of cells\n", lab[i], 100.0 * st_hist_cells[i] / st_cells);
}
#define STAT(x) x
#else
#define STAT(x)
#endif

static void ss_reserve(size_t n)
{
        if(n <= sb.cap) return;
        free(sb.A); free(sb.GA); free(sb.GB); free(sb.QP); free(sb.BON); free(sb.lin);
        n = (n + 63) & ~(size_t)63;
        sb.A = aligned_alloc(64, n * sizeof(float));
        sb.GA = aligned_alloc(64, n * sizeof(float));
        sb.GB = aligned_alloc(64, n * sizeof(float));
        sb.QP = aligned_alloc(64, 23 * n * sizeof(float));
        sb.BON = aligned_alloc(64, n * sizeof(float));
        sb.lin = aligned_alloc(64, n * sizeof(float));
        if(!sb.A || !sb.GA || !sb.GB || !sb.QP || !sb.BON || !sb.lin) abort();
        memset(sb.BON, 0, n * sizeof(float));
        sb.cap = n;
}

static inline float lane_get(float32x4_t v, int l)
{
        switch(l){
        case 0: return vgetq_lane_f32(v, 0);
        case 1: return vgetq_lane_f32(v, 1);
        case 2: return vgetq_lane_f32(v, 2);
        default: return vgetq_lane_f32(v, 3);
        }
}

/* dir = +1 forward (m->f, rows starta..enda-1), -1 backward (m->b,
   rows enda_2-1 .. starta_2). */
static inline __attribute__((always_inline)) int ss_sweep(struct aln_mem* m, const int dir, const int NV)
{
        struct states* s = dir > 0 ? m->f : m->b;
        const uint8_t* seq1 = m->seq1;
        const uint8_t* seq2 = m->seq2;
        const int startb = m->startb;
        const int endb = m->endb;
        const int L = 4 * NV;
        const int W = endb - startb - 1;
        const int S = (W + L - 1) / L;
        const int SL = S * L;
        const float gpo = m->ap->gpo;
        const float gpe = m->ap->gpe;
        const float tgpe = m->ap->tgpe;
        float** subm = m->ap->subm;
        const float soff = m->ap->subm_offset;
        const struct sparse_bonus* cons = m->consistency;
        const int lcol = dir > 0 ? startb : endb;
        const int rcol = dir > 0 ? endb : startb;
        const int lterm = dir > 0 ? (startb == 0) : (endb == m->len_b);
        const int rterm = dir > 0 ? (endb == m->len_b) : (startb == 0);
        const int y0 = dir > 0 ? startb : endb - 1;     /* y(t) = y0 + dir*t */
        const int yr = dir > 0 ? endb - 1 : startb;
        uint32_t built = 0;
        int r, rend, t, k, l, v, pass;

        ss_reserve((size_t)SL);
        STAT({ int hh = dir > 0 ? m->enda - m->starta : m->enda_2 - m->starta_2; int b = 0;
               while(b < 7 && hh >= (4 << b)) b++;
               st_calls++; st_cells += (long long)hh * W; st_hist_cells[b] += (long long)hh * W; st_rows += hh; });
        float* restrict A = sb.A;
        float* restrict GA = sb.GA;
        float* restrict GB = sb.GB;
        float* restrict BON = sb.BON;

        /* ---- first row (identical to the reference init) ---- */
        s[lcol].a = s[0].a;
        s[lcol].ga = s[0].ga;
        s[lcol].gb = s[0].gb;
        {
                float pga = s[lcol].ga, pa = s[lcol].a;
                for(t = 0; t < W; t++){
                        float g = lterm ? MAX(pga, pa) - tgpe : MAX(pga - gpe, pa - gpo);
                        sb.lin[t] = g;
                        pga = g;
                        pa = -FLT_MAX;
                }
        }
        for(l = 0, t = 0; l < L; l++){
                for(k = 0; k < S; k++, t++){
                        A[k * L + l] = -FLT_MAX;
                        GB[k * L + l] = -FLT_MAX;
                        GA[k * L + l] = t < W ? sb.lin[t] : -FLT_MAX;
                }
        }
        float ra = -FLT_MAX, rga = -FLT_MAX, rgb = -FLT_MAX;
        float la = s[lcol].a, lga = s[lcol].ga, lgb = s[lcol].gb;

        const float32x4_t vgpo = vdupq_n_f32(gpo);
        const float32x4_t vgpe = vdupq_n_f32(gpe);
        const float32x4_t vneg = vdupq_n_f32(-FLT_MAX);
        const int tl = W - 1, kl = tl % S, ll = tl / S;
        const int vl = ll / 4, ql = ll % 4;

        if(dir > 0){ r = m->starta; rend = m->enda; }
        else { r = m->enda_2 - 1; rend = m->starta_2 - 1; }

        for(; r != rend; r += dir){
                const int x = seq1[r];
                float* restrict qp = sb.QP + (size_t)x * SL;
                if(!(built & (1u << x))){
                        const float* subp = subm[x];
                        for(l = 0, t = 0; l < L; l++){
                                for(k = 0; k < S; k++, t++){
                                        qp[k * L + l] = t < W ? subp[seq2[y0 + dir * t]] - soff : 0.0F;
                                }
                        }
                        built |= 1u << x;
                }
                float rbon = 0.0F;
                int bbase = 0, bn = 0;
                if(cons){
                        bbase = r * cons->K;
                        for(bn = 0; bn < cons->K; bn++){
                                int c = cons->cols[bbase + bn];
                                if(c < 0) break;
                                if(c == yr){
                                        rbon += cons->vals[bbase + bn];
                                }else{
                                        int tt = (c - y0) * dir;
                                        if(tt >= 0 && tt < W){
                                                BON[(tt % S) * L + tt / S] += cons->vals[bbase + bn];
                                        }
                                }
                        }
                }

                /* left boundary column */
                const float pa0 = la, pga0 = lga, pgb0 = lgb;
                la = -FLT_MAX;
                lga = -FLT_MAX;
                lgb = lterm ? MAX(pgb0, pa0) - tgpe : MAX(pgb0 - gpe, pa0 - gpo);
                const float dleft = MAX3(pa0, pga0 - gpo, pgb0 - gpo);

                float32x4_t diag[4], pga[4], pa[4], dkl = vneg;
                {
                        float32x4_t dlast[4];
                        const int o = (S - 1) * L;
                        for(v = 0; v < NV; v++){
                                dlast[v] = vmaxq_f32(vmaxq_f32(vld1q_f32(A + o + 4 * v),
                                                               vsubq_f32(vld1q_f32(GA + o + 4 * v), vgpo)),
                                                     vsubq_f32(vld1q_f32(GB + o + 4 * v), vgpo));
                        }
                        for(v = 0; v < NV; v++){
                                diag[v] = vextq_f32(v ? dlast[v - 1] : vdupq_n_f32(dleft), dlast[v], 3);
                                pga[v] = vneg;
                                pa[v] = vneg;
                        }
                        pga[0] = vsetq_lane_f32(lga, pga[0], 0);
                        pa[0] = vsetq_lane_f32(la, pa[0], 0);
                }

                /* steps 1-3, fused */
                STAT(st_main_k += S;)
                for(k = 0; k < S; k++){
                        const int o = k * L;
                        for(v = 0; v < NV; v++){
                                float32x4_t a = vld1q_f32(A + o + 4 * v);
                                float32x4_t ga = vld1q_f32(GA + o + 4 * v);
                                float32x4_t gb = vld1q_f32(GB + o + 4 * v);
                                float32x4_t d = vmaxq_f32(vmaxq_f32(a, vsubq_f32(ga, vgpo)), vsubq_f32(gb, vgpo));
                                float32x4_t na = vaddq_f32(diag[v], vld1q_f32(qp + o + 4 * v));
                                if(cons) na = vaddq_f32(na, vld1q_f32(BON + o + 4 * v));
                                float32x4_t g = vmaxq_f32(vsubq_f32(pga[v], vgpe), vsubq_f32(pa[v], vgpo));
                                vst1q_f32(A + o + 4 * v, na);
                                vst1q_f32(GB + o + 4 * v, vmaxq_f32(vsubq_f32(gb, vgpe), vsubq_f32(a, vgpo)));
                                vst1q_f32(GA + o + 4 * v, g);
                                if(k == kl && v == vl) dkl = d;
                                diag[v] = d;
                                pga[v] = g;
                                pa[v] = na;
                        }
                }
                const float dright = lane_get(dkl, ql);

                /* lazy fix-up of ga across lanes */
                {
                        float32x4_t c[4];
                        for(v = 0; v < NV; v++){
                                c[v] = vmaxq_f32(vsubq_f32(vextq_f32(v ? pga[v - 1] : vneg, pga[v], 3), vgpe),
                                                 vsubq_f32(vextq_f32(v ? pa[v - 1] : vneg, pa[v], 3), vgpo));
                        }
                        for(pass = 0; pass < L - 1; pass++){
                                for(k = 0; k < S; k++){
                                        const int o = k * L;
                                        uint32x4_t ch = vdupq_n_u32(0);
                                        STAT(st_fix_k++;)
                                        for(v = 0; v < NV; v++){
                                                float32x4_t g = vld1q_f32(GA + o + 4 * v);
                                                float32x4_t n = vmaxq_f32(g, c[v]);
                                                ch = vorrq_u32(ch, vcgtq_f32(n, g));
                                                vst1q_f32(GA + o + 4 * v, n);
                                                c[v] = vsubq_f32(n, vgpe);
                                        }
                                        if(vmaxvq_u32(ch) == 0) goto ga_done;
                                }
                                {
                                        float32x4_t prev = vneg;
                                        for(v = 0; v < NV; v++){
                                                float32x4_t cur = c[v];
                                                c[v] = vextq_f32(prev, cur, 3);
                                                prev = cur;
                                        }
                                }
                        }
ga_done:;
                }

                /* right column (scalar, as the reference) */
                {
                        const float ca = ra;
                        float p = dright + (subm[x][seq2[yr]] - soff);
                        if(cons) p += rbon;
                        ra = p;
                        rga = -FLT_MAX;
                        rgb = rterm ? MAX(rgb, ca - gpo) - tgpe : MAX(rgb - gpe, ca - gpo);
                }

                if(cons){
                        for(int q = 0; q < bn; q++){
                                int c = cons->cols[bbase + q];
                                int tt = (c - y0) * dir;
                                if(c != yr && tt >= 0 && tt < W){
                                        BON[(tt % S) * L + tt / S] = 0.0F;
                                }
                        }
                }
        }

        /* write back */
        s[lcol].a = la;
        s[lcol].ga = lga;
        s[lcol].gb = lgb;
        for(l = 0, t = 0; l < L; l++){
                for(k = 0; k < S && t < W; k++, t++){
                        int c = lcol + dir * (t + 1);
                        s[c].a = A[k * L + l];
                        s[c].ga = GA[k * L + l];
                        s[c].gb = GB[k * L + l];
                }
        }
        s[rcol].a = ra;
        s[rcol].ga = rga;
        s[rcol].gb = rgb;
        return 0;
}

static inline int ss_small_f(struct aln_mem* m)
{
        return m->endb - m->startb - 1 < SS_MIN_W || m->enda - m->starta < SS_MIN_H;
}

static inline int ss_small_b(struct aln_mem* m)
{
        return m->endb - m->startb - 1 < SS_MIN_W || m->enda_2 - m->starta_2 < SS_MIN_H;
}

#define SS_KERNEL(NVV, LANES)                                                   \
        static int ss_fwd_##LANES(struct aln_mem* m)                            \
        {                                                                       \
                if(ss_small_f(m)) return aln_seqseq_foward(m);                  \
                return ss_sweep(m, 1, NVV);                                     \
        }                                                                       \
        static int ss_bwd_##LANES(struct aln_mem* m)                            \
        {                                                                       \
                if(ss_small_b(m)) return aln_seqseq_backward(m);                \
                return ss_sweep(m, -1, NVV);                                    \
        }                                                                       \
        static const struct kdp_kernel ss_striped_kernel_##LANES = {            \
                .name = "ss_striped" #LANES,                                    \
                .description = "NEON striped seq-seq fwd/bwd, " #LANES " lanes (exact)", \
                .fwd = { [KDP_SS] = ss_fwd_##LANES },                           \
                .bwd = { [KDP_SS] = ss_bwd_##LANES },                           \
        };                                                                      \
        KDP_REGISTER_KERNEL(ss_striped_kernel_##LANES)

SS_KERNEL(1, 4)
SS_KERNEL(2, 8)
SS_KERNEL(4, 16)

#endif
