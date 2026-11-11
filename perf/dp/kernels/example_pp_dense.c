/* example_pp_dense -- TEMPLATE prototype kernel (not an optimisation).
 *
 * Profile-profile forward/backward identical to lib/src/aln_profileprofile.c
 * except that the substitution score is a dense 23-term dot product
 * (prof1[0..22] . prof2[32..54]) instead of the sparse loop over the non-zero
 * residues of the prof1 column.  The summation order differs, so scores may
 * differ in the last float bits: this shows what `bench_dp --check` reports
 * for a numerically non-identical kernel.  The consistency bonus is added
 * exactly as in the reference (same row/column indices, also in the backward
 * pass).  Note: compiled without -ffast-math the 23-term float reduction is
 * a serial dependency chain, so this "dense" version is about 3x SLOWER than
 * the sparse reference.  Copy this file to start a new
 * prototype; all other kinds/ops fall back to "ref".
 */
#include <float.h>

#include "kdp_kernel.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MAX3(a, b, c) MAX(MAX(a, b), c)

static inline float dot23(const float* restrict p1, const float* restrict p2s)
{
        float s = 0.0F;
        int c;
        for(c = 0; c < 23; c++){
                s += p1[c] * p2s[c];
        }
        return s;
}

static int pp_dense_fwd(struct aln_mem* m)
{
        const float* prof1 = m->prof1;
        const float* prof2 = m->prof2;
        struct states* s = m->f;
        float pa, pga, pgb, ca, xa, xga;
        int i, j;

        prof1 += (m->starta) << 6;
        prof2 += (m->startb) << 6;
        s[m->startb].a = s[0].a;
        s[m->startb].ga = s[0].ga;
        s[m->startb].gb = s[0].gb;
        if(m->startb){
                for(j = m->startb + 1; j < m->endb; j++){
                        prof2 += 64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j - 1].ga + prof2[28], s[j - 1].a + prof2[27]);
                        s[j].gb = -FLT_MAX;
                }
                prof2 += 64;
        }else{
                for(j = m->startb + 1; j < m->endb; j++){
                        prof2 += 64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j - 1].ga, s[j - 1].a) + prof2[29];
                        s[j].gb = -FLT_MAX;
                }
                prof2 += 64;
        }
        prof2 -= (m->endb - m->startb) << 6;
        s[m->endb].a = -FLT_MAX;
        s[m->endb].ga = -FLT_MAX;
        s[m->endb].gb = -FLT_MAX;

        for(i = m->starta; i < m->enda; i++){
                prof1 += 64;
                pa = s[m->startb].a;
                pga = s[m->startb].ga;
                pgb = s[m->startb].gb;
                s[m->startb].a = -FLT_MAX;
                s[m->startb].ga = -FLT_MAX;
                xa = s[m->startb].a;
                xga = s[m->startb].ga;
                if(m->startb){
                        s[m->startb].gb = MAX(pgb + prof1[28], pa + prof1[27]);
                }else{
                        s[m->startb].gb = MAX(pgb, pa) + prof1[29];
                }
                for(j = m->startb + 1; j < m->endb; j++){
                        prof2 += 64;
                        ca = s[j].a;
                        pa = MAX3(pa, pga + prof2[-37], pgb + prof1[-37]);
                        pa += dot23(prof1, prof2 + 32);
                        if(m->consistency){
                                pa += sparse_bonus_lookup(m->consistency, i, j);
                        }
                        s[j].a = pa;
                        pga = s[j].ga;
                        s[j].ga = MAX(xga + prof2[28], xa + prof2[27]);
                        pgb = s[j].gb;
                        s[j].gb = MAX(pgb + prof1[28], ca + prof1[27]);
                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                prof2 += 64;
                ca = s[j].a;
                pa = MAX3(pa, pga + prof2[-37], pgb + prof1[-37]);
                pa += dot23(prof1, prof2 + 32);
                if(m->consistency){
                        pa += sparse_bonus_lookup(m->consistency, i, j);
                }
                s[j].a = pa;
                s[j].ga = -FLT_MAX;
                if(m->endb != m->len_b){
                        s[j].gb = MAX(s[j].gb + prof1[28], ca + prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb, ca) + prof1[29];
                }
                prof2 -= (m->endb - m->startb) << 6;
        }
        return 0;
}

static int pp_dense_bwd(struct aln_mem* m)
{
        struct states* s = m->b;
        const float* prof1 = m->prof1;
        const float* prof2 = m->prof2;
        float pa, pga, pgb, ca, xa, xga;
        int i, j;

        prof1 += (m->enda_2 + 1) << 6;
        prof2 += (m->endb + 1) << 6;
        s[m->endb].a = s[0].a;
        s[m->endb].ga = s[0].ga;
        s[m->endb].gb = s[0].gb;
        if(m->endb != m->len_b){
                for(j = m->endb - 1; j > m->startb; j--){
                        prof2 -= 64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j + 1].ga + prof2[28], s[j + 1].a + prof2[27]);
                        s[j].gb = -FLT_MAX;
                }
                prof2 -= 64;
        }else{
                for(j = m->endb - 1; j > m->startb; j--){
                        prof2 -= 64;
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j + 1].ga, s[j + 1].a) + prof2[29];
                        s[j].gb = -FLT_MAX;
                }
                prof2 -= 64;
        }
        s[m->startb].a = -FLT_MAX;
        s[m->startb].ga = -FLT_MAX;
        s[m->startb].gb = -FLT_MAX;

        i = m->enda_2 - m->starta_2;
        while(i--){
                prof1 -= 64;
                pa = s[m->endb].a;
                pga = s[m->endb].ga;
                pgb = s[m->endb].gb;
                s[m->endb].a = -FLT_MAX;
                s[m->endb].ga = -FLT_MAX;
                xa = s[m->endb].a;
                xga = s[m->endb].ga;
                if(m->endb != m->len_b){
                        s[m->endb].gb = MAX(pgb + prof1[28], pa + prof1[27]);
                }else{
                        s[m->endb].gb = MAX(pgb, pa) + prof1[29];
                }
                prof2 += (m->endb - m->startb) << 6;
                for(j = m->endb - 1; j > m->startb; j--){
                        prof2 -= 64;
                        ca = s[j].a;
                        pa = MAX3(pa, pga + prof2[91], pgb + prof1[91]);
                        pa += dot23(prof1, prof2 + 32);
                        if(m->consistency){
                                pa += sparse_bonus_lookup(m->consistency, m->starta_2 + i, j);
                        }
                        s[j].a = pa;
                        pga = s[j].ga;
                        s[j].ga = MAX(xga + prof2[28], xa + prof2[27]);
                        pgb = s[j].gb;
                        s[j].gb = MAX(pgb + prof1[28], ca + prof1[27]);
                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                prof2 -= 64;
                ca = s[j].a;
                pa = MAX3(pa, pga + prof2[91], pgb + prof1[91]);
                pa += dot23(prof1, prof2 + 32);
                if(m->consistency){
                        pa += sparse_bonus_lookup(m->consistency, m->starta_2 + i, j);
                }
                s[j].a = pa;
                s[j].ga = -FLT_MAX;
                if(m->startb){
                        s[j].gb = MAX(s[j].gb + prof1[28], ca + prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb, ca) + prof1[29];
                }
        }
        return 0;
}

static const struct kdp_kernel example_pp_dense = {
        .name = "example_pp_dense",
        .description = "TEMPLATE: PP fwd/bwd with dense 23-term dot product (rest = ref)",
        .fwd = { [KDP_PP] = pp_dense_fwd },
        .bwd = { [KDP_PP] = pp_dense_bwd },
};
KDP_REGISTER_KERNEL(example_pp_dense)
