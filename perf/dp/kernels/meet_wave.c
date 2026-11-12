/* meet_wave -- SIMD-filtered Hirschberg meetup (NEON or AVX2 via
 * lib/src/wave_simd.h), bit-identical to aln_{seqseq,seqprofile,
 * profileprofile}_meetup.
 *
 * The meetup scans the split row and keeps the best and second-best of six
 * candidate scores per column (transition codes 1, 2, 3, 5, 6, 7, in that
 * order, strict > comparisons).  The candidates of 8 columns are computed at
 * once with the reference expressions (same operations, same order, so the
 * same floats).  If none of them exceeds the current second-best, the scan
 * cannot change anything for those columns and skips them; otherwise the
 * stored candidates go through the reference update in the reference order.
 * Everything after the scan (last column, margins, flips) is the reference
 * code.
 *
 * Kernel: meet_wave (meetup only; passes fall back to ref).
 *
 * RESULT (negative, 2026-10-02): exact on all 216 dumps, but no measurable
 * gain on M3 (full-alignment ns/cell within noise of ref): assembling the six
 * AoS fields per column costs about what the scalar arithmetic does.  Not
 * used in the library.
 */
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "kdp_kernel.h"
#include "aln_seqseq.h"
#include "aln_seqprofile.h"
#include "aln_profileprofile.h"
#include "wave_simd.h"

#if WAVE_SIMD

/* the reference top-2 update for one candidate */
#define MEET_UPD(val, code, col)                                                \
        do {                                                                    \
                const float s_tmp_ = (val);                                     \
                if(s_tmp_ > max){                                               \
                        max2 = max; c2 = c; transition2 = transition;           \
                        max = s_tmp_;                                           \
                        transition = (code);                                    \
                        c = (col);                                              \
                }else if(s_tmp_ > max2){                                        \
                        max2 = s_tmp_; c2 = (col); transition2 = (code);        \
                }                                                               \
        } while(0)

static inline WAVE_TARGET vf8 vf8_hmax_bcast(vf8 v)
{
        float t[8];
        vf8_storeu(t, v);
        float x = t[0];
        for(int k = 1; k < 8; k++) x = t[k] > x ? t[k] : x;
        return vf8_set1(x);
}

/* kind: KDP_SS, KDP_SP, KDP_PP */
static WAVE_TARGET int meet_wave(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score, const int kind)
{
        struct states* f = m->f;
        struct states* b = m->b;
        const float gpo = m->ap->gpo;
        const float gpe = m->ap->gpe;
        const float open = m->ap->gpo * m->sip;        /* sp */
        const float* prof1 = NULL;
        const float* prof2 = NULL;
        int i;
        int c = -1;
        int c2 = -1;
        int transition = -1;
        int transition2 = -1;
        float max = -FLT_MAX;
        float max2 = -FLT_MAX;
        const float middle = (float)(old_cor[3] - old_cor[2]) / 2.0F + (float)old_cor[2];
        float sub;

        if(kind != KDP_SS){
                prof1 = m->prof1 + ((old_cor[4] + 1) << 6);
        }
        if(kind == KDP_PP){
                prof2 = m->prof2;       /* column i uses prof2 + (i+1)*64 */
        }

        const vf8 vmid = vf8_set1(middle);
        const vf8 v1000 = vf8_set1(1000.0F);
        const vf8 vgpo = vf8_set1(gpo);
        const vf8 vgpe = vf8_set1(gpe);
        const vf8 vopen = vf8_set1(open);

        float cand[6][8];
        i = old_cor[2];
        for(; i + 8 <= old_cor[3]; i += 8){
                float fa[8], fga[8], fgb[8], ba[8], bga[8], bgb[8], fi[8];
                for(int k = 0; k < 8; k++){
                        fa[k] = f[i + k].a; fga[k] = f[i + k].ga; fgb[k] = f[i + k].gb;
                        ba[k] = b[i + k].a; bga[k] = b[i + k].ga; bgb[k] = b[i + k].gb;
                        fi[k] = (float)(i + k);
                }
                const vf8 vfa = vf8_loadu(fa), vfga = vf8_loadu(fga), vfgb = vf8_loadu(fgb);
                const vf8 vba = vf8_loadu(ba), vbga = vf8_loadu(bga), vbgb = vf8_loadu(bgb);
                /* sub = fabsf(middle - i) / 1000 */
                vf8 d = vf8_sub(vmid, vf8_loadu(fi));
                {
                        float tt[8];
                        vf8_storeu(tt, d);
                        for(int k = 0; k < 8; k++) tt[k] = fabsf(tt[k]);
                        d = vf8_loadu(tt);
                }
                const vf8 vsub = vf8_div(d, v1000);
                vf8 c1, c2v, c3, c5, c6, c7;
                c1 = vf8_sub(vf8_add(vfa, vba), vsub);
                if(kind == KDP_SS){
                        c2v = vf8_sub(vf8_sub(vf8_add(vfa, vbga), vgpo), vsub);
                        c3 = vf8_sub(vf8_sub(vf8_add(vfa, vbgb), vgpo), vsub);
                        c5 = vf8_sub(vf8_sub(vf8_add(vfga, vba), vgpo), vsub);
                        c6 = vf8_sub(vf8_sub(vf8_add(vfgb, vbgb), vgpe), vsub);
                        c7 = vf8_sub(vf8_sub(vf8_add(vfgb, vba), vgpo), vsub);
                }else if(kind == KDP_SP){
                        c2v = vf8_sub(vf8_sub(vf8_add(vfa, vbga), vopen), vsub);
                        c3 = vf8_sub(vf8_add(vf8_add(vfa, vbgb), vf8_set1(prof1[27])), vsub);
                        c5 = vf8_sub(vf8_sub(vf8_add(vfga, vba), vopen), vsub);
                        c6 = vf8_sub(vf8_add(vf8_add(vfgb, vbgb), vf8_set1(prof1[28])), vsub);
                        c7 = vf8_sub(vf8_add(vf8_add(vfgb, vba), vf8_set1(prof1[-37])), vsub);
                }else{
                        float p27[8], pm37[8];
                        for(int k = 0; k < 8; k++){
                                const float* p2 = prof2 + ((size_t)(i + k + 1) << 6);
                                p27[k] = p2[27];
                                pm37[k] = p2[-37];
                        }
                        c2v = vf8_sub(vf8_add(vf8_add(vfa, vbga), vf8_loadu(p27)), vsub);
                        c3 = vf8_sub(vf8_add(vf8_add(vfa, vbgb), vf8_set1(prof1[27])), vsub);
                        c5 = vf8_sub(vf8_add(vf8_add(vfga, vba), vf8_loadu(pm37)), vsub);
                        c6 = vf8_sub(vf8_add(vf8_add(vfgb, vbgb), vf8_set1(prof1[28])), vsub);
                        c7 = vf8_sub(vf8_add(vf8_add(vfgb, vba), vf8_set1(prof1[-37])), vsub);
                }
                vf8_storeu(cand[0], c1);
                vf8_storeu(cand[1], c2v);
                vf8_storeu(cand[2], c3);
                vf8_storeu(cand[3], c5);
                vf8_storeu(cand[4], c6);
                vf8_storeu(cand[5], c7);
                if(i == 0){
                        /* gb -> gb join in column 0: no extension term */
                        cand[4][0] = f[0].gb + b[0].gb - fabsf(middle - 0.0F) / 1000.0F;
                }
                const vf8 vm = vf8_max(vf8_max(vf8_max(c1, c2v), vf8_max(c3, c5)), vf8_max(c6, c7));
                float bm[8];
                vf8_storeu(bm, vf8_hmax_bcast(vm));
                if(i == 0 && cand[4][0] > bm[0]) bm[0] = cand[4][0];
                if(!(bm[0] > max2)) continue;   /* nothing here can change max or max2 */
                for(int k = 0; k < 8; k++){
                        MEET_UPD(cand[0][k], 1, i + k);
                        MEET_UPD(cand[1][k], 2, i + k);
                        MEET_UPD(cand[2][k], 3, i + k);
                        MEET_UPD(cand[3][k], 5, i + k);
                        MEET_UPD(cand[4][k], 6, i + k);
                        MEET_UPD(cand[5][k], 7, i + k);
                }
        }
        /* remaining columns: reference code */
        for(; i < old_cor[3]; i++){
                sub = fabsf(middle - (float)i);
                sub /= 1000.0F;
                if(kind == KDP_SS){
                        MEET_UPD(f[i].a + b[i].a - sub, 1, i);
                        MEET_UPD(f[i].a + b[i].ga - gpo - sub, 2, i);
                        MEET_UPD(f[i].a + b[i].gb - gpo - sub, 3, i);
                        MEET_UPD(f[i].ga + b[i].a - gpo - sub, 5, i);
                        if(i == 0){
                                MEET_UPD(f[i].gb + b[i].gb - sub, 6, i);
                        }else{
                                MEET_UPD(f[i].gb + b[i].gb - gpe - sub, 6, i);
                        }
                        MEET_UPD(f[i].gb + b[i].a - gpo - sub, 7, i);
                }else if(kind == KDP_SP){
                        MEET_UPD(f[i].a + b[i].a - sub, 1, i);
                        MEET_UPD(f[i].a + b[i].ga - open - sub, 2, i);
                        MEET_UPD(f[i].a + b[i].gb + prof1[27] - sub, 3, i);
                        MEET_UPD(f[i].ga + b[i].a - open - sub, 5, i);
                        if(i == 0){
                                MEET_UPD(f[i].gb + b[i].gb - sub, 6, i);
                        }else{
                                MEET_UPD(f[i].gb + b[i].gb + prof1[28] - sub, 6, i);
                        }
                        MEET_UPD(f[i].gb + b[i].a + prof1[-37] - sub, 7, i);
                }else{
                        const float* p2 = prof2 + ((size_t)(i + 1) << 6);
                        MEET_UPD(f[i].a + b[i].a - sub, 1, i);
                        MEET_UPD(f[i].a + b[i].ga + p2[27] - sub, 2, i);
                        MEET_UPD(f[i].a + b[i].gb + prof1[27] - sub, 3, i);
                        MEET_UPD(f[i].ga + b[i].a + p2[-37] - sub, 5, i);
                        if(i == 0){
                                MEET_UPD(f[i].gb + b[i].gb - sub, 6, i);
                        }else{
                                MEET_UPD(f[i].gb + b[i].gb + prof1[28] - sub, 6, i);
                        }
                        MEET_UPD(f[i].gb + b[i].a + prof1[-37] - sub, 7, i);
                }
        }
        /* last column (reference) */
        i = old_cor[3];
        sub = fabsf(middle - (float)i);
        sub /= 1000.0F;
        if(kind == KDP_SS){
                MEET_UPD(f[i].a + b[i].gb - gpo - sub, 3, i);
                if(m->endb == m->len_b){
                        MEET_UPD(f[i].gb + b[i].gb - sub, 6, i);
                }else{
                        MEET_UPD(f[i].gb + b[i].gb - gpe - sub, 6, i);
                }
        }else{
                MEET_UPD(f[i].a + b[i].gb + prof1[27] - sub, 3, i);
                if(m->endb == m->len_b){
                        MEET_UPD(f[i].gb + b[i].gb - sub, 6, i);
                }else{
                        MEET_UPD(f[i].gb + b[i].gb + prof1[28] - sub, 6, i);
                }
        }

        /* margins and flips: reference code */
        if(max2 > -FLT_MAX){
                float _margin = max - max2;
                if(m->flip_margins != NULL && m->margin_count < m->flip_margin_alloc){
                        m->flip_margins[m->margin_count] = _margin;
                }
                m->margin_sum += _margin;
                m->margin_count++;
        }
        if(m->flip_threshold > 0.0F && c2 >= 0 && max2 > -FLT_MAX){
                float margin = max - max2;
                if(margin < m->flip_threshold){
                        if(m->flip_bit_map != NULL){
                                if(m->flip_counter < m->flip_n_uncertain){
                                        int bit = m->flip_bit_map[m->flip_counter];
                                        if(bit >= 0 && ((1U << bit) & m->flip_mask)){
                                                c = c2;
                                                transition = transition2;
                                        }
                                }
                        }else if(m->flip_mask != 0){
                                if((1U << (m->flip_counter % m->flip_stride)) & m->flip_mask){
                                        c = c2;
                                        transition = transition2;
                                }
                        }else if(m->flip_trial > 0){
                                if(m->flip_counter % m->flip_stride == m->flip_trial - 1){
                                        c = c2;
                                        transition = transition2;
                                }
                        }
                        m->flip_counter++;
                }
        }
        *meet = c;
        *t = transition;
        *score = max;
        return 0;
}

static WAVE_TARGET int meet_ss(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score)
{
        if(!kalign_simd_on()) return aln_seqseq_meetup(m, old_cor, meet, t, score);
        return meet_wave(m, old_cor, meet, t, score, KDP_SS);
}
static WAVE_TARGET int meet_sp(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score)
{
        if(!kalign_simd_on()) return aln_seqprofile_meetup(m, old_cor, meet, t, score);
        return meet_wave(m, old_cor, meet, t, score, KDP_SP);
}
static WAVE_TARGET int meet_pp(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score)
{
        if(!kalign_simd_on()) return aln_profileprofile_meetup(m, old_cor, meet, t, score);
        return meet_wave(m, old_cor, meet, t, score, KDP_PP);
}

static const struct kdp_kernel meet_wave_kernel = {
        .name = "meet_wave",
        .description = "SIMD-filtered Hirschberg meetup, all kinds (exact); passes = ref",
        .meet = { meet_ss, meet_sp, meet_pp },
};
KDP_REGISTER_KERNEL(meet_wave_kernel)

#endif
