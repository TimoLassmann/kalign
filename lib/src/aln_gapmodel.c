/* Per-member gap model (run option gap_model = 1; with gap_model = 0 none of this code runs).

   The model. When the DP inserts a column of gaps into profile X at boundary k (k columns of X to the left), every
   member of X gets a gap character there, and what it costs depends on that member:
       '=' opens a gap    residues on both sides, no gap next to the boundary  -> gpo to open, gpo to close, gpe for
                          every further column (Kalign's open/close convention)
       '-' extends a gap  an internal gap immediately left or right of the boundary   -> gpe per column
       '.' terminal       no residue before (not started) or after (ended)            -> tgpe per column
   with n= + n- + n. = all members of X. The cost of one gap column is multiplied, as in Kalign's sum-of-pairs
   scoring, by the number of residues in the opposite column. Written per column so that forward and backward DP agree:
       every gap column, when entered      res_opposite * W(k)   W  = (n= + n-) gpe + n. tgpe
       first column of a run, extra        res_opposite * Xo(k)  Xo = n= (gpo - gpe)
       last column of a run, extra         res_opposite * Xc(k)  Xc = n= gpo
   At the two ends of a profile every member is terminal, so Xo = Xc = 0 and W = members x tgpe, which is Kalign's
   existing terminal rule (without the extra gap-open charge Kalign adds when a terminal gap ends or begins).

   The counts per profile column (slots GM_* of the 64-float column) are kept by gm_leaf_counts and update_n_gm:
       n.(k) = lead(k) + trail(k+1)          lead(0) = trail(len+1) = all members
       n-(k) = gopen(k) + gext(k) + gopen(k+1)
       n=(k) = members - n. - n-
   The DP functions below are Kalign's profile-profile and sequence-profile kernels with identical loop structure
   (cells, states, transitions, Hirschberg split), only the gap terms replaced. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tldevel.h"
#include "kalign_simd.h"
#include "tl_scratch.h"
#include "aln_param.h"
#include "aln_struct.h"
#include "anchor_consistency.h"
#include "msa_struct.h"
#include "task.h"
#include "aln_gapmodel.h"

#define GM_MAX(a, b) ((a) > (b) ? (a) : (b))
#define GM_MAX3(a, b, c) GM_MAX(GM_MAX(a, b), c)

struct gm_side {
        int len;
        float nsip;
        float* nres;   /* [0..len+1], residues per column (0 at the two boundary columns) */
        float* w;      /* [0..len] per-column cost of a gap column at boundary k */
        float* xo;     /* [0..len] extra for the first column of a run */
        float* xc;     /* [0..len] extra for the last column of a run */
};

struct gm {
        struct gm_side A;   /* prof1 */
        struct gm_side B;   /* prof2, or the sequence seq2 */
        void* wv;           /* wavefront SIMD column tables (gm_wave_build), NULL without SIMD */
};

static void* gm_wave_build(const struct gm* g, const struct aln_mem* m);
static void gm_wave_release(void* wv);

/* ------------------------------------------------------------------ weights */

static int side_alloc(struct gm_side* s, int len)
{
        s->len = len;
        s->nres = calloc((size_t)len + 2, sizeof(float));
        s->w = calloc((size_t)len + 1, sizeof(float));
        s->xo = calloc((size_t)len + 1, sizeof(float));
        s->xc = calloc((size_t)len + 1, sizeof(float));
        return s->nres && s->w && s->xo && s->xc;
}

static void side_free(struct gm_side* s)
{
        free(s->nres); free(s->w); free(s->xo); free(s->xc);
        s->nres = s->w = s->xo = s->xc = NULL;
}

static void weights(float nR, float nG, float nT, const struct aln_param* ap, float* w, float* xo, float* xc)
{
        *w = (nR + nG) * ap->gpe + nT * ap->tgpe;
        *xo = nR * (ap->gpo - ap->gpe);
        *xc = nR * ap->gpo;
}

/* boundary k of a profile with counts: k columns to its left; P(c, slot) reads column c (1..len) */
#define P(c, slot) (prof[((c) << 6) + (slot)])

static void boundary_counts(const float* prof, int len, float nsip, int k, float* nR, float* nG, float* nT)
{
        float lead = (k == 0) ? nsip : P(k, GM_LEAD);
        float trail = (k == len) ? nsip : P(k + 1, GM_TRAIL);
        float g = (k == 0) ? 0.0f : P(k, GM_GOPEN) + P(k, GM_GEXT);
        float gop = (k == len) ? 0.0f : P(k + 1, GM_GOPEN);
        *nT = lead + trail;
        *nG = g + gop;
        *nR = nsip - *nT - *nG;
}

static float profile_nsip(const float* prof)
{
        return P(1, GM_NRES) + P(1, GM_LEAD) + P(1, GM_TRAIL) + P(1, GM_GOPEN) + P(1, GM_GEXT);
}

static int side_from_profile(struct gm_side* s, const float* prof, int len, const struct aln_param* ap)
{
        float nsip;
        if(len < 1) return 0;
        nsip = profile_nsip(prof);
        if(nsip < 0.5f) return 0;   /* profile without counts (e.g. built by aln_add.c): keep the default model */
        if(!side_alloc(s, len)) return 0;
        s->nsip = nsip;
        for(int c = 1; c <= len; c++){
                s->nres[c] = P(c, GM_NRES);
        }
        for(int k = 0; k <= len; k++){
                float nR, nG, nT;
                boundary_counts(prof, len, nsip, k, &nR, &nG, &nT);
                weights(nR, nG, nT, ap, &s->w[k], &s->xo[k], &s->xc[k]);
        }
        return 1;
}

static int side_from_sequence(struct gm_side* s, int len, const struct aln_param* ap)
{
        if(len < 1) return 0;
        if(!side_alloc(s, len)) return 0;
        s->nsip = 1.0f;
        for(int c = 1; c <= len; c++){
                s->nres[c] = 1.0f;
        }
        for(int k = 0; k <= len; k++){
                if(k == 0 || k == len){
                        weights(0.0f, 0.0f, 1.0f, ap, &s->w[k], &s->xo[k], &s->xc[k]);
                }else{
                        weights(1.0f, 0.0f, 0.0f, ap, &s->w[k], &s->xo[k], &s->xc[k]);
                }
        }
        return 1;
}

int gm_prepare(struct aln_mem* m)
{
        struct gm* g = NULL;
        if(m->gm || m->seq1 || !m->prof1 || !m->ap || !m->ap->gap_model) return 0;
        g = calloc(1, sizeof(struct gm));
        if(!g) return 0;
        if(!side_from_profile(&g->A, m->prof1, m->len_a, m->ap)){
                free(g);
                return 0;
        }
        if(m->prof2){
                if(!side_from_profile(&g->B, m->prof2, m->len_b, m->ap)){
                        side_free(&g->A); free(g);
                        return 0;
                }
        }else{
                if(!m->seq2 || !side_from_sequence(&g->B, m->len_b, m->ap)){
                        side_free(&g->A); free(g);
                        return 0;
                }
        }
        g->wv = gm_wave_build(g, m);
        m->gm = g;
        return 1;
}

void gm_free(struct aln_mem* m)
{
        struct gm* g = m->gm;
        if(!g) return;
        gm_wave_release(g->wv);
        side_free(&g->A);
        side_free(&g->B);
        free(g);
        m->gm = NULL;
}

/* ------------------------------------------------------------------ DP: profile-profile */

static int gm_pp_forward_scalar(struct aln_mem* m)
{
        const struct gm* g = m->gm;
        const float* nA = g->A.nres; const float* wA = g->A.w; const float* xoA = g->A.xo; const float* xcA = g->A.xc;
        const float* nB = g->B.nres; const float* wB = g->B.w; const float* xoB = g->B.xo; const float* xcB = g->B.xc;
        unsigned int freq[24];
        const float* prof1 = m->prof1;
        const float* prof2 = m->prof2;
        struct states* s = m->f;
        float pa, pga, pgb, ca, xa, xga;
        int i, j, c, f, col;
        const int k0 = m->starta;

        prof1 += (m->starta) << 6;
        prof2 += (m->startb) << 6;
        s[m->startb].a = s[0].a;
        s[m->startb].ga = s[0].ga;
        s[m->startb].gb = s[0].gb;
        /* the row before A column starta+1: gaps in A at boundary starta */
        for(j = m->startb + 1; j < m->endb; j++){
                s[j].a = -FLT_MAX;
                s[j].ga = GM_MAX(s[j-1].ga - nB[j] * wA[k0], s[j-1].a - nB[j] * (wA[k0] + xoA[k0]));
                s[j].gb = -FLT_MAX;
        }
        s[m->endb].a = -FLT_MAX;
        s[m->endb].ga = -FLT_MAX;
        s[m->endb].gb = -FLT_MAX;

        for(i = m->starta; i < m->enda; i++){
                prof1 += 64;
                col = i + 1;    /* A column of this row; gaps in A in this row sit at boundary col */
                f = 0;
                for(j = 0; j < 23; j++){
                        if(prof1[j]){
                                freq[f] = j;
                                f++;
                        }
                }
                f--;

                pa = s[m->startb].a;
                pga = s[m->startb].ga;
                pgb = s[m->startb].gb;
                s[m->startb].a = -FLT_MAX;
                s[m->startb].ga = -FLT_MAX;

                xa = s[m->startb].a;
                xga = s[m->startb].ga;

                s[m->startb].gb = GM_MAX(pgb - nA[col] * wB[m->startb],
                                         pa - nA[col] * (wB[m->startb] + xoB[m->startb]));

                for(j = m->startb + 1; j < m->endb; j++){
                        prof2 += 64;
                        ca = s[j].a;

                        /* match (col, j) after: match; a gap run in A (boundary col-1) ending at B column j-1; a gap
                           run in B (boundary j-1) ending at A column col-1 */
                        pa = GM_MAX3(pa, pga - nB[j-1] * xcA[col-1], pgb - nA[col-1] * xcB[j-1]);

                        prof2 += 32;
                        for(c = f; c >= 0; c--){
                                pa += prof1[freq[c]] * prof2[freq[c]];
                        }
                        prof2 -= 32;
                        if(m->consistency){
                                pa += sparse_bonus_lookup(m->consistency, i, j - 1);
                        }

                        s[j].a = pa;

                        pga = s[j].ga;
                        s[j].ga = GM_MAX(xga - nB[j] * wA[col], xa - nB[j] * (wA[col] + xoA[col]));

                        pgb = s[j].gb;
                        s[j].gb = GM_MAX(pgb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xoB[j]));

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                prof2 += 64;
                ca = s[j].a;

                pa = GM_MAX3(pa, pga - nB[j-1] * xcA[col-1], pgb - nA[col-1] * xcB[j-1]);

                prof2 += 32;
                for(c = f; c >= 0; c--){
                        pa += prof1[freq[c]] * prof2[freq[c]];
                }
                prof2 -= 32;
                if(m->consistency){
                        pa += sparse_bonus_lookup(m->consistency, i, j - 1);
                }

                s[j].a = pa;
                s[j].ga = -FLT_MAX;
                s[j].gb = GM_MAX(s[j].gb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xoB[j]));

                prof2 -= (m->endb - m->startb) << 6;
        }
        return OK;
}

static int gm_pp_backward_scalar(struct aln_mem* m)
{
        const struct gm* g = m->gm;
        const float* nA = g->A.nres; const float* wA = g->A.w; const float* xoA = g->A.xo; const float* xcA = g->A.xc;
        const float* nB = g->B.nres; const float* wB = g->B.w; const float* xoB = g->B.xo; const float* xcB = g->B.xc;
        unsigned int freq[24];
        struct states* s = m->b;
        const float* prof1 = m->prof1;
        const float* prof2 = m->prof2;
        float pa, pga, pgb, ca, xa, xga;
        int i, j, c, f, col;
        const int kE = m->enda_2;

        prof1 += (m->enda_2 + 1) << 6;
        s[m->endb].a = s[0].a;
        s[m->endb].ga = s[0].ga;
        s[m->endb].gb = s[0].gb;
        /* the row after A column enda_2: gaps in A at boundary enda_2; cell j covers B column j+1 */
        for(j = m->endb - 1; j > m->startb; j--){
                s[j].a = -FLT_MAX;
                s[j].ga = GM_MAX(s[j+1].ga - nB[j+1] * wA[kE], s[j+1].a - nB[j+1] * (wA[kE] + xcA[kE]));
                s[j].gb = -FLT_MAX;
        }
        prof2 += (m->startb + 1) << 6;

        s[m->startb].a = -FLT_MAX;
        s[m->startb].ga = -FLT_MAX;
        s[m->startb].gb = -FLT_MAX;

        i = m->enda_2 - m->starta_2;
        while(i--){
                prof1 -= 64;
                col = m->starta_2 + i + 1;   /* A column of this row; gaps in A in this row sit at boundary col-1 */
                f = 0;
                for(j = 0; j < 23; j++){
                        if(prof1[j]){
                                freq[f] = j;
                                f++;
                        }
                }
                f--;

                pa = s[m->endb].a;
                pga = s[m->endb].ga;
                pgb = s[m->endb].gb;
                s[m->endb].a = -FLT_MAX;
                s[m->endb].ga = -FLT_MAX;

                xa = s[m->endb].a;
                xga = s[m->endb].ga;

                s[m->endb].gb = GM_MAX(pgb - nA[col] * wB[m->endb], pa - nA[col] * (wB[m->endb] + xcB[m->endb]));

                prof2 += (m->endb - m->startb) << 6;
                for(j = m->endb - 1; j > m->startb; j--){
                        prof2 -= 64;
                        ca = s[j].a;

                        /* match (col, j+1) before: match; a gap run in A (boundary col) starting at B column j+2; a
                           gap run in B (boundary j+1) starting at A column col+1 */
                        pa = GM_MAX3(pa, pga - nB[j+2] * xoA[col], pgb - nA[col+1] * xoB[j+1]);

                        prof2 += 32;
                        for(c = f; c >= 0; c--){
                                pa += prof1[freq[c]] * prof2[freq[c]];
                        }
                        prof2 -= 32;
                        if(m->consistency){
                                pa += sparse_bonus_lookup(m->consistency, m->starta_2 + i, j);
                        }

                        s[j].a = pa;

                        pga = s[j].ga;
                        s[j].ga = GM_MAX(xga - nB[j+1] * wA[col-1], xa - nB[j+1] * (wA[col-1] + xcA[col-1]));

                        pgb = s[j].gb;
                        s[j].gb = GM_MAX(pgb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xcB[j]));

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                prof2 -= 64;
                ca = s[j].a;

                pa = GM_MAX3(pa, pga - nB[j+2] * xoA[col], pgb - nA[col+1] * xoB[j+1]);
                prof2 += 32;
                for(c = f; c >= 0; c--){
                        pa += prof1[freq[c]] * prof2[freq[c]];
                }
                prof2 -= 32;
                if(m->consistency){
                        pa += sparse_bonus_lookup(m->consistency, m->starta_2 + i, j);
                }
                s[j].a = pa;
                s[j].ga = -FLT_MAX;
                s[j].gb = GM_MAX(s[j].gb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xcB[j]));
        }
        return OK;
}

/* Shared by both meetups: the joins at row boundary mid. The forward half ends at A column mid, the backward half
   starts at A column mid+1; f[i] / b[i] are B cell i. */
static void meetup_core(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score)
{
        const struct gm* g = m->gm;
        const float* nA = g->A.nres; const float* xoA = g->A.xo; const float* xcA = g->A.xc;
        const float* nB = g->B.nres; const float* xoB = g->B.xo; const float* xcB = g->B.xc;
        struct states* f = m->f;
        struct states* b = m->b;
        const int mid = old_cor[4];
        int i;
        int c = -1;
        int c2 = -1;
        int transition = -1;
        int transition2 = -1;
        float s_tmp;
        float max = -FLT_MAX;
        float max2 = -FLT_MAX;
        float middle = (float)(old_cor[3] - old_cor[2]) / 2.0F + (float)old_cor[2];
        float sub = 0.0F;

#define GM_TRY(val, tr) do{ s_tmp = (val); \
                if(s_tmp > max){ max2 = max; c2 = c; transition2 = transition; max = s_tmp; transition = (tr); c = i; } \
                else if(s_tmp > max2){ max2 = s_tmp; c2 = i; transition2 = (tr); } }while(0)

        for(i = old_cor[2]; i < old_cor[3]; i++){
                sub = fabsf(middle - (float)i) / 1000.0F;
                GM_TRY(f[i].a + b[i].a - sub, 1);                                  /* a -> a   */
                GM_TRY(f[i].a + b[i].ga - nB[i+1] * xoA[mid] - sub, 2);             /* a -> ga  */
                GM_TRY(f[i].a + b[i].gb - nA[mid+1] * xoB[i] - sub, 3);             /* a -> gb  */
                GM_TRY(f[i].ga + b[i].a - nB[i] * xcA[mid] - sub, 5);               /* ga -> a  */
                GM_TRY(f[i].gb + b[i].gb - sub, 6);                                 /* gb -> gb */
                GM_TRY(f[i].gb + b[i].a - nA[mid] * xcB[i] - sub, 7);               /* gb -> a  */
        }
        i = old_cor[3];
        sub = fabsf(middle - (float)i) / 1000.0F;
        GM_TRY(f[i].a + b[i].gb - nA[mid+1] * xoB[i] - sub, 3);
        GM_TRY(f[i].gb + b[i].gb - sub, 6);
#undef GM_TRY

        /* identical to the default meetups from here: margins and the optional flips */
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
}

int gm_pp_meetup(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score)
{
        meetup_core(m, old_cor, meet, t, score);
        return OK;
}

/* ------------------------------------------------------------------ DP: sequence-profile (prof1 = profile A, seq2 = sequence B) */

static int gm_sp_forward_scalar(struct aln_mem* m)
{
        const struct gm* g = m->gm;
        const float* nA = g->A.nres; const float* wA = g->A.w; const float* xoA = g->A.xo; const float* xcA = g->A.xc;
        const float* nB = g->B.nres; const float* wB = g->B.w; const float* xoB = g->B.xo; const float* xcB = g->B.xc;
        struct states* s = m->f;
        const float* prof1 = m->prof1;
        const uint8_t* seq2 = m->seq2;
        float pa, pga, pgb, ca, xa, xga;
        int i, j, col;
        const int k0 = m->starta;

        prof1 += (m->starta) << 6;
        s[m->startb].a = s[0].a;
        s[m->startb].ga = s[0].ga;
        s[m->startb].gb = s[0].gb;
        for(j = m->startb + 1; j < m->endb; j++){
                s[j].a = -FLT_MAX;
                s[j].ga = GM_MAX(s[j-1].ga - nB[j] * wA[k0], s[j-1].a - nB[j] * (wA[k0] + xoA[k0]));
                s[j].gb = -FLT_MAX;
        }
        s[m->endb].a = -FLT_MAX;
        s[m->endb].ga = -FLT_MAX;
        s[m->endb].gb = -FLT_MAX;
        seq2--;

        for(i = m->starta; i < m->enda; i++){
                prof1 += 64;
                col = i + 1;

                pa = s[m->startb].a;
                pga = s[m->startb].ga;
                pgb = s[m->startb].gb;
                s[m->startb].a = -FLT_MAX;
                s[m->startb].ga = -FLT_MAX;

                xa = s[m->startb].a;
                xga = s[m->startb].ga;

                s[m->startb].gb = GM_MAX(pgb - nA[col] * wB[m->startb],
                                         pa - nA[col] * (wB[m->startb] + xoB[m->startb]));
                for(j = m->startb + 1; j < m->endb; j++){
                        ca = s[j].a;

                        pa = GM_MAX3(pa, pga - nB[j-1] * xcA[col-1], pgb - nA[col-1] * xcB[j-1]);

                        pa += prof1[32 + seq2[j]];
                        if(m->consistency){
                                pa += sparse_bonus_lookup(m->consistency, i, j - 1);
                        }

                        s[j].a = pa;

                        pga = s[j].ga;
                        s[j].ga = GM_MAX(xga - nB[j] * wA[col], xa - nB[j] * (wA[col] + xoA[col]));

                        pgb = s[j].gb;
                        s[j].gb = GM_MAX(pgb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xoB[j]));

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                ca = s[j].a;

                pa = GM_MAX3(pa, pga - nB[j-1] * xcA[col-1], pgb - nA[col-1] * xcB[j-1]);

                pa += prof1[32 + seq2[j]];
                if(m->consistency){
                        pa += sparse_bonus_lookup(m->consistency, i, j - 1);
                }

                s[j].a = pa;
                s[j].ga = -FLT_MAX;
                s[j].gb = GM_MAX(s[j].gb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xoB[j]));
        }
        return OK;
}

static int gm_sp_backward_scalar(struct aln_mem* m)
{
        const struct gm* g = m->gm;
        const float* nA = g->A.nres; const float* wA = g->A.w; const float* xoA = g->A.xo; const float* xcA = g->A.xc;
        const float* nB = g->B.nres; const float* wB = g->B.w; const float* xoB = g->B.xo; const float* xcB = g->B.xc;
        struct states* s = m->b;
        const float* prof1 = m->prof1;
        const uint8_t* seq2 = m->seq2;
        float pa, pga, pgb, ca, xa, xga;
        int i, j, col;
        const int kE = m->enda_2;

        prof1 += (m->enda_2 + 1) << 6;

        s[m->endb].a = s[0].a;
        s[m->endb].ga = s[0].ga;
        s[m->endb].gb = s[0].gb;
        for(j = m->endb - 1; j > m->startb; j--){
                s[j].a = -FLT_MAX;
                s[j].ga = GM_MAX(s[j+1].ga - nB[j+1] * wA[kE], s[j+1].a - nB[j+1] * (wA[kE] + xcA[kE]));
                s[j].gb = -FLT_MAX;
        }

        s[m->startb].a = -FLT_MAX;
        s[m->startb].ga = -FLT_MAX;
        s[m->startb].gb = -FLT_MAX;

        i = m->enda_2 - m->starta_2;
        while(i--){
                prof1 -= 64;
                col = m->starta_2 + i + 1;

                pa = s[m->endb].a;
                pga = s[m->endb].ga;
                pgb = s[m->endb].gb;
                s[m->endb].a = -FLT_MAX;
                s[m->endb].ga = -FLT_MAX;

                xa = s[m->endb].a;
                xga = s[m->endb].ga;

                s[m->endb].gb = GM_MAX(pgb - nA[col] * wB[m->endb], pa - nA[col] * (wB[m->endb] + xcB[m->endb]));

                for(j = m->endb - 1; j > m->startb; j--){
                        ca = s[j].a;

                        pa = GM_MAX3(pa, pga - nB[j+2] * xoA[col], pgb - nA[col+1] * xoB[j+1]);
                        pa += prof1[32 + seq2[j]];
                        if(m->consistency){
                                pa += sparse_bonus_lookup(m->consistency, m->starta_2 + i, j);
                        }

                        s[j].a = pa;

                        pga = s[j].ga;
                        s[j].ga = GM_MAX(xga - nB[j+1] * wA[col-1], xa - nB[j+1] * (wA[col-1] + xcA[col-1]));

                        pgb = s[j].gb;
                        s[j].gb = GM_MAX(pgb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xcB[j]));

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                ca = s[j].a;

                pa = GM_MAX3(pa, pga - nB[j+2] * xoA[col], pgb - nA[col+1] * xoB[j+1]);
                pa += prof1[32 + seq2[j]];
                if(m->consistency){
                        pa += sparse_bonus_lookup(m->consistency, m->starta_2 + i, j);
                }

                s[j].a = pa;
                s[j].ga = -FLT_MAX;
                s[j].gb = GM_MAX(s[j].gb - nA[col] * wB[j], ca - nA[col] * (wB[j] + xcB[j]));
        }
        return OK;
}

int gm_sp_meetup(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score)
{
        meetup_core(m, old_cor, meet, t, score);
        return OK;
}

/* ------------------------------------------------------------------ profile counts */

void gm_leaf_counts(float* prof, int len)
{
        for(int c = 1; c <= len; c++){
                P(c, GM_NRES) = 1.0f;
                P(c, GM_LEAD) = 0.0f;
                P(c, GM_TRAIL) = 0.0f;
                P(c, GM_GOPEN) = 0.0f;
                P(c, GM_GEXT) = 0.0f;
        }
}
#undef P

/* members of one side at the boundary before its next column `nxt` (pointer to that column; nxt - 64 is the column to
   the left), k = columns consumed so far, len = its length */
static void side_boundary(const float* nxt, int k, int len, float nsip, float* lead, float* trail, float* g_left,
                          float* gop_right)
{
        *lead = (k == 0) ? nsip : (nxt - 64)[GM_LEAD];
        *trail = (k == len) ? nsip : nxt[GM_TRAIL];
        *g_left = (k == 0) ? 0.0f : (nxt - 64)[GM_GOPEN] + (nxt - 64)[GM_GEXT];
        *gop_right = (k == len) ? 0.0f : nxt[GM_GOPEN];
}

/* update_n with the member gap model: same profile merge (including the sequence-weight rebalancing), the member
   counts carried along, and the charge stored in a new gap column computed per member. */
KALIGN_ALWAYS_INLINE int update_n_gm_body(const float* profa, const float* profb, float* newp, struct aln_param* ap, int* path, int sipa, int sipb)
{
        float gp;
        int i, j, c;
        int ka = 0, kb = 0, len_a = 0, len_b = 0;
        int prev_a_ins = 0, prev_b_ins = 0;
        const float fa = (float)sipa;
        const float fb = (float)sipb;

        for(c = 1; path[c] != 3; c++){
                if(!path[c]){ len_a++; len_b++; }
                else if(path[c] & 1){ len_b++; }
                else if(path[c] & 2){ len_a++; }
        }

        float scaleA = 1.0f;
        float scaleB = 1.0f;
        int do_rebalance = 0;
        if(ap->use_seq_weights > 0.0f && sipa > 0 && sipb > 0){
                float pseudo = ap->use_seq_weights;
                float total = (float)(sipa + sipb);
                float denom = total + 2.0f * pseudo;
                scaleA = total * ((float)sipa + pseudo) / (denom * (float)sipa);
                scaleB = total * ((float)sipb + pseudo) / (denom * (float)sipb);
                do_rebalance = 1;
        }

        /* first boundary row */
        if(do_rebalance){
                for(i = 0; i < 23; i++){
                        newp[i] = profa[i] * scaleA + profb[i] * scaleB;
                }
                for(i = 23; i < 64; i++){
                        newp[i] = profa[i] + profb[i];
                }
        }else{
                for(i = 64; i--;){
                        newp[i] = profa[i] + profb[i];
                }
        }

        profa += 64;
        profb += 64;
        newp += 64;

        c = 1;
        while(path[c] != 3){
                if(!path[c]){
                        /* match column: counts add */
                        if(do_rebalance){
                                for(i = 0; i < 23; i++){
                                        newp[i] = profa[i] * scaleA + profb[i] * scaleB;
                                }
                                for(i = 23; i < 64; i++){
                                        newp[i] = profa[i] + profb[i];
                                }
                                {
                                        float** subm = ap->subm;
                                        float dA = scaleA - 1.0f;
                                        float dB = scaleB - 1.0f;
                                        for(j = 0; j < 23; j++){
                                                float delta = 0.0f;
                                                int aa;
                                                for(aa = 0; aa < 23; aa++){
                                                        delta += (profa[aa]*dA + profb[aa]*dB) * subm[aa][j];
                                                }
                                                newp[32 + j] += delta;
                                        }
                                }
                        }else{
                                for(i = 64; i--;){
                                        newp[i] = profa[i] + profb[i];
                                }
                        }
                        /* a side whose previous merged column was a gap run inserted into it: its members' gaps that
                           opened in this column now continue that run */
                        if(prev_a_ins){
                                newp[GM_GOPEN] -= profa[GM_GOPEN];
                                newp[GM_GEXT] += profa[GM_GOPEN];
                        }
                        if(prev_b_ins){
                                newp[GM_GOPEN] -= profb[GM_GOPEN];
                                newp[GM_GEXT] += profb[GM_GOPEN];
                        }
                        profa += 64;
                        profb += 64;
                        ka++; kb++;
                        prev_a_ins = 0;
                        prev_b_ins = 0;
                }else if(path[c] & 1){
                        /* gap column inserted into A at boundary ka; B column kb+1 consumed */
                        float lead, trail, gl, gr, nT, nG, nR, nI, w, xo, xc;
                        int first, last;
                        for(i = 64; i--;){
                                newp[i] = profb[i];
                        }
                        if(prev_b_ins){
                                newp[GM_GOPEN] -= profb[GM_GOPEN];
                                newp[GM_GEXT] += profb[GM_GOPEN];
                        }
                        side_boundary(profa, ka, len_a, fa, &lead, &trail, &gl, &gr);
                        nT = lead + trail;
                        nG = gl + gr;
                        nR = fa - nT - nG;
                        nI = fa - nT;
                        newp[GM_LEAD] += lead;
                        newp[GM_TRAIL] += trail;
                        if(prev_a_ins){
                                newp[GM_GEXT] += nI;
                        }else{
                                newp[GM_GEXT] += gl;
                                newp[GM_GOPEN] += nI - gl;
                        }
                        first = !prev_a_ins;
                        last = !(path[c+1] != 3 && (path[c+1] & 1));
                        w = (nR + nG) * ap->gpe + nT * ap->tgpe;
                        xo = nR * (ap->gpo - ap->gpe);
                        xc = nR * ap->gpo;
                        gp = w + (first ? xo : 0.0f) + (last ? xc : 0.0f);
                        for(j = 32; j < 55; j++){
                                newp[j] -= gp;
                        }
                        profb += 64;
                        kb++;
                        prev_a_ins = 1;
                        prev_b_ins = 0;
                }else if(path[c] & 2){
                        /* gap column inserted into B at boundary kb; A column ka+1 consumed */
                        float lead, trail, gl, gr, nT, nG, nR, nI, w, xo, xc;
                        int first, last;
                        for(i = 64; i--;){
                                newp[i] = profa[i];
                        }
                        if(prev_a_ins){
                                newp[GM_GOPEN] -= profa[GM_GOPEN];
                                newp[GM_GEXT] += profa[GM_GOPEN];
                        }
                        side_boundary(profb, kb, len_b, fb, &lead, &trail, &gl, &gr);
                        nT = lead + trail;
                        nG = gl + gr;
                        nR = fb - nT - nG;
                        nI = fb - nT;
                        newp[GM_LEAD] += lead;
                        newp[GM_TRAIL] += trail;
                        if(prev_b_ins){
                                newp[GM_GEXT] += nI;
                        }else{
                                newp[GM_GEXT] += gl;
                                newp[GM_GOPEN] += nI - gl;
                        }
                        first = !prev_b_ins;
                        last = !(path[c+1] != 3 && (path[c+1] & 2));
                        w = (nR + nG) * ap->gpe + nT * ap->tgpe;
                        xo = nR * (ap->gpo - ap->gpe);
                        xc = nR * ap->gpo;
                        gp = w + (first ? xo : 0.0f) + (last ? xc : 0.0f);
                        for(j = 32; j < 55; j++){
                                newp[j] -= gp;
                        }
                        profa += 64;
                        ka++;
                        prev_b_ins = 1;
                        prev_a_ins = 0;
                }
                newp += 64;
                c++;
        }
        /* last boundary row */
        if(do_rebalance){
                for(i = 0; i < 23; i++){
                        newp[i] = profa[i] * scaleA + profb[i] * scaleB;
                }
                for(i = 23; i < 64; i++){
                        newp[i] = profa[i] + profb[i];
                }
        }else{
                for(i = 64; i--;){
                        newp[i] = profa[i] + profb[i];
                }
        }
        return OK;
}

/* update_n_gm is compiled twice: the compiler vectorises its column loops (element-wise sums; the one reduction keeps
   its order, and this file is compiled with -ffp-contract=off), so the AVX2 copy gives the same profile.  It runs
   on a CPU with AVX2, the plain copy otherwise (kalign_simd.h). */
static int update_n_gm_plain(const float* profa, const float* profb, float* newp, struct aln_param* ap, int* path, int sipa, int sipb)
{
        return update_n_gm_body(profa, profb, newp, ap, path, sipa, sipb);
}

#if KALIGN_AVX2
static KALIGN_TARGET_AVX2 int update_n_gm_avx2(const float* profa, const float* profb, float* newp, struct aln_param* ap, int* path, int sipa, int sipb)
{
        return update_n_gm_body(profa, profb, newp, ap, path, sipa, sipb);
}
#endif

int update_n_gm(const float* profa, const float* profb, float* newp, struct aln_param* ap, int* path, int sipa, int sipb)
{
#if KALIGN_AVX2
        if(kalign_simd_on()){
                return update_n_gm_avx2(profa, profb, newp, ap, path, sipa, sipb);
        }
#endif
        return update_n_gm_plain(profa, profb, newp, ap, path, sipa, sipb);
}

/* ------------------------------------------------------------------ wavefront SIMD kernels

   The member-model forward/backward passes on the wavefront skeleton of aln_wave.c (pw_sweep / spw_sweep): strips
   of up to 8 rows advance together along anti-diagonals, one row per lane; the left and right columns and leftover
   rows are computed by scalar code on the same row arrays. Every gap cost of the model is a per-row factor times a
   per-column factor, so a lane holds its row's factors (nA, and A's W / entry extra / exit extra at the row's A
   boundary) and the columns' factors (nB, and B's W / entry / exit at the cell's B boundary) are loaded per step
   from tables built once per alignment.

   Bit-identical to the scalar kernels above: both are compiled with -ffp-contract=off and every operation of the
   scalar statement is done in the same order (products rounded, then the subtraction; MAX(a, b) = a > b ? a : b,
   which is vf8_max's semantics).

   Direction-generic notation (as aln_wave.c): cell t = 0..W-1 is lcol + dir*(t+1), y(t) = y0 + dir*t is the
   0-based B residue index of the column the cell covers (B column y+1); row r is A residue r (A column r+1).
       A side, per row r:  nA[r+1];  kA = r+1 (forward) or r (backward)
       B side, per column: nB[y+1]; kB = y+1 (forward) or y (backward)
       entry extra = xo forward / xc backward;  exit extra = xc forward / xo backward */
#ifdef KALIGN_WAVE
#include "wave_simd.h"
#else
#define WAVE_SIMD 0
#endif

#if WAVE_SIMD

#define GMW_QPAD 32
#define GMW_WL 8
#define GMW_PART_MIN 3
#define GMW_MIN_W 16

#if KALIGN_NEON
VF8_INLINE vf8 gmw_mul(vf8 a, vf8 b) { vf8 r = { vmulq_f32(a.lo, b.lo), vmulq_f32(a.hi, b.hi) }; return r; }
#else
VF8_INLINE vf8 gmw_mul(vf8 a, vf8 b) { return _mm256_mul_ps(a, b); }
#endif

struct gmw_dir {
        float* T;       /* 23 x n match scores of the B columns (profile-profile only) */
        float* NB;      /* nB[y+1] */
        float* WB;      /* wB[kB] */
        float* WE;      /* wB[kB] + entry extra */
        float* XX;      /* exit extra at kB */
};

struct gmw {
        int n;          /* len_b + 2*QPAD */
        int pp;
        struct gmw_dir d[2];    /* 0 forward (reversed index u = len_b-1-y), 1 backward (index y) */
};

static void* gm_wave_build(const struct gm* g, const struct aln_mem* m)
{
        struct gmw* w = NULL;
        if(!kalign_simd_on()) return NULL;
        const int lb = m->len_b;
        if(lb < GMW_MIN_W + 1) return NULL;
        w = calloc(1, sizeof(*w));
        if(!w) return NULL;
        w->n = lb + 2 * GMW_QPAD;
        w->pp = m->prof2 != NULL;
        for(int d = 0; d < 2; d++){
                struct gmw_dir* D = &w->d[d];
                if(w->pp) D->T = calloc((size_t)23 * w->n, sizeof(float));
                D->NB = calloc((size_t)w->n, sizeof(float));
                D->WB = calloc((size_t)w->n, sizeof(float));
                D->WE = calloc((size_t)w->n, sizeof(float));
                D->XX = calloc((size_t)w->n, sizeof(float));
                if((w->pp && !D->T) || !D->NB || !D->WB || !D->WE || !D->XX){
                        gm_wave_release(w);
                        return NULL;
                }
                for(int y = 0; y < lb; y++){
                        const int ix = GMW_QPAD + (d == 0 ? lb - 1 - y : y);
                        const int kB = d == 0 ? y + 1 : y;
                        D->NB[ix] = g->B.nres[y + 1];
                        D->WB[ix] = g->B.w[kB];
                        D->WE[ix] = g->B.w[kB] + (d == 0 ? g->B.xo[kB] : g->B.xc[kB]);
                        D->XX[ix] = d == 0 ? g->B.xc[kB] : g->B.xo[kB];
                        if(w->pp){
                                const float* p = m->prof2 + (size_t)(y + 1) * 64;
                                for(int c = 0; c < 23; c++){
                                        D->T[(size_t)c * w->n + ix] = p[32 + c];
                                }
                        }
                }
        }
        return w;
}

static void gm_wave_release(void* wv)
{
        struct gmw* w = wv;
        if(!w) return;
        for(int d = 0; d < 2; d++){
                free(w->d[d].T); free(w->d[d].NB); free(w->d[d].WB); free(w->d[d].WE); free(w->d[d].XX);
        }
        free(w);
}

struct gmw_buf {
        float* ra;
        float* rga;
        float* rgb;
        float* rd;
        float* bl;
        size_t cap;
};

static _Thread_local struct gmw_buf gb = {0};

static void gmw_reserve(int W)
{
        size_t n = (size_t)W + 2 * GMW_QPAD;
        if(n > gb.cap){
                free(gb.ra); free(gb.rga); free(gb.rgb); free(gb.rd); free(gb.bl);
                gb.ra = malloc(n * sizeof(float));
                gb.rga = malloc(n * sizeof(float));
                gb.rgb = malloc(n * sizeof(float));
                gb.rd = malloc(n * sizeof(float));
                gb.bl = calloc(n * GMW_WL, sizeof(float));
                if(!gb.ra || !gb.rga || !gb.rgb || !gb.rd || !gb.bl) abort();
                gb.cap = n;
                tl_scratch_track(TL_SCRATCH_GAPMODEL, (void*[]){gb.ra, gb.rga, gb.rgb, gb.rd, gb.bl}, 5);
        }
}

/* match score of row r (A column r+1, non-zero residues nz[0..f]) against B residue y, added to pa in the scalar
   kernel's order */
static inline float gmw_match(const struct aln_mem* m, int pp, float pa, const float* p1, const int* nz, int f, int y)
{
        if(pp){
                const float* p2 = m->prof2 + (size_t)(y + 1) * 64;
                for(int c = f; c >= 0; c--){
                        pa += p1[nz[c]] * p2[32 + nz[c]];
                }
        }else{
                pa += p1[32 + m->seq2[y]];
        }
        return pa;
}

static inline __attribute__((always_inline)) WAVE_TARGET int gmw_sweep(struct aln_mem* m, const int dir, const int pp)
{
        const struct gm* g = m->gm;
        const struct gmw* wv = g->wv;
        const struct gmw_dir* D = &wv->d[dir > 0 ? 0 : 1];
        struct states* s = dir > 0 ? m->f : m->b;
        const float* P1 = m->prof1;
        const int startb = m->startb;
        const int endb = m->endb;
        const int lb = m->len_b;
        const int W = endb - startb - 1;
        const struct sparse_bonus* cons = m->consistency;
        const int lcol = dir > 0 ? startb : endb;
        const int rcol = dir > 0 ? endb : startb;
        const int y0 = dir > 0 ? startb : endb - 1;
        const int yr = dir > 0 ? endb - 1 : startb;
        const int yl = y0 - dir;                 /* the left column's y */
        const int bstride = W + 2 * GMW_QPAD;
        const int c0 = dir > 0 ? (lb - 1 - y0) : y0;
        const int n = wv->n;
        const float* T = pp ? D->T + GMW_QPAD : NULL;
        const float* NB = D->NB + GMW_QPAD;
        const float* WB = D->WB + GMW_QPAD;
        const float* WE = D->WE + GMW_QPAD;
        const float* XX = D->XX + GMW_QPAD;
        const float* nA = g->A.nres; const float* wA = g->A.w;
        const float* entA = dir > 0 ? g->A.xo : g->A.xc;
        const float* exA = dir > 0 ? g->A.xc : g->A.xo;
        const float* nBa = g->B.nres; const float* wBa = g->B.w;
        const float* entB = dir > 0 ? g->B.xo : g->B.xc;
        const float* exB = dir > 0 ? g->B.xc : g->B.xo;
        int r, rend, h, t, q, sstep;

#define KA(rr) (dir > 0 ? (rr) + 1 : (rr))
#define KB(yy) (dir > 0 ? (yy) + 1 : (yy))
#define YT(tt) (y0 + dir * (tt))
#define CIX(tt) (dir > 0 ? lb - 1 - YT(tt) : YT(tt))
#define ROWP(rr) (P1 + (size_t)((rr) + 1) * 64)

        gmw_reserve(W);
        float* restrict RA = gb.ra + GMW_QPAD;
        float* restrict RGA = gb.rga + GMW_QPAD;
        float* restrict RGB = gb.rgb + GMW_QPAD;
        float* restrict RD = gb.rd + GMW_QPAD;

        /* the edge columns' B factors */
        const int kBL = KB(yl), kBR = KB(yr);
        const float nBL = nBa[yl + 1];
        const float wBL = wBa[kBL], weBL = wBa[kBL] + entB[kBL], xxBL = exB[kBL];
        const float wBR = wBa[kBR], weBR = wBa[kBR] + entB[kBR];

        /* ---- the row before the first swept row (scalar init row of the kernels) ---- */
        const int r_init = dir > 0 ? m->starta - 1 : m->enda_2;
        {
                const int ki = KA(r_init);
                const float wi = wA[ki], wei = wA[ki] + entA[ki];
                float pga, pa;
                s[lcol].a = s[0].a;
                s[lcol].ga = s[0].ga;
                s[lcol].gb = s[0].gb;
                RA[-1] = s[lcol].a;
                RGA[-1] = s[lcol].ga;
                RGB[-1] = s[lcol].gb;
                pga = RGA[-1];
                pa = RA[-1];
                for(t = 0; t < W; t++){
                        const float nb = NB[CIX(t)];
                        const float x1 = pga - nb * wi;
                        const float x2 = pa - nb * wei;
                        const float gg = x1 > x2 ? x1 : x2;
                        RA[t] = -FLT_MAX;
                        RGA[t] = gg;
                        RGB[t] = -FLT_MAX;
                        pga = gg;
                        pa = -FLT_MAX;
                }
                /* exit values (the next row's diagonal) */
                {
                        const float na = nA[r_init + 1], xa = exA[ki];
                        float u = RGA[-1] - nBL * xa;
                        float v = RA[-1] > u ? RA[-1] : u;
                        u = RGB[-1] - na * xxBL;
                        RD[-1] = v > u ? v : u;
                        for(t = 0; t < W; t++){
                                const int ix = CIX(t);
                                u = RGA[t] - NB[ix] * xa;
                                v = RA[t] > u ? RA[t] : u;
                                u = RGB[t] - na * XX[ix];
                                RD[t] = v > u ? v : u;
                        }
                }
        }
        for(t = W; t < W + GMW_WL + 1; t++){
                RA[t] = -FLT_MAX; RGA[t] = -FLT_MAX; RGB[t] = -FLT_MAX; RD[t] = -FLT_MAX;
        }
        float ra = -FLT_MAX, rga = -FLT_MAX, rgb = -FLT_MAX;
        const vf8 vneg = vf8_set1(-FLT_MAX);

        if(dir > 0){ r = m->starta; rend = m->enda; h = m->enda - m->starta; }
        else { r = m->enda_2 - 1; rend = m->starta_2 - 1; h = m->enda_2 - m->starta_2; }

        for(int hs = 0; h >= GMW_PART_MIN; h -= hs, r += dir * hs){
                hs = h < GMW_WL ? h : GMW_WL;
                float lgb[GMW_WL], dleft[GMW_WL], rbon[GMW_WL], dw[GMW_WL];
                float dcap[GMW_WL][GMW_WL];
                float na_l[GMW_WL], wa_l[GMW_WL], wea_l[GMW_WL], xa_l[GMW_WL];
                const float* p1q[GMW_WL];
                const float* bp[GMW_WL];
                uint32_t umask = 0;
                int ulist[23], nu = 0;

                /* lane constants and the left column chain */
                {
                        float pa0 = RA[-1], pgb0 = RGB[-1];
                        for(q = 0; q < GMW_WL; q++){
                                const int rq = r + dir * (q < hs ? q : hs - 1);
                                const int ka = KA(rq);
                                float u, v;
                                p1q[q] = ROWP(rq);
                                na_l[q] = nA[rq + 1];
                                wa_l[q] = wA[ka];
                                wea_l[q] = wA[ka] + entA[ka];
                                xa_l[q] = exA[ka];
                                u = pgb0 - na_l[q] * wBL;
                                v = pa0 - na_l[q] * weBL;
                                lgb[q] = u > v ? u : v;
                                u = -FLT_MAX - nBL * xa_l[q];
                                v = -FLT_MAX > u ? -FLT_MAX : u;
                                u = lgb[q] - na_l[q] * xxBL;
                                dleft[q] = v > u ? v : u;
                                pa0 = -FLT_MAX;
                                pgb0 = lgb[q];
                                if(pp){
                                        for(int c = 0; c < 23; c++) if(p1q[q][c]) umask |= 1u << c;
                                }
                        }
                }
                for(int c = 22; c >= 0; c--) if(umask & (1u << c)) ulist[nu++] = c;
                const vf8 vna = vf8_loadu(na_l);
                const vf8 vwa = vf8_loadu(wa_l);
                const vf8 vwea = vf8_loadu(wea_l);
                const vf8 vxa = vf8_loadu(xa_l);
                vf8 vw[23];
                for(int k = 0; k < nu; k++){
                        float tmp[GMW_WL];
                        for(q = 0; q < GMW_WL; q++) tmp[q] = p1q[q][ulist[k]];
                        vw[k] = vf8_loadu(tmp);
                }
                for(q = 0; q < GMW_WL; q++){
                        const int rq = r + dir * (q < hs ? q : hs - 1);
                        rbon[q] = 0.0F;
                        bp[q] = gb.bl + (size_t)q * bstride + GMW_QPAD - q;
                        if(cons && q < hs){
                                const int bb = rq * cons->K;
                                for(int e = 0; e < cons->K; e++){
                                        int c = cons->cols[bb + e];
                                        if(c < 0) break;
                                        if(c == yr){
                                                rbon[q] += cons->vals[bb + e];
                                        }else{
                                                int tt = (c - y0) * dir;
                                                if(tt >= 0 && tt < W) gb.bl[(size_t)q * bstride + GMW_QPAD + tt] += cons->vals[bb + e];
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
                        if(pp){
                                for(int k = 0; k < nu; k++){
                                        a = vf8_add(a, gmw_mul(vw[k], vf8_loadu(T + (size_t)ulist[k] * n + ci)));
                                }
                        }else{
                                float subl[GMW_WL];
                                for(int z = 0; z < GMW_WL; z++){
                                        const int tt = sstep - z;
                                        subl[z] = (tt >= 0 && tt < W) ? p1q[z][32 + m->seq2[y0 + dir * tt]] : 0.0F;
                                }
                                a = vf8_add(a, vf8_loadu(subl));
                        }
                        if(cons){
                                a = vf8_add(a, vf8_gather(bp, sstep));
                        }
                        const vf8 upA = vf8_shift_in(pA, RA[sstep]);
                        const vf8 upGB = vf8_shift_in(pGB, RGB[sstep]);
                        const vf8 nb = vf8_loadu(NB + ci);
                        const vf8 wbv = vf8_loadu(WB + ci);
                        const vf8 webv = vf8_loadu(WE + ci);
                        const vf8 xxv = vf8_loadu(XX + ci);
                        vf8 ga = vf8_max(vf8_sub(pGA, gmw_mul(nb, vwa)), vf8_sub(pA, gmw_mul(nb, vwea)));
                        vf8 gbv = vf8_max(vf8_sub(upGB, gmw_mul(vna, wbv)), vf8_sub(upA, gmw_mul(vna, webv)));
                        vf8 d = vf8_max(vf8_max(a, vf8_sub(ga, gmw_mul(nb, vxa))), vf8_sub(gbv, gmw_mul(vna, xxv)));
                        if(sstep < GMW_WL - 1){
                                const int ql = sstep + 1;
                                a = vf8_with_lane(a, ql, -FLT_MAX);
                                ga = vf8_with_lane(ga, ql, -FLT_MAX);
                                gbv = vf8_with_lane(gbv, ql, lgb[ql]);
                                d = vf8_with_lane(d, ql, dleft[ql]);
                        }
                        if(sstep >= hs - 1){
                                const int tw = sstep - (hs - 1);
                                if(hs == GMW_WL){
                                        RA[tw] = vf8_last(a);
                                        RGA[tw] = vf8_last(ga);
                                        RGB[tw] = vf8_last(gbv);
                                        RD[tw] = vf8_last(d);
                                }else{
                                        RA[tw] = vf8_lane(a, hs - 1);
                                        RGA[tw] = vf8_lane(ga, hs - 1);
                                        RGB[tw] = vf8_lane(gbv, hs - 1);
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
                        pGB = gbv;
                }
                for(q = 0; q < hs; q++) dw[q] = dcap[q][q];

                /* right column, rows in order (scalar) */
                for(q = 0; q < hs; q++){
                        const float* p1 = p1q[q];
                        int nz[23], f = 0;
                        float u, v;
                        if(pp){
                                for(int c = 0; c < 23; c++) if(p1[c]) nz[f++] = c;
                        }
                        const float ca = ra;
                        float p = gmw_match(m, pp, q ? dw[q - 1] : rdw_prev, p1, nz, f - 1, yr);
                        if(cons) p += rbon[q];
                        ra = p;
                        rga = -FLT_MAX;
                        u = rgb - na_l[q] * wBR;
                        v = ca - na_l[q] * weBR;
                        rgb = u > v ? u : v;
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
                                        if(c != yr && tt >= 0 && tt < W) gb.bl[(size_t)q * bstride + GMW_QPAD + tt] = 0.0F;
                                }
                        }
                }
        }

        /* ---- leftover rows: scalar on the row arrays ---- */
        for(; r != rend; r += dir){
                const float* p1 = ROWP(r);
                const int ka = KA(r);
                const float na = nA[r + 1], wa = wA[ka], wea = wA[ka] + entA[ka], xa = exA[ka];
                int nz[23], f = 0;
                float u, v;
                if(pp){
                        for(int c = 0; c < 23; c++) if(p1[c]) nz[f++] = c;
                }
                f--;
                float pa0 = RA[-1], pgb0 = RGB[-1];
                float pa = RD[-1];
                float xa_ = -FLT_MAX, xga = -FLT_MAX;
                RA[-1] = -FLT_MAX;
                RGA[-1] = -FLT_MAX;
                u = pgb0 - na * wBL;
                v = pa0 - na * weBL;
                RGB[-1] = u > v ? u : v;
                u = RGA[-1] - nBL * xa;
                v = RA[-1] > u ? RA[-1] : u;
                u = RGB[-1] - na * xxBL;
                RD[-1] = v > u ? v : u;
                for(t = 0; t < W; t++){
                        const int ix = CIX(t);
                        const float ca = RA[t], cgb = RGB[t];
                        const float dnext = RD[t];
                        const float nb = NB[ix];
                        float a = gmw_match(m, pp, pa, p1, nz, f, YT(t));
                        if(cons) a += sparse_bonus_lookup(cons, r, YT(t));
                        RA[t] = a;
                        u = xga - nb * wa;
                        v = xa_ - nb * wea;
                        RGA[t] = u > v ? u : v;
                        u = cgb - na * WB[ix];
                        v = ca - na * WE[ix];
                        RGB[t] = u > v ? u : v;
                        u = RGA[t] - nb * xa;
                        v = RA[t] > u ? RA[t] : u;
                        u = RGB[t] - na * XX[ix];
                        RD[t] = v > u ? v : u;
                        xa_ = RA[t];
                        xga = RGA[t];
                        pa = dnext;
                }
                {
                        const float ca = ra;
                        float p = gmw_match(m, pp, pa, p1, nz, f, yr);
                        if(cons) p += sparse_bonus_lookup(cons, r, yr);
                        ra = p;
                        rga = -FLT_MAX;
                        u = rgb - na * wBR;
                        v = ca - na * weBR;
                        rgb = u > v ? u : v;
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
#undef KA
#undef KB
#undef YT
#undef CIX
#undef ROWP
}

static int gmw_usable(const struct aln_mem* m, int dir)
{
        const struct gm* g = m->gm;
        if(!g->wv || m->endb - m->startb - 1 < GMW_MIN_W) return 0;
        return dir > 0 ? (m->enda > m->starta) : (m->enda_2 > m->starta_2);
}

#else   /* no SIMD */

static void* gm_wave_build(const struct gm* g, const struct aln_mem* m) { (void)g; (void)m; return NULL; }
static void gm_wave_release(void* wv) { (void)wv; }

#endif

#if WAVE_SIMD
/* entry points with a constant direction (an AVX2 function cannot be inlined into the plain dispatcher) */
static WAVE_TARGET int gmw_sweep_ppf(struct aln_mem* m) { return gmw_sweep(m, 1, 1); }
static WAVE_TARGET int gmw_sweep_ppb(struct aln_mem* m) { return gmw_sweep(m, -1, 1); }
static WAVE_TARGET int gmw_sweep_spf(struct aln_mem* m) { return gmw_sweep(m, 1, 0); }
static WAVE_TARGET int gmw_sweep_spb(struct aln_mem* m) { return gmw_sweep(m, -1, 0); }
#else
#define gmw_sweep_ppf NULL
#define gmw_sweep_ppb NULL
#define gmw_sweep_spf NULL
#define gmw_sweep_spb NULL
#endif

/* the wavefront kernel when the band is wide enough and the vector kernels are on, else the scalar kernel */
static int gm_dispatch(struct aln_mem* m, int dir, int (*wave)(struct aln_mem*), int (*scalar)(struct aln_mem*))
{
#if WAVE_SIMD
        if(kalign_simd_on() && gmw_usable(m, dir)){
                return wave(m);
        }
#else
        (void)dir; (void)wave;
#endif
        return scalar(m);
}

int gm_pp_forward(struct aln_mem* m) { return gm_dispatch(m, 1, gmw_sweep_ppf, gm_pp_forward_scalar); }
int gm_pp_backward(struct aln_mem* m) { return gm_dispatch(m, -1, gmw_sweep_ppb, gm_pp_backward_scalar); }
int gm_sp_forward(struct aln_mem* m) { return gm_dispatch(m, 1, gmw_sweep_spf, gm_sp_forward_scalar); }
int gm_sp_backward(struct aln_mem* m) { return gm_dispatch(m, -1, gmw_sweep_spb, gm_sp_backward_scalar); }
