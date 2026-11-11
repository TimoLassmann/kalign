/* Independent path scorer (double precision).
 *
 * Scores a raw path (path[i] = matched column j, 1-based, or -1 for a gapped
 * row) under the affine model that Kalign's top-level forward/backward
 * passes implement:
 *   match (i,j)          subst(i,j); entering from a gap adds the "close"
 *                        cost = open cost of the previous row/column
 *   internal gap         open [27] for the first position, extend [28] after
 *   terminal gap         [29] per position, no open (row 0 / row len_a for
 *                        gaps in A, column 0 / column len_b for gaps in B)
 * Parameters per kind (A = rows, B = columns):
 *   SS  subm[seq1][seq2] - subm_offset; gaps -gpo/-gpe/-tgpe on both sides
 *   SP  prof1[32 + seq2]; A gaps prof1[27/28/29]; B gaps -gpo/-gpe/-tgpe * sip
 *   PP  sum_c prof1[c] * prof2[32 + c]; gaps from prof1 / prof2 [27/28/29]
 *   plus, when the problem carries a consistency bonus and it is in use,
 *   sparse_bonus_lookup(bonus, i - 1, j) on every match (the kernels index
 *   bonus rows 0-based and columns 1-based).
 * The scorer is used to compare a prototype's path with the reference path
 * when they differ; it is not bit-identical to Kalign's float arithmetic.
 * APPROXIMATE: on the 216 reference dumps it matches the top-level meetup
 * score within 0.05 on 72 problems only; the remainder differ at gap-type
 * switches / Hirschberg junctions, and with a consistency bonus (the
 * backward kernels look the bonus up one column left of the substituted
 * column, see README).  Use it only as a relative proto-vs-ref measure.
 */
#include <stddef.h>

#include "kdp_kernel.h"

enum { ST_START, ST_M, ST_GA, ST_GB };

struct sc {
        const struct kdp_problem* p;
        int la, lb;
};

static double sub0(const struct sc* s, int i, int j);

static double sub(const struct sc* s, int i, int j)
{
        const struct kdp_problem* p = s->p;
        double v = sub0(s, i, j);
        if(p->h.consistency && p->use_bonus){
                v += (double)sparse_bonus_lookup(&p->bonus, i - 1, j);
        }
        return v;
}

static double sub0(const struct sc* s, int i, int j)
{
        const struct kdp_problem* p = s->p;
        int c;
        double v = 0.0;
        switch(p->h.kind){
        case KDP_SS:
                return (double)p->h.subm[p->seq1[i - 1] * 23 + p->seq2[j - 1]] - (double)p->h.subm_offset;
        case KDP_SP:
                return (double)p->prof1[(size_t)i * 64 + 32 + p->seq2[j - 1]];
        default:
                for(c = 0; c < 23; c++){
                        float a = p->prof1[(size_t)i * 64 + c];
                        if(a != 0.0F){
                                v += (double)a * (double)p->prof2[(size_t)j * 64 + 32 + c];
                        }
                }
                return v;
        }
}

/* which: 0 open, 1 extend, 2 terminal */
static double gapA(const struct sc* s, int i, int which)
{
        const struct kdp_problem* p = s->p;
        if(p->h.kind == KDP_SS){
                return which == 0 ? -(double)p->h.gpo : which == 1 ? -(double)p->h.gpe : -(double)p->h.tgpe;
        }
        return (double)p->prof1[(size_t)i * 64 + 27 + which];
}

static double gapB(const struct sc* s, int j, int which)
{
        const struct kdp_problem* p = s->p;
        if(p->h.kind == KDP_SS){
                return which == 0 ? -(double)p->h.gpo : which == 1 ? -(double)p->h.gpe : -(double)p->h.tgpe;
        }
        if(p->h.kind == KDP_SP){
                /* computed in float exactly as aln_seqprofile.c does */
                float v = which == 0 ? p->h.gpo * (float)p->h.sip
                        : which == 1 ? p->h.gpe * (float)p->h.sip
                        : p->h.tgpe * (float)p->h.sip;
                return -(double)v;
        }
        return (double)p->prof2[(size_t)j * 64 + 27 + which];
}

static double step_ga(const struct sc* s, int* st, int i, int j)
{
        double c;
        if(i == 0 || i == s->la){
                c = gapB(s, j, 2);
        }else{
                c = (*st == ST_GA) ? gapB(s, j, 1) : gapB(s, j, 0);
        }
        *st = ST_GA;
        return c;
}

static double step_gb(const struct sc* s, int* st, int i, int j)
{
        double c;
        if(j == 0 || j == s->lb){
                c = gapA(s, i, 2);
        }else{
                c = (*st == ST_GB) ? gapA(s, i, 1) : gapA(s, i, 0);
        }
        *st = ST_GB;
        return c;
}

double kdp_path_score(const struct kdp_problem* p, const int* path)
{
        struct sc s;
        int st = ST_START;
        int i, j, jj, pj = 0;
        double score = 0.0;
        s.p = p;
        s.la = p->h.len_a;
        s.lb = p->h.len_b;
        for(i = 1; i <= s.la; i++){
                j = path[i];
                if(j < 0){
                        score += step_gb(&s, &st, i, pj);
                        continue;
                }
                for(jj = pj + 1; jj < j; jj++){
                        score += step_ga(&s, &st, i - 1, jj);
                }
                if(st == ST_GA){
                        score += gapB(&s, j - 1, 0);
                }else if(st == ST_GB){
                        score += gapA(&s, i - 1, 0);
                }
                score += sub(&s, i, j);
                st = ST_M;
                pj = j;
        }
        for(jj = pj + 1; jj <= s.lb; jj++){
                score += step_ga(&s, &st, s.la, jj);
        }
        return score;
}
