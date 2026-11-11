#include "tldevel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "msa_struct.h"
#include "task.h"
#include "aln_param.h"
#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_refine.h"

#define ALN_VSM_IMPORT
#include "aln_vsm.h"

/* Per-merge variable scoring matrix (KALIGN_VSM_MERGE)
 * =====================================================
 *
 * Idea (MAFFT's VSM): subtract an offset from every substitution score when
 * the two groups being merged are close, so that segments that only align
 * because a global aligner must align something are left unaligned; add
 * nothing when the groups are distant.
 *
 * 1. Divergence of the merge.  The legacy offset used the mean normalised
 *    edit distance of each member to the tree anchors, averaged over both
 *    groups.  That number describes each sequence's place in the whole input,
 *    not the pair being merged, and edit distances of unrelated proteins
 *    saturate near 0.6, so the offset was nearly constant.  Here the
 *    divergence is measured on the merge itself: one DP of group a against
 *    group b without offset, then the identity of all aligned cross-group
 *    residue pairs, computed exactly from the profile residue counts
 *    (sum_k ca[k] cb[k] / (sum ca)(sum cb) over aligned columns).  The
 *    p-distance p = 1 - identity is corrected for multiple substitutions
 *    (Kimura 1983 for proteins, Jukes-Cantor for nucleotides), so d is in
 *    substitutions per site and keeps growing where identity changes
 *    slowly (protein d = 0.80 at 50% identity, 1.11 at 40%, 1.60 at 30%,
 *    1.98 at 25%, diverging near 15%).  Identity along an optimal
 *    alignment of unrelated sequences sits near 10-15%, where the
 *    correction has diverged, so unrelated groups get no offset.
 *
 * 2. Bound.  In the DP every aligned residue pair scores s(x,y) - off and
 *    every residue opposite a gap costs gpe (profile DPs scale both by the
 *    number of residue pairs, so the per-pair comparison below holds for
 *    profiles too).  Take a block of L residues in each group.  Aligning it
 *    scores  sum(s) - L*off.  The staggered alternative (block of a against
 *    gaps, then block of b against gaps) scores  -2*gpo - 2*L*gpe.  So the
 *    block is aligned iff  mean(s) > off - 2*gpe - 2*gpo/L;  for long blocks
 *    (the worst case) the threshold on the mean pair score is  off - 2*gpe.
 *    Without offset the threshold is -2*gpe.  With mu_bg the expected score
 *    of a random pair (sum_xy fa(x) fb(y) s(x,y), from the two groups'
 *    residue compositions) we require
 *
 *            off - 2*gpe <= mu_bg      i.e.   off <= 2*gpe + mu_bg = B.
 *
 *    With off <= B the offset can at most make the aligner indifferent
 *    between aligning and staggering a block that scores like random
 *    sequence; a block with any signal above random is still cheaper to
 *    align than to stagger.  If 2*gpe <= -mu_bg the gap penalties already
 *    stagger random blocks and B = 0: no offset is allowed.  The legacy
 *    offset (~1.05 with gpe 0.684 and PFASUM60, mu_bg ~ -1.14, B ~ 0.23)
 *    violated this by ~0.8 per pair.  Terminal gaps (tgpe) are not part of
 *    the bound: leaving terminal extensions unaligned is intended.
 *
 * 3. Offset.  off = min(kappa,1) * B * max(0, 1 - d/dzero).  kappa (fraction
 *    of the bound used for identical groups) and dzero (distance at which the
 *    offset reaches zero) are the free parameters (kalign_run_config
 *    vsm_kappa, vsm_dzero).  kappa is clamped to 1 so the bound holds for
 *    every setting.
 *
 * 4. Where the offset is applied.  In legacy mode the offset is added when a
 *    leaf profile is first built, so a profile-profile merge used whatever
 *    offset its leaves had at their first merge.  Here leaf profiles carry no
 *    offset; for each merge the offset of that merge is applied to temporary
 *    copies of the two profiles (vsm_shift_profiles) and, for sequence-
 *    sequence DPs, through aln_param.subm_offset.  The profiles merged by
 *    update_n stay offset-free. */

float vsm_corrected_distance(float p, int biotype)
{
        float x;
        if(p <= 0.0f){
                return 0.0f;
        }
        if(biotype == ALN_BIOTYPE_PROTEIN){
                x = 1.0f - p - 0.2f * p * p;
                if(x <= 1e-4f){
                        return KALIGN_VSM_DMAX;
                }
                x = -logf(x);
        }else{
                x = 1.0f - 4.0f / 3.0f * p;
                if(x <= 1e-4f){
                        return KALIGN_VSM_DMAX;
                }
                x = -0.75f * logf(x);
        }
        if(x > KALIGN_VSM_DMAX){
                x = KALIGN_VSM_DMAX;
        }
        return x;
}

float vsm_offset_bound(float gpe, float mu_bg)
{
        float b = 2.0f * gpe + mu_bg;
        return (b > 0.0f) ? b : 0.0f;
}

float vsm_offset(float d, float bound, float kappa, float dzero)
{
        float r;
        if(kappa <= 0.0f || dzero <= 0.0f || bound <= 0.0f){
                return 0.0f;
        }
        if(kappa > 1.0f){
                kappa = 1.0f;
        }
        r = 1.0f - d / dzero;
        if(r <= 0.0f){
                return 0.0f;
        }
        if(r > 1.0f){
                r = 1.0f;
        }
        return kappa * bound * r;
}

/* Number of standard residue letters (identity / composition ignore B, Z, X, N). */
static int vsm_n_letters(int biotype)
{
        return (biotype == ALN_BIOTYPE_PROTEIN) ? 20 : 4;
}

/* Residue composition of a profile (normalised over the standard letters). */
static void vsm_composition(const float* prof, int len, int nl, double* f)
{
        double tot = 0.0;
        int i, k;
        for(k = 0; k < nl; k++){
                f[k] = 0.0;
        }
        for(i = 1; i <= len; i++){
                const float* col = prof + 64 * i;
                for(k = 0; k < nl; k++){
                        f[k] += col[k];
                }
        }
        for(k = 0; k < nl; k++){
                tot += f[k];
        }
        if(tot > 0.0){
                for(k = 0; k < nl; k++){
                        f[k] /= tot;
                }
        }
}

static int vsm_debug(void)
{
        static _Atomic int flag = -1;   /* read by several threads; every writer stores the same value */
        if(flag < 0){
                const char* e = getenv("KALIGN_VSM_DEBUG");
                flag = (e != NULL && e[0] != '\0' && e[0] != '0') ? 1 : 0;
        }
        return flag;
}

int vsm_merge_offset(struct msa* msa, struct aln_tasks* t, const struct aln_param* ap,
                     int a, int b, int len_a, int len_b, float* off)
{ KPERF_REGION(KP_VSM);
        struct aln_mem* ml = NULL;
        struct aln_param lap;
        double fa[20];
        double fb[20];
        double mu_bg = 0.0;
        double num = 0.0;
        double den = 0.0;
        float bound;
        float p;
        float d;
        int nl;
        int i, j, k;

        *off = 0.0f;
        if(ap->vsm_mode != KALIGN_VSM_MERGE || ap->vsm_kappa <= 0.0f || ap->vsm_dzero <= 0.0f){
                return OK;
        }
        if(len_a <= 0 || len_b <= 0){
                return OK;
        }

        nl = vsm_n_letters(msa->biotype);
        vsm_composition(t->profile[a], len_a, nl, fa);
        vsm_composition(t->profile[b], len_b, nl, fb);
        for(i = 0; i < nl; i++){
                for(j = 0; j < nl; j++){
                        mu_bg += fa[i] * fb[j] * (double)ap->subm[i][j];
                }
        }
        bound = vsm_offset_bound(ap->gpe, (float)mu_bg);
        if(bound <= 0.0f){
                return OK;
        }

        /* Offset-free DP of a against b (no consistency bonus). */
        lap = *ap;
        lap.subm_offset = 0.0f;
        RUN(alloc_aln_mem(&ml, 256));
        ml->ap = &lap;
        ml->mode = ALN_MODE_FULL;
        ml->run_parallel = 0;
        ml->len_a = len_a;
        ml->len_b = len_b;
        RUN(init_alnmem(ml));
        RUN(dispatch_alignment(msa, ml, t->profile[a], t->profile[b], a, b, len_a, len_b));

        /* ml->path[i] (i = 1..len_a) is the column of b aligned to column i of a, or -1. */
        for(i = 1; i <= len_a; i++){
                const float* ca;
                const float* cb;
                double na = 0.0;
                double nb = 0.0;
                j = ml->path[i];
                if(j < 1 || j > len_b){
                        continue;
                }
                ca = t->profile[a] + 64 * i;
                cb = t->profile[b] + 64 * j;
                for(k = 0; k < nl; k++){
                        num += (double)ca[k] * (double)cb[k];
                        na += ca[k];
                        nb += cb[k];
                }
                den += na * nb;
        }
        free_aln_mem(ml);
        ml = NULL;

        if(den <= 0.0){
                return OK;
        }
        p = (float)(1.0 - num / den);
        d = vsm_corrected_distance(p, msa->biotype);
        *off = vsm_offset(d, bound, ap->vsm_kappa, ap->vsm_dzero);

        if(vsm_debug()){
                fprintf(stderr, "VSM2 nsip_a=%d nsip_b=%d p=%.3f d=%.3f mu_bg=%.3f bound=%.3f off=%.3f\n",
                        msa->nsip[a], msa->nsip[b], p, d, mu_bg, bound, *off);
        }
        return OK;
ERROR:
        if(ml){
                free_aln_mem(ml);
        }
        return FAIL;
}

static int vsm_shift_one(float* orig, int len, float off, float** copy)
{
        float* p = NULL;
        int i, k;
        MMALLOC(p, sizeof(float) * 64 * (len + 2));
        memcpy(p, orig, sizeof(float) * 64 * (len + 2));
        for(i = 1; i <= len; i++){
                float* col = p + 64 * i;
                float n = 0.0f;
                for(k = 0; k < 23; k++){
                        n += col[k];
                }
                n *= off;
                for(k = 0; k < 23; k++){
                        col[32 + k] -= n;
                }
        }
        *copy = p;
        return OK;
ERROR:
        return FAIL;
}

int vsm_shift_profiles(struct aln_tasks* t, int a, int b, int len_a, int len_b,
                       float off, float** save_a, float** save_b)
{ KPERF_REGION(KP_VSM);
        float* ca = NULL;
        float* cb = NULL;
        *save_a = NULL;
        *save_b = NULL;
        RUN(vsm_shift_one(t->profile[a], len_a, off, &ca));
        RUN(vsm_shift_one(t->profile[b], len_b, off, &cb));
        *save_a = t->profile[a];
        *save_b = t->profile[b];
        t->profile[a] = ca;
        t->profile[b] = cb;
        return OK;
ERROR:
        if(ca){
                MFREE(ca);
        }
        return FAIL;
}

void vsm_unshift_profiles(struct aln_tasks* t, int a, int b,
                          float** save_a, float** save_b)
{ KPERF_REGION(KP_VSM);
        if(*save_a){
                MFREE(t->profile[a]);
                t->profile[a] = *save_a;
                *save_a = NULL;
        }
        if(*save_b){
                MFREE(t->profile[b]);
                t->profile[b] = *save_b;
                *save_b = NULL;
        }
}
