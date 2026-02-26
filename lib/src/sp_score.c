#include "tldevel.h"

#include <string.h>
#include <stdint.h>

#include "msa_struct.h"
#include "aln_param.h"

#define SP_SCORE_IMPORT
#include "sp_score.h"

#define SP_ALPHA 23

/* Build a frequency profile from a group of sequences.
   For each profile column, counts residues (0..SP_ALPHA-1) and gaps.

   freq: pre-allocated array of prof_len * SP_ALPHA ints (zeroed)
   n_gap: pre-allocated array of prof_len ints (zeroed)

   Sequences are expanded one at a time via gap arrays to avoid
   O(N * L) simultaneous memory allocation. */
static int build_profile(struct msa* msa, int* sip, int nsip,
                         int prof_len, int* freq, int* n_gap)
{
        int i, j, k, pos, si;
        int8_t* cols = NULL;

        for(i = 0; i < nsip; i++){
                si = sip[i];
                /* Expand gap array into column array */
                MMALLOC(cols, sizeof(int8_t) * prof_len);
                pos = 0;
                for(j = 0; j < msa->sequences[si]->len; j++){
                        for(k = 0; k < msa->sequences[si]->gaps[j]; k++){
                                cols[pos++] = -1;
                        }
                        cols[pos++] = (int8_t)msa->sequences[si]->s[j];
                }
                for(k = 0; k < msa->sequences[si]->gaps[msa->sequences[si]->len]; k++){
                        cols[pos++] = -1;
                }

                /* Accumulate into frequency profile */
                for(j = 0; j < prof_len; j++){
                        if(cols[j] >= 0 && cols[j] < SP_ALPHA){
                                freq[j * SP_ALPHA + (int)cols[j]]++;
                        }else{
                                n_gap[j]++;
                        }
                }
                MFREE(cols);
                cols = NULL;
        }

        return OK;
ERROR:
        if(cols) MFREE(cols);
        return FAIL;
}

/* Compute profile-based SP score for cross-group sequence pairs.

   Instead of iterating all O(|A|*|B|) pairs, builds residue frequency
   profiles for each group and computes the expected score analytically:
     substitution: sum_r sum_s freq_a[r] * freq_b[s] * subm[r][s]
     gap penalty:  gpo at path-level gap opens, gpe/tgpe per gap-residue pair

   The substitution component is exact (cross-group pairs are independent).
   Gap-open penalties are tracked at the path level (gap-in-A/B runs) which
   captures the dominant inter-group gap structure. Internal gap opens within
   profiles are not tracked (would require O(N^2) per-pair state).

   Complexity: O((|A|+|B|) * L) for profile construction
             + O(SP_ALPHA^2 * path_len) for scoring. */
int compute_sp_score(struct msa* msa, struct aln_param* ap,
                     int* path, int* sip_a, int nsip_a,
                     int* sip_b, int nsip_b, float* score)
{
        int* freq_a = NULL;
        int* freq_b = NULL;
        int* gap_a = NULL;
        int* gap_b = NULL;
        int path_len;
        int i, j, c, si;
        int pos_a, pos_b;
        float total = 0.0F;

        const float gpo = ap->gpo;
        const float gpe = ap->gpe;
        const float tgpe = ap->tgpe;
        float** subm = ap->subm;

        path_len = path[0];

        /* Compute profile lengths from first sequence in each group */
        si = sip_a[0];
        int prof_a_len = msa->sequences[si]->len;
        for(i = 0; i <= msa->sequences[si]->len; i++){
                prof_a_len += msa->sequences[si]->gaps[i];
        }

        si = sip_b[0];
        int prof_b_len = msa->sequences[si]->len;
        for(i = 0; i <= msa->sequences[si]->len; i++){
                prof_b_len += msa->sequences[si]->gaps[i];
        }

        /* Allocate and build frequency profiles */
        MMALLOC(freq_a, sizeof(int) * prof_a_len * SP_ALPHA);
        MMALLOC(gap_a, sizeof(int) * prof_a_len);
        memset(freq_a, 0, sizeof(int) * prof_a_len * SP_ALPHA);
        memset(gap_a, 0, sizeof(int) * prof_a_len);
        RUN(build_profile(msa, sip_a, nsip_a, prof_a_len, freq_a, gap_a));

        MMALLOC(freq_b, sizeof(int) * prof_b_len * SP_ALPHA);
        MMALLOC(gap_b, sizeof(int) * prof_b_len);
        memset(freq_b, 0, sizeof(int) * prof_b_len * SP_ALPHA);
        memset(gap_b, 0, sizeof(int) * prof_b_len);
        RUN(build_profile(msa, sip_b, nsip_b, prof_b_len, freq_b, gap_b));

        /* Walk the path and score using profiles.
           Track path-level gap runs for gap-open penalties. */
        pos_a = 0;
        pos_b = 0;
        int in_a_gap = 0;   /* currently in a gap-in-A run */
        int in_b_gap = 0;   /* currently in a gap-in-B run */
        for(c = 1; c <= path_len; c++){
                int step = path[c] & 3;
                int is_terminal = path[c] & 32;
                float pen = is_terminal ? tgpe : gpe;

                if(step == 0){
                        /* Match: both profiles advance */
                        int* fa = freq_a + pos_a * SP_ALPHA;
                        int* fb = freq_b + pos_b * SP_ALPHA;

                        /* Substitution: exact cross-group sum */
                        for(i = 0; i < SP_ALPHA; i++){
                                if(fa[i] == 0) continue;
                                for(j = 0; j < SP_ALPHA; j++){
                                        if(fb[j] == 0) continue;
                                        total += (float)(fa[i] * fb[j]) * subm[i][j];
                                }
                        }

                        /* Gap penalty: residue-gap cross pairs (gpe only,
                           internal gap opens not tracked) */
                        int n_res_a = nsip_a - gap_a[pos_a];
                        int n_gap_b = gap_b[pos_b];
                        int n_gap_a = gap_a[pos_a];
                        int n_res_b = nsip_b - gap_b[pos_b];
                        total -= (float)(n_res_a * n_gap_b + n_gap_a * n_res_b) * pen;

                        in_a_gap = 0;
                        in_b_gap = 0;
                        pos_a++;
                        pos_b++;
                }else if(step == 1){
                        /* Gap in A: only B advances.
                           All A sequences are gapped at the path level.
                           Only B sequences with residues contribute penalty. */
                        int n_res_b = nsip_b - gap_b[pos_b];
                        int n_pairs = nsip_a * n_res_b;
                        if(!in_a_gap){
                                total -= (float)n_pairs * gpo;
                        }
                        total -= (float)n_pairs * pen;
                        in_a_gap = 1;
                        in_b_gap = 0;
                        pos_b++;
                }else if(step == 2){
                        /* Gap in B: only A advances.
                           All B sequences are gapped at the path level.
                           Only A sequences with residues contribute penalty. */
                        int n_res_a = nsip_a - gap_a[pos_a];
                        int n_pairs = n_res_a * nsip_b;
                        if(!in_b_gap){
                                total -= (float)n_pairs * gpo;
                        }
                        total -= (float)n_pairs * pen;
                        in_a_gap = 0;
                        in_b_gap = 1;
                        pos_a++;
                }
        }

        *score = total;

        MFREE(freq_a);
        MFREE(gap_a);
        MFREE(freq_b);
        MFREE(gap_b);

        return OK;
ERROR:
        if(freq_a) MFREE(freq_a);
        if(gap_a) MFREE(gap_a);
        if(freq_b) MFREE(freq_b);
        if(gap_b) MFREE(gap_b);
        return FAIL;
}
