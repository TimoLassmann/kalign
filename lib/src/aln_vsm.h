#ifndef ALN_VSM_H
#define ALN_VSM_H

#ifdef ALN_VSM_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

/* Variable scoring matrix (VSM) modes (aln_param.vsm_mode /
   kalign_run_config.vsm_mode). */
#define KALIGN_VSM_LEGACY 0     /* offset = max(0, vsm_amax - mean distance-to-anchors),
                                   baked into leaf profiles (behaviour up to v3.6) */
#define KALIGN_VSM_MERGE  1     /* per-merge, divergence-dependent, bounded offset */

struct msa;
struct aln_param;
struct aln_tasks;

/* Evolutionary distance from a p-distance (1 - identity):
   protein: Kimura (1983) d = -ln(1 - p - 0.2 p^2);
   nucleotide: Jukes-Cantor d = -3/4 ln(1 - 4/3 p).
   Returns KALIGN_VSM_DMAX once the correction diverges. */
EXTERN float vsm_corrected_distance(float p, int biotype);

/* Largest admissible offset for a merge: max(0, 2*gpe + mu_bg).
   See the derivation in aln_vsm.c. */
EXTERN float vsm_offset_bound(float gpe, float mu_bg);

/* offset(d) = min(kappa,1) * bound * max(0, 1 - d/dzero). Never exceeds bound. */
EXTERN float vsm_offset(float d, float bound, float kappa, float dzero);

/* Per-merge offset for aligning groups a and b (KALIGN_VSM_MERGE only;
   returns 0 otherwise). Profiles t->profile[a], t->profile[b] must be built
   (and gap penalties set) and must not carry an offset. ap holds the gap
   penalties used for this merge's DP. Runs one offset-free DP of a vs b to
   measure the cross-group identity of the aligned residue pairs. */
EXTERN int vsm_merge_offset(struct msa* msa, struct aln_tasks* t, const struct aln_param* ap,
                            int a, int b, int len_a, int len_b, float* off);

/* Swap t->profile[a] and t->profile[b] for copies whose substitution scores
   (slots 32..54) are lowered by off * (residues in that column), so that the
   profile DPs charge off once per aligned residue pair. The originals are
   returned in save_a/save_b and put back by vsm_unshift_profiles(). */
EXTERN int vsm_shift_profiles(struct aln_tasks* t, int a, int b, int len_a, int len_b,
                              float off, float** save_a, float** save_b);
EXTERN void vsm_unshift_profiles(struct aln_tasks* t, int a, int b,
                                 float** save_a, float** save_b);

#define KALIGN_VSM_DMAX 10.0f

#undef ALN_VSM_IMPORT
#undef EXTERN

#endif
