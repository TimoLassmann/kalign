#ifndef PROBMSA_DIRICHLET_H
#define PROBMSA_DIRICHLET_H

#include "tldevel.h"

#define DIRI_NUM_COMPONENTS 20
#define DIRI_ALPHA_SIZE 20

/*
 * Sjölander 20-component Dirichlet mixture prior for protein sequences.
 * From: Sjölander, Karplus, Brown, Hughey, Krogh, Mian, Haussler (1996)
 * "Dirichlet mixtures: a method for improved detection of weak but
 *  significant protein sequence homology" CABIOS 12(4):327-345
 *
 * Parameters from the MEME suite recode3.20comp (UCSC SAM toolkit).
 * Alphabet order: ACDEFGHIKLMNPQRSTVWY (matches kalign internal encoding).
 */

struct probmsa_diri_mix {
        int n_comp;
        int alpha_size;
        double q[DIRI_NUM_COMPONENTS];                              /* mixing weights */
        double alpha[DIRI_NUM_COMPONENTS][DIRI_ALPHA_SIZE];         /* pseudocounts */
        double alpha_sum[DIRI_NUM_COMPONENTS];                      /* precomputed |alpha_k| */
        double log_gamma_alpha[DIRI_NUM_COMPONENTS][DIRI_ALPHA_SIZE]; /* precomputed lgamma(alpha_ki) */
        double log_gamma_alpha_sum[DIRI_NUM_COMPONENTS];            /* precomputed lgamma(|alpha_k|) */
};

/* Initialize the 20-component Sjölander mixture.
 * prior_scale: multiplier for all alpha pseudocounts.
 *   1.0 = full Sjölander prior (alpha sums ~5-11, suitable for HMMER/MEME)
 *   0.1 = gentle regularization (alpha sums ~0.5-1.1, for progressive MSA)
 *   0.0 = no prior (not recommended, use NULL diri instead) */
#ifdef PROBMSA_DIRICHLET_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int probmsa_diri_init(struct probmsa_diri_mix** out, double prior_scale);
EXTERN void probmsa_diri_free(struct probmsa_diri_mix* m);

/*
 * Compute posterior mean amino acid probabilities given observed counts.
 *
 * counts[i]: number of times amino acid i was observed (can be fractional)
 * n_total:   sum of counts (= |n|)
 * prob_out[i]: output posterior mean probability for amino acid i
 *
 * Formula (Sjölander 1996, Appendix A.4):
 *   p_i = sum_k P(k|n) * (n_i + alpha_ki) / (|n| + |alpha_k|)
 * where P(k|n) = q_k * Prob(n|alpha_k) / sum_j q_j * Prob(n|alpha_j)
 */
EXTERN int probmsa_diri_posterior_mean(const struct probmsa_diri_mix* m,
                                        const double* counts, double n_total,
                                        double* prob_out,
                                        double* concentration_out);

/*
 * Compute joint probability matrix from Dirichlet mixture components.
 * joint[s][t] = Σ_k w_k * (α_k_s / |α_k|) * (α_k_t / |α_k|)
 * bg[s] = Σ_k w_k * (α_k_s / |α_k|)
 *
 * Captures co-occurrence patterns: amino acids that tend to appear in
 * the same column type get higher joint probability. Provides a
 * substitution model derived purely from the Dirichlet prior (no BLOSUM).
 *
 * Output is in Dirichlet alphabet order (ACDEFGHIKLMNPQRSTVWY).
 * Caller must remap to kalign-23 if needed.
 */
EXTERN int probmsa_diri_joint_matrix(const struct probmsa_diri_mix* m,
                                      double joint[][DIRI_ALPHA_SIZE],
                                      double* background);

#undef EXTERN
#endif /* PROBMSA_DIRICHLET_H */
