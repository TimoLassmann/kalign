#ifndef PROBMSA_PROFILE_H
#define PROBMSA_PROFILE_H

#include <stdint.h>
#include "probmsa_submat.h"
#include "probmsa_dirichlet.h"

#ifdef PROBMSA_PROFILE_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct probmsa_column {
        double prob[PROBMSA_ALPHA_MAX];  /* probability distribution over alphabet */
        double quality;                   /* confidence (high = certain) */
        int depth;                        /* number of sequences contributing */
        double gap_frac;                  /* fraction of seqs with gap here */
        double concentration;             /* effective Dirichlet concentration (data + prior) */
};

struct probmsa_profile {
        struct probmsa_column* cols;
        int length;
        int n_sequences;
        int alpha_size;
};

/* Create leaf profile from encoded sequence (kalign uint8_t* s).
   If diri is non-NULL, apply Dirichlet mixture regularization. */
EXTERN int probmsa_profile_from_seq(const uint8_t* seq, int len, int alpha_size,
                                     const struct probmsa_diri_mix* diri,
                                     struct probmsa_profile** out);

/* Create merged profile from MEA path and match posteriors.
   posteriors[i][j] = P(Match at i,j | data), stored in the Y matrix after MEA.
   gap_post_a[i] = P(position i in A is gapped), 0-indexed.
   gap_post_b[j] = P(position j in B is gapped), 0-indexed.
   path format: path[0]=aln_len, path[1..aln_len]=0/1/2, path[aln_len+1]=3. */
EXTERN int probmsa_profile_merge(struct probmsa_profile* a,
                                  struct probmsa_profile* b,
                                  const int* path,
                                  double** posteriors,
                                  const double* gap_post_a,
                                  const double* gap_post_b,
                                  int alpha_size,
                                  const struct probmsa_diri_mix* diri,
                                  struct probmsa_profile** out);

EXTERN void probmsa_profile_free(struct probmsa_profile* p);

#undef PROBMSA_PROFILE_IMPORT
#undef EXTERN

#endif
