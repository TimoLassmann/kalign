#ifndef DSSIM_H
#define DSSIM_H

#include "kalign/kalign.h"

#ifdef DSSIM_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int dssim_get_fasta(struct msa **msa, int n_seq, int n_obs, int dna,int len, int seed);
/* EXTERN int dssim_get_fasta(struct msa** msa, int n_seq, int n_obs,int len,int seed); */

#undef DSSIM_IMPORT
#undef EXTERN


#endif
