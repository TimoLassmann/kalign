#ifndef IDATA_H
#define IDATA_H

#ifdef IDATA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct msa;
struct aln_param;

EXTERN int get_internal_data(struct msa* msa, struct aln_param* ap,double** out,int* n_out);

#endif
