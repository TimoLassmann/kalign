#ifndef LIBKALIGN_H
#define LIBKALIGN_H

#ifdef LIBKALIGN_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN int kalign(char **seq, int *len,int numseq, char ***aligned, int *out_aln_len);
#undef LIBKALIGN_IMPORT
#undef EXTERN


#endif
