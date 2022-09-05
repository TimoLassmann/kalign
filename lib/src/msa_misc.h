#ifndef MSA_MISC_H
#define MSA_MISC_H

#ifdef MSA_MISC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct msa;
EXTERN int GCGMultchecksum(struct msa* msa);
EXTERN int GCGchecksum(char *seq, int len);

#undef MSA_MISC_IMPORT
#undef EXTERN


#endif
