#ifndef ALN_CONTROLLER_H
#define ALN_CONTROLLER_H

#include <stdint.h>

#ifdef ALN_CONTROLLER_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int aln_runner(struct aln_mem* m,struct aln_param* ap,int* path);

#undef ALN_CONTROLLER_IMPORT
#undef EXTERN
#endif
