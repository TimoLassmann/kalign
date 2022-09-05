#ifndef ALN_CONTROLLER_H
#define ALN_CONTROLLER_H

#include <stdint.h>

#ifdef ALN_CONTROLLER_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct aln_mem;
EXTERN int aln_runner(struct aln_mem* m);
EXTERN int aln_runner_serial(struct aln_mem* m);

#undef ALN_CONTROLLER_IMPORT
#undef EXTERN
#endif
