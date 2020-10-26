#ifndef ALN_MEM_H
#define ALN_MEM_H

#ifdef ALN_MEM_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct aln_mem;

EXTERN int alloc_aln_mem(struct aln_mem** mem, int x);
EXTERN int resize_aln_mem(struct aln_mem* m);
EXTERN void free_aln_mem(struct aln_mem* m);

#undef ALN_MEM_IMPORT
#undef EXTERN
#endif
