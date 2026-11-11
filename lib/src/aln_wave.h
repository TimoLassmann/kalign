#ifndef ALN_WAVE_H
#define ALN_WAVE_H

#ifdef ALN_WAVE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct aln_mem;

/* KALIGN_WAVE (default ON): see aln_wave.c.  prepare returns 1 if it
   attached tables to m (the caller then owns them and must call
   aln_wave_free); aln_wave_ss/pp return 0 if they ran the pass, 1 if the
   caller must run the scalar kernel. */
EXTERN int aln_wave_prepare(struct aln_mem* m);
EXTERN void aln_wave_free(struct aln_mem* m);
EXTERN int aln_wave_ss(struct aln_mem* m, int dir);
EXTERN int aln_wave_pp(struct aln_mem* m, int dir);
EXTERN int aln_wave_sp(struct aln_mem* m, int dir);

#undef ALN_WAVE_IMPORT
#undef EXTERN

#endif
