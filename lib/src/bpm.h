#ifndef BPM_H
#define BPM_H

#include <stdint.h>
#ifdef BPM_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif


/* BPM picks the SIMD kernel at run time (inside bpm_block) when the CPU has
   it; no caller needs to know. */
#define BPM(a,b,len_a,len_b) bpm_block(a,b,len_a,len_b)

/* #define LOG_MSG(...) do {                       \ */
/*                 log_message( __VA_ARGS__ );     \ */
/*         }while (0) */


/* Kept for compatibility; does nothing.  Safe to call on any CPU, any number
   of times, from any thread.  bpm_256 needs no set-up. */
EXTERN  void set_broadcast_mask(void);

EXTERN uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m);
EXTERN uint8_t bpm(const uint8_t* t,const uint8_t* p,int n,int m);

EXTERN int bpm_block(const uint8_t *t, const uint8_t *p, int n, int m);
/* Plain-C versions of bpm_block and bpm_256: same results, never use SIMD. */
EXTERN int bpm_block_plain(const uint8_t *t, const uint8_t *p, int n, int m);
EXTERN uint8_t bpm_256_plain(const uint8_t* t,const uint8_t* p,int n,int m);
EXTERN uint8_t dyn_256(const uint8_t* t,const uint8_t* p,int n,int m);
#undef BPM_IMPORT
#undef EXTERN

#endif
