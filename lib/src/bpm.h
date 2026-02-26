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


/* #ifdef HAVE_AVX2 */
/* #define BPM(a,b,len_a,len_b) bpm_256(a,b,len_a,len_b) */
/* #else */
#define BPM(a,b,len_a,len_b) bpm_block(a,b,len_a,len_b)
/* #endif */

/* #define LOG_MSG(...) do {                       \ */
/*                 log_message( __VA_ARGS__ );     \ */
/*         }while (0) */


/* Must be called before bpm_256!!!!  */
EXTERN  void set_broadcast_mask(void);

EXTERN uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m);
EXTERN uint8_t bpm(const uint8_t* t,const uint8_t* p,int n,int m);

EXTERN int bpm_block(const uint8_t *t, const uint8_t *p, int n, int m);
EXTERN uint8_t dyn_256(const uint8_t* t,const uint8_t* p,int n,int m);
#undef BPM_IMPORT
#undef EXTERN

#endif
