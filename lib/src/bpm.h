/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

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
