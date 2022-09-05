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

#ifndef SEQUENCE_DISTANCE_H
#define SEQUENCE_DISTANCE_H

#include <stdint.h>

#ifdef SEQUENCE_DISTANCE_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif





struct msa;
/* #include "alignment_parameters.h" */

EXTERN float** d_estimation(struct msa* msa, int* samples, int num_samples,int pair);
EXTERN float calc_distance(uint8_t *seq_a, uint8_t *seq_b, int len_a,
                           int len_b);

#undef SEQUENCE_DISTANCE_IMPORT
#undef EXTERN

#endif
