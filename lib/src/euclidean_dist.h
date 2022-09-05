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

#ifndef EUCLIDIAN_DIST_H
#define EUCLIDIAN_DIST_H

#ifdef EUCLIDEAN_DIST_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

EXTERN int edist_256(const float* a,const float* b, const int len, float* ret);
EXTERN int edist_serial(const float* a,const float* b,const int len, float* ret);
EXTERN int edist_serial_d(const double* a,const double* b,const int len, double* ret);


#undef EUCLIDEAN_DIST_IMPORT
#undef EXTERN

#endif
