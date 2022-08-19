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

#ifndef MISC_H
#define MISC_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>


#include <stdint.h>


struct rng_state;

extern int byg_start(char* pattern,char*text);
extern int byg_end(char* pattern,char*text);
extern int byg_count(char* pattern,char*text);

extern int shuffle_arr_r(int* arr,int n, struct rng_state* rng);

/* The following two hash functions are taken from the supplementary of: */

/* Steinegger, Martin, and Johannes Söding. "Clustering huge protein sequence sets in linear time." Nature communications 9.1 (2018): 2542. */

// (c) 2017 Johannes Soeding & Martin Steinegger, Gnu Public License version 3
//unsigned circ_hash(const int * x, unsigned length);
uint16_t circ_hash(const uint8_t* x, const uint8_t length);
//unsigned circ_hash_next(const int * x, unsigned length, int x_first, short unsigned h);
uint16_t circ_hash_next(const uint8_t * x,const uint8_t length,const uint8_t x_first, uint16_t h);

extern char* basename(const char* name);

#endif
