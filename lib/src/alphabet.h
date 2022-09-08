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

#ifndef ALPHABET_H
#define ALPHABET_H

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif

/* #include "global.h" */

#include <inttypes.h>

#define ALPHA_defPROTEIN 21
#define ALPHA_ambigiousPROTEIN 23
#define ALPHA_redPROTEIN 13
#define ALPHA_redPROTEIN2 8
#define ALPHA_defDNA 5

#define ALPHA_UNKNOWN 255
#define ALPHA_UNDEFINED -1

struct alphabet{
        int8_t to_internal[128];
        int8_t to_external[32];
        int type;
        int L;
};


kalign_extern struct alphabet* create_alphabet(int type);
kalign_extern int switch_alphabet(struct alphabet* a, int type);

#endif
