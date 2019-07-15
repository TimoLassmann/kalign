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

#include "global.h"

#define defPROTEIN 21
#define redPROTEIN 13
#define defDNA 5


struct alphabet{
        int8_t to_internal[128];
        int8_t to_external[32];
        int type;
        int L;
};

extern struct alphabet* create_alphabet(int type);
extern int switch_alphabet(struct alphabet* a, int type);



#endif
