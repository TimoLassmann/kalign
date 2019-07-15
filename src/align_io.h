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

#ifndef ALIGN_IO_H
#define ALIGN_IO_H

#include <unistd.h>
#include "parameters.h"

#define SEEK_START 0
#define SEEK_END 2




extern struct alignment* read_alignment(char* infile);
extern struct alignment* detect_and_read_sequences(struct parameters* param);
extern int make_dna(struct alignment* aln);
extern void free_aln(struct alignment* aln);

extern int convert_alignment_to_internal(struct alignment* aln, int type);

extern int dealign(struct alignment* aln);

extern int output(struct alignment* aln,struct parameters* param);

extern int make_aliged_seq(uint8_t* aligned, uint8_t* unaligned, int* gaps,int len);

#endif
