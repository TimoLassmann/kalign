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

#ifndef ALN_DATA_H
#define ALN_DATA_H

#include "tldevel.h"


struct kalign_sequence{
        int* s;
        char* seq;
        int alloc_seq_len;
        int len;
        uint8_t* name;
        int id;
};


struct kalign_alignmment{
        struct kalign_sequence** s_arr;
        int numseq;
        int alloc_numseq;
};

extern struct kalign_alignmment* kalign_aln_alloc(void);
extern int kalign_alignment_resize(struct kalign_alignmment* aln);
extern void kalign_alignmment_free(struct kalign_alignmment* aln);


extern int kalign_seq_resize(struct kalign_sequence* ks);

#endif
