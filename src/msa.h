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

#ifndef MSA_H
#define MSA_H


#define MSA_NAME_LEN 128
#define FORMAT_FA 1
#define FORMAT_MSF 2
#define FORMAT_CLU 3

#include <stdint.h>

struct msa_seq{
        char* name;
        char* seq;
        uint8_t* s;
        int* gaps;
        int len;
        int alloc_len;
};

struct msa{
        struct msa_seq** sequences;
        int** sip;
        int* nsip;
        int* plen;
        int numseq;
        int num_profiles;
        int alloc_numseq;
        int aligned;
        int letter_freq[128];
        int L;



};

/* dealign */
int dealign_msa(struct msa* msa);


/* convert */
int convert_msa_to_internal(struct msa* msa, int type);
/* rw functions */

struct msa* read_input(char* infile,struct msa* msa);
int write_msa(struct msa* msa, char* outfile, int type);
void free_msa(struct msa* msa);

#endif
