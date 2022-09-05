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

#ifndef ALIGNMENT_PARAMETERS_H
#define ALIGNMENT_PARAMETERS_H

//#include "global.h"
#include "parameters.h"
//#include "misc.h"
#include <stdint.h>
struct rng_state;

#define ALN_MODE_SCORE_ONLY 2
#define ALN_MODE_FULL 1

struct aln_param{
        /* struct rng_state* rng; */
        /* run specific settings */
        /* const float* prof1; */
        /* const float* prof2; */
        /* const uint8_t* seq1; */
        /* const uint8_t* seq2; */

        /* int sip; */
        /* int mode; */

        int nthreads;
        int chaos;
        /* actual parameters  */
        float** subm;
        /* int* tree; */
        float gpo;
        float gpe;
        float tgpe;
        float score;

};

extern int init_ap(struct aln_param** aln_param, struct parameters* param,int L);
extern void free_ap(struct aln_param* ap);
#endif
