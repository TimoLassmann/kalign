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


#ifndef PHMM_H
#define PHMM_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>


#include "tldevel.h"


#define INDEXMM 0
#define INDEXGPO 1
#define INDEXGPE 2
#define INDEXTM 3

struct phmm{

        float** fM ;
        float** fX ;
        float** fY ;

        float** bM ;
        float** bX ;
        float** bY ;

        float emit_M[23][23];
        float emit_background[23];

        float emit_M_e[23][23];
        float emit_background_e[23];

        float transition[4];
        float transition_e[4];

        float tau;
        float tau_e;
        float f_score;
        float b_score;
        float eta;
        int alloc_x;
        int alloc_y;
        int L;                  /* alphabet len */
};

int forward_phmm(struct phmm* phmm,uint8_t* seq_a,uint8_t* seq_b, int len_a,int len_b);
int backward_phmm(struct phmm* phmm,uint8_t* seq_a,uint8_t* seq_b, int len_a,int len_b);
int collect_phmm(struct phmm* phmm,uint8_t* seq_a,uint8_t* seq_b, int len_a,int len_b);
int re_estimate(struct phmm* phmm);

int phmm_transitions(struct phmm* phmm);
int print_phmm(struct phmm* phmm,int len_a,int len_b);
int add_pseudocounts(struct phmm* phmm, float w);
int clear_phmm_e(struct phmm* phmm);

int simple_init(struct phmm*phmm);
/* memory functions */
struct phmm* alloc_phmm(int size);

void free_phmm(struct phmm* phmm);



#endif
