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

#include "global.h"
#include "msa.h"

#include "alignment_parameters.h"

extern float** d_estimation(struct msa* msa, int* samples, int num_samples,int pair);

extern float** aln_distance(struct msa* msa,struct aln_param* ap);

extern float** bpm_distance_thin(struct msa* msa, int* seeds, int num_seeds);
extern float** bpm_distance_pair(struct msa* msa, int* selection, int num_sel);

extern float** bpm_distance(struct msa* msa, int* seeds, int num_seeds);


extern float** kmer_distance(struct msa* msa, int* seeds, int num_seeds, int kmer_len);

extern float** kmer_bpm_distance(struct msa* msa, int kmer_len, int num_seeds);

extern float** protein_wu_distance(struct msa* msa, float zlevel,int* seeds, int num_anchors);

extern float** dna_wu_distance(struct msa* msa, float zlevel, int* seeds, int num_anchors);




#endif
