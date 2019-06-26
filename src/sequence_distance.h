#ifndef SEQUENCE_DISTANCE_H
#define SEQUENCE_DISTANCE_H

#include "global.h"
#include "msa.h"

#include "alignment_parameters.h"

extern float** aln_distance(struct msa* msa,struct aln_param* ap);

extern float** bpm_distance_thin(struct msa* msa, int* seeds, int num_seeds);
extern float** bpm_distance_pair(struct msa* msa, int* selection, int num_sel);

extern float** bpm_distance(struct msa* msa, int* seeds, int num_seeds);


extern float** kmer_distance(struct msa* msa, int* seeds, int num_seeds, int kmer_len);

extern float** kmer_bpm_distance(struct msa* msa, int kmer_len, int num_seeds);

extern float** protein_wu_distance(struct msa* msa, float zlevel, int nj,int* seeds, int num_anchors);

extern float** dna_distance(struct msa* msa, float zlevel, int nj);




#endif
