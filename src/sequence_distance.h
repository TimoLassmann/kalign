#ifndef SEQUENCE_DISTANCE_H
#define SEQUENCE_DISTANCE_H

#include "global.h"

#include "alignment_parameters.h"

extern float** aln_distance(struct alignment* aln,struct aln_param* ap);

extern float** bpm_distance_thin(struct alignment* aln, int* seeds, int num_seeds);
extern float** bpm_distance_pair(struct alignment* aln, int* selection, int num_sel);

extern float** bpm_distance(struct alignment* aln, int* seeds, int num_seeds);


extern float** kmer_distance(struct alignment* aln, int* seeds, int num_seeds, int kmer_len);

extern float** kmer_bpm_distance(struct alignment* aln, int kmer_len, int num_seeds);

extern float** protein_wu_distance(struct alignment* aln, float zlevel, int nj,int* seeds, int num_anchors);

extern float** dna_distance(struct alignment* aln, float zlevel, int nj);




#endif
