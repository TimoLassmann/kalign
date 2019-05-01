#ifndef SEQUENCE_DISTANCE_H
#define SEQUENCE_DISTANCE_H

#include "global.h"

extern float** kmer_distance(struct alignment* aln, int* seeds, int num_seeds, int kmer_len);

extern float** protein_wu_distance(struct alignment* aln, float zlevel, int nj,int* seeds, int num_anchors);

extern float** dna_distance(struct alignment* aln, float zlevel, int nj);




#endif
