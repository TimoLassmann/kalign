#ifndef COUNTS_FROM_RANDOM_TREES_H
#define COUNTS_FROM_RANDOM_TREES_H




#include "global.h"
#include "alignment_parameters.h"
#include "bisectingKmeans.h"
#include "alignment.h"
#include "weave_alignment.h"
#include "align_io.h"


extern int counts_from_random_trees(struct alignment* aln, struct aln_param* ap, int num_iter);
#endif
