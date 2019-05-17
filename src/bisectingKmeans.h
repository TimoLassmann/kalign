#ifndef BISECTINGKMEANS_H
#define BISECTINGKMEANS_H

#include "global.h"

#include "rng.h"
#include "alignment_parameters.h"


extern int random_tree(struct aln_param* ap, int numseq);
extern int build_tree_kmeans(struct alignment* aln,struct parameters* param, struct aln_param* ap);
#endif
