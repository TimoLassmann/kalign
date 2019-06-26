#ifndef WEAVE_ALIGNMENT_H
#define WEAVE_ALIGNMENT_H

#include "global.h"
#include "msa.h"
extern int weave(struct msa* msa, int** map, int* tree);

extern int clean_aln(struct msa* msa);

#endif
