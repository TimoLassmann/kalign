#ifndef WEAVE_ALIGNMENT_H
#define WEAVE_ALIGNMENT_H

#include "global.h"

extern int weave(struct alignment* aln, int** map, int* tree);

extern int clean_aln(struct alignment* aln);

#endif
