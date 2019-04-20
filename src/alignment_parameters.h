#ifndef ALIGNMENT_PARAMETERS_H
#define ALIGNMENT_PARAMETERS_H

#include "global.h"
#include "parameters.h"
#include "misc.h"
struct aln_param{
        float** subm;
        float gpo;
        float gpe;
        float tgpe;

        int* tree;
};


extern struct aln_param* init_ap(struct parameters* param,int numseq);
extern void free_ap(struct aln_param* ap);
#endif
