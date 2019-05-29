#ifndef GLOBAL_H
#define GLOBAL_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>


#include "tldevel.h"

struct sequence_info{
        struct sequence_info* next;
        char* name;
        char* value;
};

struct feature{
        struct feature *next;
        char* type;
        char* note;
        int start;
        int end;
        int color;
};

struct alignment{
        struct feature** ft;
        struct sequence_info** si;
        unsigned int** sip;
        unsigned int* nsip;
        unsigned int* sl;
        unsigned int* lsn;
        uint8_t** s;
        int** gaps;
        char** seq;
        char** sn;
        int aln_len;
        int max_len;
        int numseq;
        int num_profiles;
        int dna;
        int L;
};
#endif
