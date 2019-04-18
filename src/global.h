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
        int** s;
        char**seq;
        char** sn;
        int numseq;
        int num_profiles;
};

#endif
