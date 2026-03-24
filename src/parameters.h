#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdint.h>


#ifdef PARAMETERS_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct parameters{
        char **infile;
        char *input;
        char *outfile;
        char* format;
        int type;
        float gpo;
        float gpe;
        float tgpe;
        int num_infiles;
        int out_format;
        int nthreads;
        int min_support;
        char* load_poar;
        char* mode;  /* "fast", "default", "recall", "accurate" (NULL = default) */
        int help_flag;
        int quiet;
        int dump_internal;
        int reformat;
        int rename;
        int clean;
        int unalign;
};

EXTERN struct parameters* init_param(void);
EXTERN void free_parameters(struct parameters *param);

EXTERN int check_msa_format_string(char* format);

#undef PARAMETERS_IMPORT
#undef EXTERN

#endif
