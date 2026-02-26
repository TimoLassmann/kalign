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


#define KALIGNDIST_ALN 0
#define KALIGNDIST_BPM 1
#define KALIGNDIST_WU 2

struct parameters{
        char **infile;
        char *input;
        char *outfile;
        char* format;
        char* aln_param_file;
        int type;
        float gpo;
        float gpe;
        float tgpe;
        float matadd;
        int chaos;
        int out_format;
        int param_set;
        int dist_method;
        int num_infiles;
        int reformat;
        int rename;             /* rename sequences - to make bali_score swallow the alignments */
        int dump_internal;
        int nthreads;
        int clean;
        int unalign;
        int refine;
        int adaptive_budget;
        int ensemble;
        uint64_t ensemble_seed;
        int min_support;
        char* save_poar;
        char* load_poar;
        int consistency_anchors;
        float consistency_weight;
        int realign;
        float vsm_amax;
        int mode;  /* 0=default, 1=fast, 2=precise */
        int help_flag;
        int quiet;
};

EXTERN struct parameters* init_param(void);
EXTERN void free_parameters(struct parameters *param);

EXTERN int check_msa_format_string(char* format);

#undef PARAMETERS_IMPORT
#undef EXTERN

#endif
