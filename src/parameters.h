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
        int vsm_mode;      /* -1 = preset; 0 = legacy VSM; 1 = per-merge VSM */
        float vsm_kappa;   /* -1 = preset */
        float vsm_dzero;   /* -1 = preset */
        int guide_dist;    /* -1 = preset; 0 = BPM (shipped); 1 = pairwise-alignment distance */
        int num_infiles;
        int out_format;
        int nthreads;
        int fill_in;                /* -1 = preset, 0 = off, 1 = on */
        float fill_anchor_frac;     /* < 0 = preset */
        int fill_support;           /* < 1 = preset */
        char* mode;  /* "fast", "default", "recall", "accurate" (NULL = default) */
        char* add_file;             /* new sequences to add to existing alignment */
        char* existing_file;        /* existing alignment to add sequences to */
        float confidence_threshold;  /* mask columns below this confidence (0=off) */
        int confidence_style;       /* KALIGN_MASK_LOWERCASE or KALIGN_MASK_REMOVE */
        char* confidence_output;    /* write per-column confidence to file (NULL=off) */
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
