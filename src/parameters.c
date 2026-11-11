#include "tldevel.h"
#include "kalign/kalign.h"

#include <string.h>
#include <strings.h>
#include <unistd.h>
#define PARAMETERS_IMPORT
#include "parameters.h"

static int get_default_thread_count(void)
{
        int cores = 1;

#if defined(_SC_NPROCESSORS_ONLN)
        cores = sysconf(_SC_NPROCESSORS_ONLN);
        if (cores <= 0) cores = 1;
#endif

        if (cores > 1) cores = cores - 1;
        if (cores > 16) cores = 16;
        if (cores < 1) cores = 1;

        return cores;
}

struct parameters*init_param(void)
{
        struct parameters* param = NULL;
        MMALLOC(param, sizeof(struct parameters));
        param->infile = NULL;
        param->num_infiles = 0;
        param->input = NULL;
        param->outfile = NULL;
        param->format = NULL;
        param->help_flag = 0;
        param->dump_internal = 0;
        param->type = -1;
        param->gpo = -1.0;
        param->gpe = -1.0;
        param->tgpe = -1.0;
        param->vsm_mode = -1;
        param->vsm_kappa = -1.0f;
        param->vsm_dzero = -1.0f;
        param->guide_dist = -1;
        param->nthreads = get_default_thread_count();
        param->fill_in = -1;
        param->fill_anchor_frac = -1.0f;
        param->fill_support = 0;
        param->mode = NULL;
        param->quiet = 0;
        param->out_format = 0;
        param->reformat = 0;
        param->rename = 0;
        param->clean = 0;
        param->unalign = 0;
        param->add_file = NULL;
        param->existing_file = NULL;
        param->confidence_threshold = 0.0f;
        param->confidence_style = 0;  /* KALIGN_MASK_LOWERCASE */
        param->confidence_output = NULL;
        return param;
ERROR:
        free_parameters(param);
        return NULL;
}

void free_parameters(struct parameters* param)
{
        if(param){
                if(param->infile){
                        MFREE(param->infile);
                }
                MFREE(param);
        }
}

int check_msa_format_string(char* format)
{
        if(format){
                if(strcasecmp(format,"fasta") != 0 &&
                   strcasecmp(format,"fa") != 0 &&
                   strcasecmp(format,"afa") != 0 &&
                   strcasecmp(format,"msf") != 0 &&
                   strcasecmp(format,"clu") != 0 &&
                   strcasecmp(format,"clustal") != 0){
                        ERROR_MSG("Format %s not recognized.",format);
                }
        }
        return OK;
ERROR:
        return FAIL;
}
