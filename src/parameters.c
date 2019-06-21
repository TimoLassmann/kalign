#include "global.h"
#include "parameters.h"





struct parameters*init_param(void)
{
        struct parameters* param = NULL;


        MMALLOC(param, sizeof(struct parameters));
        param->dist_method = KALIGNDIST_BPM;
        param->aln_param_file = NULL;
        param->param_set = -1;
        param->infile = NULL;
        param->num_infiles = 0;
        param->input = NULL;
        param->outfile = NULL;
        param->format = NULL;
        param->reformat = 0;
        param->rename = 1;
        param->help_flag = 0;
        param->quiet = 0;
        return param;
ERROR:
        free_parameters(param);
        return NULL;

}


void free_parameters(struct parameters* param)
{
        if(param){
                if(param->num_infiles){
                        MFREE(param->infile);
                }
                MFREE(param);
        }
}
