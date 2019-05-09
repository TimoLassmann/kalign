#include "global.h"
#include "parameters.h"





struct parameters* init_param(void)
{
        struct parameters* param = NULL;


        MMALLOC(param, sizeof(struct parameters));
        param->dist_method = KALIGNDIST_BPM;
        param->gpo = -1.0;
        param->gpe = -1.0;
        param->tgpe = -1.0;
        param->secret = -1.0;
        param->zlevel = 58.8;
        param->sub_matrix = 0;
        param->aa = 0;


        param->infile = NULL;
        param->num_infiles = 0;
        param->input = NULL;
        param->outfile = NULL;
        param->format = NULL;
        param->help_flag = 0;
        param->quiet = 0;
        param->id = -1;
        param->distance = "wu";
        param->reformat = 0;
        param->sort = 0;

        param->print_svg_tree = 0;

        param->dna = -1;

        param->feature_type = NULL;
        param->alignment_type = "default";
        param->tree = "upgma";
        param->ntree = 2;
        param->print_tree = NULL;
        param->alter_gaps = 0;
        param->alter_range = 0.5;
        param->alter_weight = 100;

        param->internal_gap_weight = 0;
        param->smooth_window = 1;
        param->gap_inc = 0.0;
        param->same_feature_score = 75;
        param->diff_feature_score = -5;

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
