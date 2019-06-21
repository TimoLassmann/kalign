#ifndef PARAMETERS_H
#define PARAMETERS_H


#define KALIGNDIST_ALN 0
#define KALIGNDIST_BPM 1
#define KALIGNDIST_WU 2

struct parameters{
        char **infile;
        char *input;
        char *outfile;
        char* format;
        char* aln_param_file;
        //int reformat;
        /* char* feature_type; */
        /* char* alignment_type; */
        /* char* feature_mode; */
        /* char* distance; */
        /* char* tree; */
        /* char* sort; */
        /* char* sub_matrix; */
        /* char* print_tree; */
        /* char* print_svg_tree; */
        /* float gpo; */
        /* float gpe; */
        /* float tgpe; */
        /* float secret; */
        /* float zlevel; */
        /* float same_feature_score; */
        /* float diff_feature_score; */

        int param_set;
        int dist_method;
        int num_infiles;
        int reformat;
        int rename;             /* rename sequences - to make bali_score swallow the alignments */
        /* int id; */
        /* int aa; */
        /* int alter_gaps; */
        /* int ntree; */
        int help_flag;
        int quiet;

        /* int dna; */
        /* float alter_range; */
        /* int alter_weight; */
        /* float internal_gap_weight; */
        /* int smooth_window; */
        /* float gap_inc; */
};

extern struct parameters* init_param(void);
extern void free_parameters(struct parameters* param);
#endif
