#ifndef PARAMETERS_H
#define PARAMETERS_H

static char  usage[] = "\n\
        Usage: kalign2   [INFILE] [OUTFILE] [OPTIONS]\n\
        \n\
	Options:\n\n\
        -s,	-gapopen          Gap open penalty\n\
        	-gap_open\n\
        	-gpo\n\
        	\n\
        -e,	-gapextension     Gap extension penalty\n\
        	-gap_ext\n\
        	-gpe\n\
        	\n\
        -t,	-terminal_gap_extension_penalty	Terminal gap penalties\n\
        	-tgpe\n\
        	\n\
        -m,	-matrix_bonus     A constant added to the substitution matrix.\n\
        	-bonus\n\
        	\n\
        -c,	-sort            The order in which the sequences appear in the output alignment.\n\
		                   <input, tree, gaps.>\n\
		\n\
        -g,	-feature          Selects feature mode and specifies which features are to be used:\n\
		                   e.g. all, maxplp, STRUCT, PFAM-A....\n\
           	-same_feature_score          Score for aligning same features\n\
		-diff_feature_score          Penalty for aligning different features\n\
        	\n\
        -d,	-distance         Distance method.\n\
		                   <wu,pair>\n\
		\n\
        -b,	-guide-tree       Guide tree method.\n\
		-tree             <nj,upgma>\n\
		\n\
	-z,	-zcutoff         Parameter used in the wu-manber based distance calculation\n\
		\n\
        -i,	-input            The input file.\n\
        	-infile\n\
        	-in\n\
        	\n\
        -o,	-output           The output file.\n\
        	-outfile\n\
        	-out\n\
        	\n\
        -a,	-gap_inc           Parameter increases gap penalties depending on the number of existing gaps\n\
        	\n\
        -f,	-format           The output format:\n\
		                   <fasta, msf, aln, clu, macsim>\n\
		\n\
	-q,	-quiet            Print nothing to STDERR.\n\
		                  Read nothing from STDIN\n\
	\n\
	Examples:\n\n\
	Using pipes:\n\
		kalign2 [OPTIONS] < [INFILE]   > [OUTFILE]\n\
		more [INFILE] |  kalign2 [OPTIONS] > [OUTFILE]\n\
         \n\
	Relaxed gap penalties:\n\
		kalign2 -gpo 60 -gpe 9 -tgpe 0 -bonus 0 < [INFILE]   > [OUTFILE]\n\
         \n\
        Feature alignment with pairwise alignment based distance method and NJ guide tree:\n\
        	kalign2 -in test.xml -distance pair -tree nj -sort gaps -feature STRUCT -format macsim -out test.macsim\n\
        ";



struct parameters{
        char **infile;
        char *input;
        char *outfile;
        char* format;
        //int reformat;
        char* feature_type;
        char* alignment_type;
        char* feature_mode;
        char* distance;
        char* tree;
        char* sort;
        char* sub_matrix;
        char* print_tree;
        char* print_svg_tree;
        float gpo;
        float gpe;
        float tgpe;
        float secret;
        float zlevel;
        float same_feature_score;
        float diff_feature_score;

        int num_infiles;
        int reformat;
        int id;
        int aa;
        int alter_gaps;
        int ntree;
        int help_flag;
        int quiet;

        int dna;
        float alter_range;
        int alter_weight;
        float internal_gap_weight;
        int smooth_window;
        float gap_inc;
};

extern struct parameters* init_param(void);
extern void free_parameters(struct parameters* param);
#endif
