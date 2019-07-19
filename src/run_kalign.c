/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "global.h"
#include "msa.h"
#include "parameters.h"
#include "align_io.h"
#include "alignment_parameters.h"
// #include "estimate_aln_param.h"
// #include "tree_building.h"
#include "bisectingKmeans.h"
#include "alignment.h"
#include "weave_alignment.h"

//#include "counts_from_random_trees.h"
#include "misc.h"
#include <getopt.h>
#include "alphabet.h"

#define OPT_SET 1
#define OPT_ALNPARAM 2
#define OPT_RENAME 3
#define OPT_REFORMAT 4
#define OPT_SHOWW 5

int run_kalign(struct parameters* param);


int print_kalign_header(void);
int print_kalign_help(int argc, char * argv[]);
int print_kalign_warranty(void);

int print_kalign_help(int argc, char * argv[])
{
        const char usage[] = " -i <seq file> -o <out aln> ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--format","Output format." ,"[Fasta]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--reformat","Reformat existing alignment." ,"[NA]"  );
        fprintf(stdout,"\n");
        return OK;
}

int print_kalign_warranty(void)
{
        fprintf(stdout,"Here is the Disclaimer of Warranty section of the GNU General Public License (GPL):\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"15. Disclaimer of Warranty.\n");
        fprintf(stdout,"THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\n");
        fprintf(stdout,"APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n");
        fprintf(stdout,"HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\n");
        fprintf(stdout,"OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n");
        fprintf(stdout,"THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\n");
        fprintf(stdout,"PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n");
        fprintf(stdout,"IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\n");
        fprintf(stdout,"ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"A complete copy of the GPL can be found in the \"COPYING\" file.\n");
        return OK;
}

int print_kalign_header(void)
{
        fprintf(stdout,"\n");
        fprintf(stdout,"Kalign (%s)\n", PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2006,2019 Timo Lassmann\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY; for details type:\n");
        fprintf(stdout,"`kalign -showw'.\n");
        fprintf(stdout,"This is free software, and you are welcome to redistribute it\n");
        fprintf(stdout,"under certain conditions; consult the COPYING file for details.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Please cite:\n");


        fprintf(stdout,"  Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer.\n");
        fprintf(stdout,"  \"Kalign2: high-performance multiple alignment of protein and\n");
        fprintf(stdout,"  nucleotide sequences allowing external features.\"\n");
        fprintf(stdout,"  Nucleic acids research 37.3 (2008): 858-865.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"  Lassmann, Timo, and Erik LL Sonnhammer. \"Kalign–an accurate and\n");
        fprintf(stdout,"  fast multiple sequence alignment algorithm.\"\n  BMC bioinformatics 6.1 (2005): 298.\n");
        fprintf(stdout,"\n");

        return OK;
}



int main(int argc, char *argv[])
{
        int c;
        int showw = 0;
        struct parameters* param = NULL;

        RUNP(param = init_param());

        while (1){
                static struct option long_options[] ={
                        {"showw", 0,0,OPT_SHOWW },
                        {"alnp", required_argument,0,OPT_ALNPARAM},
                        {"set", required_argument,0,OPT_SET},
                        {"format",  required_argument, 0, 'f'},
                        {"reformat",  0, 0, OPT_REFORMAT},
                        {"rename",  0, 0, OPT_RENAME},
                        {"input",  required_argument, 0, 'i'},
                        {"infile",  required_argument, 0, 'i'},
                        {"in",  required_argument, 0, 'i'},
                        {"output",  required_argument, 0, 'o'},
                        {"outfile",  required_argument, 0, 'o'},
                        {"out",  required_argument, 0, 'o'},
                        {"help",   no_argument,0,'h'},
                        {"quiet",  0, 0, 'q'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;

                c = getopt_long_only (argc, argv,"i:o:f:hq",long_options, &option_index);
                //c = getopt (argc, argv, "hi:o:");
                /* Detect the end of the options. */

                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_SHOWW:
                        showw = 1;
                        break;
                case OPT_ALNPARAM:
                        param->aln_param_file = optarg;
                        break;
                case OPT_SET:
                        param->param_set = atoi(optarg);
                        break;
                case OPT_RENAME:
                        param->rename = 1;
                        break;
                case 'f':
                        param->format = optarg;
                        break;
                case OPT_REFORMAT:
                        param->reformat = 1;
                        break;
                case 'h':
                        param->help_flag = 1;
                        break;
                case 'i':
                        param->num_infiles =1;
                        MMALLOC(param->infile, sizeof(char*));
                        param->infile[0] = optarg;

                        break;
                case 'o':
                        param->outfile = optarg;
                        break;
                case '?':
                        free_parameters(param);
                        exit(1);
                        break;
                default:

                        abort ();
                }
        }

        print_kalign_header();

        if(showw){
                print_kalign_warranty();
                free_parameters(param);
                return EXIT_SUCCESS;

        }
        if(param->help_flag){
                RUN(print_kalign_help(argc, argv));
                free_parameters(param);
                return EXIT_SUCCESS;
        }

        if (optind < argc){
                c = 0;
                //fprintf(stderr,"EXTRA :%d\n",argc - optind);
                c = param->num_infiles;
                param->num_infiles += argc-optind;
                MREALLOC(param->infile, sizeof(char*) * param->num_infiles);
                while (optind < argc){
                        param->infile[c] =  argv[optind++];
                        c++;
                }
        }

        /*for (i = 0; i< param->num_infiles; i++){
                fprintf(stdout, "%d %s\n", i, param->infile[i]);

                }*/

        //exit(0);

        if(param->quiet){
                fclose(stderr);
        }

        //fprintf(stderr,"%s", license);
        if (param->help_flag){
               //fprintf(stderr,"%s\n", usage);
                exit(1);
        }
        //exit(0);

        if(!param->format){
                param->out_format = FORMAT_FA;
        }else{
                if(strstr(param->format,"msf")){
                        param->out_format = FORMAT_MSF;
                }else if(strstr(param->format,"clu")){
                        param->out_format = FORMAT_CLU;
                }else if(strstr(param->format,"fasta")){
                        param->out_format  =FORMAT_FA;
                }else if(strstr(param->format,"fa")){
                        param->out_format  =FORMAT_FA;
                }else{
                        ERROR_MSG("Format %s not recognized.",param->format);
                }
        }
        //LOG_MSG("SET: %d",param->param_set);

        if (param->num_infiles == 0){
                LOG_MSG("No infiles");
//fprintf(stderr,"%s\n", usage);
                return EXIT_SUCCESS;
        }


        log_command_line(argc, argv);
        RUN(run_kalign(param));
        free_parameters(param);
        return EXIT_SUCCESS;
ERROR:
        free_parameters(param);
        return EXIT_FAILURE;
}


int run_kalign(struct parameters* param)
{
        struct msa* msa = NULL;
        //struct alignment* aln = NULL;
        struct aln_param* ap = NULL;
        //float** dm;
        int** map = NULL;       /* holds all alignment paths  */
        int i;

        DECLARE_TIMER(t1);
        /* Step 1: read all input sequences & figure out output  */
        //LOG_MSG("Reading input.");

        for(i = 0; i < param->num_infiles;i++){
                RUNP(msa = read_input(param->infile[i],msa));
        }

        //RUNP(aln = detect_and_read_sequences(param));
        /* copy dna parameter to alignment */
        //aln->dna = param->dna;

        LOG_MSG("Detected: %d sequences.", msa->numseq);
        //LOG_MSG("Output is %s in format %s.", param->outfile,param->format);
        //LOG_MSG("Is DNA: %d", aln->dna);
        //param->dna = aln->dna;
        /* If we just want to reformat end here */
        if(param->reformat){
                //LOG_MSG("%s reformat",param->reformat);
                for (i = 0 ;i < msa->numseq;i++){
                        msa->nsip[i] = i;
                        // No idea what this was supposed to do */
                        /*for (j = 0; j < aln->sl[i];j++){
                                aln->s[i][j] = 0;
                                }*/
                }
                if(param->rename){
                        for (i = 0 ;i < msa->numseq;i++){

                                snprintf(msa->sequences[i]->name, 128, "SEQ%d", i+1);
                                //aln->lsn[i] = strnlen(aln->sn[i], 128);
                        }
                }
                //sparam->format = param->reformat;

                if (byg_start(param->format,"fastaFASTAfaFA") != -1){
                        RUN(dealign_msa(msa));
                }
                RUN(write_msa(msa, param->outfile, param->out_format));
                //RUN(output(aln, param));
                free_msa(msa);
                return OK;
        }
        if(msa->aligned){
                RUN(dealign_msa(msa));
        }

        //param->param_set = 4;
        /* allocate aln parameters  */

        RUNP(ap = init_ap(msa->numseq,msa->L ));

        //counts_from_random_trees(aln, ap, 10);

        //LOG_MSG("Building guide tree.");
        //START_TIMER(t1);
        //random_tree(ap, aln->numseq);
        //RUN(convert_alignment_to_internal(aln,defPROTEIN ));
        //param->dist_method = KALIGNDIST_WU;
        //RUN(build_tree(aln,param,ap));
        //RUN(convert_alignment_to_internal(aln,redPROTEIN));
        RUN(build_tree_kmeans(msa,ap));

        if(msa->L == redPROTEIN){
                RUN(convert_msa_to_internal(msa, defPROTEIN));
        }


        //STOP_TIMER(t1);
        //LOG_MSG("Took %f sec.", GET_TIMING(t1));

        //LOG_MSG("Building guide tree.");
        //START_TIMER(t1);
        //STOP_TIMER(t1);
        //LOG_MSG("Took %f sec.", GET_TIMING(t1));

        /* Start alignment stuff */
        LOG_MSG("Aligning");
        START_TIMER(t1);
        RUNP(map = hirschberg_alignment(msa, ap));
        STOP_TIMER(t1);
        LOG_MSG("Took %f sec.", GET_TIMING(t1));

        /* set to aligned */
        msa->aligned = 1;

        RUN(weave(msa , map, ap->tree));

        /* clean up map */
        for(i = 0; i < msa->num_profiles ;i++){
               if(map[i]){
                       MFREE(map[i]);
               }
        }
        MFREE(map);
        map = NULL;

        RUN(write_msa(msa, param->outfile, param->out_format));
        //RUN(output(aln, param));
        free_msa(msa);
        free_ap(ap);
        return OK;
ERROR:
        return FAIL;
}
