/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019, 2020, 2021 Timo Lassmann

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

#include "tldevel.h"
#include "tlmisc.h"
#include "kalign/kalign.h"
#include "version.h"
#include "parameters.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#define OPT_SET 1
#define OPT_SHOWW 5
#define OPT_GPO 6
#define OPT_GPE 7
#define OPT_TGPE 8

#define OPT_NTHREADS 10

#define OPT_ALN_TYPE 13

static int set_aln_type(char* in, int* type );

static int run_kalign(struct parameters* param);

static int print_kalign_header(void);
static int print_kalign_help(char * argv[]);
static int print_kalign_warranty(void);

int print_kalign_help(char * argv[])
{
        const char usage[] = " -i <seq file> -o <out aln> ";
        char* basename = NULL;


        RUN(tlfilename(argv[0], &basename));

        fprintf(stdout,"\nUsage: %s %s\n\n",basename ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--format","Output format." ,"[Fasta]"  );

        /* fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--reformat","Reformat existing alignment." ,"[NA]"  ); */
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--type","Alignment type (rna, dna, internal)." ,"[rna]"  );
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",MESSAGE_MARGIN-3,"","Options: protein, divergent (protein)" ,""  );
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",MESSAGE_MARGIN-3,"","         rna, dna, internal (nuc)." ,""  );

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open penalty." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extension penalty." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Terminal gap extension penalty." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-n/--nthreads","Number of threads." ,"[4]");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--version (-V/-v)","Prints version." ,"[NA]"  );

        fprintf(stdout,"\nExamples:\n\n");

        fprintf(stdout,"Passing sequences via stdin:\n\n   cat input.fa | kalign -f fasta > out.afa\n\n");
        fprintf(stdout,"Combining multiple input files:\n\n   kalign seqsA.fa seqsB.fa seqsC.fa -f fasta > combined.afa\n\n");

        if(basename){
                MFREE(basename);
        }
        return OK;
ERROR:
        if(basename){
                MFREE(basename);
        }
        return FAIL;
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
        fprintf(stdout,"Kalign (%s)\n", KALIGN_PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2006,2019,2020,2021,2023 Timo Lassmann\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY; for details type:\n");
        fprintf(stdout,"`kalign -showw'.\n");
        fprintf(stdout,"This is free software, and you are welcome to redistribute it\n");
        fprintf(stdout,"under certain conditions; consult the COPYING file for details.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Please cite:\n");

        fprintf(stdout,"  Lassmann, Timo.\n");
        fprintf(stdout,"  \"Kalign 3: multiple sequence alignment of large data sets.\"\n");
        fprintf(stdout,"  Bioinformatics (2019) \n");
        fprintf(stdout,"  https://doi.org/10.1093/bioinformatics/btz795\n");
        fprintf(stdout,"\n");
        return OK;
}

int main(int argc, char *argv[])
{
        int version = 0;
        int c;
        int showw = 0;
        struct parameters* param = NULL;
        char* in = NULL;
        char* in_type = NULL;
        RUNP(param = init_param());

        param->num_infiles = 0;

        while (1){
                static struct option long_options[] ={
                        {"showw", 0,0,OPT_SHOWW },
                        {"set", required_argument,0,OPT_SET},
                        {"format",  required_argument, 0, 'f'},
                        {"type",  required_argument, 0, OPT_ALN_TYPE},
                        {"gpo",  required_argument, 0, OPT_GPO},
                        {"gpe",  required_argument, 0, OPT_GPE},
                        {"tgpe",  required_argument, 0, OPT_TGPE},
                        {"nthreads",  required_argument, 0, 'n'},
                        {"input",  required_argument, 0, 'i'},
                        {"infile",  required_argument, 0, 'i'},
                        {"in",  required_argument, 0, 'i'},
                        {"output",  required_argument, 0, 'o'},
                        {"outfile",  required_argument, 0, 'o'},
                        {"out",  required_argument, 0, 'o'},
                        {"help",   no_argument,0,'h'},
                        {"version",   no_argument,0,'v'},
                        {"quiet",  0, 0, 'q'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;

                c = getopt_long_only (argc, argv,"i:o:f:n:hqvV",long_options, &option_index);

                /* Detect the end of the options. */
                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_ALN_TYPE:
                        in_type = optarg;
                        break;
                case OPT_SHOWW:
                        showw = 1;
                        break;
                case OPT_SET:
                        param->param_set = atoi(optarg);
                        break;
                case 'f':
                        param->format = optarg;
                        break;
                case 'n':
                        param->nthreads = atoi(optarg);
                        break;
                case 'q':
                        param->quiet = 1;
                        break;
                case OPT_GPO:
                        param->gpo = atof(optarg);
                        break;
                case OPT_GPE:
                        param->gpe = atof(optarg);
                        break;
                case OPT_TGPE:
                        param->tgpe = atof(optarg);
                        break;
                case 'h':
                        param->help_flag = 1;
                        break;
                case 'v':
                case 'V':
                        version = 1;
                        break;
                case 'i':
                        in = optarg;
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

        if(version){
                fprintf(stdout,"%s %s\n",KALIGN_PACKAGE_NAME, KALIGN_PACKAGE_VERSION);
                free_parameters(param);
                return EXIT_SUCCESS;
        }

        if(!param->dump_internal){
                if(!param->quiet){
                print_kalign_header();
                }
        }

        if(showw){
                print_kalign_warranty();
                free_parameters(param);
                return EXIT_SUCCESS;
        }

        if(param->help_flag){
                RUN(print_kalign_help(argv));
                free_parameters(param);
                return EXIT_SUCCESS;
        }

        if(param->nthreads < 1){
                RUN(print_kalign_help(argv));
                LOG_MSG("Number of threads has to be >= 1.");
                free_parameters(param);
                return EXIT_FAILURE;
        }

        param->num_infiles = 0;

        /* Not sure if this is the best way to detect stdin input */
        /* I fixed this by attempting to read and checking for valid
           sequences */

        if (!isatty(fileno(stdin))){
                param->num_infiles++;
        }

        if(in){
                param->num_infiles++;
        }

        if (optind < argc){
                param->num_infiles += argc-optind;
        }

        if(param->num_infiles == 0){
                RUN(print_kalign_help(argv));
                LOG_MSG("No input files");
                free_parameters(param);
                return EXIT_SUCCESS;
        }
        //fprintf(stdout,"%d fgiles\n",param->num_infiles);
        MMALLOC(param->infile, sizeof(char*) * param->num_infiles);

        c = 0;
        if (!isatty(fileno(stdin))){
                param->infile[c] = NULL;
                c++;
        }

        if(in){
                param->infile[c] = in;
                c++;
        }

        if (optind < argc){
                while (optind < argc){
                        param->infile[c] =  argv[optind++];
                        c++;
                }
        }

        RUN(check_msa_format_string(param->format));
        RUN(set_aln_type(in_type, &param->type));

        if(param->num_infiles == 0){
                if (!isatty(fileno(stdin))){
                        LOG_MSG("Attempting stdin");
                        param->num_infiles =1;
                        MMALLOC(param->infile, sizeof(char*));
                        param->infile[0] = NULL;
                }else{
                        LOG_MSG("No infiles");
                        free_parameters(param);
                        return EXIT_SUCCESS;
                }
        }

        /* if(param->chaos){ */
        /*         if(param->chaos == 1){ */
        /*                 ERROR_MSG("Param chaos need to be bigger than 1 (currently %d)", param->chaos); */
        /*         } */
        /*         if(param->chaos > 10){ */
        /*                 ERROR_MSG("Param chaos bigger than 10 (currently %d)",param->chaos); */
        /*         } */
        /* } */

        RUN(run_kalign(param));

        /* if(devtest){ */
        /*         for(c = 0; c < param->num_infiles;c++){ */
        /*                 MFREE(param->infile[c]); */
        /*         } */
        /* } */

        free_parameters(param);
        return EXIT_SUCCESS;
ERROR:
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_kalign(struct parameters* param)
{
        struct msa* msa = NULL;

        if(param->num_infiles == 1){
                RUN(kalign_read_input(param->infile[0], &msa,param->quiet));
        }else{
                for(int i = 0; i < param->num_infiles;i++){
                        RUN(kalign_read_input(param->infile[i], &msa,param->quiet));
                }
        }



        RUN(kalign_run(msa,
                       param->nthreads,
                       param->type,
                       param->gpo,
                       param->gpe,
                       param->tgpe));


        RUN(kalign_write_msa(msa, param->outfile, param->format));
        kalign_free_msa(msa);
        return OK;
ERROR:
        kalign_free_msa(msa);
        return FAIL;
}


int set_aln_type(char* in, int* type )
{
        int t = 0;
        if(in){
                if(strstr(in,"rna")){
                        t = KALIGN_TYPE_RNA;
                }else if(strstr(in,"dna")){
                        t = KALIGN_TYPE_DNA;
                }else if(strstr(in,"internal")){
                        t = KALIGN_TYPE_DNA_INTERNAL;
                }else if(strstr(in,"protein")){
                        t = KALIGN_TYPE_PROTEIN;
                }else if(strstr(in,"divergent")){
                        t = KALIGN_TYPE_PROTEIN_DIVERGENT;
                }else{
                        ERROR_MSG("In %s not recognized.",in);
                }
        }else{
                t = KALIGN_TYPE_UNDEFINED;
        }
        *type = t;
        return OK;
ERROR:
        return FAIL;
}

