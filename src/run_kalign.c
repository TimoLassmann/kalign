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

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "global.h"

#include "tlmisc.h"

#include "msa.h"
#include "parameters.h"
//#include "align_io.h"
#include "alignment_parameters.h"

#include "bisectingKmeans.h"
/* #include "alignment.h" */
#include "weave_alignment.h"

#include "esl_stopwatch.h"
#include "misc.h"
#include <getopt.h>
#include <string.h>
#include "idata.h"
#include "alphabet.h"

#include "aln_run.h"
#include "aln_task.h"

#define OPT_SET 1
#define OPT_ALNPARAM 2
#define OPT_RENAME 3
#define OPT_REFORMAT 4
#define OPT_SHOWW 5
#define OPT_GPO 6
#define OPT_GPE 7
#define OPT_TGPE 8
#define OPT_MATADD 9

#define OPT_DEVTEST 10
#define OPT_CHAOS 11
#define OPT_DUMP_INTERNAL 12
#define OPT_NTHREADS 14

#define OPT_CLEAN 15
#define OPT_UNALIGN 16

static int run_kalign(struct parameters* param);
static int check_for_sequences(struct msa* msa);

static int print_kalign_header(void);
static int print_kalign_help(char * argv[]);
static int print_kalign_warranty(void);
static int print_AVX_warning(void);

int print_kalign_help(char * argv[])
{
        const char usage[] = " -i <seq file> -o <out aln> ";
        char* basename = NULL;


        RUN(tlfilename(argv[0], &basename));

        fprintf(stdout,"\nUsage: %s %s\n\n",basename ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--format","Output format." ,"[Fasta]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--reformat","Reformat existing alignment." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open penalty." ,"[5.5]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extension penalty." ,"[2.0]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Terminal gap extension penalty." ,"[1.0]"  );

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
        fprintf(stdout,"Kalign (%s)\n", PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2006,2019,2020,2021 Timo Lassmann\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY; for details type:\n");
        fprintf(stdout,"`kalign -showw'.\n");
        fprintf(stdout,"This is free software, and you are welcome to redistribute it\n");
        fprintf(stdout,"under certain conditions; consult the COPYING file for details.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Please cite:\n");

        /*        fprintf(stdout,"  Kalign 3: multiple sequence alignment of large data sets
Timo Lassmann
Bioinformatics, btz795, https://doi.org/10.1093/bioinformatics/btz795
        */
        fprintf(stdout,"  Lassmann, Timo.\n");
        fprintf(stdout,"  \"Kalign 3: multiple sequence alignment of large data sets.\"\n");
        fprintf(stdout,"  Bioinformatics (2019) \n");
        fprintf(stdout,"  https://doi.org/10.1093/bioinformatics/btz795\n");
        fprintf(stdout,"\n");

        /*fprintf(stdout,"  Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer.\n");
        fprintf(stdout,"  \"Kalign2: high-performance multiple alignment of protein and\n");
        fprintf(stdout,"  nucleotide sequences allowing external features.\"\n");
        fprintf(stdout,"  Nucleic acids research 37.3 (2008): 858-865.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"  Lassmann, Timo, and Erik LL Sonnhammer. \"Kalign–an accurate and\n");
        fprintf(stdout,"  fast multiple sequence alignment algorithm.\"\n  BMC bioinformatics 6.1 (2005): 298.\n");
        fprintf(stdout,"\n");*/

        return OK;
}

int print_AVX_warning(void)
{
        fprintf(stdout,"\n");
        fprintf(stdout,"WARNING: AVX2 instruction set not found!\n");
        fprintf(stdout,"         Kalign will not run optimally.\n");
        fprintf(stdout,"\n");
        return OK;
}


int main(int argc, char *argv[])
{
        int version = 0;
        int c;
        int showw = 0;
        int devtest = 0;
        struct parameters* param = NULL;
        char* in = NULL;
        RUNP(param = init_param());

        param->num_infiles = 0;

        while (1){
                static struct option long_options[] ={
                        {"showw", 0,0,OPT_SHOWW },
                        {"alnp", required_argument,0,OPT_ALNPARAM},
                        {"set", required_argument,0,OPT_SET},
                        {"format",  required_argument, 0, 'f'},
                        {"reformat",  required_argument, 0, OPT_REFORMAT},
                        {"devtest",  no_argument, 0, OPT_DEVTEST},
                        {"dumpinternal",no_argument, 0, OPT_DUMP_INTERNAL},
                        {"changename",  0, 0, OPT_RENAME},
                        {"clean",  0, 0, OPT_CLEAN},
                        {"unalign",  0, 0, OPT_UNALIGN},
                        {"gpo",  required_argument, 0, OPT_GPO},
                        {"gpe",  required_argument, 0, OPT_GPE},
                        {"tgpe",  required_argument, 0, OPT_TGPE},
                        {"matadd",  required_argument, 0, OPT_MATADD},
                        {"chaos",   required_argument, 0, OPT_CHAOS},
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
                case OPT_CLEAN:
                        param->clean = 1;
                        break;
                case OPT_UNALIGN:
                        param->unalign = 1;
                        break;
                case OPT_CHAOS:
                        param->chaos = atoi(optarg);
                        break;
                case OPT_DUMP_INTERNAL:
                        param->dump_internal = 1;
                        break;
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
                case 'n':
                        param->nthreads = atoi(optarg);
                        break;
                case OPT_REFORMAT:
                        param->format = optarg;
                        param->reformat = 1;
                        break;
                case OPT_DEVTEST:
                        devtest =1 ;
                        break;
                case OPT_GPO:
                        param->gpo = atof(optarg);
                        break;
                case OPT_GPE:
                        param->gpe = atof(optarg);
                        break;
                case OPT_TGPE :
                        param->tgpe = atof(optarg);
                        break;
                case OPT_MATADD :
                        param->matadd = atof(optarg);
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
                fprintf(stdout,"%s %s\n",PACKAGE_NAME, PACKAGE_VERSION);
                free_parameters(param);
                return EXIT_SUCCESS;
        }

        if(!param->dump_internal){
                print_kalign_header();
        }

#ifndef HAVE_AVX2
        RUN(print_AVX_warning());
#endif

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

        /* Just for internal development & testing */
        if(devtest){
                char* datadir = NULL;
                char* tmp = NULL;
                int f;
                f = 0;
                for(c = 0; c < param->num_infiles;c++){
                        param->infile[f] = param->infile[c];
                        if(param->infile[c] != NULL){
                                f++;
                        }
                }
                param->num_infiles = f;
                LOG_MSG("Start aln tests.");
                datadir = getenv("testdatafiledir");
                for(c = 0; c < param->num_infiles;c++){
                        tmp = NULL;
                        MMALLOC(tmp, sizeof(char) * 1024);
                        if(!datadir){
                                snprintf(tmp,1024, "%s",param->infile[c]);
                        }else{
                                snprintf(tmp,1024, "%s/%s", datadir, param->infile[c]);
                        }
                        param->infile[c] = tmp;
                }
        }

        if(!param->format){
                param->out_format = FORMAT_FA;
        }else{
                if(strstr(param->format,"msf")){
                        param->out_format = FORMAT_MSF;
                }else if(strstr(param->format,"clu")){
                        param->out_format = FORMAT_CLU;
                }else if(strstr(param->format,"fasta")){
                        param->out_format  = FORMAT_FA;
                }else if(strstr(param->format,"fa")){
                        param->out_format  = FORMAT_FA;
                }else{
                        ERROR_MSG("Format %s not recognized.",param->format);
                }
        }

        if (param->num_infiles == 0){
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

        if(param->chaos){
                if(param->chaos == 1){
                        ERROR_MSG("Param chaos need to be bigger than 1 (currently %d)", param->chaos);
                }
                if(param->chaos > 10){
                        ERROR_MSG("Param chaos bigger than 10 (currently %d)",param->chaos);
                }
        }

        RUN(run_kalign(param));

        if(devtest){
                for(c = 0; c < param->num_infiles;c++){
                        MFREE(param->infile[c]);
                }
        }

        free_parameters(param);
        return EXIT_SUCCESS;
ERROR:
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_kalign(struct parameters* param)
{
        struct msa* msa = NULL;
        struct msa* tmp_msa = NULL;
        struct aln_param* ap = NULL;
        struct aln_tasks* tasks = NULL;

        int i;

#ifdef HAVE_OPENMP
        omp_set_nested(1);
        omp_set_num_threads(param->nthreads);
#endif

        if(param->num_infiles == 1){
                RUN(read_input(param->infile[0],&msa));
        }else{
                for(i = 0; i < param->num_infiles;i++){
                        RUN(read_input(param->infile[i],&tmp_msa));
                        if(tmp_msa){
                                RUN(merge_msa(&msa, tmp_msa));
                                free_msa(tmp_msa);
                                tmp_msa = NULL;
                        }
                }
        }
        /* check if we have sequences  */
        RUN(check_for_sequences(msa));

        /* extra checks for input alignments */
        if(param->clean){
                RUN(run_extra_checks_on_msa(msa));
        }

        LOG_MSG("Detected: %d sequences.", msa->numseq);
        /* If we just want to reformat end here */
        if(param->reformat){
                //LOG_MSG("%s reformat",param->reformat);
                for (i = 0 ;i < msa->numseq;i++){
                        msa->nsip[i] = i;
                }

                if(param->rename){
                        for (i = 0 ;i < msa->numseq;i++){
                                snprintf(msa->sequences[i]->name, 128, "SEQ%d", i+1);
                        }
                }
                if(param->unalign){
                        /* if(param->out_format == FORMAT_FA){ */
                        RUN(dealign_msa(msa));
                }

                if(param->out_format != FORMAT_FA && msa->aligned != ALN_STATUS_ALIGNED){
                        ERROR_MSG("Input sequences are not aligned - cannot write to MSA format: %s", param->format);
                }

                RUN(write_msa(msa, param->outfile, param->out_format));

                free_msa(msa);
                return OK;
        }

        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }

        /* allocate aln parameters  */
        RUN(init_ap(&ap,param,msa->L ));

        if(param->dump_internal){
                double* s;
                int s_len;
                RUN(get_internal_data(msa,ap,&s, &s_len));
                for(i = 0; i < s_len;i++){
                        fprintf(stdout,"%0.2f ", s[i]);
                }
                fprintf(stdout,"\n");
                MFREE(s);
                free_msa(msa);
                return OK;
        }

        /* Allocate tasks  */
        RUN(alloc_tasks(&tasks, msa->numseq));

        /* Start bi-secting K-means sequence clustering */
        if(!param->chaos){
#ifdef HAVE_OPENMP
                i = floor(log((double) param->nthreads) / log(2.0)) + 4;
                i = MACRO_MIN(i, 10);
                /* LOG_MSG("Set %d level (%d)", i, param->nthreads); */
                omp_set_max_active_levels(i);
#endif
                RUN(build_tree_kmeans(msa,ap,&tasks));
        }
        /* by default all protein sequences are converted into a reduced alphabet
           when read from file. Here we turn them back into the default representation. */
        if(msa->L == ALPHA_redPROTEIN){
                //RUN(convert_msa_to_internal(msa, defPROTEIN));
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }
        /* allocate aln parameters  */
        RUN(init_ap(&ap,param,msa->L ));

        /* Start alignment stuff */
        DECLARE_TIMER(t1);
        LOG_MSG("Aligning");
        START_TIMER(t1);

        /* testing  */
#ifdef HAVE_OPENMP
        i = floor(log((double) param->nthreads) / log(2.0)) + 2;
        i = MACRO_MIN(i, 10);
        /* LOG_MSG("Set %d level (%d)", i, param->nthreads); */
        omp_set_max_active_levels(i);
#endif
        if(param->chaos){

                RUN(create_chaos_msa_openMP(msa, ap,tasks));

        }else{
                RUN(create_msa_tree(msa, ap, tasks,param->nthreads));
        }
        /* RUN(create_msa_openMP(msa,ap, tasks)); */

/*         if(param->chaos){ */
/* #ifdef HAVE_OPENMP */
/*                 RUN(create_chaos_msa_openMP(msa, ap,tasks)); */
/* #else */
/*                 RUN(create_chaos_msa_serial(msa, ap,tasks)); */
/* #endif */
/*         }else{ */
/* #ifdef HAVE_OPENMP */
/*                 RUN(create_msa_openMP(msa,ap, tasks)); */
/* #else */
/*                 RUN(create_msa_serial(msa,ap, tasks)); */
/* #endif */
/*         } */
        //RUNP(map = hirschberg_alignment(msa, ap));
        STOP_TIMER(t1);
        GET_TIMING(t1);

        /* set to aligned */
        msa->aligned = ALN_STATUS_ALIGNED;

        /* clean up map */
        free_tasks(tasks);

        /* We are done. */
        RUN(write_msa(msa, param->outfile, param->out_format));

        free_msa(msa);
        free_ap(ap);
        DESTROY_TIMER(t1);
        return OK;
ERROR:
        free_msa(tmp_msa);
        free_msa(msa);
        return FAIL;
}

int check_for_sequences(struct msa* msa)

{
        if(!msa){
                ERROR_MSG("No sequences were found in the input files or standard input.");
        }
        if(msa->numseq < 2){
                if(msa->numseq == 0){
                        ERROR_MSG("No sequences were found in the input files or standard input.");
                }else if (msa->numseq == 1){
                        ERROR_MSG("Only 1 sequence was found in the input files or standard input");
                }
        }
        return OK;
ERROR:
        return FAIL;
}
