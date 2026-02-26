#include "tldevel.h"
#include "tlmisc.h"
#include "kalign/kalign.h"
/* #include "version.h" */
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
#define OPT_REFINE 14
#define OPT_ADAPTIVE_BUDGET 15
#define OPT_ENSEMBLE 16
#define OPT_ENSEMBLE_SEED 17
#define OPT_MIN_SUPPORT 18
#define OPT_SAVE_POAR 19
#define OPT_LOAD_POAR 20
#define OPT_CONSISTENCY 21
#define OPT_FAST 22
#define OPT_PRECISE 23
#define OPT_REALIGN 24
#define OPT_VSM_AMAX 25
#define OPT_CONSISTENCY_WEIGHT 26

static int set_aln_type(char* in, int* type );
static int set_refine_mode(char* in, int* refine);

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

        fprintf(stdout,"Modes:\n\n");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"(default)","Consistency anchors + VSM (best general-purpose)");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"--fast","VSM only, no consistency (fastest)");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"--precise","Ensemble(3) + VSM + realign (highest precision)");

        fprintf(stdout,"\nOptions:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--format","Output format." ,"[Fasta]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--type","Alignment type (rna, dna, internal)." ,"[rna]"  );
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",MESSAGE_MARGIN-3,"","Options: protein, divergent (protein)" ,""  );
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",MESSAGE_MARGIN-3,"","         rna, dna, internal (nuc)." ,""  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open penalty." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extension penalty." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Terminal gap extension penalty." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--refine","Refinement mode." ,"[none]");
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",MESSAGE_MARGIN-3,"","Options: none, all, confident" ,""  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-n/--nthreads","Number of threads." ,"[auto: N-1, max 16]");

        fprintf(stdout,"\nEnsemble options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--ensemble","Number of ensemble runs." ,"[off; 5 if no value given]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--ensemble-seed","RNG seed for ensemble." ,"[42]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--min-support","Explicit consensus threshold." ,"[auto]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--save-poar","Save POAR table to file." ,"[off]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--load-poar","Load POAR table for re-threshold." ,"[off]");

        fprintf(stdout,"\nAdvanced (usually managed by modes):\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--consistency","Anchor consistency (K anchors)." ,"[5]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--consistency-weight","Consistency anchor weight." ,"[2.0]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--realign","Alignment-guided tree rebuild iters." ,"[0]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--vsm-amax","VSM amplitude (0 to disable)." ,"[auto]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--adaptive-budget","Scale refinement trials by uncertainty." ,"[off]");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--version (-V/-v)","Prints version." ,"[NA]"  );

        fprintf(stdout,"\nExamples:\n\n");

        fprintf(stdout,"Passing sequences via stdin:\n\n   cat input.fa | kalign -i - -f fasta > out.afa\n\n");
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
        fprintf(stdout,"Disclaimer of Warranty (Apache License, Version 2.0, Section 7):\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Unless required by applicable law or agreed to in writing, Licensor\n");
        fprintf(stdout,"provides the Work (and each Contributor provides its Contributions)\n");
        fprintf(stdout,"on an \"AS IS\" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,\n");
        fprintf(stdout,"either express or implied, including, without limitation, any\n");
        fprintf(stdout,"warranties or conditions of TITLE, NON-INFRINGEMENT,\n");
        fprintf(stdout,"MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"See the COPYING file for the full Apache License, Version 2.0.\n");
        return OK;
}

int print_kalign_header(void)
{
        fprintf(stderr,"\n");
        fprintf(stderr,"Kalign (%s)\n", KALIGN_PACKAGE_VERSION);
        fprintf(stderr,"\n");
        fprintf(stderr,"Copyright (C) 2006-2026 Timo Lassmann\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"Licensed under the Apache License, Version 2.0.\n");
        fprintf(stderr,"See the COPYING file or http://www.apache.org/licenses/LICENSE-2.0\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"Please cite:\n");

        fprintf(stderr,"  Lassmann, Timo.\n");
        fprintf(stderr,"  \"Kalign 3: multiple sequence alignment of large data sets.\"\n");
        fprintf(stderr,"  Bioinformatics (2019) \n");
        fprintf(stderr,"  https://doi.org/10.1093/bioinformatics/btz795\n");
        fprintf(stderr,"\n");
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
        char* in_refine = NULL;
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
                        {"refine",  required_argument, 0, OPT_REFINE},
                        {"adaptive-budget",  no_argument, 0, OPT_ADAPTIVE_BUDGET},
                        {"ensemble",  optional_argument, 0, OPT_ENSEMBLE},
                        {"ensemble-seed",  required_argument, 0, OPT_ENSEMBLE_SEED},
                        {"min-support",  required_argument, 0, OPT_MIN_SUPPORT},
                        {"save-poar",  required_argument, 0, OPT_SAVE_POAR},
                        {"load-poar",  required_argument, 0, OPT_LOAD_POAR},
                        {"consistency",  required_argument, 0, OPT_CONSISTENCY},
                        {"consistency-weight",  required_argument, 0, OPT_CONSISTENCY_WEIGHT},
                        {"fast",  no_argument, 0, OPT_FAST},
                        {"precise",  no_argument, 0, OPT_PRECISE},
                        {"realign",  required_argument, 0, OPT_REALIGN},
                        {"vsm-amax",  required_argument, 0, OPT_VSM_AMAX},
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
                case OPT_REFINE:
                        in_refine = optarg;
                        break;
                case OPT_ADAPTIVE_BUDGET:
                        param->adaptive_budget = 1;
                        break;
                case OPT_ENSEMBLE:
                        if(optarg){
                                param->ensemble = atoi(optarg);
                        }else if(optind < argc && argv[optind][0] != '-'){
                                param->ensemble = atoi(argv[optind]);
                                optind++;
                        }else{
                                param->ensemble = 5;
                        }
                        break;
                case OPT_ENSEMBLE_SEED:
                        param->ensemble_seed = (uint64_t)strtoull(optarg, NULL, 10);
                        break;
                case OPT_MIN_SUPPORT:
                        param->min_support = atoi(optarg);
                        break;
                case OPT_SAVE_POAR:
                        param->save_poar = optarg;
                        break;
                case OPT_LOAD_POAR:
                        param->load_poar = optarg;
                        break;
                case OPT_CONSISTENCY:
                        param->consistency_anchors = atoi(optarg);
                        break;
                case OPT_CONSISTENCY_WEIGHT:
                        param->consistency_weight = atof(optarg);
                        break;
                case OPT_FAST:
                        param->mode = 1;
                        break;
                case OPT_PRECISE:
                        param->mode = 2;
                        break;
                case OPT_REALIGN:
                        param->realign = atoi(optarg);
                        break;
                case OPT_VSM_AMAX:
                        param->vsm_amax = atof(optarg);
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

        /* Use "-" to indicate stdin (like samtools/bcftools).
         * NULL in the infile array signals read_file_stdin() to use stdin. */
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
        MMALLOC(param->infile, sizeof(char*) * param->num_infiles);

        c = 0;
        if(in){
                param->infile[c] = (strcmp(in, "-") == 0) ? NULL : in;
                c++;
        }

        if (optind < argc){
                while (optind < argc){
                        if(strcmp(argv[optind], "-") == 0){
                                param->infile[c] = NULL; /* stdin */
                        }else{
                                param->infile[c] = argv[optind];
                        }
                        c++;
                        optind++;
                }
        }

        RUN(check_msa_format_string(param->format));
        RUN(set_aln_type(in_type, &param->type));
        RUN(set_refine_mode(in_refine, &param->refine));

        /* Apply mode presets. Explicit params (already set above) will have
         * overridden their init_param() defaults, so we only fill in mode
         * values for fields that are still at their init_param() defaults. */
        if(param->mode == 1){
                /* fast: no consistency anchors (unless user explicitly set --consistency) */
                if(param->consistency_anchors == 5){  /* still at default */
                        param->consistency_anchors = 0;
                }
        }else if(param->mode == 2){
                /* precise: ensemble + realign */
                if(param->ensemble == 0){
                        param->ensemble = 3;
                }
                if(param->realign == 0){
                        param->realign = 1;
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



        if(param->load_poar != NULL){
                RUN(kalign_consensus_from_poar(msa,
                                               param->load_poar,
                                               param->min_support > 0 ? param->min_support : 2));
        }else if(param->ensemble > 0){
                RUN(kalign_ensemble(msa,
                                    param->nthreads,
                                    param->type,
                                    param->ensemble,
                                    param->gpo,
                                    param->gpe,
                                    param->tgpe,
                                    param->ensemble_seed,
                                    param->min_support,
                                    param->save_poar,
                                    param->refine, 0.0f, param->vsm_amax,
                                    param->realign, -1.0f,
                                    param->consistency_anchors, param->consistency_weight));
        }else if(param->realign > 0){
                RUN(kalign_run_realign(msa,
                                       param->nthreads,
                                       param->type,
                                       param->gpo,
                                       param->gpe,
                                       param->tgpe,
                                       param->refine,
                                       param->adaptive_budget,
                                       0.0f, param->vsm_amax,
                                       param->realign, -1.0f,
                                       param->consistency_anchors, param->consistency_weight));
        }else{
                RUN(kalign_run_seeded(msa,
                                      param->nthreads,
                                      param->type,
                                      param->gpo,
                                      param->gpe,
                                      param->tgpe,
                                      param->refine,
                                      param->adaptive_budget,
                                      0, 0.0f, 0.0f, param->vsm_amax, -1.0f,
                                      param->consistency_anchors, param->consistency_weight));
        }


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
                }else if(strstr(in,"pfasum43")){
                        t = KALIGN_TYPE_PROTEIN_PFASUM43;
                }else if(strstr(in,"pfasum60")){
                        t = KALIGN_TYPE_PROTEIN_PFASUM60;
                }else if(strstr(in,"pfasum")){
                        t = KALIGN_TYPE_PROTEIN_PFASUM_AUTO;
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

int set_refine_mode(char* in, int* refine)
{
        if(in){
                if(strstr(in,"all")){
                        *refine = KALIGN_REFINE_ALL;
                }else if(strstr(in,"confident")){
                        *refine = KALIGN_REFINE_CONFIDENT;
                }else if(strstr(in,"none")){
                        *refine = KALIGN_REFINE_NONE;
                }else{
                        ERROR_MSG("Refine mode '%s' not recognized. Use: none, all, confident.", in);
                }
        }
        /* When in is NULL, keep the default from init_param() */
        return OK;
ERROR:
        return FAIL;
}
