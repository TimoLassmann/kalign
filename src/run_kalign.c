#include "tldevel.h"
#include "tlmisc.h"
#include "kalign/kalign.h"
#include "parameters.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#define OPT_SHOWW 5
#define OPT_GPO 6
#define OPT_GPE 7
#define OPT_TGPE 8
#define OPT_ALN_TYPE 13
#define OPT_MODE 22
#define OPT_LOAD_POAR 20
#define OPT_MIN_SUPPORT 18

static int set_aln_type(char* in, int* type);

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

        fprintf(stdout,"Modes (--mode <name>):\n\n");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"fast","Single run, fastest");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"default","Single run with consistency anchors (default)");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"recall","Ensemble, optimized for recall");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"accurate","Ensemble, highest precision");

        fprintf(stdout,"\nOptions:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--format","Output format." ,"[Fasta]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--type","Sequence type." ,"[auto]"  );
        fprintf(stdout,"%*s%-*s  %s %s\n",3,"",MESSAGE_MARGIN-3,"","Options: protein, dna, rna." ,""  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open penalty (overrides preset)." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extension penalty (overrides preset)." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Terminal gap extension penalty (overrides preset)." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-n/--nthreads","Number of threads." ,"[auto]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--load-poar","Load POAR table for re-threshold." ,"[off]");

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
        RUNP(param = init_param());

        param->num_infiles = 0;

        while (1){
                static struct option long_options[] ={
                        {"showw", 0,0,OPT_SHOWW },
                        {"format",  required_argument, 0, 'f'},
                        {"type",  required_argument, 0, OPT_ALN_TYPE},
                        {"gpo",  required_argument, 0, OPT_GPO},
                        {"gpe",  required_argument, 0, OPT_GPE},
                        {"tgpe",  required_argument, 0, OPT_TGPE},
                        {"mode",  required_argument, 0, OPT_MODE},
                        {"load-poar",  required_argument, 0, OPT_LOAD_POAR},
                        {"min-support",  required_argument, 0, OPT_MIN_SUPPORT},
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
                case OPT_MODE:
                        param->mode = optarg;
                        break;
                case OPT_LOAD_POAR:
                        param->load_poar = optarg;
                        break;
                case OPT_MIN_SUPPORT:
                        param->min_support = atoi(optarg);
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
        struct kalign_run_config runs[KALIGN_MAX_PRESET_RUNS];
        struct kalign_ensemble_config ens = kalign_ensemble_config_defaults();
        int n_runs = 0;

        if(param->num_infiles == 1){
                RUN(kalign_read_input(param->infile[0], &msa, param->quiet));
        }else{
                for(int i = 0; i < param->num_infiles; i++){
                        RUN(kalign_read_input(param->infile[i], &msa, param->quiet));
                }
        }

        if(param->load_poar != NULL){
                RUN(kalign_consensus_from_poar(msa,
                                               param->load_poar,
                                               param->min_support > 0 ? param->min_support : 2));
        }else{
                /* Use mode preset (fast/default/recall/accurate).
                   kalign_get_mode_preset auto-selects protein vs nucleotide
                   based on the detected biotype. */
                const char* mode = param->mode ? param->mode : "default";
                int ret = kalign_get_mode_preset(mode,
                                                  kalign_msa_get_biotype(msa),
                                                  runs, &n_runs, &ens);
                if(ret != 0){
                        ERROR_MSG("Unknown mode: '%s'. Use: fast, default, recall, accurate.", mode);
                }

                /* Override preset values with explicit user parameters.
                   Sentinel -1.0 means "use preset default". */
                for(int k = 0; k < n_runs; k++){
                        if(param->gpo >= 0.0f)  runs[k].gpo  = param->gpo;
                        if(param->gpe >= 0.0f)  runs[k].gpe  = param->gpe;
                        if(param->tgpe >= 0.0f) runs[k].tgpe = param->tgpe;
                }

                RUN(kalign_align_full(msa, runs, n_runs, &ens, param->nthreads));
        }

        RUN(kalign_write_msa(msa, param->outfile, param->format));
        kalign_free_msa(msa);
        return OK;
ERROR:
        kalign_free_msa(msa);
        return FAIL;
}

int set_aln_type(char* in, int* type)
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
                }else{
                        ERROR_MSG("Sequence type '%s' not recognized. Use: protein, dna, rna.",in);
                }
        }else{
                t = KALIGN_TYPE_UNDEFINED;
        }
        *type = t;
        return OK;
ERROR:
        return FAIL;
}
