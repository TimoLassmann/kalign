#include "tldevel.h"
#include "tlmisc.h"
#include "kalign/kalign.h"
#include "msa_struct.h"
#include "aln_add.h"
#include "parameters.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <strings.h>

#define OPT_SHOWW 5
#define OPT_GPO 6
#define OPT_GPE 7
#define OPT_TGPE 8
#define OPT_ALN_TYPE 13
#define OPT_MODE 22
#define OPT_CONF_THRESHOLD 30
#define OPT_CONF_STYLE 31
#define OPT_CONF_OUTPUT 32
#define OPT_ADD 33
#define OPT_EXISTING 34
#define OPT_VSM_MODE 35
#define OPT_VSM_KAPPA 36
#define OPT_VSM_DZERO 37
#define OPT_FILL_IN 38
#define OPT_FILL_ANCHOR_FRAC 39
#define OPT_FILL_SUPPORT 40
#define OPT_GUIDE_DIST 41

static int set_aln_type(char* in, int* type);
static int apply_seq_type(struct msa* msa, struct parameters* param, const char* what);

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
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"fast","Single run, fastest (default). Alias: default");
        fprintf(stdout,"%*s%-*s: %s\n",3,"",MESSAGE_MARGIN-3,"accurate","5-member ensemble with per-column confidence; at most 5,000 sequences. Alias: recall");

        fprintf(stdout,"\nGeneral options (both modes):\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-i/--input FILE","Input file (or - for stdin)." ,"[]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-o/--output FILE","Output file." ,"[stdout]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-f/--format","Output format: fasta (fa, afa), msf or clu (clustal)." ,"[fasta]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--type","Sequence type: protein, dna or rna (overrides detection)." ,"[auto]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-n/--nthreads","Number of threads (at least 1; capped at the number of cores)." ,"[cores-1, max 16]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-q/--quiet","Do not print the banner or progress messages." ,"[off]");

        fprintf(stdout,"\nFast only (cannot be used with --mode accurate):\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open penalty (number >= 0)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extension penalty (number >= 0)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Terminal gap extension penalty (number >= 0)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--vsm-mode","Substitution offset: 0 from mean distance to anchor sequences, 1 from identity at each merge (bounded)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--vsm-kappa","Per-merge offset: fraction of the offset bound (number >= 0)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--vsm-dzero","Per-merge offset: distance where the offset is 0 (number >= 0)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--guide-dist","Guide-tree distance: 0 BPM, 1 pairwise alignment." ,"[fitted]");

        fprintf(stdout,"\nAccurate only (cannot be used with --mode fast):\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--fill-in 0-3","Place unaligned consensus residues (0 off, 1 realign, 2 project, 3 attach)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--fill-anchor-frac F","Fill-in anchor needs >= F*numseq residues (min 2; 0-1)." ,"[fitted]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--fill-support N","Attach fill-in: votes needed to join a column (integer >= 1)." ,"[fitted]");

        fprintf(stdout,"\nConfidence (per-column confidence exists only with --mode accurate;\nwith --mode fast these options have no effect):\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--confidence-threshold","Mask columns below this confidence (0-1)." ,"[off]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--confidence-style","Masking style: lowercase or remove." ,"[lowercase]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--confidence-output","Write per-column confidence to file." ,"[off]");

        fprintf(stdout,"\nAdding sequences (cannot be combined with --mode or the fast-only and accurate-only options):\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--add FILE","New sequences to add to existing alignment." ,"[off]");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--existing FILE","Existing alignment (new seqs added to this)." ,"[off]");

        fprintf(stdout,"\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-V/-v/--version","Prints version." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-h/--help","Prints this help." ,"[NA]"  );

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

/* Strict numeric option parsing: the whole argument must be a number. */
static int parse_int_arg(const char* s, int* out)
{
        char* end = NULL;
        long v;
        errno = 0;
        v = strtol(s, &end, 10);
        if(end == s || *end != '\0' || errno != 0 || v < INT_MIN || v > INT_MAX){
                return FAIL;
        }
        *out = (int) v;
        return OK;
}

static int parse_float_arg(const char* s, float* out)
{
        char* end = NULL;
        float v;
        errno = 0;
        v = strtof(s, &end);
        if(end == s || *end != '\0' || errno != 0){
                return FAIL;
        }
        *out = v;
        return OK;
}

/* Finite number >= 0. */
static int parse_nonneg_float_arg(const char* s, float* out)
{
        float v;
        if(parse_float_arg(s, &v) != OK || !isfinite(v) || v < 0.0f){
                return FAIL;
        }
        *out = v;
        return OK;
}

static int mode_is_accurate(const char* mode)
{
        return mode != NULL && (strcasecmp(mode, "accurate") == 0 || strcasecmp(mode, "recall") == 0);
}

static int mode_is_known(const char* mode)
{
        return mode == NULL || mode_is_accurate(mode) ||
                strcasecmp(mode, "fast") == 0 || strcasecmp(mode, "default") == 0;
}

int main(int argc, char *argv[])
{
        int version = 0;
        int c;
        int showw = 0;
        struct parameters* param = NULL;
        char* in = NULL;
        char* in_type = NULL;
        const char* fast_only_opt = NULL;      /* first fast-only option seen */
        const char* accurate_only_opt = NULL;  /* first accurate-only option seen */
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
                        {"vsm-mode",  required_argument, 0, OPT_VSM_MODE},
                        {"vsm-kappa",  required_argument, 0, OPT_VSM_KAPPA},
                        {"vsm-dzero",  required_argument, 0, OPT_VSM_DZERO},
                        {"guide-dist",  required_argument, 0, OPT_GUIDE_DIST},
                        {"mode",  required_argument, 0, OPT_MODE},
                        {"fill-in",  required_argument, 0, OPT_FILL_IN},
                        {"fill-anchor-frac",  required_argument, 0, OPT_FILL_ANCHOR_FRAC},
                        {"fill-support",  required_argument, 0, OPT_FILL_SUPPORT},
                        {"nthreads",  required_argument, 0, 'n'},
                        {"confidence-threshold", required_argument, 0, OPT_CONF_THRESHOLD},
                        {"confidence-style", required_argument, 0, OPT_CONF_STYLE},
                        {"confidence-output", required_argument, 0, OPT_CONF_OUTPUT},
                        {"add",  required_argument, 0, OPT_ADD},
                        {"existing",  required_argument, 0, OPT_EXISTING},
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
                        if(parse_int_arg(optarg, &param->nthreads) != OK || param->nthreads < 1){
                                ERROR_MSG("-n/--nthreads must be an integer >= 1, got '%s'", optarg);
                        }
                        break;
                case 'q':
                        param->quiet = 1;
                        break;
                case OPT_GPO:
                        if(parse_nonneg_float_arg(optarg, &param->gpo) != OK){
                                ERROR_MSG("--gpo must be a number >= 0, got '%s'", optarg);
                        }
                        if(!fast_only_opt) fast_only_opt = "--gpo";
                        break;
                case OPT_GPE:
                        if(parse_nonneg_float_arg(optarg, &param->gpe) != OK){
                                ERROR_MSG("--gpe must be a number >= 0, got '%s'", optarg);
                        }
                        if(!fast_only_opt) fast_only_opt = "--gpe";
                        break;
                case OPT_TGPE:
                        if(parse_nonneg_float_arg(optarg, &param->tgpe) != OK){
                                ERROR_MSG("--tgpe must be a number >= 0, got '%s'", optarg);
                        }
                        if(!fast_only_opt) fast_only_opt = "--tgpe";
                        break;
                case OPT_VSM_MODE:
                        if(parse_int_arg(optarg, &param->vsm_mode) != OK ||
                           (param->vsm_mode != 0 && param->vsm_mode != 1)){
                                ERROR_MSG("--vsm-mode must be 0 (distance to anchors) or 1 (per-merge identity), got '%s'", optarg);
                        }
                        if(!fast_only_opt) fast_only_opt = "--vsm-mode";
                        break;
                case OPT_VSM_KAPPA:
                        if(parse_nonneg_float_arg(optarg, &param->vsm_kappa) != OK){
                                ERROR_MSG("--vsm-kappa must be a number >= 0, got '%s'", optarg);
                        }
                        if(!fast_only_opt) fast_only_opt = "--vsm-kappa";
                        break;
                case OPT_VSM_DZERO:
                        if(parse_nonneg_float_arg(optarg, &param->vsm_dzero) != OK){
                                ERROR_MSG("--vsm-dzero must be a number >= 0, got '%s'", optarg);
                        }
                        if(!fast_only_opt) fast_only_opt = "--vsm-dzero";
                        break;
                case OPT_GUIDE_DIST:
                        if(parse_int_arg(optarg, &param->guide_dist) != OK ||
                           (param->guide_dist != KALIGN_GUIDE_DIST_BPM &&
                            param->guide_dist != KALIGN_GUIDE_DIST_PAIRALN)){
                                ERROR_MSG("--guide-dist must be 0 (BPM) or 1 (pairwise alignment), got '%s'", optarg);
                        }
                        if(!fast_only_opt) fast_only_opt = "--guide-dist";
                        break;
                case OPT_MODE:
                        param->mode = optarg;
                        break;
                case OPT_FILL_IN:
                        if(parse_int_arg(optarg, &param->fill_in) != OK ||
                           param->fill_in < KALIGN_FILL_OFF || param->fill_in > KALIGN_FILL_ATTACH){
                                ERROR_MSG("--fill-in must be 0 (off), 1 (realign), 2 (project) or 3 (attach), got '%s'", optarg);
                        }
                        if(!accurate_only_opt) accurate_only_opt = "--fill-in";
                        break;
                case OPT_FILL_ANCHOR_FRAC:
                        if(parse_float_arg(optarg, &param->fill_anchor_frac) != OK ||
                           !(param->fill_anchor_frac >= 0.0f && param->fill_anchor_frac <= 1.0f)){
                                ERROR_MSG("--fill-anchor-frac must be a number in [0, 1], got '%s'", optarg);
                        }
                        if(!accurate_only_opt) accurate_only_opt = "--fill-anchor-frac";
                        break;
                case OPT_FILL_SUPPORT:
                        if(parse_int_arg(optarg, &param->fill_support) != OK || param->fill_support < 1){
                                ERROR_MSG("--fill-support must be an integer >= 1, got '%s'", optarg);
                        }
                        if(!accurate_only_opt) accurate_only_opt = "--fill-support";
                        break;
                case OPT_CONF_THRESHOLD:
                        if(parse_float_arg(optarg, &param->confidence_threshold) != OK ||
                           !(param->confidence_threshold >= 0.0f && param->confidence_threshold <= 1.0f)){
                                ERROR_MSG("--confidence-threshold must be a number in [0, 1], got '%s'", optarg);
                        }
                        break;
                case OPT_CONF_STYLE:
                        if(strcmp(optarg, "remove") == 0){
                                param->confidence_style = KALIGN_MASK_REMOVE;
                        }else if(strcmp(optarg, "lowercase") == 0){
                                param->confidence_style = KALIGN_MASK_LOWERCASE;
                        }else{
                                ERROR_MSG("--confidence-style must be 'lowercase' or 'remove', got '%s'", optarg);
                        }
                        break;
                case OPT_CONF_OUTPUT:
                        param->confidence_output = optarg;
                        break;
                case OPT_ADD:
                        param->add_file = optarg;
                        break;
                case OPT_EXISTING:
                        param->existing_file = optarg;
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

        if(argc == 1){
                RUN(print_kalign_help(argv));
                free_parameters(param);
                return EXIT_SUCCESS;
        }

        if(!mode_is_known(param->mode)){
                ERROR_MSG("Unknown mode: '%s'. Use: fast or accurate.", param->mode);
        }

        if(param->add_file == NULL && param->existing_file == NULL){
                if(mode_is_accurate(param->mode) && fast_only_opt){
                        ERROR_MSG("%s cannot be used with --mode accurate: the accurate ensemble uses its fitted settings. Use --mode fast to set alignment parameters.", fast_only_opt);
                }
                if(!mode_is_accurate(param->mode) && accurate_only_opt){
                        ERROR_MSG("%s cannot be used with --mode fast: it only applies to the accurate ensemble. Use --mode accurate.", accurate_only_opt);
                }
        }

        param->num_infiles = 0;

        if(in){
                param->num_infiles++;
        }

        if (optind < argc){
                param->num_infiles += argc-optind;
        }

        /* --add mode doesn't need -i input files */
        if(param->add_file != NULL && param->existing_file != NULL){
                const char* add_conflict = NULL;
                if(param->mode) add_conflict = "--mode";
                else if(fast_only_opt) add_conflict = fast_only_opt;
                else if(accurate_only_opt) add_conflict = accurate_only_opt;
                if(add_conflict){
                        ERROR_MSG("%s cannot be used with --add: adding sequences uses its own settings.", add_conflict);
                }
                param->num_infiles = 0;  /* not using normal input */
        }else if(param->add_file != NULL || param->existing_file != NULL){
                LOG_MSG("Both --add and --existing must be specified together.");
                free_parameters(param);
                return EXIT_FAILURE;
        }else if(param->num_infiles == 0){
                ERROR_MSG("No input files.");
        }
        if(param->num_infiles > 0){
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

        /* --add mode: add new sequences to existing alignment */
        if(param->add_file != NULL && param->existing_file != NULL){
                struct msa* new_seqs = NULL;

                /* Read existing alignment */
                RUN(kalign_read_input(param->existing_file, &msa, param->quiet));

                /* Read new sequences (allows single sequence) */
                RUN(kalign_read_sequences(param->add_file, &new_seqs, param->quiet));
                RUN(apply_seq_type(msa, param, param->existing_file));
                if(new_seqs){
                        RUN(apply_seq_type(new_seqs, param, param->add_file));
                }

                if(!param->quiet){
                        LOG_MSG("Adding %d sequences to existing alignment of %d sequences",
                                new_seqs->numseq, msa->numseq);
                }

                /* Add new sequences to existing alignment */
                RUN(kalign_add_sequences(msa, new_seqs, param->nthreads));
                kalign_free_msa(new_seqs);

                /* Write result */
                RUN(kalign_write_msa(msa, param->outfile, param->format));
                kalign_free_msa(msa);
                return OK;
        }

        if(param->num_infiles == 1){
                RUN(kalign_read_input(param->infile[0], &msa, param->quiet));
        }else{
                for(int i = 0; i < param->num_infiles; i++){
                        RUN(kalign_read_input(param->infile[i], &msa, param->quiet));
                }
        }
        if(msa == NULL){
                ERROR_MSG("No sequences were found in the input files or standard input.");
        }
        RUN(apply_seq_type(msa, param, param->num_infiles == 1 && param->infile[0] ? param->infile[0] : "input"));

        /* Use mode preset (fast/default/recall/accurate).
           kalign_get_mode_preset auto-selects protein vs nucleotide
           based on the detected biotype. */
        const char* mode = param->mode ? param->mode : "fast";
        if(!mode_is_accurate(mode) &&
           (param->confidence_threshold > 0.0f || param->confidence_output != NULL) &&
           !param->quiet){
                WARNING_MSG("Per-column confidence needs --mode accurate; the confidence options have no effect with --mode fast.");
        }
        if(!kalign_mode_size_ok(mode, kalign_msa_get_biotype(msa), msa->numseq)){
                WARNING_MSG("Mode '%s' is not feasible above %d sequences (input: %d). Use --mode fast.",
                            mode, KALIGN_ACCURATE_MAX_SEQ, msa->numseq);
                goto ERROR;
        }
        int ret = kalign_get_mode_preset(mode,
                                          kalign_msa_get_biotype(msa),
                                          runs, &n_runs, &ens);
        if(ret != 0){
                ERROR_MSG("Unknown mode: '%s'. Use: fast or accurate.", mode);
        }

        /* Override preset values with explicit user parameters.
           Sentinel -1.0 means "use preset default". */
        for(int k = 0; k < n_runs; k++){
                if(param->gpo >= 0.0f)  runs[k].gpo  = param->gpo;
                if(param->gpe >= 0.0f)  runs[k].gpe  = param->gpe;
                if(param->tgpe >= 0.0f) runs[k].tgpe = param->tgpe;
                if(param->vsm_mode >= 0) runs[k].vsm_mode = param->vsm_mode;
                if(param->vsm_kappa >= 0.0f) runs[k].vsm_kappa = param->vsm_kappa;
                if(param->vsm_dzero >= 0.0f) runs[k].vsm_dzero = param->vsm_dzero;
                if(param->guide_dist >= 0) runs[k].guide_dist = param->guide_dist;
        }
        if(param->fill_in >= 0) ens.fill_in = param->fill_in;
        if(param->fill_anchor_frac >= 0.0f) ens.fill_anchor_frac = param->fill_anchor_frac;
        if(param->fill_support >= 1) ens.fill_support = param->fill_support;

        /* Without confidence options the alignment is written straight from
           the gap counts (kalign_write_msa), never held as numseq x alnlen */
        if(param->confidence_threshold <= 0.0f && param->confidence_output == NULL){
                msa->defer_finalise = 1;
        }
        RUN(kalign_align_full(msa, runs, n_runs, &ens, param->nthreads));

        /* Apply confidence masking if requested */
        if(param->confidence_threshold > 0.0f){
                RUN(kalign_mask_by_confidence(msa, param->confidence_threshold,
                                              param->confidence_style));
        }

        /* Write per-column confidence scores if requested */
        if(param->confidence_output != NULL){
                RUN(kalign_write_confidence(msa, param->confidence_output));
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
        if(kalign_parse_seq_type(in, type) != 0){
                ERROR_MSG("Sequence type '%s' not recognized. Use: auto, protein, dna, rna.",in);
        }
        return OK;
ERROR:
        return FAIL;
}

/* Apply --type to msa; warn (unless quiet) when it overrides the detected type. */
static int apply_seq_type(struct msa* msa, struct parameters* param, const char* what)
{
        int overridden = 0;
        const char* detected = kalign_seq_type_name(msa);
        if(kalign_set_seq_type(msa, param->type, &overridden) != 0){
                ERROR_MSG("Invalid sequence type.");
        }
        if(overridden && !param->quiet){
                WARNING_MSG("%s: sequences look like %s, aligning as %s as requested (--type).",
                            what, detected, kalign_seq_type_name(msa));
        }
        return OK;
ERROR:
        return FAIL;
}
