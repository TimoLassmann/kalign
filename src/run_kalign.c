#include "global.h"
#include "parameters.h"
#include "align_io.h"
#include "alignment_parameters.h"
#include "estimate_aln_param.h"
#include "tree_building.h"
#include "bisectingKmeans.h"
#include "alignment.h"
#include "weave_alignment.h"
#include "misc.h"
#include <getopt.h>


int run_kalign(struct parameters* param);
int detect_dna(struct alignment* aln);


int main(int argc, char *argv[])
{
        int c;
        struct parameters* param = NULL;

        RUNP(param = init_param());

        while (1){
                static struct option long_options[] ={
                        {"gapopen",  required_argument, 0,'s'},
                        {"gpo",  required_argument, 0, 's'},
                        {"gapextension",  required_argument, 0, 'e'},
                        {"gpe",  required_argument, 0, 'e'},
                        {"secret",  required_argument, 0, 'm'},
                        {"bonus",  required_argument, 0, 'm'},
                        {"terminalgapextension",  required_argument, 0, 't'},
                        {"tgpe",  required_argument, 0, 't'},
                        {"zcutoff",  required_argument, 0, 0},
                        {"distance",  required_argument, 0, 'd'},
                        {"ntree",  required_argument, 0, 0},
                        {"tree",  required_argument, 0, 0},
                        {"format",  required_argument, 0, 'f'},
                        {"reformat",  0, 0, 'r'},
                        {"sort",required_argument,0,'c'},
                        {"feature",  required_argument, 0, 0},
                        {"type",  required_argument, 0, 0},
                        {"alter_gaps",  required_argument, 0, 0},
                        {"altergaps",  required_argument, 0, 0},
                        {"alter_range",  required_argument, 0, 0},
                        {"alter_weight",  required_argument, 0, 0},
                        {"internal_gap_weight",  required_argument, 0, 0},
                        {"smooth_window",  required_argument, 0, 0},
                        {"gap_inc",  required_argument, 0, 'a'},
                        {"matrix",  required_argument, 0, 0},
                        {"mmbonus",  required_argument, 0, 0},
                        {"nuc",  0, 0, 0},
                        {"dna",  0, 0, 0},
                        {"rna",  0, 0, 0},
                        {"protein",  0, 0, 0},
                        {"profile", 0, 0, 0},
                        {"prof", 0, 0, 0},
                        {"id", required_argument, 0, 0},
                        {"printtree", required_argument, 0, 0},
                        {"svgtree", required_argument, 0, 0},
                        {"svg_tree", required_argument, 0, 0},
                        {"pairwise", 0, 0, 0},
                        {"same_feature_score", required_argument, 0, 0},
                        {"diff_feature_score", required_argument, 0, 0},

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
                c = getopt_long_only (argc, argv,"i:o:hqs:e:m:t:z:c:f:d:b:a:r",long_options, &option_index);
                //c = getopt (argc, argv, "hi:o:");
                /* Detect the end of the options. */
                if (c == -1){
                        break;
                }
                switch (c){
                case 0:
                        if (long_options[option_index].flag != 0){
                                break;
                        }
                        switch (option_index){
                        case 0:
                        case 1:
                                fprintf(stderr,"GAGAGA\n");
                                param->gpo = atof(optarg);
                                break;
                        case 2:
                        case 3:
                                param->gpe = atof(optarg);
                                break;
                        case 4:
                        case 5:
                                param->secret = atof(optarg);
                                break;
                        case 6:
                        case 7:
                                param->tgpe = atof(optarg);
                                break;
                        case 8:
                                param->zlevel = atof(optarg);
                                break;
                        case 9:
                                param->distance = optarg;
                                break;
                        case 10:
                                param->ntree = atoi(optarg);
                                break;
                        case 11:
                                param->tree = optarg;
                                break;
                        case 12:
                                param->format = optarg;
                                break;
                        case 13:
                                param->reformat = 1;
                                break;
                        case 14:
                                param->sort = optarg;
                                break;
                        case 15:
                                param->feature_type = optarg;
                                break;
                        case 16:
                                param->alignment_type = optarg;
                                break;
                        case 17:
                        case 18:
                                param->alter_gaps = atoi(optarg);
                                break;
                        case 19:
                                param->alter_range = atof(optarg);
                                break;
                        case 20:
                                param->alter_weight = atoi(optarg);
                                break;
                        case 21:
                                param->internal_gap_weight = atof(optarg);
                                break;
                        case 22:
                                param->smooth_window = atoi(optarg);
                                break;
                        case 23:
                                param->gap_inc = atof(optarg);
                                break;
                        case 24:
                                param->sub_matrix = optarg;
                                break;
                        case 25:
                                param->aa = atoi(optarg);
                                break;
                        case 26:
                        case 27:
                        case 28:
                                param->dna = 1;
                                break;
                        case 29:
                                param->dna = 0;
                                break;
                        case 30:
                        case 31:
                                param->alignment_type = "profile";
                                break;
                        case 32:
                                param->id = atoi(optarg);
                                break;
                        case 33:
                                param->print_tree = optarg;
                                break;
                        case 34:
                        case 35:
                                param->print_svg_tree = optarg;
                                break;
                        case 36:
                                param->alignment_type = "pairwise";
                                break;
                        case 37:
                                param->same_feature_score = atof(optarg);//"pairwise";
                                break;
                        case 38:
                                param->diff_feature_score = atof(optarg);//lignment_type = "pairwise";
                                break;


                        default:
                                break;
                        }
                        //printf ("option%d %s",option_index,long_options[option_index].name);
                        //if (optarg){
                        //	printf (" with arg %s\n", optarg);
                        //}
                        break;
                case 's':
                        param->gpo = atof(optarg);
                        //param->help_flag = 1;
                        break;
                case 'e':
                        param->gpe = atof(optarg);
                        break;
                case 'm':
                        param->secret = atof(optarg);
                        break;
                case 't':
                        param->tgpe = atof(optarg);
                        break;
                case 'z':
                        param->zlevel = atof(optarg);
                        break;

                case 'c':
                        param->sort = optarg;
                        break;
                case 'f':
                        param->format = optarg;
                        break;
                case 'r':
                        param->reformat = 1;
                        break;

                case 'd':
                        param->distance = optarg;
                        break;
                case 'b':
                        param->tree = optarg;
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
                case 'a':
                        param->gap_inc = atof(optarg);
                        break;
                case 'q':
                        param->quiet = 1;
                        break;
                case '?':
                        free_parameters(param);
                        exit(1);
                        break;
                default:
                        abort ();
                }
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
        if(!strncmp("pairwise", param->alignment_type, 8)){
                param->ntree = 1;
                if(param->tgpe == -1.0){
                        param->tgpe =  0.0f;
                }
        }
        if(param->gap_inc < 0.0){
                //fprintf(stderr,"%s\n", usage);
                fprintf(stderr,"Invalid parameter setting: gap_inc needs to be > 0 \n");
                exit(1);
        }

        if(param->quiet){
                fclose(stderr);
        }
        //fprintf(stderr,"%s", license);
        if (param->help_flag){
                //fprintf(stderr,"%s\n", usage);
                exit(1);
        }
        //exit(0);




        if (param->num_infiles == 0){
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
        struct alignment* aln = NULL;
        struct aln_param* ap = NULL;
        //float** dm;
        int** map = NULL;       /* holds all alignment paths  */
        int i,j;

        DECLARE_TIMER(t1);
        /* Step 1: read all input sequences & figure out output  */
        //LOG_MSG("Reading input.");
        RUNP(aln = detect_and_read_sequences(param));
        /* copy dna parameter to alignment */
        aln->dna = param->dna;
        if(param->ntree > aln->numseq){
                param->ntree = aln->numseq;
        }
        RUN(detect_dna(aln));
        LOG_MSG("Detected: %d sequences.", aln->numseq);
        //LOG_MSG("Output is %s in format %s.", param->outfile,param->format);
        //LOG_MSG("Is DNA: %d", aln->dna);
        param->dna = aln->dna;
        /* If we just want to reformat end here */
        if(param->reformat){
                for (i = 0 ;i < aln->numseq;i++){
                        aln->nsip[i] = i;
                        for (j = 0; j < aln->sl[i];j++){
                                aln->s[i][j] = 0;
                        }
                }
                param->format = "fasta";//param->reformat;
                RUN(output(aln, param));
                free_aln(aln);
                return OK;
        }
        /* allocate aln parameters  */
        RUNP(ap = init_ap(param,aln->numseq));
        //fprintf(stderr,"        %0.8f	gap open penalty\n",ap->gpo);
        //fprintf(stderr,"        %0.8f	gap extension\n",(float)gpe/10);
        //fprintf(stderr,"        %0.8f	gap extension\n",ap->gpe);
        //fprintf(stderr,"        %0.8f	terminal gap penalty\n",(float)tgpe/10);
        //fprintf(stderr,"        %0.8f	terminal gap penalty\n",ap->tgpe);
        //fprintf(stderr,"        %0.8f	bonus\n",param->secret/10);
        //fprintf(stderr,"        %0.8f	bonus\n",param->secret);

        //RUN(estimate_aln_param(aln, ap));

        LOG_MSG("Building guide tree.");
        START_TIMER(t1);

        //RUN(build_tree(aln,param,ap));
        RUN(build_tree_kmeans(aln,param,ap));
        STOP_TIMER(t1);
        LOG_MSG("Took %f sec.", GET_TIMING(t1));

        //LOG_MSG("Building guide tree.");
        //START_TIMER(t1);

        //STOP_TIMER(t1);
        //LOG_MSG("Took %f sec.", GET_TIMING(t1));

        /* Start alignment stuff */
        LOG_MSG("Aligning");
        START_TIMER(t1);
        RUNP(map = hirschberg_alignment(aln, ap));
        STOP_TIMER(t1);
        LOG_MSG("Took %f sec.", GET_TIMING(t1));

        RUN(weave(aln , map, ap->tree));


        RUN(output(aln, param));
        /* clean up map */
        for(i = 0; i < aln->num_profiles ;i++){
               if(map[i]){
                        MFREE(map[i]);
                }
        }
        MFREE(map);
        free_aln(aln);
        free_ap(ap);
        return OK;
ERROR:
        return FAIL;
}


int detect_dna(struct alignment* aln)
{
        int i;
        ASSERT(aln != NULL, "No alignment.");

        if(aln->dna){
                for(i = 0; i < aln->numseq;i++){
                        aln->dna = byg_detect(aln->s[i],aln->sl[i]);
                        if(aln->dna){
                                break;
                        }
                }
        }

        if(aln->dna == 1){
                //brief sanity check...
                for(i = 0; i < aln->numseq;i++){
                        if(aln->sl[i] < 6){
                                ERROR_MSG("Dna/Rna alignments are only supported for sequences longer than 6.");
                        }
                }
                RUN(make_dna(aln));
        }
        return OK;
ERROR:
        return FAIL;
}
