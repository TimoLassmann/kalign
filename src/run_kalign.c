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

#define OPT_SET 1

int run_kalign(struct parameters* param);
int detect_dna(struct alignment* aln);


int main(int argc, char *argv[])
{
        int c;
        struct parameters* param = NULL;

        RUNP(param = init_param());

        while (1){
                static struct option long_options[] ={
                        {"set", required_argument,0,OPT_SET},
                        {"format",  required_argument, 0, 'f'},
                        {"reformat",  0, 0, 'r'},
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
                c = getopt_long_only (argc, argv,"i:o:r:f:hq",long_options, &option_index);
                //c = getopt (argc, argv, "hi:o:");
                /* Detect the end of the options. */

                if (c == -1){

                        break;
                }
                switch(c) {
                case OPT_SET:
                        param->param_set = atoi(optarg);
                        break;
                case 'f':
                        param->format = optarg;
                        break;
                case 'r':
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
        //aln->dna = param->dna;
        RUN(detect_dna(aln));
        LOG_MSG("Detected: %d sequences.", aln->numseq);
        //LOG_MSG("Output is %s in format %s.", param->outfile,param->format);
        //LOG_MSG("Is DNA: %d", aln->dna);
        //param->dna = aln->dna;
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



        //param->param_set = 4;
        /* allocate aln parameters  */
        RUNP(ap = init_ap(param,aln->numseq,aln->L));
        //fprintf(stderr,"        %0.8f	gap open penalty\n",ap->gpo);
        //fprintf(stderr,"        %0.8f	gap extension\n",(float)gpe/10);
        //fprintf(stderr,"        %0.8f	gap extension\n",ap->gpe);
        //fprintf(stderr,"        %0.8f	terminal gap penalty\n",(float)tgpe/10);
        //fprintf(stderr,"        %0.8f	terminal gap penalty\n",ap->tgpe);
        //fprintf(stderr,"        %0.8f	bonus\n",param->secret/10);
        //fprintf(stderr,"        %0.8f	bonus\n",param->secret);
        /*LOG_MSG("Estimating alignment param.");
        START_TIMER(t1);

        RUN(estimate_aln_param(aln, ap));
        STOP_TIMER(t1);
        LOG_MSG("Took %f sec.", GET_TIMING(t1));*/




        LOG_MSG("Building guide tree.");
        START_TIMER(t1);

        //random_tree(ap, aln->numseq);


        //param->dist_method = KALIGNDIST_WU;
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

        int iter;
        for(iter = 0; iter < 0;iter++){
                RUN(weave(aln , map, ap->tree));
                param->dist_method = KALIGNDIST_ALN;
                build_tree(aln, param,ap);
                clean_aln(aln);


                for(i = 0; i < aln->num_profiles ;i++){
                        if(map[i]){
                                MFREE(map[i]);
                        }
                }
                MFREE(map);
                map = NULL;


                LOG_MSG("Aligning");
                START_TIMER(t1);
                RUNP(map = hirschberg_alignment(aln, ap));
                STOP_TIMER(t1);
                LOG_MSG("Took %f sec.", GET_TIMING(t1));
        }

                int max_len = 0;
        for(i = 0; i < aln->numseq;i++){
                if(aln->sl[i] > max_len){
                        max_len = aln->sl[i];
                }

        }
        aln->gaps = galloc(aln->gaps,aln->numseq,max_len+1,0);

        RUN(weave(aln , map, ap->tree));

        /* clean up map */
        for(i = 0; i < aln->num_profiles ;i++){
               if(map[i]){
                       MFREE(map[i]);
               }
        }
        MFREE(map);
        map = NULL;

        /* alignment output order .... */
        for (i = 0; i < aln->numseq;i++){
                aln->nsip[i] = i;
        }

        RUN(output(aln, param));

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
