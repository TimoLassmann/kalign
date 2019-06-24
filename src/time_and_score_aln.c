#include "global.h"

#include <sys/time.h>
#include <getopt.h>

#include <limits.h>
#include <stdlib.h>

#include "align_io.h"
#include "misc.h"

#define OPT_TESTALN 1
#define OPT_REFALN 2
#define OPT_PROGRAM 3
#define OPT_TMPDIR 4


int timescore(char* test,char* ref, char* program,char* scratch);

int print_help_score_and_align(char **argv);



int main(int argc, char *argv[])
{
        char* testseq = NULL;
        char* refseq = NULL;
        char* program = NULL;
        char* scratch = NULL;

        int help = 0;
        int c;

        while (1){
                static struct option long_options[] ={
                        {"test", required_argument,0,OPT_TESTALN},
                        {"ref", required_argument,0,OPT_REFALN},
                        {"program",  required_argument, 0, OPT_PROGRAM},
                        {"scratch",  required_argument, 0, OPT_TMPDIR},
                        {"help",   no_argument,0,'h'},
                        {"quiet",  0, 0, 'q'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"h",long_options, &option_index);
                //c = getopt (argc, argv, "hi:o:");
                /* Detect the end of the options. */

                if (c == -1){
                        break;
                }
                switch(c) {
                case  OPT_TESTALN:
                        testseq = optarg;
                        break;

                case OPT_REFALN:
                        refseq = optarg;
                        break;
                case OPT_PROGRAM:
                        program = optarg;
                        break;
                case OPT_TMPDIR:
                        scratch  = optarg;
                        break;
                case 'h':
                        help = 1;
                        break;
                case '?':
                        help = 1;

                        break;
                default:

                        abort ();
                }
        }

        if(help){
                RUN(print_help_score_and_align(argv));
        }
        /* check options  */
        if(!testseq){

                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No test file");


        }else{
                if(!my_file_exists(testseq)){
                        ERROR_MSG("test file: %s does not exist",testseq);
                }
        }

        if(!refseq){

                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No ref file");


        }else{
                if(!my_file_exists(refseq)){
                        ERROR_MSG("ref file: %s does not exist",testseq);
                }
        }

        if(!program){
                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No program specified");
        }

        if(scratch){
                struct stat sb;
                if(stat(scratch, &sb) == 0 && S_ISDIR(sb.st_mode)) {
                        LOG_MSG("use tmpdir: %s", scratch);
                }else{
                        ERROR_MSG("tmpdir %s does not exist.", scratch);
                }
        }

        RUN(timescore(testseq, refseq, program, scratch));
        return EXIT_SUCCESS;
ERROR:
        WARNING_MSG("Something went wrong. Use this program like this:\n\n");
        RUN(print_help_score_and_align(argv));
        return EXIT_FAILURE;
}

int timescore(char* test,char* ref, char* program,char* scratch)
{

        struct alignment* ref_aln = NULL;
        struct alignment* test_aln = NULL;
        struct parameters* param = NULL;

        char path_buffer[BUFFER_LEN+BUFFER_LEN];
        char cmd[BUFFER_LEN];
        char ret[BUFFER_LEN];
        char* path = NULL;
        FILE* pipe = NULL;
        double time;
        double SP,TC;

        struct timeval  tv1, tv2;
        int pos_test;
        int i,j,c;

        path = realpath(scratch, path_buffer);
        if (path) {
                LOG_MSG("scratch is: %s", path);
        }else{
                snprintf(path_buffer, BUFFER_LEN+BUFFER_LEN, "%s",".");
                path= path_buffer;
        }



        if(!strcmp(program , "kalign")){
                snprintf(cmd, BUFFER_LEN, "kalign %s -o %s/test.msf",test,path);
        }

        LOG_MSG("Running: %s", cmd);

        gettimeofday(&tv1, NULL);

        /* run program */
        RUNP(pipe = popen(cmd,"r"));

        while (fgets(ret, BUFFER_LEN, pipe)){
                fprintf(stderr,"%s", ret);
        }
        pclose(pipe);

        gettimeofday(&tv2, NULL);

        time = (double) (tv2.tv_usec - tv1.tv_usec) / (double) CLOCKS_PER_SEC  + (double) (tv2.tv_sec - tv1.tv_sec);

        printf ("Total time = %f seconds\n",time);


        /* extract all sequences from the test alignment that are found in the reference...  */
        RUNP(param = init_param());

        param->infile = NULL;
        param->num_infiles = 1;
        MMALLOC(param->infile, sizeof(char*) * param->num_infiles);

        param->infile[0] = NULL;
        MMALLOC(param->infile[0], sizeof(char) * BUFFER_LEN);

        snprintf(param->infile[0], BUFFER_LEN, "%s",ref);

        RUNP(ref_aln = detect_and_read_sequences(param));

        snprintf(param->infile[0], BUFFER_LEN, "%s/test.msf",path);

        RUNP(test_aln = detect_and_read_sequences(param));

        /* find all reference sequences in the test alignment  */
        /* delete all other sequences  */
        /* print out the ref sequences now aligned with the test MSA program */


        pos_test = 0;
        for(i = 0; i < ref_aln->numseq;i++){
                LOG_MSG("Searching for: %s", ref_aln->sn[i]);
                for(j = 0; j < test_aln->numseq;j++){
                        if(  ref_aln->lsn[i] == test_aln->lsn[j]){
                        if(!strcmp( ref_aln->sn[i],  test_aln->sn[j])){
                                if(pos_test != j){
                                LOG_MSG("Found:: %s (%d)", test_aln->sn[j], pos_test);

                                /* copy name */
                                if(test_aln->lsn[pos_test] < test_aln->lsn[j]){
                                        MREALLOC(test_aln->sn[pos_test], test_aln->lsn[j]);


                                        test_aln->lsn[pos_test] = test_aln->lsn[j];
                                }
                                snprintf(test_aln->sn[pos_test] ,test_aln->lsn[pos_test],"%s",test_aln->sn[j]);
                                /* copy sequence */

                                if(test_aln->sl[pos_test] < test_aln->sl[j]){
                                        MREALLOC(test_aln->seq[pos_test]  , test_aln->sl[j]);


                                        test_aln->sl[pos_test] = test_aln->sl[j];
                                }
                                snprintf(test_aln->seq[pos_test],test_aln->sl[pos_test],"%s",test_aln->seq[j]);

                                /* copy gaps  */

                                for(c = 0; c < test_aln->sl[pos_test]+1;c++){
                                        test_aln->gaps[pos_test][c] = test_aln->gaps[j][c];
                                }
                                }
                                pos_test++;
                                break;
                        }
                        }

                }
        }

        int org_len;


        org_len = test_aln->numseq;
        test_aln->numseq = pos_test;

        MMALLOC(param->outfile, sizeof(char) * BUFFER_LEN);

        snprintf(param->outfile, BUFFER_LEN, "%s/evalaln.msf",path);
        param->format = "msf";
        for (i = 0; i < test_aln->numseq;i++){
                test_aln->nsip[i] = i;
        }

        output(test_aln,param);

        test_aln->numseq = org_len;


        /* score */


        snprintf(cmd, BUFFER_LEN, "bali_score %s %s",  ref, param->outfile);
        // Execute a process listing



        // Setup our pipe for reading and execute our command.
        RUNP(pipe = popen(cmd,"r"));


        // Grab data from process execution


        while (fgets(ret, BUFFER_LEN, pipe)){
                //fprintf(stdout,"%s", ret);
        }

        fprintf(stdout,"-%s-", ret);

        sscanf(ret, "%*s %*s %lf %lf", &SP,&TC);
        pclose(pipe);

        LOG_MSG("SP:%f TC:%f", SP,TC);

        free_aln(ref_aln);
        free_aln(test_aln);
        MFREE(param->outfile);
        MFREE(param->infile[0]);
        free_parameters(param);

        return OK;
ERROR:
        return FAIL;
}



int print_help_score_and_align(char **argv)
{
        const char usage[] = "  -test <test sequences> -ref <reference alignment> -program <kalign|clustalo|muscle>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--scratch","Scratch directory." ,"[NA]"  );

        return OK;
}
