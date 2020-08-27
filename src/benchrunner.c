
#include "global.h"

#include "tldevel.h"
#include "tlmisc.h"

#include "msa.h"
#include <sys/time.h>
#include <sys/stat.h>
#include <getopt.h>

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define OPT_TESTALN 1
#define OPT_REFALN 2
#define OPT_PROGRAM 3
#define OPT_TMPDIR 4
#define OPT_OUTPUT 5

#define OPT_GPO 6
#define OPT_GPE 7
#define OPT_TGPE 8

//#define BUFFER_LEN 512

struct parameters_br{
        char* testseq;
        char* refseq;
        char* program;
        char* scratch;
        char* output;

        float gpo;
        float gpe;
        float tgpe;
};



int print_help_score_and_align(char **argv);

int run_and_score(struct parameters_br* param);

int main(int argc, char *argv[])
{

        struct parameters_br* param = NULL;
        MMALLOC(param, sizeof(struct parameters_br));
        param->testseq = NULL;
        param->refseq = NULL;
        param->program = NULL;
        param->scratch = NULL;
        param->output = NULL;

        param->gpo = FLT_MAX;
        param->gpe = FLT_MAX;
        param->tgpe = FLT_MAX;
        
        int help = 0;
        int c;

        while (1){
                static struct option long_options[] ={
                        {"test", required_argument,0,OPT_TESTALN},
                        {"ref", required_argument,0,OPT_REFALN},
                        {"program",  required_argument, 0, OPT_PROGRAM},
                        {"scratch",  required_argument, 0, OPT_TMPDIR},
                        {"output",  required_argument, 0, OPT_OUTPUT},
                        {"gpo",  required_argument, 0, OPT_GPO},
                        {"gpe",  required_argument, 0, OPT_GPE},
                        {"tgpe",  required_argument, 0, OPT_TGPE},
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
                        param->testseq = optarg;
                        break;

                case OPT_REFALN:
                        param->refseq = optarg;
                        break;
                case OPT_PROGRAM:
                        param->program = optarg;
                        break;
                case OPT_TMPDIR:
                        param->scratch  = optarg;
                        break;
                case OPT_OUTPUT:
                        param->output  = optarg;
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
        if(!param->testseq){

                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No test file\n");


        }else{
                if(!my_file_exists(param->testseq)){
                        ERROR_MSG("test file: %s does not exist\n",param->testseq);
                }
        }

        if(!param->refseq){

                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No ref file\n");


        }else{
                if(!my_file_exists(param->refseq)){
                        ERROR_MSG("ref file: %s does not exist",param->testseq);
                }
        }

        if(!param->program){
                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No program specified\n");
        }
        if(!param->output){
                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No output specified\n");

        }

        if(param->scratch){
                struct stat sb;
                if(stat(param->scratch, &sb) == 0 && S_ISDIR(sb.st_mode)) {
                        LOG_MSG("use tmpdir: %s", param->scratch);
                }else{
                        ERROR_MSG("tmpdir %s does not exist.", param->scratch);
                }
        }


        RUN(run_and_score(param));
        MFREE(param);
        //RUN(timescore(testseq, refseq, program, scratch,output));
        return EXIT_SUCCESS;
ERROR:
        WARNING_MSG("Something went wrong. Use this program like this:\n\n");
        RUN(print_help_score_and_align(argv));
        return EXIT_FAILURE;
}

int run_and_score(struct parameters_br* param)
{

        FILE* out_ptr = NULL;
        struct msa* ref_aln = NULL;
        struct msa* test_aln = NULL;

        struct msa_seq* seq_ptr = NULL;
        //struct parameters* param = NULL;
        auto int rc;
        //char path_buffer[BUFFER_LEN];
        char cmd[BUFFER_LEN*2];
        char ret[BUFFER_LEN*2];
        char* path = NULL;
        FILE* pipe = NULL;
        double time;
        double SP,TC;
        double average_seq_len = 0.0;

        struct timeval  tv1, tv2;
        int pos_test;
        int i,j;

        char *envpaths = getenv("PATH");

        path = realpath(param->scratch, NULL);
        if (path) {
                LOG_MSG("scratch is: %s", path);
        }else{
                MMALLOC(path, sizeof(char) * 10);
                rc = snprintf(path, 10, "%s",".");
        }



        if(!strcmp(param->program , "kalign")){
                if (system("which kalign")){
                        ERROR_MSG("kalign is not found in your path:\n%s\n",envpaths);
                }

                if(param->gpo != FLT_MAX && param->gpe != FLT_MAX && param->tgpe != FLT_MAX){
                        rc = snprintf(cmd, BUFFER_LEN*2, "kalign -gpo %f -gpe %f -tgpe %f %s -f msf -o %s/test.msf",param->gpo,param->gpe, param->tgpe, param->testseq,path);
                }else{
                rc = snprintf(cmd, BUFFER_LEN*2, "kalign %s -f msf -o %s/test.msf", param->testseq,path);
                }
        }else if(!strcmp(param->program,"kalign2")){
                if (system("which kalign2")){
                            ERROR_MSG("kalign2 is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFFER_LEN*2, "kalign2 -i  %s -f msf -o %s/test.msf",param->testseq,path);
        }else if(!strcmp(param->program,"muscle")){
                if (system("which muscle3.8.31_i86linux32")){
                        ERROR_MSG("muscle3.8.31_i86linux32 is not found in your path:\n%s\n",envpaths);
                }
                rc = snprintf(cmd, BUFFER_LEN*2, "muscle3.8.31_i86linux32  -maxiters 2 -msf -in %s -out %s/test.msf",param->testseq,path);
        }else if(!strcmp(param->program,"clustal")){
                if (system("which clustalo-1.2.4-Ubuntu-x86_64")){
                        ERROR_MSG("clustalo-1.2.4-Ubuntu-x86_64 is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFFER_LEN*2, "clustalo-1.2.4-Ubuntu-x86_64 --threads=1 --iterations=2 --outfmt=msf  --in %s --force --out %s/test.msf --verbose ",param->testseq,path);

        }else{
                ERROR_MSG("Program %s not recognize\n", param->program);
        }

        LOG_MSG("Running: %s", cmd);

        gettimeofday(&tv1, NULL);


        RUNP(pipe = popen(cmd,"r"));
        //if (system("which gnuplot"))
        while (fgets(ret, BUFFER_LEN, pipe)){
                //fprintf(stderr,"%s", ret);
        }

        rc = pclose(pipe);

        if (rc == EXIT_SUCCESS) { // == 0
                LOG_MSG("Running:\n%s was a success", cmd);
                gettimeofday(&tv2, NULL);

                time = (double) (tv2.tv_usec - tv1.tv_usec) / (double) CLOCKS_PER_SEC  + (double) (tv2.tv_sec - tv1.tv_sec);

                printf ("Total time = %f seconds\n",time);
                // extract all sequences from the test alignment that are found in the reference...  
                LOG_MSG("read reference alignment");

                RUN(read_input( param->refseq , &ref_aln));
                //RUN(read_input(param->infile[0],&ref_aln));//  detect_and_read_sequences(param));

                snprintf( ret, BUFFER_LEN*2, "%s/test.msf",path);
                LOG_MSG("read test alignment");
                RUN(read_input(ret,&test_aln));//  detect_and_read_sequences(param));


                average_seq_len = 0.0;
                for(j = 0; j < test_aln->numseq;j++){
                        average_seq_len += test_aln->sequences[j]->len;//   test_aln->sl[j];
                }

                average_seq_len = average_seq_len / (double) test_aln->numseq;
                pos_test = 0;
                LOG_MSG("Search for matching sequences ");
                for(i = 0; i < ref_aln->numseq;i++){
                        //LOG_MSG("Searching for: %s", ref_aln->sequences[i]->name);// sn[i]);
                        for(j = 0; j < test_aln->numseq;j++){
                                if(strnlen(ref_aln->sequences[i]->name, MSA_NAME_LEN) == strnlen(test_aln->sequences[j]->name, MSA_NAME_LEN)){
                                        //if(  ref_aln->lsn[i] == test_aln->lsn[j]){
                                        if(!strcmp( ref_aln->sequences[i]->name, test_aln->sequences[j]->name)){

                                                //LOG_MSG("Found: %s",ref_aln->sequences[i]->name);
                                                seq_ptr = test_aln->sequences[pos_test];
                                                test_aln->sequences[pos_test] = test_aln->sequences[j];
                                                test_aln->sequences[j] = seq_ptr;
                                                pos_test++;
                                        }
                                }

                        }
                }

                int org_len;


                org_len = test_aln->numseq;
                test_aln->numseq = pos_test;

                

                snprintf(ret, BUFFER_LEN, "%s/evalaln.msf",path);
                write_msa(test_aln, ret, FORMAT_MSF);

                
                //output(test_aln,param);

                test_aln->numseq = org_len;



                snprintf(cmd, BUFFER_LEN, "bali_score %s %s",  param->refseq, ret);
                // Execute a process listing



                // Setup our pipe for reading and execute our command.
                RUNP(pipe = popen(cmd,"r"));


                // Grab data from process execution


                while (fgets(ret, BUFFER_LEN, pipe)){
                        //fprintf(stdout,"%s", ret);
                }


                sscanf(ret, "%*s %*s %lf %lf", &SP,&TC);
                pclose(pipe);
                if(!strcmp(param->output, "stdout")){
                        out_ptr = stdout;
                }else if(!my_file_exists(param->output)){
                        RUNP(out_ptr = fopen(param->output ,"w"));
                        fprintf(out_ptr,"Program,Alignment,AVGLEN,NUMSEQ,SP,TC,Time\n");
                }else{
                        RUNP(out_ptr = fopen(param->output,"a"));
                }

                fprintf(out_ptr,"%s,%s,%f,%d,%f,%f,%f\n",param->program,basename(param->refseq), average_seq_len,test_aln->numseq,SP,TC,time);
                fclose(out_ptr);

                LOG_MSG("SP:%f TC:%f", SP,TC);
                free_msa(ref_aln);
                free_msa(test_aln );
                //free_aln(test_aln);
                //MFREE(param->outfile);
                //MFREE(param->infile[0]);
                //free_parameters(param);


        } else {  // EXIT_FAILURE is not used by all programs, maybe needs some adaptation.
                WARNING_MSG("Running:\n%s failed.", cmd);
                WARNING_MSG("Error code was: %d", rc);
        }





        MFREE(path);

        return OK;
ERROR:
        return FAIL;
}



int print_help_score_and_align(char **argv)
{
        const char usage[] = "  -test <test sequences> -ref <reference alignment> -program <kalign|clustalo|muscle> -o <outfile>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"NOTE: the program appends results to the output file.\n\n");
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--scratch","Scratch directory." ,"[NA]"  );

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open ." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extend." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Gap extend terminal." ,"[NA]"  );

        return OK;
}
