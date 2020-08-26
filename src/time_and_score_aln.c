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

#include "tlmisc.h"

#include "msa.h"
#include <sys/time.h>
#include <sys/stat.h>
#include <getopt.h>

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "align_io.h"
#include "misc.h"

#define OPT_TESTALN 1
#define OPT_REFALN 2
#define OPT_PROGRAM 3
#define OPT_TMPDIR 4
#define OPT_OUTPUT 5

int timescore(char* test,char* ref, char* program,char* scratch,char* out_file_name);
int print_help_score_and_align(char **argv);

int main(int argc, char *argv[])
{
        char* testseq = NULL;
        char* refseq = NULL;
        char* program = NULL;
        char* scratch = NULL;
        char* output = NULL;
        int help = 0;
        int c;

        while (1){
                static struct option long_options[] ={
                        {"test", required_argument,0,OPT_TESTALN},
                        {"ref", required_argument,0,OPT_REFALN},
                        {"program",  required_argument, 0, OPT_PROGRAM},
                        {"scratch",  required_argument, 0, OPT_TMPDIR},
                        {"output",  required_argument, 0, OPT_OUTPUT},
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
                case OPT_OUTPUT:
                        output  = optarg;
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
                ERROR_MSG("No test file\n");


        }else{
                if(!my_file_exists(testseq)){
                        ERROR_MSG("test file: %s does not exist\n",testseq);
                }
        }

        if(!refseq){

                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No ref file\n");


        }else{
                if(!my_file_exists(refseq)){
                        ERROR_MSG("ref file: %s does not exist",testseq);
                }
        }

        if(!program){
                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No program specified\n");
        }
        if(!output){
                RUN(print_help_score_and_align(argv));
                ERROR_MSG("No output specified\n");

        }

        if(scratch){
                struct stat sb;
                if(stat(scratch, &sb) == 0 && S_ISDIR(sb.st_mode)) {
                        LOG_MSG("use tmpdir: %s", scratch);
                }else{
                        ERROR_MSG("tmpdir %s does not exist.", scratch);
                }
        }

        RUN(timescore(testseq, refseq, program, scratch,output));
        return EXIT_SUCCESS;
ERROR:
        WARNING_MSG("Something went wrong. Use this program like this:\n\n");
        RUN(print_help_score_and_align(argv));
        return EXIT_FAILURE;
}

int timescore(char* test,char* ref, char* program,char* scratch,char* out_file_name)
{

        FILE* out_ptr = NULL;
        struct msa* ref_aln = NULL;
        struct msa* test_aln = NULL;

        struct msa_seq* seq_ptr = NULL;
        struct parameters* param = NULL;
        auto int rc;
        //char path_buffer[BUFFER_LEN];
        char cmd[BUFFER_LEN*2];
        char ret[BUFFER_LEN];
        char* path = NULL;
        FILE* pipe = NULL;
        double time;
        double SP,TC;
        double average_seq_len = 0.0;

        struct timeval  tv1, tv2;
        int pos_test;
        int i,j;

        char *envpaths = getenv("PATH");

        path = realpath(scratch, NULL);
        if (path) {
                LOG_MSG("scratch is: %s", path);
        }else{
                MMALLOC(path, sizeof(char) * 10);
                rc = snprintf(path, 10, "%s",".");
        }



        if(!strcmp(program , "kalign")){
                if (system("which kalign")){
                        ERROR_MSG("kalign is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFFER_LEN*2, "kalign %s -f msf -o %s/test.msf",test,path);
        }else if(!strcmp(program,"kalign2")){
                if (system("which kalign2")){
                            ERROR_MSG("kalign2 is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFFER_LEN*2, "kalign2 -i  %s -f msf -o %s/test.msf",test,path);
        }else if(!strcmp(program,"muscle")){
                if (system("which muscle3.8.31_i86linux32")){
                        ERROR_MSG("muscle3.8.31_i86linux32 is not found in your path:\n%s\n",envpaths);
                }
                rc = snprintf(cmd, BUFFER_LEN*2, "muscle3.8.31_i86linux32  -maxiters 2 -msf -in %s -out %s/test.msf",test,path);
        }else if(!strcmp(program,"clustal")){
                if (system("which clustalo-1.2.4-Ubuntu-x86_64")){
                        ERROR_MSG("clustalo-1.2.4-Ubuntu-x86_64 is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFFER_LEN*2, "clustalo-1.2.4-Ubuntu-x86_64 --threads=1 --iterations=2 --outfmt=msf  --in %s --force --out %s/test.msf --verbose ",test,path);

        }else{
                ERROR_MSG("Program %s not recognize\n", program);
        }

        LOG_MSG("Running: %s", cmd);

        gettimeofday(&tv1, NULL);

        /* run program */
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
                /* extract all sequences from the test alignment that are found in the reference...  */
                RUNP(param = init_param());

                param->infile = NULL;
                param->num_infiles = 1;
                MMALLOC(param->infile, sizeof(char*) * param->num_infiles);

                param->infile[0] = NULL;
                MMALLOC(param->infile[0], sizeof(char) * BUFFER_LEN*2);

                snprintf(param->infile[0], BUFFER_LEN*2, "%s",ref);
                LOG_MSG("read reference alignment");
                RUN(read_input(param->infile[0],&ref_aln));//  detect_and_read_sequences(param));

                snprintf(param->infile[0], BUFFER_LEN*2, "%s/test.msf",path);
                LOG_MSG("read test alignment");
                RUN(read_input(param->infile[0],&test_aln));//  detect_and_read_sequences(param));

                /* find all reference sequences in the test alignment  */
                /* delete all other sequences  */
                /* print out the ref sequences now aligned with the test MSA program */


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
                                                /* swap test sequences to beginning.  */
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

                MMALLOC(param->outfile, sizeof(char) * BUFFER_LEN*2);

                snprintf(param->outfile, BUFFER_LEN*2, "%s/evalaln.msf",path);
                param->format = "msf";
                for (i = 0; i < test_aln->numseq;i++){
                        test_aln->nsip[i] = i;
                }

                write_msa(test_aln, param->outfile, FORMAT_MSF);
                //output(test_aln,param);

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


                sscanf(ret, "%*s %*s %lf %lf", &SP,&TC);
                pclose(pipe);
                if(!my_file_exists(out_file_name)){
                        RUNP(out_ptr = fopen(out_file_name,"w"));
                        fprintf(out_ptr,"Program,Alignment,AVGLEN,NUMSEQ,SP,TC,Time\n");
                }else{
                        RUNP(out_ptr = fopen(out_file_name,"a"));
                }

                fprintf(out_ptr,"%s,%s,%f,%d,%f,%f,%f\n",program,basename(ref), average_seq_len,test_aln->numseq,SP,TC,time);
                fclose(out_ptr);

                LOG_MSG("SP:%f TC:%f", SP,TC);
                free_msa(ref_aln);
                free_msa(test_aln );
                //free_aln(test_aln);
                MFREE(param->outfile);
                MFREE(param->infile[0]);
                free_parameters(param);


        } else {  // EXIT_FAILURE is not used by all programs, maybe needs some adaptation.
                WARNING_MSG("Running:\n%s failed.", cmd);
                WARNING_MSG("Error code was: %d", rc);
        }
        /*return result;
          pclose(pipe);*/






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

        return OK;
}
