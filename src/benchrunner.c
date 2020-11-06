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
#define OPT_MATADD 9

#define OPT_RUNNAME 10

//#define BUFFER_LEN 512

struct parameters_br{
        char* testseq;
        char* refseq;
        char* program;
        char* scratch;
        char* output;
        char* run_name;
        char** options;
        int uniq;
        int n_options;
};

static int read_clean_alignments(char* ref_filename, char* test_filename, struct msa** ref, struct msa** test);
static int sort_msa_by_name(const void *a, const void *b);

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
        param->run_name = NULL;
        param->options = NULL;
        param->uniq = 0;
        param->n_options = 0;
        int help = 0;
        int c;

        while (1){
                static struct option long_options[] ={
                        {"test", required_argument,0,OPT_TESTALN},
                        {"ref", required_argument,0,OPT_REFALN},
                        {"program",  required_argument, 0, OPT_PROGRAM},
                        {"scratch",  required_argument, 0, OPT_TMPDIR},
                        {"runname",  required_argument, 0, OPT_RUNNAME },
                        {"output",  required_argument, 0, OPT_OUTPUT},
                        {"uniq",  required_argument, 0, 'u'},
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
                case OPT_RUNNAME:
                        param->run_name  = optarg;
                        break;
                case 'u':
                        param->uniq = atoi(optarg);
                        break;
                case 'h':
                        help = 1;
                        break;
                case '?':
                        LOG_MSG("Print %c", optopt);
                        help = 1;
                        break;
                default:
                        abort ();
                }
        }

        if (optind < argc){
                param->n_options = argc-optind;
                if(param->n_options != 1){
                        while (optind < argc){
                                fprintf(stdout,"%s",argv[optind++]);
                        }
                        ERROR_MSG("benrunner can only accept one additional argument");
                }
                MMALLOC(param->options, sizeof(char*) * param->n_options);
                c = 0;
                while (optind < argc){
                        param->options[c] =  argv[optind++];
                        c++;
                }
                param->n_options = c;
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
        if(param->options){
                MFREE(param->options);
        }
        MFREE(param);
        //RUN(timescore(testseq, refseq, program, scratch,output));
        return EXIT_SUCCESS;
ERROR:
        if(param->options){
                MFREE(param->options);
        }
        MFREE(param);
        /* RUN(print_help_score_and_align(argv)); */
        return EXIT_FAILURE;
}

int run_and_score(struct parameters_br* param)
{

        FILE* out_ptr = NULL;
        struct msa* ref_aln = NULL;
        struct msa* test_aln = NULL;

        char* basename = NULL;
        //struct parameters* param = NULL;
        auto int rc;
        //char path_buffer[BUFFER_LEN];
        char cmd[BUFSIZ*2];
        char ret[BUFSIZ];
        char progline[BUFSIZ];
        char* options = NULL;
        char* path = NULL;
        FILE* pipe = NULL;
        double time;
        double SP,TC;
        double average_seq_len = 0.0;

        struct timeval  tv1, tv2;

        int i,j,l;

        char *envpaths = getenv("PATH");

        path = realpath(param->scratch, NULL);
        if (path) {
                LOG_MSG("scratch is: %s", path);
        } else{
                MMALLOC(path, sizeof(char) * 10);
                rc = snprintf(path, 10, "%s",".");
        }

        if(param->n_options){

                l = strlen(param->options[0]);
                MMALLOC(options, sizeof(char) * ( l+1));

                for(i = 0; i < l;i++){
                        if(param->options[0][i] == '@' || param->options[0][i] == ':'){
                                options[i] = ' ';
                        }else{
                                options[i] = param->options[0][i];
                        }
                }
                options[l] = 0;

        }else{
                MMALLOC(options, sizeof(char) * 1);

                options[0] = 0;

        }



        if(!strcmp(param->program , "kalign")){
                if (system("which kalign")){
                        ERROR_MSG("kalign is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFSIZ, "kalign %s %s -f msf -o %s/test_%d.msf",options, param->testseq,path,param->uniq);

        }else if(!strcmp(param->program,"kalign2")){
                if (system("which kalign2")){
                        ERROR_MSG("kalign2 is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFFER_LEN*2, "kalign2 -i  %s -f msf -o %s/test.msf",param->testseq,path);
        }else if(!strcmp(param->program,"muscle")){
                if (system("which muscle3.8.31_i86linux32")){
                        ERROR_MSG("muscle3.8.31_i86linux32 is not found in your path:\n%s\n",envpaths);
                }
                rc = snprintf(cmd, BUFFER_LEN*2, "muscle3.8.31_i86linux32  -maxiters 2 -msf -in %s -out %s/test_%d.msf",param->testseq,path,param->uniq);
        }else if(!strcmp(param->program,"clustal")){
                if (system("which clustalo-1.2.4-Ubuntu-x86_64")){
                        ERROR_MSG("clustalo-1.2.4-Ubuntu-x86_64 is not found in your path:\n%s\n",envpaths);
                }

                rc = snprintf(cmd, BUFFER_LEN*2, "clustalo-1.2.4-Ubuntu-x86_64 --threads=1 --iterations=2 --outfmt=msf  --in %s --force --out %s/test.msf --verbose ",param->testseq,path);

        }else{
                ERROR_MSG("Program %s not recognize\n", param->program);
        }

        //LOG_MSG("Running: %s", cmd);

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
                snprintf( ret, BUFFER_LEN*2, "%s/test_%d.msf",path,param->uniq);
                RUN(read_clean_alignments(param->refseq, ret,  &ref_aln,&test_aln));


                average_seq_len = 0.0;
                for(j = 0; j < test_aln->numseq;j++){
                        average_seq_len += test_aln->sequences[j]->len;//   test_aln->sl[j];
                }

                average_seq_len = average_seq_len / (double) test_aln->numseq;


                snprintf(ret, BUFSIZ, "%s/refaln_%d.msf",path,param->uniq);
                write_msa(ref_aln, ret, FORMAT_MSF);

                snprintf(ret, BUFSIZ, "%s/evalaln_%d.msf",path,param->uniq);
                write_msa(test_aln, ret, FORMAT_MSF);

                /* snprintf(cmd,BUFSIZ*2 , "bali_score %s %s",  param->refseq, ret); */
                snprintf(cmd,BUFSIZ*2 , "mumsa %s/refaln_%d.msf %s/evalaln_%d.msf",path,param->uniq,path,param->uniq);

                // Setup our pipe for reading and execute our command.
                RUNP(pipe = popen(cmd,"r"));

                // Grab data from process execution
                while (fgets(ret, BUFSIZ, pipe)){
                        //fprintf(stdout,"%s", ret);
                }

                TC = 0.0;
                sscanf(ret, "%*s %*s %lf", &SP);
                rc = pclose(pipe);
                if(rc == EXIT_SUCCESS){
                        /* construct program name  */
                        l = strlen(options);

                        if(l > 1){
                                for(i = 0; i < l;i++){
                                        if(options[i] == ' '){
                                                options[i] = '_';
                                        }

                                }
                                snprintf(progline, BUFSIZ, "%s%s", param->program,options);
                        }else{
                                snprintf(progline, BUFSIZ, "%s", param->program);
                        }

                        if(!strcmp(param->output, "server")){

                                fprintf(stdout,"%f\n",SP);


                        }else if(!my_file_exists(param->output)){
                                RUN(tlfilename(param->refseq, &basename));
                                RUNP(out_ptr = fopen(param->output ,"w"));

                                fprintf(out_ptr,"Runname,Program,Alignment,AVGLEN,NUMSEQ,SP,TC,Time\n");
                                if(param->run_name){
                                        fprintf(out_ptr,"%s,%s,%s,%f,%d,%f,%f,%f\n",
                                                param->run_name,
                                                progline,
                                                basename,
                                                average_seq_len,
                                                test_aln->numseq,
                                                SP,
                                                TC,
                                                time);
                                }else{
                                        fprintf(out_ptr,"%s,%s,%s,%f,%d,%f,%f,%f\n",
                                                param->program,
                                                progline,
                                                basename,
                                                average_seq_len,
                                                test_aln->numseq,
                                                SP,
                                                TC,
                                                time);
                                }
                                fclose(out_ptr);
                                if(basename){
                                        MFREE(basename);
                                }

                        }else{
                                RUN(tlfilename(param->refseq, &basename));
                                RUNP(out_ptr = fopen(param->output,"a"));
                                if(param->run_name){
                                        fprintf(out_ptr,"%s,%s,%s,%f,%d,%f,%f,%f\n",
                                                param->run_name,
                                                progline,
                                                basename,
                                                average_seq_len,
                                                test_aln->numseq,
                                                SP,
                                                TC,
                                                time);
                                }else{
                                        fprintf(out_ptr,"%s,%s,%s,%f,%d,%f,%f,%f\n",
                                                param->program,
                                                progline,
                                                basename,
                                                average_seq_len,
                                                test_aln->numseq,
                                                SP,
                                                TC,
                                                time);
                                }
                                fclose(out_ptr);
                                if(basename){
                                        MFREE(basename);
                                }
                        }

                }else{
                        ERROR_MSG("Running %s failed", cmd);
                }


                //LOG_MSG("SP:%f TC:%f", SP,TC);

                //free_aln(test_aln);
                //MFREE(param->outfile);
                //MFREE(param->infile[0]);
                //free_parameters(param);
                free_msa(ref_aln);
                free_msa(test_aln );

        } else {  // EXIT_FAILURE is not used by all programs, maybe needs some adaptation.
                ERROR_MSG("Running:\n%s failed.", cmd);
                /* WARNING_MSG("Error code was: %d", rc); */
        }

        snprintf(ret, BUFFER_LEN, "%s/refaln_%d.msf",path,param->uniq);
        if(my_file_exists(ret)){
                remove(ret);
        }
        snprintf(ret, BUFFER_LEN, "%s/evalaln_%d.msf",path,param->uniq);
        if(my_file_exists(ret)){
                remove(ret);
        }

        snprintf( ret, BUFFER_LEN*2, "%s/test_%d.msf",path,param->uniq);
        if(my_file_exists(ret)){
                remove(ret);
        }

        if(path){
                MFREE(path);
        }
        MFREE(options);
        //MFREE(path);

        return OK;
ERROR:
        WARNING_MSG("Failed");
        WARNING_MSG("ref: %s", param->refseq);

        snprintf(ret, BUFFER_LEN, "%s/refaln_%d.msf",path,param->uniq);
        if(my_file_exists(ret)){
                remove(ret);
        }
        snprintf(ret, BUFFER_LEN, "%s/evalaln_%d.msf",path,param->uniq);
        if(my_file_exists(ret)){
                /* remove(ret); */
        }

        snprintf( ret, BUFFER_LEN*2, "%s/test_%d.msf",path,param->uniq);
        if(my_file_exists(ret)){
                /* remove(ret); */
        }
        if(options){
                MFREE(options);
        }
        if(ref_aln){
                free_msa(ref_aln);
        }
        if(test_aln){
                free_msa(test_aln);
        }
        if(path){
                MFREE(path);
        }
        return FAIL;
}



int print_help_score_and_align(char **argv)
{
        const char usage[] = "  -test <test sequences> -ref <reference alignment> -program <kalign|clustalo|muscle> -o <outfile>";

        char* basename = NULL;


        RUN(tlfilename(argv[0], &basename));
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename ,usage);
        fprintf(stdout,"NOTE: the program appends results to the output file.\n\n");
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--scratch","Scratch directory." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--runname","Name of run." ,"[NA]"  );


        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open ." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extend." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Gap extend terminal." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--matadd","Addition to sub matrix." ,"[NA]"  );

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

int read_clean_alignments(char* ref_filename, char* test_filename, struct msa** ref, struct msa** test)
{
        struct msa* ref_aln;
        struct msa* tmp_aln;
        struct msa* test_aln;

        int i,j,c,n;

        RUN(read_input( ref_filename, &ref_aln));
        RUN(read_input(test_filename, &test_aln));//  detect_and_read_sequences(param));

        /* sort test and ref alignment  */
        qsort(test_aln->sequences , test_aln->numseq, sizeof(struct msa_seq *), sort_msa_by_name);
        qsort(ref_aln->sequences , ref_aln->numseq, sizeof(struct msa_seq *), sort_msa_by_name);

        /* if test contains more sequences than ref
           copy ref; free sequences; copy sequences from test over
        */
        /* Case one -> all good */
        if(ref_aln->numseq == test_aln->numseq){
                *ref = ref_aln;
                *test = test_aln;

        }else if(ref_aln->numseq < test_aln->numseq){ /* Case 2 more test then ref  */

                /* read ref in again */
                RUN(read_input( ref_filename, &tmp_aln));

                /* free sequences  */
                for(i = 0; i < tmp_aln->numseq;i++){
                        struct msa_seq* seq;
                        seq = tmp_aln->sequences[i];
                        MFREE(seq->name);
                        MFREE(seq->seq);
                        MFREE(seq->s);
                        MFREE(seq->gaps);
                        MFREE(seq);
                }
                /* copy test sequences matching ref names over. */
                for(i = 0; i < ref_aln->numseq;i++){
                        c = -1;
                        n = 0;
                        for(j = 0; j < test_aln->numseq;j++){
                                //if(  ref_aln->lsn[i] == test_aln->lsn[j]){
                                if(!strncmp( ref_aln->sequences[i]->name, test_aln->sequences[j]->name, MSA_NAME_LEN)){
                                        n++;
                                        c = j;

                                }

                        }
                        if(c == -1){
                                ERROR_MSG("Sequence %s not found in test alignment!", ref_aln->sequences[i]->name);
                        }else if( n > 1){
                                ERROR_MSG("Multiple sequences with name %s found",ref_aln->sequences[i]->name);
                        }else{
                                tmp_aln->sequences[i] = test_aln->sequences[c];
                                test_aln->sequences[c] = NULL;
                        }
                }
                free_msa(test_aln);
                *ref = ref_aln;
                *test = tmp_aln;
        }else{                  /* case 3: more ref than test??? should never happen */
                ERROR_MSG("Reference contains more sequences than test aln");
        }
        return OK;
ERROR:
        if(tmp_aln){
                free_msa(tmp_aln);
        }
        if(test_aln){
                free_msa(test_aln);
        }
        if(ref_aln){
                free_msa(ref_aln);
        }
        return FAIL;
}

int sort_msa_by_name(const void *a, const void *b)
{
        struct msa_seq* const *one = a;

        struct msa_seq* const *two = b;


        if(strncmp((*one)->name, (*two)->name, MSA_NAME_LEN) < 0){
                return -1;
        }else{
                return 1;
        }
}
