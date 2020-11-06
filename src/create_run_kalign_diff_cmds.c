#include "tldevel.h"
#include "tlmisc.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

struct parameters{
        char* inputdir;
        char* scratchdir;
        char* outfile;
};

static int help(char * argv[]);

int help(char * argv[])
{
        const char usage[] = " -i <input dir> -o <out file> -s <scratchdir>";
        fprintf(stdout,"\nUsage: %s %s\n\n",argv[0] ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"\nTo plot results :\n\n");
        fprintf(stdout,"dat = read.csv(\"results_DNA.csv\",header = T)\n");
        fprintf(stdout,"l = str_split(dat$Alignment, \"_\")\n");
        fprintf(stdout,"c = data.frame(sapply(l, \"[[\", 1))\n");
        fprintf(stdout,"colnames(c) = c(\"Case\");\n");
        fprintf(stdout,"dat = cbind(dat,c)\n");
        fprintf(stdout,"ggplot(dat, aes(Case, SP))  +  geom_boxplot(aes(colour = Program)) + scale_color_discrete(name = \"Program\")\n");

        return OK;
}


int main(int argc, char *argv[])
{
        /* idea here is call find to get input alignment names
           and produce benchrunner command with varying gpo. gpe. tge
        */
        struct parameters* param = NULL;
        FILE* pipe = NULL;

        char cmd[BUFSIZ];
        char ret[BUFSIZ];
        char out[BUFSIZ];
        char filename[512];
        char runname[512];
        char* hit;
        int rc;

        float gpo;
        float gpe;
        float tgpe;
        float matadd;
        int l;
        int c;
        int i,j;
        int go;

        MMALLOC(param, sizeof(struct parameters));
        param->outfile = NULL;
        param->inputdir = NULL;
        param->scratchdir = NULL;
        while (1){
                static struct option long_options[] ={
                        {"input",  required_argument, 0, 'i'},
                        {"infile",  required_argument, 0, 'i'},
                        {"in",  required_argument, 0, 'i'},
                        {"output",  required_argument, 0, 'o'},
                        {"outfile",  required_argument, 0, 'o'},
                        {"out",  required_argument, 0, 'o'},
                        {"scratch",  required_argument, 0, 's'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;

                c = getopt_long_only (argc, argv,"i:o:s:",long_options, &option_index);

                /* Detect the end of the options. */
                if (c == -1){
                        break;
                }
                switch(c) {
                case 'i':
                        param->inputdir = optarg;
                        break;
                case 's':
                        param->scratchdir = optarg;
                        break;
                case 'o':
                        param->outfile = optarg;
                        break;
                case '?':
                        MFREE(param);
                        exit(1);
                        break;
                default:
                        abort ();
                }
        }


        if(!param->inputdir){
                RUN(help(argv));
                ERROR_MSG("No input dir! use -i <> ");
        }else{
                if(!my_file_exists(param->inputdir)){
                        ERROR_MSG("Input dir: %s not found.",param->inputdir);
                }
        }

        if(!param->scratchdir){
                RUN(help(argv));
                ERROR_MSG("No scratch dir! use -i <> ");
        }else{
                if(!my_file_exists(param->scratchdir)){
                        ERROR_MSG("Scratch dir: %s not found.",param->scratchdir);
                }
        }

        if(!param->outfile){
                RUN(help(argv));
                ERROR_MSG("No output file! use -o <> ");
        }else{
                if(my_file_exists(param->outfile)){
                        WARNING_MSG("output file %s exists - will append results!", param->outfile);
                }
        }

        l = strlen(param->inputdir);
        j = 0;
        for(i = 0; i < l;i++){
                if(param->inputdir[i] == '/'){
                        if(i != l-1){
                                j= 0;
                        }else{
                                runname[j] = 0;

                        }
                }else{
                        runname[j] = param->inputdir[i];
                        j++;
                }
        }
        runname[j] = 0;


        snprintf(cmd, BUFSIZ, "find %s  -name \"*.fa\"", param->inputdir);
        //snprintf(cmd, BUFSIZ, "find ~/kalignbenchmark/data/data-set1  ~/kalignbenchmark/data/data-set2 -name \"*.msf\" | grep structural");
        RUNP(pipe = popen(cmd,"r"));
        //if (system("which gnuplot"))
        c = 1;
        while (fgets(ret, BUFSIZ, pipe)){
                //fprintf(stderr,"%s", ret);
                ret[strcspn(ret, "\r\n")] = 0;
                go = 1;
                snprintf(filename, 512, "%s", ret);

                l = strlen(filename);

                filename[l-3] = '.';
                filename[l-2] = 'a';
                filename[l-1] = 'f';
                filename[l] = 'a';
                filename[l+1] = 0;
                //fprintf(stderr,"%s\n", filename);
                if(!my_file_exists(filename)){
                        //LOG_MSG("Doesn't exits: %s",filename);
                        filename[l-3] = '.';
                        filename[l-2] = 'm';
                        filename[l-1] = 's';
                        filename[l] = 'f';
                        filename[l+1] = 0;
                        if(!my_file_exists(filename)){
                                go = 0;
                                WARNING_MSG("Can't find reference alignment for: %s",ret);
                        }
                }
                if(go){
                        rc = snprintf(out, BUFSIZ, "benchrunner --scratch %s --runname %s -program muscle -test %s -ref %s -o %s -u %d \n", param->scratchdir, runname,ret,filename,param->outfile, c);
                        fprintf(stdout,"%s",out);
                        c++;

                        rc = snprintf(out, BUFSIZ, "benchrunner --scratch %s --runname %s -program kalign -test %s -ref %s -o %s -u %d \n", param->scratchdir, runname,ret, filename,param->outfile, c);
                        fprintf(stdout,"%s",out);
                        c++;

                        rc = snprintf(out, BUFSIZ, "benchrunner --scratch %s @-n:1@ --runname %s -program kalign -test %s -ref %s -o %s -u %d \n", param->scratchdir, runname,ret, filename,param->outfile, c);
                        fprintf(stdout,"%s",out);
                        c++;

                        rc = snprintf(out, BUFSIZ, "benchrunner --scratch %s  @-n:4@  --runname %s -program kalign -test %s -ref %s -o %s -u %d \n", param->scratchdir, runname,ret, filename,param->outfile, c);
                        fprintf(stdout,"%s",out);
                        c++;

                        /* rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @--tgpe:0.5@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c); */
                        /* fprintf(stdout,"%s",out); */
                        /* c++; */

                        /* rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @--tgpe:1.0@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c); */
                        /* fprintf(stdout,"%s",out); */
                        /* c++; */

                        /* rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @--tgpe:1.5@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c); */
                        /* fprintf(stdout,"%s",out); */
                        /* c++; */


                        /* rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @-chaos:2@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c); */
                        /* fprintf(stdout,"%s",out); */
                        /* c++; */

                        /* rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @-chaos:4@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c); */
                        /* fprintf(stdout,"%s",out); */
                        /* c++; */

                        /* rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @-chaos:6@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c); */
                        /* fprintf(stdout,"%s",out); */
                        /* c++; */

                        /* rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @-chaos:10@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c); */
                        /* fprintf(stdout,"%s",out); */
                        /* c++; */
                }else{
                        WARNING_MSG("Can't find reference alignment for: %s",ret);
                }


        }

        rc = pclose(pipe);

        if (rc != EXIT_SUCCESS) { // == 0
                ERROR_MSG("COMMAND FAILED: %s" , cmd);
        }

        snprintf(cmd, BUFSIZ, "%s_plotting.R", param->outfile);

        RUNP(pipe = fopen(cmd, "w"));
        fprintf(pipe,"library(tidyverse)\n");
        fprintf(pipe,"dat = read.csv(\"%s\",header = T)\n",param->outfile);
        fprintf(pipe,"l =  str_split_fixed(dat$Alignment, \"_\",2)\n");
        fprintf(pipe,"c = data.frame(matrix(unlist(l), nrow=length(l), byrow=T))\n");
        fprintf(pipe,"colnames(l) = c(\"Case\",\"Aln\");\n");
        fprintf(pipe,"dat = cbind(dat,l)\n");
        fprintf(pipe,"ggplot(dat, aes(Case, SP))  +  geom_boxplot(aes(colour = Program)) + scale_color_discrete(name = \"Program\")\n");
        snprintf(cmd, BUFSIZ, "%s.pdf", param->outfile);
        fprintf(pipe,"ggsave(\"%s\",p)\n", cmd);

        fclose(pipe);

        if(param){
                MFREE(param);
        }
        return EXIT_SUCCESS;
ERROR:
        if(param){
                MFREE(param);
        }
        return EXIT_FAILURE;
}
