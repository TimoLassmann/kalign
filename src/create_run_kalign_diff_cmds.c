#include "tldevel.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
        /* idea here is call find to get input alignment names
           and produce benchrunner command with varying gpo. gpe. tge
        */
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
        snprintf(cmd, BUFSIZ, "find ~/kalignbenchmark/data/bb3_release -name \"*.msf\"");
        RUNP(pipe = popen(cmd,"r"));
        //if (system("which gnuplot"))
        c = 1;
        while (fgets(ret, BUFSIZ, pipe)){
                //fprintf(stderr,"%s", ret);
                ret[strcspn(ret, "\r\n")] = 0;

                hit = strstr(ret, "bb3_release/RV");
                if(hit){
                        hit += 12;
                        l = strlen(hit);
                        for(i = 0; i < l;i++){
                                if(hit[i] == '/'){
                                        runname[i] = 0;
                                        break;
                                }
                                runname[i] = hit[i];

                        }
                }else{
                        snprintf(runname, 512, "%s", "default");
                }
                snprintf(filename, 512, "%s", ret);
                l = strlen(filename);
                filename[l-4] = '.';
                filename[l-3] = 't';
                filename[l-2] = 'f';
                filename[l-1] = 'a';
                //fprintf(stderr,"%s\n", filename);
                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program muscle -test %s -ref %s -o kalign_diff_out.csv -u %d \n", runname,filename,ret,c);
                fprintf(stdout,"%s",out);
                c++;

                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c);
                fprintf(stdout,"%s",out);
                c++;

                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @-chaos:2@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c);
                fprintf(stdout,"%s",out);
                c++;

                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @-chaos:4@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c);
                fprintf(stdout,"%s",out);
                c++;

                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname %s -program kalign @-chaos:6@ -test %s -ref %s -o kalign_diff_out.csv -u %d \n",runname, filename,ret,c);
                fprintf(stdout,"%s",out);
                c++;

        }

        rc = pclose(pipe);

        if (rc == EXIT_SUCCESS) { // == 0
        }

/* dat = read.csv("big_run_out.csv",header = T) */
/* dat %>% group_by(Alignment) %>% filter(SP == max(SP)) %>%  print(n=500)  */
/* dat %>% group_by(Alignment) %>% filter(SP == max(SP)) %>% summarise_at(vars(SP,GPO,GPE,TGPE), c(max,mean)) %>%  print(n=500) */

/* summary(dat %>% group_by(Alignment) %>% filter(GPO == 5.5) %>% filter(GPE == 2.0) %>% filter(TGPE == 1.0)) */
/* summary(dat %>% group_by(Alignment) %>% filter(SP == max(SP)) %>% summarise_at(vars(SP,GPO,GPE,TGPE), c(max,mean))) */
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
