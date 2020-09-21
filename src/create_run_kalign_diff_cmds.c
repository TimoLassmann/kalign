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
        int rc;

        float gpo;
        float gpe;
        float tgpe;
        float matadd;
        int l;
        int c;

        snprintf(cmd, BUFSIZ, "find ~/kalignbenchmark/data/bb3_release -name \"*.msf\"");
        RUNP(pipe = popen(cmd,"r"));
        //if (system("which gnuplot"))
        c = 1;
        while (fgets(ret, BUFSIZ, pipe)){
                //fprintf(stderr,"%s", ret);
                ret[strcspn(ret, "\r\n")] = 0;

                snprintf(filename, 512, "%s", ret);
                l = strlen(filename);
                filename[l-4] = '.';
                filename[l-3] = 't';
                filename[l-2] = 'f';
                filename[l-1] = 'a';
                //fprintf(stderr,"%s\n", filename);
                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname muscle -program muscle -test %s -ref %s -o kalign_diff_out.csv -u %d \n", filename,ret,c);
                fprintf(stdout,"%s",out);
                c++;

                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname defaultset -program kalign -test %s -ref %s -o kalign_diff_out.csv -u %d --gpo 5.5 -gpe 2.0 -tgpe 1.0 -matadd 0.0 \n", ret,ret,c);
                fprintf(stdout,"%s",out);
                c++;

                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname default -program kalign -test %s -ref %s -o kalign_diff_out.csv -u %d \n", ret,ret,c);
                fprintf(stdout,"%s",out);
                c++;

                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname coarse --gpo 5.545 -gpe 2.013  -tgpe 1.001 -matadd 0.5196 -program kalign -test %s -ref %s -o kalign_diff_out.csv -u %d \n", ret,ret,c);
                fprintf(stdout,"%s",out);
                c++;
                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --runname coarse_round --gpo 5.5 -gpe 2.013  -tgpe 1.0 -matadd 0.5 -program kalign -test %s -ref %s -o kalign_diff_out.csv -u %d \n", ret,ret,c);
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
