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

        int rc;

        float gpo_start = 5.0F;
        float gpo_step = 0.1F;

        float gpe_start = 1.5F;
        float gpe_step = 0.1F;

        float tgpe_start = 0.5F;
        float tgpe_step = 0.1F;

        float matadd_start = 0.0F;
        float matadd_step = 0.1F;

        float gpo;
        float gpe;
        float tgpe;
        float matadd;

        int i,j,c,f;

        snprintf(cmd, BUFSIZ, "find ~/kalignbenchmark/data/bb3_release -name \"*.msf\"");
        RUNP(pipe = popen(cmd,"r"));
        //if (system("which gnuplot"))
        while (fgets(ret, BUFSIZ, pipe)){
                //fprintf(stderr,"%s", ret);
                ret[strcspn(ret, "\r\n")] = 0;

                gpo = gpo_start;
                for(i = 0;i < 15;i++){
                        gpe = gpe_start;
                        for(j = 0;j < 15;j++){
                                tgpe = tgpe_start;
                                for(c = 0;c <15;c++){
                                        matadd = matadd_start;
                                        for(f = 0; f < 15; f++){
                                                rc = snprintf(out, BUFSIZ, "benchrunner --scratch ~/tmp --gpo %f -gpe %f -tgpe %f -matadd %f -program kalign -test %s -ref %s -o big_run_out.csv -u %d \n", gpo,gpe,tgpe,matadd, ret,ret,  i << 12| j << 8 | c << 4 | f);
                                                fprintf(stdout,"%s",out);
                                                matadd += matadd_step;
                                        }
                                        tgpe+= tgpe_step;
                                }
                                gpe += gpe_step;
                        }
                        gpo += gpo_step;
                }



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
