
#include "libkalign.h"
#include "tldevel.h"
#include <string.h>

int main(void)
{
        char * array[95] = {
                "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE",
                "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK",
                "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK",
                "MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK",
        };
        int numseq = 4;
        int* len = NULL;
        MMALLOC(len, sizeof(int) * numseq);

        for(int i = 0 ; i < numseq; i++){
                len[i] = strnlen(array[i], 95);
        }
        char** aln = NULL;
        int aln_len = 0;
        /* Call kalign  */
        RUN(kalign(array,len, numseq, &aln, &aln_len));

        /* Print out alignment */
        fprintf(stdout,"Aligned:\n");
        for(int i = 0; i < numseq;i++){
                fprintf(stdout,"%s\n", aln[i]);
        }
        /* Free alignment  */
        for(int i = 0; i < numseq;i++){
                free(aln[i]);
        }
        free(aln);

        free(len);

        return OK;
ERROR:
        return FAIL;
}