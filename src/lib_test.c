#include "tldevel.h"
#include "libkalign.h"
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


        fprintf(stdout,"Aligned:\n");
        for(int i = 0; i < numseq;i++){
                fprintf(stdout,"%s\n", aln[i]);
        }
        /* Free alignmenr  */
        for(int i = 0; i < numseq;i++){
                free(aln[i]);
        }
        free(aln);





        MFREE(len);


        return OK;
ERROR:
        return FAIL;
}
