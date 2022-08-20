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

        int* len = NULL;
        MMALLOC(len, sizeof(int) * 4);

        for(int i = 0 ; i < 4; i++){
                len[i] = strnlen(array[i], 95);
        }
        char** aln = NULL;
        int aln_len = 0;

        RUN(kalign(array,len, 4, &aln, &aln_len));

        MFREE(len);


        return OK;
ERROR:
        return FAIL;
}
