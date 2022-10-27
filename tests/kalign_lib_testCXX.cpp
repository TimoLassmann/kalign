
#include <array>
#include <iostream>
#include <string>
#include "kalign/kalign.h"

int main() {
        // Initialize array
        char * inseq[95] = {
                "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE",
                "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK",
                "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK",
                "MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK",
        };
        int numseq = 4;

        int* L = new int[numseq];        
        L[0] = 83;
        L[1] = 85;
        L[2] = 91;
        L[3] = 86;
        
        char** aln = NULL;
        int aln_len = 0;
                
        kalign(inseq,L, numseq,4 , KALIGN_TYPE_PROTEIN, -1, -1 , -1, &aln, &aln_len);
        
        std::cout << "Aligned:\n";
        for(int i = 0; i < numseq;i++){
                std::cout << aln[i] << "\n";
        
        }
        /* Free alignment  */
        for(int i = 0; i < numseq;i++){
                // delete[] aln[i];
                free(aln[i]);
        }
        free(aln);        
        delete[] L;
}
