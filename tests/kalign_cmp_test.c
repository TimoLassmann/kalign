#include "stdio.h"
#include "stdlib.h"
#include "kalign/kalign.h"

int main(int argc, char *argv[])
{
        struct msa* r = NULL;
        struct msa* t = NULL;

        fprintf(stdout,"inputs: %d\n", argc);
        if(argc == 3){
                kalign_read_input(argv[1], &r,1);
                kalign_read_input(argv[2], &t,1);
                float score = 0.0;

                kalign_msa_compare(r, t, &score);

                kalign_free_msa(r);
                kalign_free_msa(t);
                fprintf(stdout,"Aln score : %f\n", score);
        }
        return 0;
}
