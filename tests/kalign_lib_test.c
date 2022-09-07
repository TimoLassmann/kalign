#include "stdio.h"

#include "kalign/kalign.h"

int main(int argc, char *argv[])
{
        struct msa* msa = NULL;

        /* fprintf(stdout,"inputs: %d\n", argc); */
        for(int i = 1; i < argc;i++){
                /* fprintf(stdout,"reading from %s\n", argv[i]); */
                kalign_read_input(argv[i], &msa,1);
        }

        kalign_run(msa,1 , 0, -1, -1 , -1);

        /* kalign_write_msa(msa, "test.clu", "clu"); */

        /* kalign_write_msa(msa, "test.fa", "fasta"); */
        struct msa* ref = NULL;
        for(int i = 1; i < argc;i++){
                /* fprintf(stdout,"reading from %s\n", argv[i]); */
                kalign_read_input(argv[i], &ref,1);
        }


        float score = 0.0;
        kalign_msa_compare(ref, msa, 0 , &score);

        kalign_write_msa(msa, "test.clu", "clu");
        kalign_write_msa(ref, "ref.clu", "clu");

        kalign_free_msa(msa);
        return 0;
}
