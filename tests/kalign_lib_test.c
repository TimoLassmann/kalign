#include <stdio.h>
#include <stdlib.h>

#include "kalign/kalign.h"

int main(int argc, char *argv[])
{
        struct msa* msa = NULL;
        if(argc <= 1){
                fprintf(stdout,"no inputs\n");
                return EXIT_SUCCESS;
        }

        for(int i = 1; i < argc;i++){
                fprintf(stdout,"reading from %s\n", argv[i]);
                kalign_read_input(argv[i], &msa,1);
        }
        /* Align seqences */
        kalign_run(msa,1 , -1, -1, -1 , -1);
        /* write alignment in clustal format */
        kalign_write_msa(msa, "test.clu", "clu");
        /* write alignment in aligned fasta format */
        kalign_write_msa(msa, "test.fa", "fasta");
        /* cleaning up */
        kalign_free_msa(msa);
        return EXIT_SUCCESS;
}
