#include "stdio.h"

#include "kalign/kalign.h"

int main(int argc, char *argv[])
{
        struct msa* msa = NULL;

        fprintf(stdout,"inputs: %d\n", argc);
        for(int i = 1; i < argc;i++){
                fprintf(stdout,"reading from %s\n", argv[i]);
                kalign_read_input(argv[i], &msa,1);
        }

        kalign_write_msa(msa, "test.clu", "clu");

        kalign_write_msa(msa, "test.fa", "fasta");

        kalign_free_msa(msa);
        return 0;
}
