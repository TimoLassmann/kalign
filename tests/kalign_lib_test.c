#include <stdio.h>
#include <stdlib.h>

#include "kalign/kalign.h"
int test_2(void);
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
        test_2();
        return EXIT_SUCCESS;
}

int test_2(void){

        // Initialize array
        char* inseq[] = {
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCATGCATGCA",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCATGCATGCA",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCATGCATGCA",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        };

        int numseq = 4;
        int L[] = {76, 76, 76, 48};

        char** aln = NULL;
        int aln_len = 0;

        kalign(inseq, L, numseq, 4 , KALIGN_TYPE_DNA_INTERNAL, -1, -1 , -1, &aln, &aln_len);

        printf("Aligned:\n");
        for(int i = 0; i < numseq; i++){
                printf("%s\n", aln[i]);
                /* printf("%s\n", inseq[i]); */
        }

        for (int i = 0; i < numseq; ++i) {
                int num_aligned = 0;
                for (int j = 0; j < aln_len; ++j) {
                        if (aln[i][j] != '-' && aln[i][j] != '\0') {
                                num_aligned++;
                        }
                }
                if (num_aligned != L[i]) {
                        fprintf(stderr, "error: sequence %d is not fully aligned\n", i);
                        fprintf(stderr, "%s\n", inseq[i]);
                        fprintf(stderr, "%s\n", aln[i]);
                }
        }

        /* Free alignment  */
        for(int i = 0; i < numseq;i++){
                free(aln[i]);
        }
        free(aln);
        return 0;
}
