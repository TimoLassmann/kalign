#include "kalign/kalign.h"

#include "msa_struct.h"
#include "dssim.h"
#include "tldevel.h"

int main(int argc, char *argv[])
{
        struct msa* m = NULL;

        dssim_get_fasta(&m, 20, 10, 250, 0);
        /* for(int i= 0; i < m->numseq;i++){ */
        /*         fprintf(stdout,">Seq_%d\n%s\n",i, m->sequences[i]->seq); */
        /* } */

        kalign_run(m, 1, KALIGN_TYPE_PROTEIN,0.0,0.0,0.0);

        RUN(kalign_write_msa(m, NULL, "msf"));
        kalign_free_msa(m);

        return EXIT_SUCCESS;
ERROR:
        kalign_free_msa(m);
        return EXIT_FAILURE;
}
