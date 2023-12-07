#include "tldevel.h"
#include "tlrng.h"

#include "kalign/kalign.h"

#include "msa_struct.h"
#include "msa_sort.h"
#include "dssim.h"


int main(int argc, char *argv[])
{
        struct msa* m = NULL;
        struct rng_state* rng = NULL;

        rng = init_rng(0);

        dssim_get_fasta(&m, 5, 10, 45, 0);
        /* for(int i= 0; i < m->numseq;i++){ */
        /*         fprintf(stdout,">Seq_%d\n%s\n",i, m->sequences[i]->seq); */
        /* } */

        kalign_run(m, 1, KALIGN_TYPE_PROTEIN,0.0,0.0,0.0);

        RUN(kalign_write_msa(m, NULL, "msf"));

        msa_shuffle_seq(m, rng);
        RUN(kalign_write_msa(m, NULL, "msf"));

        kalign_free_msa(m);

        free_rng(rng);
        return EXIT_SUCCESS;
ERROR:
        kalign_free_msa(m);
        return EXIT_FAILURE;
}
