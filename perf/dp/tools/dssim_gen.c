/* dssim_gen -- write a FASTA file simulated with Kalign's DSSim generator
 * (tests/dssim.c, the same code behind kalign.generate_test_sequences()).
 *
 * usage: dssim_gen <n_seq> <n_obs> <dna:0|1> <model_len> <seed> <out.fa> [max_len]
 *
 * max_len (optional): sequences longer than max_len are truncated to
 * max_len residues (hard limit of the perf datasets: <= 500 residues).
 * Prints: numseq mean_len max_len n_truncated
 *
 * DSSim samples a profile HMM of <model_len> match states: each match/insert
 * emission is built from <n_obs> "observed" residues (95% the column's
 * consensus residue, 5% random) plus a background prior, so n_obs controls
 * divergence: large n_obs -> peaked emissions -> similar sequences; small
 * n_obs -> flat emissions -> divergent sequences.  Indel rate is fixed by
 * DSSim (0.04, or 0.02 when n_seq > 100).
 */
#include <stdio.h>
#include <stdlib.h>

#include "msa_struct.h"
#include "kalign/kalign.h"
#include "dssim.h"

int main(int argc, char** argv)
{
        struct msa* m = NULL;
        FILE* f = NULL;
        int n_seq, n_obs, dna, len, seed, i;
        int max_len = 0, n_trunc = 0, longest = 0;
        long total = 0;
        if(argc != 7 && argc != 8){
                fprintf(stderr, "usage: dssim_gen <n_seq> <n_obs> <dna:0|1> <model_len> <seed> <out.fa> [max_len]\n");
                return 1;
        }
        n_seq = atoi(argv[1]);
        n_obs = atoi(argv[2]);
        dna = atoi(argv[3]);
        len = atoi(argv[4]);
        seed = atoi(argv[5]);
        if(argc == 8) max_len = atoi(argv[7]);
        if(dssim_get_fasta(&m, n_seq, n_obs, dna, len, seed) != 0){
                fprintf(stderr, "dssim_get_fasta failed\n");
                return 1;
        }
        f = fopen(argv[6], "w");
        if(!f){
                fprintf(stderr, "cannot write %s\n", argv[6]);
                return 1;
        }
        for(i = 0; i < m->numseq; i++){
                int l = m->sequences[i]->len;
                if(max_len > 0 && l > max_len){
                        l = max_len;
                        n_trunc++;
                }
                fprintf(f, ">seq%d\n%.*s\n", i + 1, l, m->sequences[i]->seq);
                total += l;
                if(l > longest) longest = l;
        }
        fclose(f);
        printf("%d %.2f %d %d\n", m->numseq, (double)total / (double)m->numseq, longest, n_trunc);
        kalign_free_msa(m);
        return 0;
}
