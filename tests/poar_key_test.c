/* POAR keys for sequences longer than 4096 residues.
 *
 * Kalign up to 4.0.0 packed (pos_i << 20 | pos_j) into 32 bits, so residues
 * x and x + 4096 of the lower sequence shared a key and the support of two
 * different residue pairs was merged.  Here run 0 aligns residue 10 of
 * sequence 0 with sequence 1, and run 1 aligns residue 4106 with it: each
 * pair must carry only its own run's bit, in memory and after a file round
 * trip. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tldevel.h"
#include "poar.h"

#define L0 5000

static int check(const struct poar_table* t, const char* what)
{
        const struct poar_entry* e = NULL;
        size_t n = 0;
        struct poar_scratch sc;
        uint32_t s10 = poar_support(t, 0, 10, 1, 0);
        uint32_t s4106 = poar_support(t, 0, 4106, 1, 0);

        if(s10 != 1u || s4106 != 2u){
                fprintf(stderr, "%s: support(10,0)=%u support(4106,0)=%u, expected 1 and 2\n",
                        what, s10, s4106);
                return 1;
        }
        poar_scratch_init(&sc);
        if(poar_pair_entries(t, 0, 1, &sc, &e, &n) != OK || n != 2 ||
           e[0].key != POAR_KEY(10, 0) || e[0].support != 1u ||
           e[1].key != POAR_KEY(4106, 0) || e[1].support != 2u){
                fprintf(stderr, "%s: pair entries wrong (n=%zu)\n", what, n);
                poar_scratch_free(&sc);
                return 1;
        }
        poar_scratch_free(&sc);
        return 0;
}

int main(int argc, char* argv[])
{
        char* seqs[2];
        struct poar_table* t = NULL;
        struct poar_table* r = NULL;
        struct pos_matrix* pm = NULL;
        const char* path = argc > 1 ? argv[1] : "poar_key_test.poar";
        int fail = 0;

        seqs[0] = malloc(L0 + 1);
        seqs[1] = malloc(L0 + 1);
        memset(seqs[0], 'A', L0);
        seqs[0][L0] = 0;

        if(poar_table_alloc(&t, 2) != OK) return 1;
        for(int run = 0; run < 2; run++){
                memset(seqs[1], '-', L0);
                seqs[1][L0] = 0;
                seqs[1][run == 0 ? 10 : 4106] = 'A';
                if(pos_matrix_from_msa(&pm, seqs, 2, L0) != OK) return 1;
                if(extract_poars(t, pm, run) != OK) return 1;
                pos_matrix_free(pm);
                pm = NULL;
        }

        fail |= check(t, "implicit table");
        if(poar_table_write(t, path) != OK || poar_table_read(&r, path) != OK){
                fprintf(stderr, "POAR file round trip failed\n");
                fail = 1;
        }else{
                fail |= check(r, "stored table");
        }
        remove(path);

        poar_table_free(t);
        poar_table_free(r);
        free(seqs[0]);
        free(seqs[1]);
        fprintf(stdout, "poar_key_test: %s\n", fail ? "FAIL" : "PASS");
        return fail;
}
