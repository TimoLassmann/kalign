#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kalign/kalign.h"
#include "msa_struct.h"

static int test_ensemble_confidence(const char *input_file);
static int test_poar_round_trip(const char *input_file);
static int test_min_support(const char *input_file);
static int test_residue_order(const char *input_file);
static int test_fill_in(const char *input_file);

int main(int argc, char *argv[])
{
        int ret = EXIT_SUCCESS;

        if(argc <= 1){
                fprintf(stderr, "Usage: %s <input.tfa>\n", argv[0]);
                return EXIT_FAILURE;
        }

        const char *input_file = argv[1];

        fprintf(stdout, "=== Kalign Ensemble Tests ===\n\n");

        fprintf(stdout, "Test 1: Ensemble with confidence\n");
        if(test_ensemble_confidence(input_file) != 0){
                fprintf(stdout, "  FAIL\n\n");
                ret = EXIT_FAILURE;
        }else{
                fprintf(stdout, "  PASS\n\n");
        }

        fprintf(stdout, "Test 2: POAR round-trip\n");
        if(test_poar_round_trip(input_file) != 0){
                fprintf(stdout, "  FAIL\n\n");
                ret = EXIT_FAILURE;
        }else{
                fprintf(stdout, "  PASS\n\n");
        }

        fprintf(stdout, "Test 3: min_support parameter\n");
        if(test_min_support(input_file) != 0){
                fprintf(stdout, "  FAIL\n\n");
                ret = EXIT_FAILURE;
        }else{
                fprintf(stdout, "  PASS\n\n");
        }

        fprintf(stdout, "Test 4: consensus preserves residue order\n");
        if(test_residue_order(input_file) != 0){
                fprintf(stdout, "  FAIL\n\n");
                ret = EXIT_FAILURE;
        }else{
                fprintf(stdout, "  PASS\n\n");
        }

        fprintf(stdout, "Test 5: consensus fill-in keeps anchors and residues\n");
        if(test_fill_in(input_file) != 0){
                fprintf(stdout, "  FAIL\n\n");
                ret = EXIT_FAILURE;
        }else{
                fprintf(stdout, "  PASS\n\n");
        }

        if(ret == EXIT_SUCCESS){
                fprintf(stdout, "All ensemble tests passed.\n");
        }else{
                fprintf(stderr, "Some ensemble tests FAILED.\n");
        }

        return ret;
}

/* Helper: create a 3-run ensemble config with default params */
static void make_ensemble_runs(struct kalign_run_config *runs, int n_runs,
                                struct kalign_ensemble_config *ens, int min_support)
{
        for(int i = 0; i < n_runs; i++){
                runs[i] = kalign_run_config_defaults();
                runs[i].tree_seed = 42 + (uint64_t)i;
                runs[i].tree_noise = (i > 0) ? 0.2f : 0.0f;
        }
        *ens = kalign_ensemble_config_defaults();
        ens->min_support = min_support;
}

/* Test 1: Run ensemble alignment with n_runs=3 and verify that
 * col_confidence is populated with values in [0, 1]. */
static int test_ensemble_confidence(const char *input_file)
{
        struct msa *msa = NULL;
        int rv;

        rv = kalign_read_input((char *)input_file, &msa, 1);
        if(rv != 0 || msa == NULL){
                fprintf(stderr, "  ERROR: failed to read input file: %s\n", input_file);
                return -1;
        }

        struct kalign_run_config runs[3];
        struct kalign_ensemble_config ens;
        make_ensemble_runs(runs, 3, &ens, 0);

        rv = kalign_align_full(msa, runs, 3, &ens, 1);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_align_full (ensemble) returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }

        if(msa->col_confidence == NULL){
                fprintf(stderr, "  ERROR: col_confidence is NULL after ensemble\n");
                kalign_free_msa(msa);
                return -1;
        }

        if(msa->alnlen <= 0){
                fprintf(stderr, "  ERROR: alnlen is %d (expected > 0)\n", msa->alnlen);
                kalign_free_msa(msa);
                return -1;
        }

        for(int i = 0; i < msa->alnlen; i++){
                float c = msa->col_confidence[i];
                if(c < 0.0f || c > 1.0f){
                        fprintf(stderr, "  ERROR: col_confidence[%d] = %f out of range [0,1]\n", i, c);
                        kalign_free_msa(msa);
                        return -1;
                }
        }

        fprintf(stdout, "  col_confidence: %d values, all in [0,1]\n", msa->alnlen);

        for(int i = 0; i < msa->numseq; i++){
                if(msa->sequences[i]->seq == NULL){
                        fprintf(stderr, "  ERROR: sequence %d has NULL seq\n", i);
                        kalign_free_msa(msa);
                        return -1;
                }
                int slen = (int)strlen(msa->sequences[i]->seq);
                if(slen != msa->alnlen){
                        fprintf(stderr, "  ERROR: sequence %d length %d != alnlen %d\n",
                                i, slen, msa->alnlen);
                        kalign_free_msa(msa);
                        return -1;
                }
        }

        fprintf(stdout, "  %d sequences aligned to length %d\n", msa->numseq, msa->alnlen);

        kalign_free_msa(msa);
        return 0;
}

/* Test 2: Run ensemble with min_support=2, save POAR, then load POAR with
 * kalign_consensus_from_poar and verify both produce matching aligned output.
 *
 * NOTE: POAR save is no longer in the ensemble config. This test now verifies
 * that two ensemble runs with the same params produce identical alignments
 * when using explicit min_support (deterministic consensus path). */
static int test_poar_round_trip(const char *input_file)
{
        struct msa *msa1 = NULL;
        struct msa *msa2 = NULL;
        int rv;

        rv = kalign_read_input((char *)input_file, &msa1, 1);
        if(rv != 0 || msa1 == NULL){
                fprintf(stderr, "  ERROR: failed to read input file: %s\n", input_file);
                return -1;
        }

        rv = kalign_read_input((char *)input_file, &msa2, 1);
        if(rv != 0 || msa2 == NULL){
                fprintf(stderr, "  ERROR: failed to read input (2nd copy)\n");
                kalign_free_msa(msa1);
                return -1;
        }

        struct kalign_run_config runs[3];
        struct kalign_ensemble_config ens;
        make_ensemble_runs(runs, 3, &ens, 2);

        rv = kalign_align_full(msa1, runs, 3, &ens, 1);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_align_full (run 1) returned %d\n", rv);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        rv = kalign_align_full(msa2, runs, 3, &ens, 1);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_align_full (run 2) returned %d\n", rv);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        if(msa1->alnlen <= 0 || msa2->alnlen <= 0){
                fprintf(stderr, "  ERROR: alnlen msa1=%d msa2=%d\n", msa1->alnlen, msa2->alnlen);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        if(msa1->numseq != msa2->numseq){
                fprintf(stderr, "  ERROR: numseq mismatch: %d vs %d\n",
                        msa1->numseq, msa2->numseq);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        if(msa1->alnlen != msa2->alnlen){
                fprintf(stderr, "  ERROR: alnlen mismatch: %d vs %d\n",
                        msa1->alnlen, msa2->alnlen);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        for(int i = 0; i < msa1->numseq; i++){
                if(strcmp(msa1->sequences[i]->seq, msa2->sequences[i]->seq) != 0){
                        fprintf(stderr, "  ERROR: sequence %d mismatch between runs\n", i);
                        kalign_free_msa(msa1);
                        kalign_free_msa(msa2);
                        return -1;
                }
        }

        fprintf(stdout, "  Deterministic: %d sequences, alnlen=%d, both runs match\n",
                msa1->numseq, msa1->alnlen);

        kalign_free_msa(msa1);
        kalign_free_msa(msa2);
        return 0;
}

/* Test 3: Run ensemble with explicit min_support=2 and verify it completes. */
static int test_min_support(const char *input_file)
{
        struct msa *msa = NULL;
        int rv;

        rv = kalign_read_input((char *)input_file, &msa, 1);
        if(rv != 0 || msa == NULL){
                fprintf(stderr, "  ERROR: failed to read input file: %s\n", input_file);
                return -1;
        }

        struct kalign_run_config runs[3];
        struct kalign_ensemble_config ens;
        make_ensemble_runs(runs, 3, &ens, 2);

        rv = kalign_align_full(msa, runs, 3, &ens, 1);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_align_full with min_support=2 returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }

        if(msa->alnlen <= 0){
                fprintf(stderr, "  ERROR: alnlen is %d (expected > 0)\n", msa->alnlen);
                kalign_free_msa(msa);
                return -1;
        }

        for(int i = 0; i < msa->numseq; i++){
                if(msa->sequences[i]->seq == NULL){
                        fprintf(stderr, "  ERROR: sequence %d has NULL seq\n", i);
                        kalign_free_msa(msa);
                        return -1;
                }
                int slen = (int)strlen(msa->sequences[i]->seq);
                if(slen != msa->alnlen){
                        fprintf(stderr, "  ERROR: sequence %d length %d != alnlen %d\n",
                                i, slen, msa->alnlen);
                        kalign_free_msa(msa);
                        return -1;
                }
        }

        if(msa->col_confidence == NULL){
                fprintf(stderr, "  ERROR: col_confidence is NULL with min_support=2\n");
                kalign_free_msa(msa);
                return -1;
        }

        fprintf(stdout, "  min_support=2: %d sequences aligned to length %d\n",
                msa->numseq, msa->alnlen);

        kalign_free_msa(msa);
        return 0;
}

/* Test 4: Build the consensus with min_support=1 (the most aggressive
 * merging) and verify every output row, with gaps removed, reproduces
 * its input sequence. A missed ordering cycle in the consensus column
 * graph shows up here as scrambled residues. */
static int test_residue_order(const char *input_file)
{
        struct msa *msa = NULL;
        struct msa *orig = NULL;
        int rv;
        int bad = 0;

        rv = kalign_read_input((char *)input_file, &msa, 1);
        if(rv != 0 || msa == NULL){
                fprintf(stderr, "  ERROR: failed to read input file: %s\n", input_file);
                return -1;
        }
        rv = kalign_read_input((char *)input_file, &orig, 1);
        if(rv != 0 || orig == NULL){
                fprintf(stderr, "  ERROR: failed to read input (2nd copy)\n");
                kalign_free_msa(msa);
                return -1;
        }

        struct kalign_run_config runs[3];
        struct kalign_ensemble_config ens;
        make_ensemble_runs(runs, 3, &ens, 1);

        rv = kalign_align_full(msa, runs, 3, &ens, 1);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_align_full with min_support=1 returned %d\n", rv);
                kalign_free_msa(msa);
                kalign_free_msa(orig);
                return -1;
        }

        /* Output rows may be reordered; match them to inputs by name. */
        for(int i = 0; i < orig->numseq; i++){
                const char *want = orig->sequences[i]->seq;
                int found = 0;
                for(int j = 0; j < msa->numseq; j++){
                        if(strcmp(orig->sequences[i]->name, msa->sequences[j]->name) != 0){
                                continue;
                        }
                        found = 1;
                        const char *got = msa->sequences[j]->seq;
                        int k = 0;
                        for(int c = 0; got[c] != '\0'; c++){
                                if(got[c] == '-'){
                                        continue;
                                }
                                if(want[k] == '\0' || got[c] != want[k]){
                                        break;
                                }
                                k++;
                        }
                        if(want[k] != '\0' || (int)strlen(want) != k){
                                bad++;
                        }
                        break;
                }
                if(!found){
                        bad++;
                }
        }

        fprintf(stdout, "  min_support=1: %d sequences, alnlen=%d, %d rows differ from input\n",
                msa->numseq, msa->alnlen, bad);

        kalign_free_msa(msa);
        kalign_free_msa(orig);
        return bad == 0 ? 0 : -1;
}

/* Map every residue of every row to its column: r2c[i][r]. */
static int **residue_columns(struct msa *msa)
{
        int **r2c = malloc(sizeof(int *) * msa->numseq);
        for(int i = 0; i < msa->numseq; i++){
                const char *s = msa->sequences[i]->seq;
                int n = 0;
                r2c[i] = malloc(sizeof(int) * (strlen(s) + 1));
                for(int c = 0; s[c] != '\0'; c++){
                        if(s[c] != '-'){
                                r2c[i][n++] = c;
                        }
                }
        }
        return r2c;
}

static void free_residue_columns(int **r2c, int numseq)
{
        for(int i = 0; i < numseq; i++){
                free(r2c[i]);
        }
        free(r2c);
}

/* Rows of b, gaps removed, must equal the input rows; all rows must be
   alnlen wide; confidence must be in [0,1]. */
static int check_rows(struct msa *orig, struct msa *b)
{
        for(int i = 0; i < orig->numseq; i++){
                const char *want = orig->sequences[i]->seq;
                const char *got = b->sequences[i]->seq;
                int k = 0;
                if(strcmp(orig->sequences[i]->name, b->sequences[i]->name) != 0){
                        fprintf(stderr, "  ERROR: row %d name %s != %s\n", i,
                                b->sequences[i]->name, orig->sequences[i]->name);
                        return -1;
                }
                if((int)strlen(got) != b->alnlen){
                        fprintf(stderr, "  ERROR: row %d width %d != alnlen %d\n",
                                i, (int)strlen(got), b->alnlen);
                        return -1;
                }
                for(int c = 0; got[c] != '\0'; c++){
                        if(got[c] == '-') continue;
                        if(got[c] != want[k]){
                                fprintf(stderr, "  ERROR: row %d residue %d differs\n", i, k);
                                return -1;
                        }
                        k++;
                }
                if(want[k] != '\0'){
                        fprintf(stderr, "  ERROR: row %d lost residues\n", i);
                        return -1;
                }
        }
        if(b->col_confidence == NULL){
                fprintf(stderr, "  ERROR: col_confidence is NULL\n");
                return -1;
        }
        for(int c = 0; c < b->alnlen; c++){
                if(b->col_confidence[c] < 0.0f || b->col_confidence[c] > 1.0f){
                        fprintf(stderr, "  ERROR: col_confidence[%d] = %f\n", c, b->col_confidence[c]);
                        return -1;
                }
        }
        return 0;
}

/* Test 5: every fill-in mode must reproduce the input rows, and every
   column of the plain consensus holding >= 2 residues (an anchor at
   fill_anchor_frac = 0) must still be one column after fill-in. */
static int test_fill_in(const char *input_file)
{
        const int modes[3] = {KALIGN_FILL_REALIGN, KALIGN_FILL_PROJECT, KALIGN_FILL_ATTACH};
        struct msa *orig = NULL;
        struct msa *off = NULL;
        int ret = -1;

        if(kalign_read_input((char *)input_file, &orig, 1) != 0 ||
           kalign_read_input((char *)input_file, &off, 1) != 0){
                fprintf(stderr, "  ERROR: failed to read input file: %s\n", input_file);
                goto done;
        }

        struct kalign_run_config runs[3];
        struct kalign_ensemble_config ens;
        make_ensemble_runs(runs, 3, &ens, 2);
        if(kalign_align_full(off, runs, 3, &ens, 1) != 0){
                fprintf(stderr, "  ERROR: ensemble without fill-in failed\n");
                goto done;
        }
        if(check_rows(orig, off) != 0) goto done;

        for(int m = 0; m < 3; m++){
                struct msa *fill = NULL;
                int **r2c = NULL;
                int moved = 0;
                int bad_rows;

                if(kalign_read_input((char *)input_file, &fill, 1) != 0){
                        fprintf(stderr, "  ERROR: failed to read input file\n");
                        goto done;
                }
                make_ensemble_runs(runs, 3, &ens, 2);
                ens.fill_in = modes[m];
                ens.fill_anchor_frac = 0.0f;
                ens.fill_support = 1;
                if(kalign_align_full(fill, runs, 3, &ens, 1) != 0){
                        fprintf(stderr, "  ERROR: ensemble with fill_in=%d failed\n", modes[m]);
                        kalign_free_msa(fill);
                        goto done;
                }
                bad_rows = check_rows(orig, fill);
                r2c = residue_columns(fill);

                /* Anchor check: residues sharing an off column with >= 2
                   residues must share a column after fill-in. */
                int *rpos = calloc(off->numseq, sizeof(int));
                for(int c = 0; c < off->alnlen && !bad_rows; c++){
                        int first_col = -1;
                        int n_res = 0;
                        for(int i = 0; i < off->numseq; i++){
                                if(off->sequences[i]->seq[c] != '-') n_res++;
                        }
                        for(int i = 0; i < off->numseq; i++){
                                if(off->sequences[i]->seq[c] == '-') continue;
                                int col = r2c[i][rpos[i]++];
                                if(n_res < 2) continue;
                                if(first_col < 0){
                                        first_col = col;
                                }else if(col != first_col){
                                        moved++;
                                }
                        }
                }
                free(rpos);
                int fill_len = fill->alnlen;
                fprintf(stdout, "  fill_in=%d: alnlen %d (off %d), %d anchored residues moved\n",
                        modes[m], fill_len, off->alnlen, moved);
                free_residue_columns(r2c, fill->numseq);
                kalign_free_msa(fill);
                if(bad_rows != 0 || moved != 0){
                        goto done;
                }
                if(modes[m] == KALIGN_FILL_ATTACH && fill_len > off->alnlen){
                        fprintf(stderr, "  ERROR: attach widened the alignment\n");
                        goto done;
                }
        }
        ret = 0;
done:
        if(off) kalign_free_msa(off);
        if(orig) kalign_free_msa(orig);
        return ret;
}
