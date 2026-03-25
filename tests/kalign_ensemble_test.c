#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kalign/kalign.h"
#include "msa_struct.h"

static int test_ensemble_confidence(const char *input_file);
static int test_poar_round_trip(const char *input_file);
static int test_min_support(const char *input_file);

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
