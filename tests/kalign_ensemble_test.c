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

        /* n_runs=3, default gap penalties, seed=42, min_support=0 (auto), no POAR save */
        rv = kalign_ensemble(msa, 1, -1, 3, -1.0f, -1.0f, -1.0f, 42, 0, NULL, 0, 0.0f, -1.0f, 0, 0.0f, 0, 2.0f);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_ensemble returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }

        /* Check that col_confidence was allocated */
        if(msa->col_confidence == NULL){
                fprintf(stderr, "  ERROR: col_confidence is NULL after ensemble\n");
                kalign_free_msa(msa);
                return -1;
        }

        /* Verify alignment length is positive */
        if(msa->alnlen <= 0){
                fprintf(stderr, "  ERROR: alnlen is %d (expected > 0)\n", msa->alnlen);
                kalign_free_msa(msa);
                return -1;
        }

        /* Check that all col_confidence values are in [0, 1] */
        for(int i = 0; i < msa->alnlen; i++){
                float c = msa->col_confidence[i];
                if(c < 0.0f || c > 1.0f){
                        fprintf(stderr, "  ERROR: col_confidence[%d] = %f out of range [0,1]\n", i, c);
                        kalign_free_msa(msa);
                        return -1;
                }
        }

        fprintf(stdout, "  col_confidence: %d values, all in [0,1]\n", msa->alnlen);

        /* Also verify that sequences are aligned (have equal length = alnlen) */
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

/* Test 2: Run ensemble with save_poar, then load POAR with
 * kalign_consensus_from_poar and verify both produce aligned output. */
static int test_poar_round_trip(const char *input_file)
{
        struct msa *msa1 = NULL;
        struct msa *msa2 = NULL;
        int rv;
        const char *poar_path = "test_ensemble_poar.bin";

        /* First run: ensemble with POAR save */
        rv = kalign_read_input((char *)input_file, &msa1, 1);
        if(rv != 0 || msa1 == NULL){
                fprintf(stderr, "  ERROR: failed to read input file: %s\n", input_file);
                return -1;
        }

        /* Use explicit min_support=2 so the ensemble always takes the consensus
         * path.  kalign_consensus_from_poar() also requires min_support >= 1,
         * and both must use the same threshold for the output to match. */
        rv = kalign_ensemble(msa1, 1, -1, 3, -1.0f, -1.0f, -1.0f, 42, 2, poar_path, 0, 0.0f, -1.0f, 0, 0.0f, 0, 2.0f);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_ensemble (save) returned %d\n", rv);
                kalign_free_msa(msa1);
                return -1;
        }

        /* Verify the POAR file was created */
        FILE *fp = fopen(poar_path, "rb");
        if(fp == NULL){
                fprintf(stderr, "  ERROR: POAR file was not created: %s\n", poar_path);
                kalign_free_msa(msa1);
                return -1;
        }
        fclose(fp);

        fprintf(stdout, "  POAR saved to %s\n", poar_path);

        /* Second run: load POAR and derive consensus */
        rv = kalign_read_input((char *)input_file, &msa2, 1);
        if(rv != 0 || msa2 == NULL){
                fprintf(stderr, "  ERROR: failed to read input for POAR load\n");
                kalign_free_msa(msa1);
                return -1;
        }

        rv = kalign_consensus_from_poar(msa2, poar_path, 2);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_consensus_from_poar returned %d\n", rv);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        /* Verify both MSAs have valid aligned sequences */
        if(msa1->alnlen <= 0){
                fprintf(stderr, "  ERROR: msa1 alnlen = %d\n", msa1->alnlen);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        if(msa2->alnlen <= 0){
                fprintf(stderr, "  ERROR: msa2 alnlen = %d\n", msa2->alnlen);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        /* Both should have the same number of sequences */
        if(msa1->numseq != msa2->numseq){
                fprintf(stderr, "  ERROR: numseq mismatch: %d vs %d\n",
                        msa1->numseq, msa2->numseq);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        /* Both should have the same alignment length (same POAR, same consensus) */
        if(msa1->alnlen != msa2->alnlen){
                fprintf(stderr, "  ERROR: alnlen mismatch: %d vs %d\n",
                        msa1->alnlen, msa2->alnlen);
                kalign_free_msa(msa1);
                kalign_free_msa(msa2);
                return -1;
        }

        /* Verify aligned sequences match between direct ensemble and POAR-loaded consensus */
        for(int i = 0; i < msa1->numseq; i++){
                if(msa1->sequences[i]->seq == NULL || msa2->sequences[i]->seq == NULL){
                        fprintf(stderr, "  ERROR: NULL seq at index %d\n", i);
                        kalign_free_msa(msa1);
                        kalign_free_msa(msa2);
                        return -1;
                }
                if(strcmp(msa1->sequences[i]->seq, msa2->sequences[i]->seq) != 0){
                        fprintf(stderr, "  ERROR: sequence %d mismatch between ensemble and POAR load\n", i);
                        fprintf(stderr, "    direct:  %.60s...\n", msa1->sequences[i]->seq);
                        fprintf(stderr, "    loaded:  %.60s...\n", msa2->sequences[i]->seq);
                        kalign_free_msa(msa1);
                        kalign_free_msa(msa2);
                        return -1;
                }
        }

        fprintf(stdout, "  Round-trip: %d sequences, alnlen=%d, consensus matches\n",
                msa1->numseq, msa1->alnlen);

        kalign_free_msa(msa1);
        kalign_free_msa(msa2);

        /* Clean up temp file */
        remove(poar_path);

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

        /* min_support=2 (explicit), n_runs=3 */
        rv = kalign_ensemble(msa, 1, -1, 3, -1.0f, -1.0f, -1.0f, 42, 2, NULL, 0, 0.0f, -1.0f, 0, 0.0f, 0, 2.0f);
        if(rv != 0){
                fprintf(stderr, "  ERROR: kalign_ensemble with min_support=2 returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }

        /* Verify alignment was produced */
        if(msa->alnlen <= 0){
                fprintf(stderr, "  ERROR: alnlen is %d (expected > 0)\n", msa->alnlen);
                kalign_free_msa(msa);
                return -1;
        }

        /* Verify sequences are aligned */
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

        /* Verify col_confidence is populated */
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
