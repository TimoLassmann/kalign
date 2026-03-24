/* memcheck_stress.c — Comprehensive memory stress test for kalign.
 *
 * Exercises all major code paths repeatedly to surface memory bugs:
 * - In-memory alignment (kalign_arr_to_msa path)
 * - File-based alignment (kalign_read_input path)
 * - Realignment iterations
 * - Refinement (confident, inline)
 * - Ensemble alignment with consensus
 * - MSA comparison (simple + detailed + with mask)
 * - Consistency anchors
 * - VSM + seq_weights
 * - Align + write + read-back + compare (full benchmark loop)
 *
 * Compile with ASAN:
 *   cc -fsanitize=address -O0 -g -DDEBUG \
 *     -I../lib/include -I../lib/src \
 *     memcheck_stress.c \
 *     -L../build-asan/lib -lkalign_static -ltldevel \
 *     -fopenmp -lm -o memcheck_stress
 *
 * Compile for Valgrind:
 *   cc -O0 -g -DDEBUG \
 *     -I../lib/include -I../lib/src \
 *     memcheck_stress.c \
 *     -L../build-debug/lib -lkalign_static -ltldevel \
 *     -fopenmp -lm -o memcheck_stress
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <kalign/kalign.h>
#include <kalign/kalign_config.h>

#include "msa_struct.h"
#include "msa_alloc.h"
#include "msa_op.h"
#include "msa_cmp.h"
#include "msa_io.h"

static int n_passed = 0;
static int n_failed = 0;

#define RUN_TEST(fn, ...) do { \
    fprintf(stderr, "--- %s ---\n", #fn); \
    if (fn(__VA_ARGS__) == 0) { n_passed++; fprintf(stderr, "  PASSED\n"); } \
    else { n_failed++; fprintf(stderr, "  *** FAILED ***\n"); } \
} while(0)

/* ======== Test helpers ======== */

static char* test_protein_seqs[] = {
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAATGIQIRGIVKWFNRRKEMISAYDLLAK",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAATGIQIRGIVKWFNRRKEMISAYDLLAK",
    "MKQAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAATGIQIRGIVKWFNRRKEMISAYDLLAK",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPAAQFEVVHSLAKWKRQQIAATGIQIRGIVKWFNRRKEMISAYDLLAK",
    "MKTVYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAATGIQIRGIVKWFNRRKEMISAYDLLAK",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAATGIQIRGIV",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQD",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAATGIQIRGIVKWFNRRKEMISAYDLLAKMKTAY",
};
static int n_protein_seqs = 8;

static int get_lens(char** seqs, int n, int* lens)
{
    for (int i = 0; i < n; i++)
        lens[i] = (int)strlen(seqs[i]);
    return 0;
}

/* ======== Tests ======== */

/* Test 1: Repeated in-memory alignment */
static int test_inmem_align(int n)
{
    int lens[8];
    get_lens(test_protein_seqs, n_protein_seqs, lens);

    for (int i = 0; i < n; i++) {
        char** aligned = NULL;
        int aln_len = 0;
        int ret = kalign(test_protein_seqs, lens, n_protein_seqs, 1,
                         KALIGN_TYPE_UNDEFINED, -1.0f, -1.0f, -1.0f,
                         &aligned, &aln_len);
        if (ret != 0) { fprintf(stderr, "  failed iter %d\n", i); return 1; }
        for (int j = 0; j < n_protein_seqs; j++) free(aligned[j]);
        free(aligned);
    }
    return 0;
}

/* Test 2: Repeated file-based alignment */
static int test_file_align(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) { fprintf(stderr, "  read failed iter %d\n", i); return 1; }
        msa->quiet = 1;
        ret = kalign_run(msa, 1, KALIGN_TYPE_UNDEFINED, -1.0f, -1.0f, -1.0f, KALIGN_REFINE_NONE, 0);
        if (ret != 0) { fprintf(stderr, "  align failed iter %d\n", i); kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 3: Repeated alignment with realign */
static int test_realign(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_run_realign(msa, 1, KALIGN_TYPE_UNDEFINED,
                                 -1.0f, -1.0f, -1.0f,
                                 KALIGN_REFINE_NONE, 0,
                                 0.0f, -1.0f, 1, -1.0f, 0, 2.0f);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 4: Repeated alignment with refinement */
static int test_refine(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_run(msa, 1, KALIGN_TYPE_UNDEFINED, -1.0f, -1.0f, -1.0f,
                         KALIGN_REFINE_CONFIDENT, 0);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 5: Repeated MSA comparison (simple) */
static int test_compare(const char* ref_file, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* ref = NULL;
        struct msa* test = NULL;
        float score = 0.0f;
        int ret = kalign_read_input((char*)ref_file, &ref, 1);
        if (ret != 0 || !ref) return 1;
        ret = kalign_read_input((char*)ref_file, &test, 1);
        if (ret != 0 || !test) { kalign_free_msa(ref); return 1; }
        ret = kalign_msa_compare(ref, test, &score);
        kalign_free_msa(ref);
        kalign_free_msa(test);
        if (ret != 0) return 1;
    }
    return 0;
}

/* Test 6: Repeated detailed MSA comparison */
static int test_compare_detailed(const char* ref_file, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* ref = NULL;
        struct msa* test = NULL;
        struct poar_score score;
        int ret = kalign_read_input((char*)ref_file, &ref, 1);
        if (ret != 0 || !ref) return 1;
        ret = kalign_read_input((char*)ref_file, &test, 1);
        if (ret != 0 || !test) { kalign_free_msa(ref); return 1; }
        ret = kalign_msa_compare_detailed(ref, test, 0.2f, &score);
        kalign_free_msa(ref);
        kalign_free_msa(test);
        if (ret != 0) return 1;
    }
    return 0;
}

/* Test 7: Repeated ensemble alignment */
static int test_ensemble(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_ensemble(msa, 1, KALIGN_TYPE_UNDEFINED,
                              3, -1.0f, -1.0f, -1.0f,
                              42, 0, NULL,
                              KALIGN_REFINE_NONE, 0.0f, -1.0f,
                              0, -1.0f, 0, 2.0f);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 8: Full benchmark loop: align + write + read-back + compare */
static int test_benchmark_loop(const char* input, const char* ref_file, int n)
{
    char tmpfile[] = "/tmp/kalign_memcheck_output.fa";
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_run(msa, 1, KALIGN_TYPE_UNDEFINED, -1.0f, -1.0f, -1.0f,
                         KALIGN_REFINE_NONE, 0);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        ret = kalign_write_msa(msa, tmpfile, "fasta");
        kalign_free_msa(msa);
        if (ret != 0) return 1;

        /* Compare */
        struct msa* ref = NULL;
        struct msa* test = NULL;
        float sp = 0.0f;
        struct poar_score detailed;
        ret = kalign_read_input((char*)ref_file, &ref, 1);
        if (ret != 0 || !ref) return 1;
        ret = kalign_read_input(tmpfile, &test, 1);
        if (ret != 0 || !test) { kalign_free_msa(ref); return 1; }
        ret = kalign_msa_compare(ref, test, &sp);
        kalign_free_msa(ref); ref = NULL;
        kalign_free_msa(test); test = NULL;
        if (ret != 0) return 1;

        /* Detailed compare */
        ret = kalign_read_input((char*)ref_file, &ref, 1);
        if (ret != 0 || !ref) return 1;
        ret = kalign_read_input(tmpfile, &test, 1);
        if (ret != 0 || !test) { kalign_free_msa(ref); return 1; }
        ret = kalign_msa_compare_detailed(ref, test, 0.2f, &detailed);
        kalign_free_msa(ref); ref = NULL;
        kalign_free_msa(test); test = NULL;
        if (ret != 0) return 1;
    }
    remove(tmpfile);
    return 0;
}

/* Test 9: Consistency anchors */
static int test_consistency(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_run_seeded(msa, 1, KALIGN_TYPE_UNDEFINED,
                                -1.0f, -1.0f, -1.0f,
                                KALIGN_REFINE_NONE, 0,
                                0, 0.0f, 0.0f, -1.0f, -1.0f,
                                3, 2.0f);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 10: Ensemble + realign */
static int test_ensemble_realign(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_ensemble(msa, 1, KALIGN_TYPE_UNDEFINED,
                              3, -1.0f, -1.0f, -1.0f,
                              42, 0, NULL,
                              KALIGN_REFINE_CONFIDENT, 0.0f, -1.0f,
                              1, -1.0f, 0, 2.0f);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 11: Ensemble + VSM + seq_weights */
static int test_ensemble_vsm_sw(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_ensemble(msa, 1, KALIGN_TYPE_UNDEFINED,
                              3, -1.0f, -1.0f, -1.0f,
                              42, 0, NULL,
                              KALIGN_REFINE_CONFIDENT, 0.0f, 2.0f,
                              1, 1.0f, 0, 2.0f);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 12: Inline refinement */
static int test_inline_refine(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_run(msa, 1, KALIGN_TYPE_UNDEFINED, -1.0f, -1.0f, -1.0f,
                         KALIGN_REFINE_INLINE, 0);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 13: Mixed parameter variations (like optimizer does) */
static int test_param_sweep(const char* input, int n)
{
    float vsm_vals[] = {0.0f, 1.0f, 2.0f, 3.0f};
    float sw_vals[] = {0.0f, 1.0f};
    int cons_vals[] = {0, 3, 5};
    int idx = 0;

    for (int i = 0; i < n; i++) {
        float vsm = vsm_vals[idx % 4];
        float sw = sw_vals[(idx / 4) % 2];
        int cons = cons_vals[(idx / 8) % 3];
        idx++;

        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;
        ret = kalign_run_seeded(msa, 1, KALIGN_TYPE_UNDEFINED,
                                -1.0f, -1.0f, -1.0f,
                                KALIGN_REFINE_NONE, 0,
                                0, 0.0f, 0.0f, vsm, sw,
                                cons, 2.0f);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 14: kalign_align_full with single run configs */
static int test_align_full_single(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;

        struct kalign_run_config cfg = kalign_run_config_defaults();
        cfg.vsm_amax = 2.0f;
        cfg.consistency_anchors = 3;
        cfg.realign = 1;
        cfg.refine = KALIGN_REFINE_CONFIDENT;

        ret = kalign_align_full(msa, &cfg, 1, NULL, 1);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

/* Test 15: kalign_align_full with multi-run configs (ensemble) */
static int test_align_full_multi(const char* input, int n)
{
    for (int i = 0; i < n; i++) {
        struct msa* msa = NULL;
        int ret = kalign_read_input((char*)input, &msa, 1);
        if (ret != 0 || !msa) return 1;
        msa->quiet = 1;

        struct kalign_run_config runs[3];
        for (int k = 0; k < 3; k++) {
            runs[k] = kalign_run_config_defaults();
            runs[k].vsm_amax = 2.0f;
            runs[k].realign = 1;
            runs[k].refine = KALIGN_REFINE_CONFIDENT;
        }
        runs[1].gpo = 3.0f; runs[1].gpe = 1.5f; runs[1].tgpe = 0.5f;
        runs[1].tree_seed = 43; runs[1].tree_noise = 0.2f;
        runs[2].gpo = 8.0f; runs[2].gpe = 3.0f; runs[2].tgpe = 1.5f;
        runs[2].tree_seed = 44; runs[2].tree_noise = 0.3f;

        struct kalign_ensemble_config ens = kalign_ensemble_config_defaults();

        ret = kalign_align_full(msa, runs, 3, &ens, 1);
        if (ret != 0) { kalign_free_msa(msa); return 1; }
        kalign_free_msa(msa);
    }
    return 0;
}

int main(int argc, char* argv[])
{
    int n = 10;
    const char* input_file = NULL;
    const char* ref_file = NULL;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <unaligned.tfa> <reference.msf> [n_iters]\n", argv[0]);
        return 1;
    }
    input_file = argv[1];
    ref_file = argv[2];
    if (argc > 3) n = atoi(argv[3]);

    fprintf(stderr, "=== Kalign Memory Stress Test (%d iterations per test) ===\n\n", n);

    RUN_TEST(test_inmem_align, n);
    RUN_TEST(test_file_align, input_file, n);
    RUN_TEST(test_realign, input_file, n);
    RUN_TEST(test_refine, input_file, n);
    RUN_TEST(test_compare, ref_file, n);
    RUN_TEST(test_compare_detailed, ref_file, n);
    RUN_TEST(test_ensemble, input_file, n);
    RUN_TEST(test_benchmark_loop, input_file, ref_file, n);
    RUN_TEST(test_consistency, input_file, n);
    RUN_TEST(test_ensemble_realign, input_file, n);
    RUN_TEST(test_ensemble_vsm_sw, input_file, n);
    RUN_TEST(test_inline_refine, input_file, n);
    RUN_TEST(test_param_sweep, input_file, n * 3);
    RUN_TEST(test_align_full_single, input_file, n);
    RUN_TEST(test_align_full_multi, input_file, n);

    fprintf(stderr, "\n=== Results: %d passed, %d failed ===\n", n_passed, n_failed);
    return n_failed > 0 ? 1 : 0;
}
