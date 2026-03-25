/*
 * kalign_api_test.c — comprehensive tests for the kalign public C API.
 *
 * Tests the unified kalign_align_full entry point with various configs:
 *   - Single run with refine modes
 *   - VSM + seq_weights via run config
 *   - Seeded tree + consistency anchors
 *   - Realign iterations
 *   - Ensemble alignment
 *   - kalign_msa_compare_detailed()
 *   - kalign_msa_compare_with_mask()
 *   - kalign_check_msa()
 *   - reformat_settings_msa()
 *   - kalign_write_msa() round-trip fasta
 *
 * Each test:
 *   1. Reads input from the file passed as argv[1]
 *   2. Calls the API function
 *   3. Asserts structural correctness (alnlen>0, all seqs same length,
 *      non-gap chars preserved, scores in range)
 *   4. Returns 0 on success, -1 on failure
 *
 * Exit code: EXIT_SUCCESS if all pass, EXIT_FAILURE otherwise.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "kalign/kalign.h"
#include "msa_struct.h"
#include "msa_cmp.h"

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static int count_residues(const char *seq)
{
        int n = 0;
        for (int i = 0; seq[i]; i++) {
                if (seq[i] != '-') n++;
        }
        return n;
}

static int *snapshot_lengths(struct msa *m)
{
        int *lens = malloc(sizeof(int) * m->numseq);
        if (!lens) return NULL;
        for (int i = 0; i < m->numseq; i++) {
                lens[i] = m->sequences[i]->len;
        }
        return lens;
}

static int verify_alignment(struct msa *m, int *orig_lens, const char *label)
{
        if (m->alnlen <= 0) {
                fprintf(stderr, "  [%s] FAIL: alnlen=%d\n", label, m->alnlen);
                return -1;
        }
        for (int i = 0; i < m->numseq; i++) {
                if (m->sequences[i]->seq == NULL) {
                        fprintf(stderr, "  [%s] FAIL: seq %d is NULL\n", label, i);
                        return -1;
                }
                int slen = (int)strlen(m->sequences[i]->seq);
                if (slen != m->alnlen) {
                        fprintf(stderr, "  [%s] FAIL: seq %d strlen=%d != alnlen=%d\n",
                                label, i, slen, m->alnlen);
                        return -1;
                }
                if (orig_lens) {
                        int res = count_residues(m->sequences[i]->seq);
                        if (res != orig_lens[i]) {
                                fprintf(stderr, "  [%s] FAIL: seq %d residues=%d != orig=%d\n",
                                        label, i, res, orig_lens[i]);
                                return -1;
                        }
                }
        }
        return 0;
}

static struct msa *read_input(const char *path)
{
        struct msa *m = NULL;
        if (kalign_read_input((char *)path, &m, 1) != 0 || m == NULL) {
                fprintf(stderr, "  ERROR: cannot read %s\n", path);
                return NULL;
        }
        return m;
}

/* Helper: run a single alignment with default config */
static int align_default(struct msa *msa)
{
        struct kalign_run_config cfg = kalign_run_config_defaults();
        return kalign_align_full(msa, &cfg, 1, NULL, 1);
}

/* ------------------------------------------------------------------ */
/* individual test functions                                           */
/* ------------------------------------------------------------------ */

static int test_run_with_refine(const char *input)
{
        int modes[] = {KALIGN_REFINE_ALL, KALIGN_REFINE_CONFIDENT};
        const char *names[] = {"REFINE_ALL", "REFINE_CONFIDENT"};

        for (int m = 0; m < 2; m++) {
                struct msa *msa = read_input(input);
                if (!msa) return -1;
                int *lens = snapshot_lengths(msa);

                struct kalign_run_config cfg = kalign_run_config_defaults();
                cfg.refine = modes[m];
                int rv = kalign_align_full(msa, &cfg, 1, NULL, 1);
                if (rv != 0) {
                        fprintf(stderr, "  kalign_align_full(%s) returned %d\n", names[m], rv);
                        free(lens);
                        kalign_free_msa(msa);
                        return -1;
                }
                if (verify_alignment(msa, lens, names[m]) != 0) {
                        free(lens);
                        kalign_free_msa(msa);
                        return -1;
                }
                fprintf(stdout, "  %s: OK (alnlen=%d)\n", names[m], msa->alnlen);
                free(lens);
                kalign_free_msa(msa);
        }
        return 0;
}

static int test_run_dist_scale(const char *input)
{
        struct msa *msa = read_input(input);
        if (!msa) return -1;
        int *lens = snapshot_lengths(msa);

        struct kalign_run_config cfg = kalign_run_config_defaults();
        cfg.vsm_amax = 2.0f;
        cfg.seq_weights = 1.0f;
        int rv = kalign_align_full(msa, &cfg, 1, NULL, 1);
        if (rv != 0) {
                fprintf(stderr, "  align_full(vsm+sw) returned %d\n", rv);
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        if (verify_alignment(msa, lens, "dist_scale") != 0) {
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        fprintf(stdout, "  vsm_amax=2.0, seq_weights=1.0: OK (alnlen=%d)\n", msa->alnlen);
        free(lens);
        kalign_free_msa(msa);
        return 0;
}

static int test_run_seeded(const char *input)
{
        struct msa *msa = read_input(input);
        if (!msa) return -1;
        int *lens = snapshot_lengths(msa);

        struct kalign_run_config cfg = kalign_run_config_defaults();
        cfg.tree_seed = 42;
        cfg.vsm_amax = 2.0f;
        cfg.consistency_anchors = 5;
        cfg.consistency_weight = 2.0f;
        int rv = kalign_align_full(msa, &cfg, 1, NULL, 1);
        if (rv != 0) {
                fprintf(stderr, "  align_full(seeded+consistency) returned %d\n", rv);
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        if (verify_alignment(msa, lens, "seeded") != 0) {
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        fprintf(stdout, "  seeded + consistency: OK (alnlen=%d)\n", msa->alnlen);
        free(lens);
        kalign_free_msa(msa);
        return 0;
}

static int test_run_realign(const char *input)
{
        struct msa *msa = read_input(input);
        if (!msa) return -1;
        int *lens = snapshot_lengths(msa);

        struct kalign_run_config cfg = kalign_run_config_defaults();
        cfg.vsm_amax = 2.0f;
        cfg.realign = 1;
        int rv = kalign_align_full(msa, &cfg, 1, NULL, 1);
        if (rv != 0) {
                fprintf(stderr, "  align_full(realign=1) returned %d\n", rv);
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        if (verify_alignment(msa, lens, "realign") != 0) {
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        fprintf(stdout, "  realign=1, vsm=2.0: OK (alnlen=%d)\n", msa->alnlen);
        free(lens);
        kalign_free_msa(msa);
        return 0;
}

static int test_ensemble_with_realign(const char *input)
{
        struct msa *msa = read_input(input);
        if (!msa) return -1;
        int *lens = snapshot_lengths(msa);

        /* 3-run ensemble with realign and refine */
        struct kalign_run_config runs[3];
        for (int i = 0; i < 3; i++) {
                runs[i] = kalign_run_config_defaults();
                runs[i].vsm_amax = 2.0f;
                runs[i].realign = 1;
                runs[i].refine = KALIGN_REFINE_CONFIDENT;
                runs[i].tree_seed = 42 + (uint64_t)i;
                runs[i].tree_noise = (i > 0) ? 0.2f : 0.0f;
        }
        struct kalign_ensemble_config ens = kalign_ensemble_config_defaults();
        int rv = kalign_align_full(msa, runs, 3, &ens, 1);
        if (rv != 0) {
                fprintf(stderr, "  align_full(ensemble+realign) returned %d\n", rv);
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        if (verify_alignment(msa, lens, "ens+realign") != 0) {
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        /* Check col_confidence */
        if (msa->col_confidence == NULL) {
                fprintf(stderr, "  [ens+realign] FAIL: col_confidence is NULL\n");
                free(lens);
                kalign_free_msa(msa);
                return -1;
        }
        for (int i = 0; i < msa->alnlen; i++) {
                if (msa->col_confidence[i] < 0.0f || msa->col_confidence[i] > 1.0f) {
                        fprintf(stderr, "  [ens+realign] FAIL: col_confidence[%d]=%.3f\n",
                                i, msa->col_confidence[i]);
                        free(lens);
                        kalign_free_msa(msa);
                        return -1;
                }
        }
        fprintf(stdout, "  ensemble+realign+refine: OK (alnlen=%d)\n", msa->alnlen);
        free(lens);
        kalign_free_msa(msa);
        return 0;
}

static int test_compare_detailed(const char *input)
{
        struct msa *ref = read_input(input);
        struct msa *test_msa = read_input(input);
        if (!ref || !test_msa) return -1;

        align_default(ref);
        align_default(test_msa);

        struct poar_score out;
        memset(&out, 0, sizeof(out));
        int rv = kalign_msa_compare_detailed(ref, test_msa, -1.0f, &out);
        if (rv != 0) {
                fprintf(stderr, "  kalign_msa_compare_detailed returned %d\n", rv);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }

        if (fabs(out.recall - 1.0) > 0.001) {
                fprintf(stderr, "  FAIL: recall=%.4f (expected 1.0)\n", out.recall);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }
        if (fabs(out.precision - 1.0) > 0.001) {
                fprintf(stderr, "  FAIL: precision=%.4f (expected 1.0)\n", out.precision);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }
        if (fabs(out.f1 - 1.0) > 0.001) {
                fprintf(stderr, "  FAIL: f1=%.4f (expected 1.0)\n", out.f1);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }
        if (fabs(out.tc - 1.0) > 0.001) {
                fprintf(stderr, "  FAIL: tc=%.4f (expected 1.0)\n", out.tc);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }
        if (out.ref_pairs <= 0 || out.test_pairs <= 0 || out.common <= 0) {
                fprintf(stderr, "  FAIL: pair counts ref=%lld test=%lld common=%lld\n",
                        (long long)out.ref_pairs, (long long)out.test_pairs, (long long)out.common);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }
        if (out.ref_pairs != out.common || out.test_pairs != out.common) {
                fprintf(stderr, "  FAIL: self-compare pairs mismatch ref=%lld test=%lld common=%lld\n",
                        (long long)out.ref_pairs, (long long)out.test_pairs, (long long)out.common);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }

        fprintf(stdout, "  compare_detailed self: recall=%.3f prec=%.3f f1=%.3f tc=%.3f pairs=%lld OK\n",
                out.recall, out.precision, out.f1, out.tc, (long long)out.common);
        kalign_free_msa(ref);
        kalign_free_msa(test_msa);
        return 0;
}

static int test_compare_with_mask(const char *input)
{
        struct msa *ref = read_input(input);
        struct msa *test_msa = read_input(input);
        if (!ref || !test_msa) return -1;

        align_default(ref);
        align_default(test_msa);

        int n_cols = ref->alnlen;
        int *mask = malloc(sizeof(int) * n_cols);
        if (!mask) {
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }
        for (int i = 0; i < n_cols; i++) mask[i] = 1;

        struct poar_score out;
        memset(&out, 0, sizeof(out));
        int rv = kalign_msa_compare_with_mask(ref, test_msa, mask, n_cols, &out);
        if (rv != 0) {
                fprintf(stderr, "  kalign_msa_compare_with_mask returned %d\n", rv);
                free(mask);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }

        if (fabs(out.recall - 1.0) > 0.001 || fabs(out.precision - 1.0) > 0.001) {
                fprintf(stderr, "  FAIL: mask-all recall=%.4f prec=%.4f\n", out.recall, out.precision);
                free(mask);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }

        /* Test with partial mask */
        for (int i = n_cols / 2; i < n_cols; i++) mask[i] = 0;
        memset(&out, 0, sizeof(out));
        rv = kalign_msa_compare_with_mask(ref, test_msa, mask, n_cols, &out);
        if (rv != 0) {
                fprintf(stderr, "  kalign_msa_compare_with_mask (partial) returned %d\n", rv);
                free(mask);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }
        if (fabs(out.recall - 1.0) > 0.001) {
                fprintf(stderr, "  FAIL: partial mask recall=%.4f\n", out.recall);
                free(mask);
                kalign_free_msa(ref);
                kalign_free_msa(test_msa);
                return -1;
        }

        fprintf(stdout, "  compare_with_mask: full-mask OK, partial-mask OK\n");
        free(mask);
        kalign_free_msa(ref);
        kalign_free_msa(test_msa);
        return 0;
}

static int test_check_msa(const char *input)
{
        struct msa *msa = read_input(input);
        if (!msa) return -1;

        msa->quiet = 1;
        int rv = kalign_check_msa(msa, 0);
        if (rv != 0) {
                fprintf(stderr, "  kalign_check_msa returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }
        fprintf(stdout, "  check_msa: OK (%d sequences validated)\n", msa->numseq);
        kalign_free_msa(msa);
        return 0;
}

static int test_reformat_settings(const char *input)
{
        struct msa *msa = read_input(input);
        if (!msa) return -1;

        int rv = align_default(msa);
        if (rv != 0) {
                fprintf(stderr, "  initial align failed\n");
                kalign_free_msa(msa);
                return -1;
        }

        rv = reformat_settings_msa(msa, 1, 0);
        if (rv != 0) {
                fprintf(stderr, "  reformat_settings_msa(rename) returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }
        for (int i = 0; i < msa->numseq; i++) {
                char expected[32];
                snprintf(expected, sizeof(expected), "SEQ%d", i + 1);
                if (strncmp(msa->sequences[i]->name, expected, strlen(expected)) != 0) {
                        fprintf(stderr, "  FAIL: seq %d name='%s' expected='%s'\n",
                                i, msa->sequences[i]->name, expected);
                        kalign_free_msa(msa);
                        return -1;
                }
        }
        fprintf(stdout, "  reformat rename: OK\n");

        rv = reformat_settings_msa(msa, 0, 1);
        if (rv != 0) {
                fprintf(stderr, "  reformat_settings_msa(unalign) returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }
        for (int i = 0; i < msa->numseq; i++) {
                for (int j = 0; j <= msa->sequences[i]->len; j++) {
                        if (msa->sequences[i]->gaps[j] != 0) {
                                fprintf(stderr, "  FAIL: seq %d gaps[%d]=%d after unalign\n",
                                        i, j, msa->sequences[i]->gaps[j]);
                                kalign_free_msa(msa);
                                return -1;
                        }
                }
        }
        fprintf(stdout, "  reformat unalign: OK (gaps[] zeroed)\n");

        kalign_free_msa(msa);
        return 0;
}

static int test_write_roundtrip(const char *input)
{
        struct msa *msa = read_input(input);
        if (!msa) return -1;

        int rv = align_default(msa);
        if (rv != 0) {
                kalign_free_msa(msa);
                return -1;
        }

        const char *tmpfile = "test_roundtrip.fa";
        rv = kalign_write_msa(msa, (char *)tmpfile, "fasta");
        if (rv != 0) {
                fprintf(stderr, "  kalign_write_msa(fasta) returned %d\n", rv);
                kalign_free_msa(msa);
                return -1;
        }

        struct msa *msa2 = NULL;
        rv = kalign_read_input((char *)tmpfile, &msa2, 1);
        if (rv != 0 || msa2 == NULL) {
                fprintf(stderr, "  failed to read back written fasta\n");
                kalign_free_msa(msa);
                remove(tmpfile);
                return -1;
        }

        if (msa->numseq != msa2->numseq) {
                fprintf(stderr, "  FAIL: numseq %d vs %d\n", msa->numseq, msa2->numseq);
                kalign_free_msa(msa);
                kalign_free_msa(msa2);
                remove(tmpfile);
                return -1;
        }

        float score = 0;
        rv = kalign_msa_compare(msa, msa2, &score);
        if (rv != 0) {
                fprintf(stderr, "  kalign_msa_compare on roundtrip returned %d\n", rv);
                kalign_free_msa(msa);
                kalign_free_msa(msa2);
                remove(tmpfile);
                return -1;
        }
        if (fabsf(score - 100.0f) > 0.01f) {
                fprintf(stderr, "  FAIL: roundtrip score=%.2f (expected 100.0)\n", score);
                kalign_free_msa(msa);
                kalign_free_msa(msa2);
                remove(tmpfile);
                return -1;
        }

        fprintf(stdout, "  write fasta roundtrip: score=%.1f OK\n", score);
        kalign_free_msa(msa);
        kalign_free_msa(msa2);
        remove(tmpfile);
        return 0;
}

/* ------------------------------------------------------------------ */
/* main: run all tests                                                 */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
        if (argc <= 1) {
                fprintf(stderr, "Usage: %s <input.tfa>\n", argv[0]);
                return EXIT_FAILURE;
        }

        const char *input = argv[1];
        int ret = EXIT_SUCCESS;
        int pass = 0, fail = 0;

        struct {
                const char *name;
                int (*fn)(const char *);
        } tests[] = {
                {"kalign_align_full + refine",      test_run_with_refine},
                {"VSM + seq_weights",               test_run_dist_scale},
                {"seeded tree + consistency",        test_run_seeded},
                {"realign iterations",              test_run_realign},
                {"ensemble + realign",              test_ensemble_with_realign},
                {"kalign_msa_compare_detailed",     test_compare_detailed},
                {"kalign_msa_compare_with_mask",    test_compare_with_mask},
                {"kalign_check_msa",                test_check_msa},
                {"reformat_settings_msa",           test_reformat_settings},
                {"write fasta roundtrip",           test_write_roundtrip},
        };

        int n_tests = (int)(sizeof(tests) / sizeof(tests[0]));

        fprintf(stdout, "=== Kalign API Tests (%d tests) ===\n\n", n_tests);

        for (int i = 0; i < n_tests; i++) {
                fprintf(stdout, "Test %d/%d: %s\n", i + 1, n_tests, tests[i].name);
                if (tests[i].fn(input) != 0) {
                        fprintf(stdout, "  FAIL\n\n");
                        fail++;
                        ret = EXIT_FAILURE;
                } else {
                        fprintf(stdout, "  PASS\n\n");
                        pass++;
                }
        }

        fprintf(stdout, "=== Results: %d passed, %d failed (of %d) ===\n",
                pass, fail, n_tests);

        return ret;
}
