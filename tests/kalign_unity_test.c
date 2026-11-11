/*
 * kalign_unity_test.c - the C library and the command line give the same
 * alignment.
 *
 * For a protein and a DNA input:
 *   (a) kalign() on arrays        == kalign --mode fast
 *   (b) kalign_align_mode(msa, m) == kalign --mode m   (fast and accurate)
 *   (c) kalign() with gpo/gpe/tgpe == kalign --gpo --gpe --tgpe
 *   (d) kalign() with an explicit type == the same with --type
 *   (e) kalign_align_mode refuses an unknown mode
 *   (f) mode aliases and case: NULL and "default" == "fast"; "recall" and
 *       "Accurate" == "accurate"; presets have 1 (fast) and 5 (accurate) runs
 *   (g) kalign() refuses an unknown type; accurate refuses more than
 *       KALIGN_ACCURATE_MAX_SEQ sequences
 * Alignments are compared as FASTA (aligned sequences only, not names).
 *
 * Usage: kalign_unity_test <kalign binary> <protein file> <dna file> <work dir>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <kalign/kalign.h>
#include "msa_struct.h"
#include "msa_op.h"

static int n_failed = 0;

#define CHECK(cond, ...) do { if(!(cond)){ \
        fprintf(stderr, "  FAIL: " __VA_ARGS__); fprintf(stderr, "\n"); n_failed++; return -1; } } while(0)

struct seqs{
        char **s;
        int n;
};

static void seqs_free(struct seqs *q)
{
        for(int i = 0; i < q->n; i++) free(q->s[i]);
        free(q->s);
        q->s = NULL;
        q->n = 0;
}

/* Read a FASTA file into strings, one per record (lines joined, names dropped). */
static int read_fasta(const char *path, struct seqs *q)
{
        FILE *f = fopen(path, "r");
        char line[65536];
        size_t cap = 0;
        q->s = NULL;
        q->n = 0;
        if(!f) return -1;
        while(fgets(line, sizeof(line), f)){
                size_t l = strlen(line);
                while(l && (line[l-1] == '\n' || line[l-1] == '\r')) line[--l] = 0;
                if(line[0] == '>'){
                        q->s = realloc(q->s, sizeof(char*) * (size_t)(q->n + 1));
                        q->s[q->n] = calloc(1, 1);
                        q->n++;
                        cap = 0;
                }else if(q->n > 0){
                        char *t = realloc(q->s[q->n - 1], cap + l + 1);
                        memcpy(t + cap, line, l + 1);
                        cap += l;
                        q->s[q->n - 1] = t;
                }
        }
        fclose(f);
        return 0;
}

static int same(const struct seqs *a, const struct seqs *b, const char *what)
{
        CHECK(a->n == b->n && a->n > 0, "%s: %d vs %d sequences", what, a->n, b->n);
        for(int i = 0; i < a->n; i++){
                CHECK(strcmp(a->s[i], b->s[i]) == 0, "%s: sequence %d differs\n%s\n%s",
                      what, i, a->s[i], b->s[i]);
        }
        return 0;
}

static int run_cli(const char *bin, const char *in, const char *out, const char *extra)
{
        char cmd[4096];
        snprintf(cmd, sizeof(cmd), "\"%s\" -i \"%s\" -o \"%s\" -f fasta -q -n 1 %s >/dev/null 2>&1",
                 bin, in, out, extra);
        return system(cmd);
}

/* Unaligned sequences of a file, as arrays for kalign(). */
static int load_arrays(const char *path, char ***seq, int **len, int *n)
{
        struct msa *msa = NULL;
        if(kalign_read_input((char*)path, &msa, 1) != 0 || !msa) return -1;
        *n = msa->numseq;
        *seq = malloc(sizeof(char*) * (size_t)*n);
        *len = malloc(sizeof(int) * (size_t)*n);
        for(int i = 0; i < *n; i++){
                const char *src = msa->sequences[i]->seq;
                char *d = malloc(strlen(src) + 1);
                int w = 0;
                for(int j = 0; src[j]; j++) if(src[j] != '-') d[w++] = src[j];
                d[w] = 0;
                (*seq)[i] = d;
                (*len)[i] = w;
        }
        kalign_free_msa(msa);
        return 0;
}

static int lib_arrays(char **seq, int *len, int n, int type, float gpo, float gpe, float tgpe,
                      struct seqs *out)
{
        char **aln = NULL;
        int alen = 0;
        CHECK(kalign(seq, len, n, 1, type, gpo, gpe, tgpe, &aln, &alen) == 0, "kalign() failed");
        out->s = aln;
        out->n = n;
        return 0;
}

static int cli_result(const char *bin, const char *in, const char *work, const char *extra,
                      struct seqs *out)
{
        char outp[2048];
        snprintf(outp, sizeof(outp), "%s/cli.fa", work);
        CHECK(run_cli(bin, in, outp, extra) == 0, "kalign %s failed", extra);
        CHECK(read_fasta(outp, out) == 0, "cannot read %s", outp);
        return 0;
}

static int test_file(const char *bin, const char *in, const char *work, int is_dna)
{
        char **seq = NULL;
        int *len = NULL, n = 0;
        struct seqs lib, cli;
        const int auto_type = KALIGN_TYPE_UNDEFINED;
        const int type = is_dna ? KALIGN_TYPE_DNA : KALIGN_TYPE_PROTEIN;
        const char *tname = is_dna ? "dna" : "protein";

        fprintf(stderr, "--- %s ---\n", in);
        CHECK(load_arrays(in, &seq, &len, &n) == 0, "cannot read %s", in);

        /* (a) kalign() == --mode fast, with detected and with explicit type */
        CHECK(lib_arrays(seq, len, n, auto_type, -1, -1, -1, &lib) == 0, "lib");
        CHECK(cli_result(bin, in, work, "--mode fast", &cli) == 0, "cli");
        CHECK(same(&lib, &cli, "(a) kalign() auto vs --mode fast") == 0, "(a)");
        seqs_free(&lib);
        seqs_free(&cli);

        /* (d) explicit type == --type */
        CHECK(lib_arrays(seq, len, n, type, -1, -1, -1, &lib) == 0, "lib");
        {
                char extra[128];
                snprintf(extra, sizeof(extra), "--mode fast --type %s", tname);
                CHECK(cli_result(bin, in, work, extra, &cli) == 0, "cli");
        }
        CHECK(same(&lib, &cli, "(d) kalign() with type vs --type") == 0, "(d)");
        seqs_free(&lib);
        seqs_free(&cli);

        /* (c) gap penalties */
        CHECK(lib_arrays(seq, len, n, auto_type, 9.5f, 1.5f, 0.75f, &lib) == 0, "lib");
        CHECK(cli_result(bin, in, work, "--mode fast --gpo 9.5 --gpe 1.5 --tgpe 0.75", &cli) == 0, "cli");
        CHECK(same(&lib, &cli, "(c) kalign() gaps vs --gpo --gpe --tgpe") == 0, "(c)");
        seqs_free(&lib);
        seqs_free(&cli);

        /* the gap penalties must matter, else (c) proves nothing */
        {
                struct seqs plain;
                CHECK(lib_arrays(seq, len, n, auto_type, -1, -1, -1, &plain) == 0, "lib");
                CHECK(lib_arrays(seq, len, n, auto_type, 30.0f, 0.1f, 0.1f, &lib) == 0, "lib");
                int differ = 0;
                for(int i = 0; i < n; i++) if(strcmp(plain.s[i], lib.s[i]) != 0) differ = 1;
                CHECK(differ, "(c) changing the gap penalties did not change the alignment");
                seqs_free(&plain);
                seqs_free(&lib);
        }

        for(int i = 0; i < n; i++) free(seq[i]);
        free(seq);
        free(len);
        return 0;
}

static int test_mode(const char *bin, const char *in, const char *work, const char *mode)
{
        struct msa *msa = NULL;
        struct seqs lib, cli;
        char outp[2048], extra[128];

        CHECK(kalign_read_input((char*)in, &msa, 1) == 0 && msa, "cannot read %s", in);
        msa->quiet = 1;
        CHECK(kalign_align_mode(msa, mode, KALIGN_SEQ_AUTO, 1) == 0, "kalign_align_mode %s failed", mode);
        snprintf(outp, sizeof(outp), "%s/lib.fa", work);
        CHECK(kalign_write_msa(msa, outp, "fasta") == 0, "write failed");
        kalign_free_msa(msa);
        CHECK(read_fasta(outp, &lib) == 0, "cannot read %s", outp);

        snprintf(extra, sizeof(extra), "--mode %s", mode);
        CHECK(cli_result(bin, in, work, extra, &cli) == 0, "cli");
        {
                char what[256];
                snprintf(what, sizeof(what), "(b) kalign_align_mode(%s) vs --mode %s on %s", mode, mode, in);
                CHECK(same(&lib, &cli, what) == 0, "(b)");
        }
        seqs_free(&lib);
        seqs_free(&cli);
        return 0;
}

static int test_bad_mode(const char *in)
{
        struct msa *msa = NULL;
        CHECK(kalign_read_input((char*)in, &msa, 1) == 0 && msa, "cannot read %s", in);
        msa->quiet = 1;
        CHECK(kalign_align_mode(msa, "nonsense", KALIGN_SEQ_AUTO, 1) != 0, "(e) unknown mode accepted");
        kalign_free_msa(msa);
        return 0;
}

static int lib_mode(const char *in, const char *work, const char *mode, struct seqs *out)
{
        struct msa *msa = NULL;
        char outp[2048];
        CHECK(kalign_read_input((char*)in, &msa, 1) == 0 && msa, "cannot read %s", in);
        msa->quiet = 1;
        CHECK(kalign_align_mode(msa, mode, KALIGN_SEQ_AUTO, 1) == 0, "kalign_align_mode %s failed",
              mode ? mode : "NULL");
        snprintf(outp, sizeof(outp), "%s/alias.fa", work);
        CHECK(kalign_write_msa(msa, outp, "fasta") == 0, "write failed");
        kalign_free_msa(msa);
        CHECK(read_fasta(outp, out) == 0, "cannot read %s", outp);
        return 0;
}

static int test_aliases(const char *in, const char *work)
{
        const char *fast_alias[] = {NULL, "default", "FAST"};
        const char *acc_alias[] = {"recall", "Accurate"};
        struct seqs ref, other;
        CHECK(lib_mode(in, work, "fast", &ref) == 0, "lib");
        for(int i = 0; i < 3; i++){
                CHECK(lib_mode(in, work, fast_alias[i], &other) == 0, "lib");
                CHECK(same(&ref, &other, "(f) fast alias") == 0, "(f)");
                seqs_free(&other);
        }
        seqs_free(&ref);
        CHECK(lib_mode(in, work, "accurate", &ref) == 0, "lib");
        for(int i = 0; i < 2; i++){
                CHECK(lib_mode(in, work, acc_alias[i], &other) == 0, "lib");
                CHECK(same(&ref, &other, "(f) accurate alias") == 0, "(f)");
                seqs_free(&other);
        }
        seqs_free(&ref);
        return 0;
}

static int test_presets(void)
{
        struct kalign_run_config runs[KALIGN_MAX_PRESET_RUNS];
        struct kalign_ensemble_config ens;
        int n_runs;
        const int types[] = {ALN_BIOTYPE_PROTEIN, ALN_BIOTYPE_DNA};
        for(int t = 0; t < 2; t++){
                n_runs = 0;
                CHECK(kalign_get_mode_preset("fast", types[t], runs, &n_runs, &ens) == 0 && n_runs == 1,
                      "(f) fast preset has %d runs", n_runs);
                n_runs = 0;
                CHECK(kalign_get_mode_preset("accurate", types[t], runs, &n_runs, &ens) == 0 && n_runs == 5,
                      "(f) accurate preset has %d runs", n_runs);
        }
        return 0;
}

static int test_errors(void)
{
        char *seq[2] = {"ACDEFGHIK", "ACDFGHIK"};
        int len[2] = {9, 8};
        char **aln = NULL;
        int alen = 0;
        struct msa *msa = NULL;
        char **many = malloc(sizeof(char*) * 5001);
        int *mlen = malloc(sizeof(int) * 5001);

        CHECK(kalign(seq, len, 2, 1, 99, -1, -1, -1, &aln, &alen) != 0, "(g) type 99 accepted");
        CHECK(kalign(seq, len, 2, 1, -7, -1, -1, -1, &aln, &alen) != 0, "(g) type -7 accepted");

        for(int i = 0; i < 5001; i++){
                many[i] = "ACDEFGHIKL";
                mlen[i] = 10;
        }
        CHECK(kalign_arr_to_msa(many, mlen, 5001, &msa) == 0 && msa, "cannot build msa");
        msa->quiet = 1;
        CHECK(kalign_align_mode(msa, "accurate", KALIGN_SEQ_AUTO, 1) != 0, "(g) accurate with 5001 sequences accepted");
        CHECK(kalign_align_mode(msa, "recall", KALIGN_SEQ_AUTO, 1) != 0, "(g) recall with 5001 sequences accepted");
        kalign_free_msa(msa);
        free(many);
        free(mlen);
        CHECK(kalign_mode_size_ok("accurate", ALN_BIOTYPE_PROTEIN, KALIGN_ACCURATE_MAX_SEQ), "(g) 5000 refused");
        CHECK(!kalign_mode_size_ok("accurate", ALN_BIOTYPE_PROTEIN, KALIGN_ACCURATE_MAX_SEQ + 1), "(g) 5001 allowed");
        return 0;
}

int main(int argc, char *argv[])
{
        const char *modes[] = {"fast", "accurate"};
        if(argc != 5){
                fprintf(stderr, "usage: %s <kalign> <protein file> <dna file> <work dir>\n", argv[0]);
                return EXIT_FAILURE;
        }
        {
                char cmd[2048];
                snprintf(cmd, sizeof(cmd), "mkdir -p \"%s\"", argv[4]);
                if(system(cmd) != 0) return EXIT_FAILURE;
        }
        test_file(argv[1], argv[2], argv[4], 0);
        test_file(argv[1], argv[3], argv[4], 1);
        for(int m = 0; m < 2; m++){
                test_mode(argv[1], argv[2], argv[4], modes[m]);
                test_mode(argv[1], argv[3], argv[4], modes[m]);
        }
        test_bad_mode(argv[2]);
        test_aliases(argv[2], argv[4]);
        test_aliases(argv[3], argv[4]);
        test_presets();
        test_errors();
        if(n_failed){
                fprintf(stderr, "%d check(s) failed\n", n_failed);
                return EXIT_FAILURE;
        }
        fprintf(stderr, "all unity checks passed\n");
        return EXIT_SUCCESS;
}
