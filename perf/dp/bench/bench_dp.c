/* bench_dp -- microbenchmark and oracle for Kalign's DP kernels.
 *
 * Loads DP problems dumped from real Kalign runs (.kdp, see kdp_format.h),
 * runs a kernel (default: "ref" = current Kalign code) and reports timing
 * per problem and aggregated per kernel type and size; records and checks
 * reference outputs so prototype kernels can be validated.
 *
 * usage: bench_dp [options] <file.kdp | dir> ...
 *   --list                 list registered kernels
 *   --kernel NAME          kernel to run (default ref)
 *   --ops fwd,bwd,align    operations to time (default all three; "none" = no timing)
 *   --min-time S           min seconds per (problem, op) (default 0.05)
 *   --min-reps N           min repetitions (default 3)
 *   --max-reps N           max repetitions (default 100000)
 *   --kind ss|sp|pp        only problems of this kind
 *   --max-files N          only the first N problems (sorted by name)
 *   --record DIR           write oracle outputs of the kernel to DIR
 *   --check DIR            compare kernel outputs with oracle in DIR
 *   --no-bonus             ignore the dumped consistency bonus (m->consistency
 *                          = NULL); default: use it, as the real run did
 *   --selfcheck            compare the kernel's path with the path recorded
 *                          in the dump by the real run (skips calls that used
 *                          midpoint flips, and bonus calls under --no-bonus)
 *   --json FILE            write per-problem results and aggregates as JSON
 *   -q                     no per-problem lines
 *
 * Timed operations (cells = len_a * len_b for all three):
 *   fwd    one forward pass over the full matrix (rows 0..len_a)
 *   bwd    one backward pass over the full matrix
 *   align  full Hirschberg alignment (fwd/bwd/meetup recursion + path), as
 *          aln_runner does in Kalign; about 2x the cells of one pass.
 *
 * Oracle comparison (--check), per problem:
 *   fwd/bwd rows   last state row of a full forward / backward pass, bitwise
 *   top            top-level meetup score + meet column + transition, bitwise
 *                  (the global DP score of Kalign's model); its difference
 *                  is reported as "top-level score diff"
 *   path           identical raw Hirschberg path (rows that differ)
 *   path score     kdp_path_score(proto path) - kdp_path_score(ref path):
 *                  an independent double-precision rescoring.  APPROXIMATE:
 *                  it agrees with the top-level score to 0.05 only on about a
 *                  third of the problems (gap-type switches at Hirschberg
 *                  junctions and the consistency bonus are not modelled
 *                  exactly), so use it only as a relative measure when paths
 *                  differ; identity of path and top score are the hard checks.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>

#include "kdp_kernel.h"

#define REF_MAGIC "KDPREF2"

struct ref_header {
        char magic[8];
        int32_t kind, len_a, len_b, meet, trans, pad;
        float score;            /* top-level Hirschberg meetup score */
        float fwd_score;        /* max state at (len_a, len_b) after a full forward pass */
        double path_score;      /* kdp_path_score() of the alignment path */
};

static double now_ns(void)
{
#ifdef __APPLE__
        return (double)clock_gettime_nsec_np(CLOCK_UPTIME_RAW);
#else
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return 1e9 * (double)ts.tv_sec + (double)ts.tv_nsec;
#endif
}

static const char* kind_name(int k)
{
        return k == KDP_SS ? "ss" : k == KDP_SP ? "sp" : "pp";
}

static const char* phase_name(int ph)
{
        static const char* n[] = { "other", "tree", "consbuild", "progressive", "refine", "realign", "ensemble" };
        return (ph >= 0 && ph < 7) ? n[ph] : "?";
}

/* ---------------- options ---------------- */
struct opts {
        const char* kernel;
        int op_fwd, op_bwd, op_align;
        double min_time;
        int min_reps, max_reps;
        int kind;
        int max_files;
        const char* record;
        const char* check;
        int selfcheck;
        int no_bonus;
        const char* json;
        int quiet;
};

/* ---------------- file collection ---------------- */
static int cmp_str(const void* a, const void* b)
{
        return strcmp(*(char* const*)a, *(char* const*)b);
}

static const char* base_of(const char* p)
{
        const char* b = strrchr(p, '/');
        return b ? b + 1 : p;
}

static int cmp_base(const void* a, const void* b)
{
        return strcmp(base_of(*(char* const*)a), base_of(*(char* const*)b));
}

/* Oracle files are keyed by the .kdp basename: refuse duplicates. */
static int check_unique_basenames(char** files, int n)
{
        char** tmp = malloc(sizeof(char*) * (size_t)n);
        int i, dup = 0;
        memcpy(tmp, files, sizeof(char*) * (size_t)n);
        qsort(tmp, (size_t)n, sizeof(char*), cmp_base);
        for(i = 1; i < n; i++){
                if(strcmp(base_of(tmp[i - 1]), base_of(tmp[i])) == 0){
                        fprintf(stderr, "bench_dp: duplicate problem name %s (%s, %s)\n",
                                base_of(tmp[i]), tmp[i - 1], tmp[i]);
                        dup = 1;
                }
        }
        free(tmp);
        return dup;
}

static int has_suffix(const char* s, const char* suf)
{
        size_t n = strlen(s), m = strlen(suf);
        return n >= m && strcmp(s + n - m, suf) == 0;
}

static void add_path(char*** list, int* n, int* alloc, const char* p)
{
        if(*n == *alloc){
                *alloc = *alloc ? *alloc * 2 : 256;
                *list = realloc(*list, sizeof(char*) * (size_t)*alloc);
        }
        (*list)[(*n)++] = strdup(p);
}

static void collect(const char* arg, char*** list, int* n, int* alloc)
{
        struct stat st;
        if(stat(arg, &st) != 0){
                fprintf(stderr, "bench_dp: cannot stat %s\n", arg);
                return;
        }
        if(S_ISDIR(st.st_mode)){
                DIR* d = opendir(arg);
                struct dirent* e;
                char buf[4096];
                if(!d) return;
                while((e = readdir(d)) != NULL){
                        if(e->d_name[0] == '.') continue;
                        snprintf(buf, sizeof(buf), "%s/%s", arg, e->d_name);
                        if(has_suffix(e->d_name, ".kdp")){
                                add_path(list, n, alloc, buf);
                        }else if(stat(buf, &st) == 0 && S_ISDIR(st.st_mode)){
                                collect(buf, list, n, alloc);
                        }
                }
                closedir(d);
        }else{
                add_path(list, n, alloc, arg);
        }
}

/* ---------------- computations ---------------- */

static void run_fwd_full(const struct kdp_resolved* r, const struct kdp_problem* p, struct aln_mem* m)
{
        kdp_mem_reset(p, m);
        r->fwd[p->h.kind](m);
}

static void run_bwd_full(const struct kdp_resolved* r, const struct kdp_problem* p, struct aln_mem* m)
{
        kdp_mem_reset(p, m);
        m->starta_2 = 0;
        m->enda_2 = p->h.len_a;
        r->bwd[p->h.kind](m);
}

/* The first (top-level) Hirschberg step: forward over [0,mid), backward
   over [mid,len_a), meetup.  Its score is the global alignment score. */
static void run_top(const struct kdp_resolved* r, const struct kdp_problem* p, struct aln_mem* m,
                    int* meet, int* trans, float* score)
{
        int old_cor[5];
        int mid;
        kdp_mem_reset(p, m);
        mid = p->h.len_a / 2;
        old_cor[0] = 0;
        old_cor[1] = p->h.len_a;
        old_cor[2] = 0;
        old_cor[3] = p->h.len_b;
        old_cor[4] = mid;
        m->enda = mid;
        m->starta_2 = mid;
        m->enda_2 = p->h.len_a;
        r->fwd[p->h.kind](m);
        r->bwd[p->h.kind](m);
        r->meet[p->h.kind](m, old_cor, meet, trans, score);
}

struct timing {
        double med_ns;
        double min_ns;
        int reps;
};

static int cmp_dbl(const void* a, const void* b)
{
        double x = *(const double*)a, y = *(const double*)b;
        return (x > y) - (x < y);
}

/* op: 0 fwd, 1 bwd, 2 align */
static struct timing time_op(const struct opts* o, const struct kdp_resolved* r,
                             const struct kdp_problem* p, struct aln_mem* m, int op)
{
        struct timing t;
        double* s = NULL;
        int cap = 64, n = 0;
        double total = 0.0;
        s = malloc(sizeof(double) * (size_t)cap);
        /* warm-up */
        if(op == 0) run_fwd_full(r, p, m);
        else if(op == 1) run_bwd_full(r, p, m);
        else kdp_align(r, p, m);
        while((n < o->min_reps || total < o->min_time * 1e9) && n < o->max_reps){
                double t0, t1;
                if(op == 0){
                        kdp_mem_reset(p, m);
                        t0 = now_ns();
                        r->fwd[p->h.kind](m);
                        t1 = now_ns();
                }else if(op == 1){
                        kdp_mem_reset(p, m);
                        m->starta_2 = 0;
                        m->enda_2 = p->h.len_a;
                        t0 = now_ns();
                        r->bwd[p->h.kind](m);
                        t1 = now_ns();
                }else{
                        t0 = now_ns();
                        kdp_align(r, p, m);
                        t1 = now_ns();
                }
                if(n == cap){
                        cap *= 2;
                        s = realloc(s, sizeof(double) * (size_t)cap);
                }
                s[n++] = t1 - t0;
                total += t1 - t0;
        }
        qsort(s, (size_t)n, sizeof(double), cmp_dbl);
        t.reps = n;
        t.min_ns = s[0];
        t.med_ns = s[n / 2];
        free(s);
        return t;
}

static double avg_nonzero(const float* prof, int len)
{
        long nz = 0;
        int i, c;
        if(!prof) return 0.0;
        for(i = 1; i <= len; i++){
                for(c = 0; c < 23; c++){
                        if(prof[(size_t)i * 64 + c] != 0.0F) nz++;
                }
        }
        return len > 0 ? (double)nz / (double)len : 0.0;
}

/* ---------------- oracle ---------------- */

static void ref_name(char* buf, size_t n, const char* dir, const char* kdp_path)
{
        const char* b = strrchr(kdp_path, '/');
        b = b ? b + 1 : kdp_path;
        snprintf(buf, n, "%s/%s.ref", dir, b);
}

struct outputs {
        struct states* fwd;
        struct states* bwd;
        int* path;
        int meet, trans;
        float score;
        float fwd_score;
        double path_score;
};

static void compute_outputs(const struct kdp_resolved* r, const struct kdp_problem* p,
                            struct aln_mem* m, struct outputs* out)
{
        int lb = p->h.len_b, la = p->h.len_a;
        out->fwd = malloc(sizeof(struct states) * (size_t)(lb + 1));
        out->bwd = malloc(sizeof(struct states) * (size_t)(lb + 1));
        out->path = malloc(sizeof(int) * (size_t)(la + 2));
        run_fwd_full(r, p, m);
        memcpy(out->fwd, m->f, sizeof(struct states) * (size_t)(lb + 1));
        out->fwd_score = out->fwd[lb].a;
        if(out->fwd[lb].ga > out->fwd_score) out->fwd_score = out->fwd[lb].ga;
        if(out->fwd[lb].gb > out->fwd_score) out->fwd_score = out->fwd[lb].gb;
        run_bwd_full(r, p, m);
        memcpy(out->bwd, m->b, sizeof(struct states) * (size_t)(lb + 1));
        run_top(r, p, m, &out->meet, &out->trans, &out->score);
        kdp_align(r, p, m);
        memcpy(out->path, m->path, sizeof(int) * (size_t)(la + 2));
        out->path_score = kdp_path_score(p, out->path);
}

static void free_outputs(struct outputs* o)
{
        free(o->fwd);
        free(o->bwd);
        free(o->path);
}

static int write_ref(const char* fn, const struct kdp_problem* p, const struct outputs* o)
{
        struct ref_header h;
        FILE* f = fopen(fn, "wb");
        int i;
        if(!f){
                fprintf(stderr, "bench_dp: cannot write %s\n", fn);
                return 1;
        }
        memset(&h, 0, sizeof(h));
        memcpy(h.magic, REF_MAGIC, 8);
        h.kind = p->h.kind;
        h.len_a = p->h.len_a;
        h.len_b = p->h.len_b;
        h.meet = o->meet;
        h.trans = o->trans;
        h.score = o->score;
        h.fwd_score = o->fwd_score;
        h.path_score = o->path_score;
        fwrite(&h, sizeof(h), 1, f);
        fwrite(o->fwd, sizeof(struct states), (size_t)(p->h.len_b + 1), f);
        fwrite(o->bwd, sizeof(struct states), (size_t)(p->h.len_b + 1), f);
        for(i = 0; i < p->h.len_a + 2; i++){
                int32_t v = o->path[i];
                fwrite(&v, sizeof(v), 1, f);
        }
        fclose(f);
        return 0;
}

static int read_ref(const char* fn, const struct kdp_problem* p, struct ref_header* h, struct outputs* o)
{
        FILE* f = fopen(fn, "rb");
        int i, ok = 1;
        if(!f) return 1;
        if(fread(h, sizeof(*h), 1, f) != 1 || memcmp(h->magic, REF_MAGIC, 8) != 0 ||
           h->len_a != p->h.len_a || h->len_b != p->h.len_b){
                fclose(f);
                return 1;
        }
        o->fwd = malloc(sizeof(struct states) * (size_t)(h->len_b + 1));
        o->bwd = malloc(sizeof(struct states) * (size_t)(h->len_b + 1));
        o->path = malloc(sizeof(int) * (size_t)(h->len_a + 2));
        ok &= fread(o->fwd, sizeof(struct states), (size_t)(h->len_b + 1), f) == (size_t)(h->len_b + 1);
        ok &= fread(o->bwd, sizeof(struct states), (size_t)(h->len_b + 1), f) == (size_t)(h->len_b + 1);
        for(i = 0; i < h->len_a + 2; i++){
                int32_t v = 0;
                ok &= fread(&v, sizeof(v), 1, f) == 1;
                o->path[i] = v;
        }
        fclose(f);
        o->meet = h->meet;
        o->trans = h->trans;
        o->score = h->score;
        o->fwd_score = h->fwd_score;
        o->path_score = h->path_score;
        return ok ? 0 : 1;
}

/* max |a-b| over finite entries; counts entries that are not bit-identical */
static double row_diff(const struct states* a, const struct states* b, int n, int* n_diff)
{
        double mx = 0.0;
        int i, k;
        *n_diff = 0;
        for(i = 0; i < n; i++){
                const float* x = &a[i].a;
                const float* y = &b[i].a;
                for(k = 0; k < 3; k++){
                        if(memcmp(&x[k], &y[k], sizeof(float)) != 0){
                                double d;
                                (*n_diff)++;
                                if(x[k] <= -FLT_MAX / 2 || y[k] <= -FLT_MAX / 2){
                                        /* one side "-inf": only a difference if the other is finite */
                                        if(!(x[k] <= -FLT_MAX / 2 && y[k] <= -FLT_MAX / 2)) d = INFINITY;
                                        else d = 0.0;
                                }else{
                                        d = fabs((double)x[k] - (double)y[k]);
                                }
                                if(d > mx) mx = d;
                        }
                }
        }
        return mx;
}

static int path_diff(const int* a, const int* b, int la)
{
        int i, n = 0;
        for(i = 1; i <= la; i++){
                if(a[i] != b[i]) n++;
        }
        return n;
}

/* ---------------- aggregation ---------------- */

#define N_BINS 5
static const double bin_edges[N_BINS] = { 1e4, 1e5, 1e6, 1e7, 1e18 };
static const char* bin_names[N_BINS] = { "<1e4", "1e4-1e5", "1e5-1e6", "1e6-1e7", ">=1e7" };

struct agg {
        double ns[3];
        double cells[3];
        int n[3];
};

static int bin_of(double cells)
{
        int b;
        for(b = 0; b < N_BINS; b++){
                if(cells < bin_edges[b]) return b;
        }
        return N_BINS - 1;
}

static void usage(void)
{
        fprintf(stderr,
                "usage: bench_dp [--list] [--kernel NAME] [--ops fwd,bwd,align|none] [--min-time S]\n"
                "                [--min-reps N] [--max-reps N] [--kind ss|sp|pp] [--max-files N]\n"
                "                [--record DIR] [--check DIR] [--selfcheck] [--no-bonus] [--json FILE] [-q]\n"
                "                <file.kdp | dir> ...\n");
}

int main(int argc, char** argv)
{
        struct opts o;
        char** files = NULL;
        int n_files = 0, alloc_files = 0;
        const struct kdp_kernel* k;
        struct kdp_resolved r;
        struct agg tot[3][N_BINS];
        struct agg all[3];
        FILE* js = NULL;
        int i, fi, n_done = 0;
        int n_check = 0, n_path_ident = 0, n_fwd_ident = 0, n_bwd_ident = 0, n_score_ident = 0;
        double worst_fwd_score_diff = 0.0, worst_score_diff = 0.0, worst_path_score_diff = 0.0, sum_path_score_diff = 0.0;
        int n_self = 0, n_self_ok = 0;
        int n_rec = 0, n_subopt = 0;
        double worst_model_gap = 0.0;
        int first_json = 1;

        memset(&o, 0, sizeof(o));
        o.kernel = "ref";
        o.op_fwd = o.op_bwd = o.op_align = 1;
        o.min_time = 0.05;
        o.min_reps = 3;
        o.max_reps = 100000;
        o.kind = -1;
        memset(tot, 0, sizeof(tot));
        memset(all, 0, sizeof(all));

        for(i = 1; i < argc; i++){
                const char* a = argv[i];
                if(!strcmp(a, "--list")){
                        int j;
                        for(j = 0; j < kdp_n_kernels(); j++){
                                const struct kdp_kernel* kk = kdp_kernel_at(j);
                                printf("%-20s %s\n", kk->name, kk->description ? kk->description : "");
                        }
                        return 0;
                }else if(!strcmp(a, "--kernel") && i + 1 < argc){ o.kernel = argv[++i];
                }else if(!strcmp(a, "--ops") && i + 1 < argc){
                        const char* s = argv[++i];
                        o.op_fwd = strstr(s, "fwd") != NULL;
                        o.op_bwd = strstr(s, "bwd") != NULL;
                        o.op_align = strstr(s, "align") != NULL;
                }else if(!strcmp(a, "--min-time") && i + 1 < argc){ o.min_time = atof(argv[++i]);
                }else if(!strcmp(a, "--min-reps") && i + 1 < argc){ o.min_reps = atoi(argv[++i]);
                }else if(!strcmp(a, "--max-reps") && i + 1 < argc){ o.max_reps = atoi(argv[++i]);
                }else if(!strcmp(a, "--kind") && i + 1 < argc){
                        const char* s = argv[++i];
                        o.kind = !strcmp(s, "ss") ? KDP_SS : !strcmp(s, "sp") ? KDP_SP : KDP_PP;
                }else if(!strcmp(a, "--max-files") && i + 1 < argc){ o.max_files = atoi(argv[++i]);
                }else if(!strcmp(a, "--record") && i + 1 < argc){ o.record = argv[++i];
                }else if(!strcmp(a, "--check") && i + 1 < argc){ o.check = argv[++i];
                }else if(!strcmp(a, "--selfcheck")){ o.selfcheck = 1;
                }else if(!strcmp(a, "--no-bonus")){ o.no_bonus = 1;
                }else if(!strcmp(a, "--json") && i + 1 < argc){ o.json = argv[++i];
                }else if(!strcmp(a, "-q")){ o.quiet = 1;
                }else if(a[0] == '-'){
                        usage();
                        return 1;
                }else{
                        collect(a, &files, &n_files, &alloc_files);
                }
        }
        if(n_files == 0){
                usage();
                return 1;
        }
        qsort(files, (size_t)n_files, sizeof(char*), cmp_str);
        if((o.record || o.check) && check_unique_basenames(files, n_files)){
                fprintf(stderr, "bench_dp: problem basenames must be unique for --record/--check "
                        "(dump with a distinct KALIGN_DP_DUMP_TAG per run)\n");
                return 1;
        }

        k = kdp_find_kernel(o.kernel);
        if(!k){
                fprintf(stderr, "bench_dp: unknown kernel '%s' (try --list)\n", o.kernel);
                return 1;
        }
        if(kdp_resolve(k, &r)) return 1;
        if(o.record){
                mkdir(o.record, 0755);
        }
        if(o.json){
                js = fopen(o.json, "w");
                if(!js){
                        fprintf(stderr, "bench_dp: cannot write %s\n", o.json);
                        return 1;
                }
                fprintf(js, "{\"kernel\": \"%s\", \"problems\": [\n", k->name);
        }
        if(!o.quiet){
                printf("%-44s %-4s %-11s %5s %5s %5s | %8s %8s %8s  ns/cell\n",
                       "problem", "kind", "phase", "len_a", "len_b", "nzA", "fwd", "bwd", "align");
        }

        for(fi = 0; fi < n_files; fi++){
                struct kdp_problem* p = NULL;
                struct aln_mem* m = NULL;
                struct timing tf = {0, 0, 0}, tb = {0, 0, 0}, ta = {0, 0, 0};
                double cells, prep_ns = 0.0;
                void* state = NULL;
                const char* base;
                int kind;
                /* check results */
                int have_check = 0, fwd_nd = 0, bwd_nd = 0, pdiff = 0, meet_same = 0;
                double fwd_md = 0.0, bwd_md = 0.0, sdiff = 0.0, psdiff = 0.0, fsdiff = 0.0;
                int self_pdiff = -1;
                double model_gap = NAN;

                if(o.max_files && n_done >= o.max_files) break;
                if(kdp_load(files[fi], &p)) continue;
                kind = p->h.kind;
                if(o.no_bonus) p->use_bonus = 0;
                if(o.kind >= 0 && kind != o.kind){
                        kdp_free(p);
                        continue;
                }
                if(p->h.len_a < 1 || p->h.len_b < 1){
                        kdp_free(p);
                        continue;
                }
                if(kdp_mem_create(p, &m)){
                        kdp_free(p);
                        continue;
                }
                n_done++;
                base = strrchr(p->path, '/');
                base = base ? base + 1 : p->path;
                cells = (double)p->h.len_a * (double)p->h.len_b;

                if(k->prepare){
                        double t0 = now_ns();
                        k->prepare(m, kind, &state);
                        prep_ns = now_ns() - t0;
                }
                kdp_kernel_state = state;

                if(o.op_fwd) tf = time_op(&o, &r, p, m, 0);
                if(o.op_bwd) tb = time_op(&o, &r, p, m, 1);
                if(o.op_align) ta = time_op(&o, &r, p, m, 2);

                if(o.record || o.check || o.selfcheck){
                        struct outputs out;
                        compute_outputs(&r, p, m, &out);
                        if(o.record){
                                char fn[4096];
                                ref_name(fn, sizeof(fn), o.record, p->path);
                                write_ref(fn, p, &out);
                                /* Scorer sanity: independent model score of the Hirschberg
                                   path vs. the top-level meetup score (forward over the top
                                   half + backward over the bottom half = Kalign's global
                                   optimum under its model; a full forward pass alone is not,
                                   because it does not treat gaps in the last row as terminal).
                                   Equal up to float rounding when scorer and kernels agree. */
                                n_rec++;
                                {
                                        double g = out.path_score - (double)out.score;
                                        double tol = 1e-4 * fabs((double)out.score) + 0.05;
                                        model_gap = g;
                                        if(fabs(g) > fabs(worst_model_gap)) worst_model_gap = g;
                                        if(fabs(g) > tol) n_subopt++;
                                }
                        }
                        if(o.check){
                                char fn[4096];
                                struct ref_header rh;
                                struct outputs ref;
                                ref_name(fn, sizeof(fn), o.check, p->path);
                                if(read_ref(fn, p, &rh, &ref) == 0){
                                        have_check = 1;
                                        n_check++;
                                        fwd_md = row_diff(out.fwd, ref.fwd, p->h.len_b + 1, &fwd_nd);
                                        bwd_md = row_diff(out.bwd, ref.bwd, p->h.len_b + 1, &bwd_nd);
                                        pdiff = path_diff(out.path, ref.path, p->h.len_a);
                                        meet_same = (out.meet == ref.meet && out.trans == ref.trans);
                                        sdiff = (double)out.score - (double)ref.score;
                                        fsdiff = (double)out.fwd_score - (double)ref.fwd_score;
                                        if(fabs(fsdiff) > fabs(worst_fwd_score_diff)) worst_fwd_score_diff = fsdiff;
                                        psdiff = out.path_score - ref.path_score;
                                        if(fwd_nd == 0) n_fwd_ident++;
                                        if(bwd_nd == 0) n_bwd_ident++;
                                        if(pdiff == 0) n_path_ident++;
                                        if(sdiff == 0.0 && meet_same) n_score_ident++;
                                        if(fabs(sdiff) > fabs(worst_score_diff)) worst_score_diff = sdiff;
                                        if(psdiff < worst_path_score_diff) worst_path_score_diff = psdiff;
                                        sum_path_score_diff += psdiff;
                                        free_outputs(&ref);
                                }else{
                                        fprintf(stderr, "bench_dp: no oracle for %s in %s\n", base, o.check);
                                }
                        }
                        if(o.selfcheck && p->ref_path && !p->h.flip && (p->use_bonus || !p->h.consistency)){
                                self_pdiff = path_diff(out.path, p->ref_path, p->h.len_a);
                                n_self++;
                                if(self_pdiff == 0) n_self_ok++;
                        }
                        free_outputs(&out);
                }

                if(k->release && state){
                        k->release(state);
                }
                kdp_kernel_state = NULL;

                {
                        int b = bin_of(cells);
                        if(o.op_fwd){ tot[kind][b].ns[0] += tf.med_ns; tot[kind][b].cells[0] += cells; tot[kind][b].n[0]++;
                                all[kind].ns[0] += tf.med_ns; all[kind].cells[0] += cells; all[kind].n[0]++; }
                        if(o.op_bwd){ tot[kind][b].ns[1] += tb.med_ns; tot[kind][b].cells[1] += cells; tot[kind][b].n[1]++;
                                all[kind].ns[1] += tb.med_ns; all[kind].cells[1] += cells; all[kind].n[1]++; }
                        if(o.op_align){ tot[kind][b].ns[2] += ta.med_ns; tot[kind][b].cells[2] += cells; tot[kind][b].n[2]++;
                                all[kind].ns[2] += ta.med_ns; all[kind].cells[2] += cells; all[kind].n[2]++; }
                }

                if(!o.quiet){
                        printf("%-44.44s %-4s %-11s %5d %5d %5.1f | %8.3f %8.3f %8.3f",
                               base, kind_name(kind), phase_name(p->h.phase), p->h.len_a, p->h.len_b,
                               avg_nonzero(p->prof1, p->h.len_a),
                               o.op_fwd ? tf.med_ns / cells : 0.0,
                               o.op_bwd ? tb.med_ns / cells : 0.0,
                               o.op_align ? ta.med_ns / cells : 0.0);
                        if(have_check){
                                printf("  chk opt:%s fwd:%s bwd:%s top:%s path:%s",
                                       fsdiff == 0.0 ? "same" : "DIFF",
                                       fwd_nd ? "DIFF" : "same", bwd_nd ? "DIFF" : "same",
                                       (sdiff == 0.0 && meet_same) ? "same" : "DIFF",
                                       pdiff ? "DIFF" : "same");
                                if(pdiff) printf(" (%d rows, dscore %.4g)", pdiff, psdiff);
                        }
                        if(self_pdiff >= 0){
                                printf("  self:%s", self_pdiff ? "DIFF" : "ok");
                        }
                        printf("\n");
                }
                if(js){
                        fprintf(js, "%s  {\"file\": \"%s\", \"kind\": \"%s\", \"phase\": \"%s\", \"len_a\": %d, \"len_b\": %d,"
                                " \"cells\": %.0f, \"nz_a\": %.3f, \"nz_b\": %.3f, \"consistency\": %d, \"flip\": %d,"
                                " \"prepare_ns\": %.0f,"
                                " \"fwd_ns\": %.1f, \"bwd_ns\": %.1f, \"align_ns\": %.1f,"
                                " \"fwd_min_ns\": %.1f, \"bwd_min_ns\": %.1f, \"align_min_ns\": %.1f,"
                                " \"reps\": [%d, %d, %d]",
                                first_json ? "" : ",\n", base, kind_name(kind), phase_name(p->h.phase),
                                p->h.len_a, p->h.len_b, cells,
                                avg_nonzero(p->prof1, p->h.len_a), avg_nonzero(p->prof2, p->h.len_b),
                                p->h.consistency, p->h.flip, prep_ns,
                                tf.med_ns, tb.med_ns, ta.med_ns, tf.min_ns, tb.min_ns, ta.min_ns,
                                tf.reps, tb.reps, ta.reps);
                        if(have_check){
                                fprintf(js, ", \"check\": {\"fwd_ndiff\": %d, \"fwd_maxabs\": %.6g, \"bwd_ndiff\": %d,"
                                        " \"bwd_maxabs\": %.6g, \"top_same\": %d, \"top_score_diff\": %.6g,"
                                        " \"path_rows_diff\": %d, \"path_score_diff\": %.6g, \"opt_score_diff\": %.6g}",
                                        fwd_nd, isinf(fwd_md) ? 1e30 : fwd_md, bwd_nd, isinf(bwd_md) ? 1e30 : bwd_md,
                                        (sdiff == 0.0 && meet_same), sdiff, pdiff, psdiff, fsdiff);
                        }
                        if(self_pdiff >= 0){
                                fprintf(js, ", \"selfcheck_rows_diff\": %d", self_pdiff);
                        }
                        if(!isnan(model_gap)){
                                fprintf(js, ", \"path_score_minus_top_score\": %.6g", model_gap);
                        }
                        fprintf(js, "}");
                        first_json = 0;
                }
                kdp_mem_free(m);
                kdp_free(p);
        }

        /* summary */
        printf("\nkernel: %s   problems: %d\n", k->name, n_done);
        printf("%-4s %-9s %6s | %10s %10s %10s | %9s %9s %9s\n", "kind", "cells", "n",
               "fwd ns/c", "bwd ns/c", "align ns/c", "fwd Gc/s", "bwd Gc/s", "aln Gc/s");
        for(i = 0; i < 3; i++){
                int b;
                for(b = 0; b <= N_BINS; b++){
                        const struct agg* a = (b < N_BINS) ? &tot[i][b] : &all[i];
                        int nn = a->n[0] > a->n[1] ? a->n[0] : a->n[1];
                        if(a->n[2] > nn) nn = a->n[2];
                        if(nn == 0) continue;
                        printf("%-4s %-9s %6d | %10.3f %10.3f %10.3f | %9.3f %9.3f %9.3f\n",
                               kind_name(i), b < N_BINS ? bin_names[b] : "ALL", nn,
                               a->cells[0] > 0 ? a->ns[0] / a->cells[0] : 0.0,
                               a->cells[1] > 0 ? a->ns[1] / a->cells[1] : 0.0,
                               a->cells[2] > 0 ? a->ns[2] / a->cells[2] : 0.0,
                               a->ns[0] > 0 ? a->cells[0] / a->ns[0] : 0.0,
                               a->ns[1] > 0 ? a->cells[1] / a->ns[1] : 0.0,
                               a->ns[2] > 0 ? a->cells[2] / a->ns[2] : 0.0);
                }
        }
        if(o.check){
                printf("\ncheck vs %s: %d problems | fwd rows identical %d | bwd rows identical %d | "
                       "top score+meet identical %d | path identical %d\n",
                       o.check, n_check, n_fwd_ident, n_bwd_ident, n_score_ident, n_path_ident);
                printf("worst optimum (full forward) score diff %.6g | worst top-level score diff %.6g |"
                       " worst path score diff (proto - ref) %.6g | sum path score diff %.6g\n",
                       worst_fwd_score_diff, worst_score_diff, worst_path_score_diff, sum_path_score_diff);
        }
        if(o.record){
                printf("recorded %d oracle files in %s; scorer(path) - top-level score: worst %.4g,"
                       " %d problems outside tolerance\n", n_rec, o.record, worst_model_gap, n_subopt);
        }
        if(o.selfcheck){
                printf("selfcheck: %d/%d problems reproduce the real-run path exactly\n", n_self_ok, n_self);
        }
        if(js){
                int b;
                int firsta = 1;
                fprintf(js, "\n], \"aggregate\": [\n");
                for(i = 0; i < 3; i++){
                        for(b = 0; b <= N_BINS; b++){
                                const struct agg* a = (b < N_BINS) ? &tot[i][b] : &all[i];
                                if(a->n[0] + a->n[1] + a->n[2] == 0) continue;
                                fprintf(js, "%s  {\"kind\": \"%s\", \"bin\": \"%s\", \"n\": [%d, %d, %d],"
                                        " \"fwd_ns_per_cell\": %.4f, \"bwd_ns_per_cell\": %.4f, \"align_ns_per_cell\": %.4f}",
                                        firsta ? "" : ",\n", kind_name(i), b < N_BINS ? bin_names[b] : "ALL",
                                        a->n[0], a->n[1], a->n[2],
                                        a->cells[0] > 0 ? a->ns[0] / a->cells[0] : 0.0,
                                        a->cells[1] > 0 ? a->ns[1] / a->cells[1] : 0.0,
                                        a->cells[2] > 0 ? a->ns[2] / a->cells[2] : 0.0);
                                firsta = 0;
                        }
                }
                fprintf(js, "\n], \"check\": {\"n\": %d, \"fwd_identical\": %d, \"bwd_identical\": %d,"
                         " \"top_identical\": %d, \"path_identical\": %d, \"worst_opt_score_diff\": %.6g, \"worst_top_score_diff\": %.6g,"
                        " \"worst_path_score_diff\": %.6g, \"sum_path_score_diff\": %.6g},"
                        " \"selfcheck\": {\"n\": %d, \"ok\": %d}}\n",
                        n_check, n_fwd_ident, n_bwd_ident, n_score_ident, n_path_ident,
                        worst_fwd_score_diff, worst_score_diff, worst_path_score_diff, sum_path_score_diff, n_self, n_self_ok);
                fclose(js);
        }
        for(i = 0; i < n_files; i++) free(files[i]);
        free(files);
        return 0;
}
