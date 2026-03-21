/*
 * bench_threadpool.c — Benchmarks for the three kalign parallelism patterns.
 *
 *   1. Distance matrix  (parallel for)
 *   2. Recursive tree    (task spawning)
 *   3. Fork-join 2 tasks (Hirschberg-like)
 *
 * Each benchmark compares serial, threadpool (1..N threads),
 * and optionally OpenMP.
 */

#include "threadpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/* ── Timing ───────────────────────────────────────────────────── */

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

#define NREPS 7

static double median(double *arr, int n)
{
    qsort(arr, (size_t)n, sizeof(double), cmp_double);
    return arr[n / 2];
}

static void fmt_time(double secs, char *buf, int len)
{
    if (secs < 1e-3)
        snprintf(buf, (size_t)len, "%7.1f us", secs * 1e6);
    else if (secs < 1.0)
        snprintf(buf, (size_t)len, "%7.2f ms", secs * 1e3);
    else
        snprintf(buf, (size_t)len, "%7.3f s ", secs);
}

/* ── Workload helpers ─────────────────────────────────────────── */

static uint8_t **create_sequences(int n, int len)
{
    uint8_t **seqs = (uint8_t **)malloc((size_t)n * sizeof(uint8_t *));
    for (int i = 0; i < n; i++) {
        seqs[i] = (uint8_t *)malloc((size_t)len);
        unsigned x = (unsigned)i * 2654435761u;
        for (int j = 0; j < len; j++) {
            x = x * 1103515245u + 12345u;
            seqs[i][j] = (uint8_t)((x >> 16) % 20);
        }
    }
    return seqs;
}

static void free_sequences(uint8_t **seqs, int n)
{
    for (int i = 0; i < n; i++) free(seqs[i]);
    free(seqs);
}

/* ══════════════════════════════════════════════════════════════════
 * 1. DISTANCE MATRIX (parallel for)
 * ══════════════════════════════════════════════════════════════════ */

struct dm_ctx {
    float   **dm;
    uint8_t **seqs;
    int       seq_len;
    int       num_anchors;
};

static void dm_chunk(int start, int end, void *arg)
{
    struct dm_ctx *ctx = (struct dm_ctx *)arg;
    for (int i = start; i < end; i++) {
        for (int j = 0; j < ctx->num_anchors; j++) {
            int diff = 0;
            for (int k = 0; k < ctx->seq_len; k++)
                diff += (ctx->seqs[i][k] != ctx->seqs[j][k]);
            ctx->dm[i][j] = (float)diff / (float)ctx->seq_len;
        }
    }
}

static double bench_dm_serial(uint8_t **seqs, float **dm,
                              int N, int M, int L)
{
    struct dm_ctx ctx = { dm, seqs, L, M };
    double t0 = now();
    dm_chunk(0, N, &ctx);
    return now() - t0;
}

static double bench_dm_tp(uint8_t **seqs, float **dm,
                          int N, int M, int L, int nthreads)
{
    struct dm_ctx ctx = { dm, seqs, L, M };
    threadpool_t *pool = tp_create(nthreads);
    double t0 = now();
    tp_parallel_for(pool, 0, N, dm_chunk, &ctx);
    double elapsed = now() - t0;
    tp_destroy(pool);
    return elapsed;
}

#ifdef HAVE_OPENMP
static double bench_dm_omp(uint8_t **seqs, float **dm,
                           int N, int M, int L)
{
    double t0 = now();
    int i, j, k;
    #pragma omp parallel for shared(dm, seqs) private(i, j, k) schedule(static)
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            int diff = 0;
            for (k = 0; k < L; k++)
                diff += (seqs[i][k] != seqs[j][k]);
            dm[i][j] = (float)diff / (float)L;
        }
    }
    return now() - t0;
}
#endif

/* ══════════════════════════════════════════════════════════════════
 * 2. RECURSIVE TREE (task spawning)
 * ══════════════════════════════════════════════════════════════════ */

struct tree_node {
    struct tree_node *left;
    struct tree_node *right;
    int   leaf_id;
    int   depth;
    float result;
};

static struct tree_node *build_tree(int depth, int *next_id)
{
    struct tree_node *n = (struct tree_node *)calloc(1, sizeof(*n));
    n->depth = depth;
    if (depth == 0) {
        n->leaf_id = (*next_id)++;
        return n;
    }
    n->leaf_id = -1;
    n->left  = build_tree(depth - 1, next_id);
    n->right = build_tree(depth - 1, next_id);
    return n;
}

static void free_tree(struct tree_node *n)
{
    if (!n) return;
    free_tree(n->left);
    free_tree(n->right);
    free(n);
}

#define LEAF_ITERS 500

static float leaf_work(int id)
{
    float sum = 0.0f;
    unsigned x = (unsigned)id * 2654435761u;
    for (int i = 0; i < LEAF_ITERS; i++) {
        x = x * 1103515245u + 12345u;
        sum += (float)(x >> 16) / 65536.0f;
    }
    return sum;
}

/* Serial */
static void tree_serial(struct tree_node *n)
{
    if (n->leaf_id >= 0) { n->result = leaf_work(n->leaf_id); return; }
    tree_serial(n->left);
    tree_serial(n->right);
    n->result = n->left->result + n->right->result;
}

/* Threadpool */
struct tree_tp_arg {
    threadpool_t     *pool;
    struct tree_node *node;
};

/* Serial cutoff: avoid task creation overhead for tiny subtrees.
 * Not needed for correctness (V2 Chase-Lev bounds stack to O(depth)),
 * but worthwhile for performance — each task group costs a malloc. */
#define TREE_TASK_CUTOFF 4

static void tree_tp_task(void *arg)
{
    struct tree_tp_arg *ta = (struct tree_tp_arg *)arg;
    struct tree_node *n = ta->node;

    if (n->leaf_id >= 0) { n->result = leaf_work(n->leaf_id); return; }

    if (n->depth <= TREE_TASK_CUTOFF) {
        tree_serial(n);
        return;
    }

    struct tree_tp_arg left  = { ta->pool, n->left };
    struct tree_tp_arg right = { ta->pool, n->right };

    tp_group_t *g = tp_group_create(ta->pool);
    tp_group_submit(g, tree_tp_task, &left);
    tp_group_submit(g, tree_tp_task, &right);
    tp_group_wait(g);
    tp_group_destroy(g);

    n->result = n->left->result + n->right->result;
}

#ifdef HAVE_OPENMP
static void tree_omp_task(struct tree_node *n)
{
    if (n->leaf_id >= 0) { n->result = leaf_work(n->leaf_id); return; }

    #pragma omp task shared(n)
    tree_omp_task(n->left);
    #pragma omp task shared(n)
    tree_omp_task(n->right);
    #pragma omp taskwait

    n->result = n->left->result + n->right->result;
}
#endif

/* ══════════════════════════════════════════════════════════════════
 * 3. FORK-JOIN (Hirschberg-style: 2 independent tasks)
 * ══════════════════════════════════════════════════════════════════ */

struct fj_sum_arg {
    const float *arr;
    int   start;
    int   end;
    double result;
};

static void fj_sum_task(void *arg)
{
    struct fj_sum_arg *a = (struct fj_sum_arg *)arg;
    double s = 0.0;
    for (int i = a->start; i < a->end; i++)
        s += (double)a->arr[i];
    a->result = s;
}

static double bench_fj_serial(const float *arr, int N)
{
    double s = 0.0;
    double t0 = now();
    for (int i = 0; i < N; i++) s += (double)arr[i];
    double elapsed = now() - t0;
    (void)s;
    return elapsed;
}

static double bench_fj_tp(const float *arr, int N, int nthreads)
{
    threadpool_t *pool = tp_create(nthreads);
    int mid = N / 2;
    struct fj_sum_arg left  = { arr, 0,   mid, 0.0 };
    struct fj_sum_arg right = { arr, mid, N,   0.0 };

    double t0 = now();
    tp_group_t *g = tp_group_create(pool);
    tp_group_submit(g, fj_sum_task, &left);
    tp_group_submit(g, fj_sum_task, &right);
    tp_group_wait(g);
    tp_group_destroy(g);
    double elapsed = now() - t0;

    (void)(left.result + right.result);
    tp_destroy(pool);
    return elapsed;
}

/* ══════════════════════════════════════════════════════════════════
 * 4. OVERHEAD MEASUREMENT
 * ══════════════════════════════════════════════════════════════════ */

static void noop_task(void *arg) { (void)arg; }

static double bench_overhead(int nthreads)
{
    threadpool_t *pool = tp_create(nthreads);
    int N = 10000;

    double t0 = now();
    for (int i = 0; i < N; i++) {
        tp_group_t *g = tp_group_create(pool);
        tp_group_submit(g, noop_task, NULL);
        tp_group_wait(g);
        tp_group_destroy(g);
    }
    double elapsed = now() - t0;

    tp_destroy(pool);
    return elapsed / N;
}

/* ══════════════════════════════════════════════════════════════════
 * MAIN
 * ══════════════════════════════════════════════════════════════════ */

int main(void)
{
    int ncpus = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpus <= 0) ncpus = 4;

    /* Thread counts to test */
    int tc[] = { 1, 2, 4, 8, 16 };
    int ntc = 0;
    for (int i = 0; i < 5; i++) {
        if (tc[i] <= ncpus * 2) ntc++;
    }

    printf("threadpool benchmark\n");
    printf("====================\n");
    printf("CPUs: %d\n\n", ncpus);

    /* ── Overhead ─────────────────────────────────────────────── */
    {
        char buf[32];
        double per_task = bench_overhead(ncpus);
        fmt_time(per_task, buf, sizeof(buf));
        printf("Per-task overhead (submit+wait+destroy): %s\n\n", buf);
    }

    /* ── Distance Matrix ──────────────────────────────────────── */
    {
        struct { int N, M, L; } sizes[] = {
            {  500,  50, 200 },
            { 2000, 100, 200 },
            { 5000, 100, 200 },
        };
        int nsizes = 3;

        printf("--- Distance Matrix (parallel for) ---\n");
        printf("  %-14s %10s", "Size", "serial");
        for (int t = 0; t < ntc; t++)
            printf("  tp(%2d)  ", tc[t]);
#ifdef HAVE_OPENMP
        printf("  omp     ");
#endif
        printf("\n");

        for (int s = 0; s < nsizes; s++) {
            int N = sizes[s].N, M = sizes[s].M, L = sizes[s].L;

            uint8_t **seqs = create_sequences(N, L);
            float **dm = (float **)malloc((size_t)N * sizeof(float *));
            for (int i = 0; i < N; i++)
                dm[i] = (float *)calloc((size_t)M, sizeof(float));

            char label[32];
            snprintf(label, sizeof(label), "%dx%d", N, M);
            printf("  %-14s", label);

            /* Serial */
            double times[NREPS];
            for (int r = 0; r < NREPS; r++)
                times[r] = bench_dm_serial(seqs, dm, N, M, L);
            double serial_t = median(times, NREPS);
            char buf[32];
            fmt_time(serial_t, buf, sizeof(buf));
            printf(" %s", buf);

            /* Threadpool */
            for (int t = 0; t < ntc; t++) {
                for (int r = 0; r < NREPS; r++)
                    times[r] = bench_dm_tp(seqs, dm, N, M, L, tc[t]);
                double tp_t = median(times, NREPS);
                fmt_time(tp_t, buf, sizeof(buf));
                printf(" %s", buf);
            }

#ifdef HAVE_OPENMP
            for (int r = 0; r < NREPS; r++)
                times[r] = bench_dm_omp(seqs, dm, N, M, L);
            double omp_t = median(times, NREPS);
            fmt_time(omp_t, buf, sizeof(buf));
            printf(" %s", buf);
#endif
            printf("\n");

            for (int i = 0; i < N; i++) free(dm[i]);
            free(dm);
            free_sequences(seqs, N);
        }
        printf("\n");
    }

    /* ── Recursive Tree ───────────────────────────────────────── */
    {
        int depths[] = { 10, 14, 17 };
        int ndepths = 3;

        printf("--- Recursive Tree (task spawning) ---\n");
        printf("  %-14s %10s", "Depth/Leaves", "serial");
        for (int t = 0; t < ntc; t++)
            printf("  tp(%2d)  ", tc[t]);
#ifdef HAVE_OPENMP
        printf("  omp     ");
#endif
        printf("\n");

        for (int d = 0; d < ndepths; d++) {
            int depth = depths[d];
            int nleaves = 1 << depth;

            char label[32];
            snprintf(label, sizeof(label), "d=%d (%d)", depth, nleaves);
            printf("  %-14s", label);

            double times[NREPS];
            char buf[32];

            /* Serial */
            for (int r = 0; r < NREPS; r++) {
                int id = 0;
                struct tree_node *root = build_tree(depth, &id);
                double t0 = now();
                tree_serial(root);
                times[r] = now() - t0;
                free_tree(root);
            }
            fmt_time(median(times, NREPS), buf, sizeof(buf));
            printf(" %s", buf);

            /* Threadpool */
            for (int t = 0; t < ntc; t++) {
                for (int r = 0; r < NREPS; r++) {
                    int id = 0;
                    struct tree_node *root = build_tree(depth, &id);
                    threadpool_t *pool = tp_create(tc[t]);
                    struct tree_tp_arg targ = { pool, root };
                    double t0 = now();
                    tree_tp_task(&targ);
                    times[r] = now() - t0;
                    tp_destroy(pool);
                    free_tree(root);
                }
                fmt_time(median(times, NREPS), buf, sizeof(buf));
                printf(" %s", buf);
            }

#ifdef HAVE_OPENMP
            for (int r = 0; r < NREPS; r++) {
                int id = 0;
                struct tree_node *root = build_tree(depth, &id);
                double t0 = now();
                #pragma omp parallel
                #pragma omp single nowait
                tree_omp_task(root);
                times[r] = now() - t0;
                free_tree(root);
            }
            fmt_time(median(times, NREPS), buf, sizeof(buf));
            printf(" %s", buf);
#endif
            printf("\n");
        }
        printf("\n");
    }

    /* ── Fork-Join ────────────────────────────────────────────── */
    {
        int fsizes[] = { 10000, 100000, 1000000, 10000000 };
        int nfsizes = 4;

        printf("--- Fork-Join (2 tasks, array sum) ---\n");
        printf("  %-14s %10s", "N", "serial");
        for (int t = 0; t < ntc; t++)
            printf("  tp(%2d)  ", tc[t]);
        printf("\n");

        for (int s = 0; s < nfsizes; s++) {
            int N = fsizes[s];
            float *arr = (float *)malloc((size_t)N * sizeof(float));
            unsigned x = 12345;
            for (int i = 0; i < N; i++) {
                x = x * 1103515245u + 12345u;
                arr[i] = (float)(x >> 16) / 65536.0f;
            }

            char label[32];
            if (N >= 1000000)
                snprintf(label, sizeof(label), "%dM", N / 1000000);
            else
                snprintf(label, sizeof(label), "%dK", N / 1000);
            printf("  %-14s", label);

            double times[NREPS];
            char buf[32];

            for (int r = 0; r < NREPS; r++)
                times[r] = bench_fj_serial(arr, N);
            fmt_time(median(times, NREPS), buf, sizeof(buf));
            printf(" %s", buf);

            for (int t = 0; t < ntc; t++) {
                for (int r = 0; r < NREPS; r++)
                    times[r] = bench_fj_tp(arr, N, tc[t]);
                fmt_time(median(times, NREPS), buf, sizeof(buf));
                printf(" %s", buf);
            }
            printf("\n");

            free(arr);
        }
        printf("\n");
    }

    return 0;
}
