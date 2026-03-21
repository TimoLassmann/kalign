/*
 * test_threadpool.c — Tests for all three parallelism patterns.
 *
 *   1. Parallel for  (distance-matrix style)
 *   2. Fork-join      (Hirschberg forward/backward)
 *   3. Recursive tasks (guide-tree traversal)
 */

#include "threadpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdatomic.h>

/* ── Test harness ─────────────────────────────────────────────── */

static int g_pass, g_fail, g_total;

#define CHECK(cond) do {                                          \
    if (!(cond)) {                                                \
        fprintf(stderr, "    FAIL: %s  (line %d)\n",              \
                #cond, __LINE__);                                 \
        return -1;                                                \
    }                                                             \
} while (0)

#define RUN(name) do {                                            \
    printf("  %-45s ", #name);                                    \
    fflush(stdout);                                               \
    g_total++;                                                    \
    if (test_##name() == 0) { printf("PASS\n"); g_pass++; }      \
    else                    { printf("FAIL\n"); g_fail++; }       \
} while (0)

/* ── 1. create / destroy ─────────────────────────────────────── */

static int test_create_destroy(void)
{
    threadpool_t *p = tp_create(2);
    CHECK(p != NULL);
    CHECK(tp_get_nthreads(p) == 2);
    tp_destroy(p);
    return 0;
}

static int test_create_auto(void)
{
    threadpool_t *p = tp_create(0);
    CHECK(p != NULL);
    CHECK(tp_get_nthreads(p) >= 1);
    tp_destroy(p);
    return 0;
}

/* ── 2. parallel for ─────────────────────────────────────────── */

static void fill_squares(int start, int end, void *arg)
{
    int *arr = (int *)arg;
    for (int i = start; i < end; i++)
        arr[i] = i * i;
}

static int test_parallel_for_correctness(void)
{
    int N = 10000;
    int *arr = (int *)calloc((size_t)N, sizeof(int));
    CHECK(arr != NULL);

    threadpool_t *pool = tp_create(4);
    CHECK(pool != NULL);

    tp_parallel_for(pool, 0, N, fill_squares, arr);

    for (int i = 0; i < N; i++)
        CHECK(arr[i] == i * i);

    tp_destroy(pool);
    free(arr);
    return 0;
}

static int test_parallel_for_single_thread(void)
{
    int N = 5000;
    int *arr = (int *)calloc((size_t)N, sizeof(int));
    CHECK(arr != NULL);

    threadpool_t *pool = tp_create(1);
    tp_parallel_for(pool, 0, N, fill_squares, arr);

    for (int i = 0; i < N; i++)
        CHECK(arr[i] == i * i);

    tp_destroy(pool);
    free(arr);
    return 0;
}

static int test_parallel_for_empty_range(void)
{
    threadpool_t *pool = tp_create(2);
    /* Should be a no-op, not crash. */
    tp_parallel_for(pool, 5, 5, fill_squares, NULL);
    tp_parallel_for(pool, 10, 3, fill_squares, NULL);
    tp_destroy(pool);
    return 0;
}

/* Distance-matrix style: outer loop parallel, inner loop serial. */

struct dm_ctx {
    float **dm;
    uint8_t **seqs;
    int seq_len;
    int num_anchors;
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

static int test_parallel_for_distance_matrix(void)
{
    int N = 200, M = 20, L = 50;

    /* Allocate sequences */
    uint8_t **seqs = (uint8_t **)malloc((size_t)N * sizeof(uint8_t *));
    CHECK(seqs != NULL);
    for (int i = 0; i < N; i++) {
        seqs[i] = (uint8_t *)malloc((size_t)L);
        for (int j = 0; j < L; j++)
            seqs[i][j] = (uint8_t)((i * 7 + j * 13) % 20);
    }

    /* Serial reference */
    float **ref = (float **)malloc((size_t)N * sizeof(float *));
    for (int i = 0; i < N; i++) {
        ref[i] = (float *)calloc((size_t)M, sizeof(float));
        for (int j = 0; j < M; j++) {
            int diff = 0;
            for (int k = 0; k < L; k++)
                diff += (seqs[i][k] != seqs[j][k]);
            ref[i][j] = (float)diff / (float)L;
        }
    }

    /* Parallel computation */
    float **dm = (float **)malloc((size_t)N * sizeof(float *));
    for (int i = 0; i < N; i++)
        dm[i] = (float *)calloc((size_t)M, sizeof(float));

    struct dm_ctx ctx = { dm, seqs, L, M };
    threadpool_t *pool = tp_create(4);
    tp_parallel_for(pool, 0, N, dm_chunk, &ctx);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            CHECK(dm[i][j] == ref[i][j]);

    tp_destroy(pool);
    for (int i = 0; i < N; i++) { free(seqs[i]); free(ref[i]); free(dm[i]); }
    free(seqs); free(ref); free(dm);
    return 0;
}

/* ── 3. fork-join ─────────────────────────────────────────────── */

struct fj_arg {
    int *out;
    int  val;
};

static void fj_write(void *arg)
{
    struct fj_arg *a = (struct fj_arg *)arg;
    *a->out = a->val;
}

static int test_fork_join_two(void)
{
    threadpool_t *pool = tp_create(2);
    int a = 0, b = 0;
    struct fj_arg aa = { &a, 42 };
    struct fj_arg bb = { &b, 99 };

    tp_group_t *g = tp_group_create(pool);
    tp_group_submit(g, fj_write, &aa);
    tp_group_submit(g, fj_write, &bb);
    tp_group_wait(g);
    tp_group_destroy(g);

    CHECK(a == 42);
    CHECK(b == 99);
    tp_destroy(pool);
    return 0;
}

static int test_fork_join_four(void)
{
    threadpool_t *pool = tp_create(4);
    int results[4] = {0};
    struct fj_arg args[4];
    for (int i = 0; i < 4; i++) {
        args[i].out = &results[i];
        args[i].val = (i + 1) * 10;
    }

    tp_group_t *g = tp_group_create(pool);
    for (int i = 0; i < 4; i++)
        tp_group_submit(g, fj_write, &args[i]);
    tp_group_wait(g);
    tp_group_destroy(g);

    for (int i = 0; i < 4; i++)
        CHECK(results[i] == (i + 1) * 10);

    tp_destroy(pool);
    return 0;
}

/* ── 4. recursive tasks (parallel fibonacci) ──────────────────── */

struct fib_arg {
    threadpool_t *pool;
    int n;
    int result;
};

static void fib_task(void *arg)
{
    struct fib_arg *f = (struct fib_arg *)arg;
    if (f->n <= 1) {
        f->result = f->n;
        return;
    }
    struct fib_arg left  = { .pool = f->pool, .n = f->n - 1, .result = 0 };
    struct fib_arg right = { .pool = f->pool, .n = f->n - 2, .result = 0 };

    tp_group_t *g = tp_group_create(f->pool);
    tp_group_submit(g, fib_task, &left);
    tp_group_submit(g, fib_task, &right);
    tp_group_wait(g);
    tp_group_destroy(g);

    f->result = left.result + right.result;
}

static int serial_fib(int n)
{
    if (n <= 1) return n;
    return serial_fib(n - 1) + serial_fib(n - 2);
}

static int test_recursive_fib(void)
{
    threadpool_t *pool = tp_create(4);

    for (int n = 0; n <= 15; n++) {
        struct fib_arg f = { .pool = pool, .n = n, .result = -1 };
        fib_task(&f);
        CHECK(f.result == serial_fib(n));
    }

    tp_destroy(pool);
    return 0;
}

/* ── 5. recursive tree (simulates guide-tree traversal) ───────── */

struct tree_node {
    struct tree_node *left;
    struct tree_node *right;
    int               leaf_id;   /* >= 0 for leaves, -1 for internal */
    float             result;
};

static struct tree_node *build_tree(int depth, int *next_id)
{
    struct tree_node *n = (struct tree_node *)calloc(1, sizeof(*n));
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

struct tree_task_arg {
    threadpool_t    *pool;
    struct tree_node *node;
};

static float leaf_work(int id)
{
    float sum = 0.0f;
    unsigned x = (unsigned)id * 2654435761u;
    for (int i = 0; i < 100; i++) {
        x = x * 1103515245u + 12345u;
        sum += (float)(x >> 16) / 65536.0f;
    }
    return sum;
}

static void tree_task(void *arg)
{
    struct tree_task_arg *ta = (struct tree_task_arg *)arg;
    struct tree_node *n = ta->node;

    if (n->leaf_id >= 0) {
        n->result = leaf_work(n->leaf_id);
        return;
    }

    struct tree_task_arg left_arg  = { ta->pool, n->left };
    struct tree_task_arg right_arg = { ta->pool, n->right };

    tp_group_t *g = tp_group_create(ta->pool);
    tp_group_submit(g, tree_task, &left_arg);
    tp_group_submit(g, tree_task, &right_arg);
    tp_group_wait(g);
    tp_group_destroy(g);

    n->result = n->left->result + n->right->result;
}

static float serial_tree(struct tree_node *n)
{
    if (n->leaf_id >= 0)
        return leaf_work(n->leaf_id);
    return serial_tree(n->left) + serial_tree(n->right);
}

static int test_recursive_tree(void)
{
    int id = 0;
    struct tree_node *root = build_tree(10, &id);  /* 1024 leaves */

    float serial_result = serial_tree(root);

    /* Reset results */
    /* (serial_tree wrote into the same nodes — rebuild) */
    free_tree(root);
    id = 0;
    root = build_tree(10, &id);

    threadpool_t *pool = tp_create(4);
    struct tree_task_arg targ = { pool, root };
    tree_task(&targ);

    /* Compare within float tolerance */
    float diff = root->result - serial_result;
    if (diff < 0) diff = -diff;
    CHECK(diff < 0.001f);

    tp_destroy(pool);
    free_tree(root);
    return 0;
}

/* ── 6. many small tasks (stress test) ────────────────────────── */

static void increment(void *arg)
{
    atomic_int *counter = (atomic_int *)arg;
    atomic_fetch_add(counter, 1);
}

static int test_many_tasks(void)
{
    int N = 10000;
    threadpool_t *pool = tp_create(4);
    atomic_int counter;
    atomic_init(&counter, 0);

    tp_group_t *g = tp_group_create(pool);
    for (int i = 0; i < N; i++)
        tp_group_submit(g, increment, &counter);
    tp_group_wait(g);
    tp_group_destroy(g);

    CHECK(atomic_load(&counter) == N);
    tp_destroy(pool);
    return 0;
}

/* ── 7. nested groups (non-recursive) ─────────────────────────── */

struct nested_arg {
    threadpool_t *pool;
    atomic_int   *counter;
};

static void inner_task(void *arg)
{
    atomic_int *c = (atomic_int *)arg;
    atomic_fetch_add(c, 1);
}

static void outer_task(void *arg)
{
    struct nested_arg *na = (struct nested_arg *)arg;
    tp_group_t *g = tp_group_create(na->pool);
    for (int i = 0; i < 5; i++)
        tp_group_submit(g, inner_task, na->counter);
    tp_group_wait(g);
    tp_group_destroy(g);
}

static int test_nested_groups(void)
{
    threadpool_t *pool = tp_create(4);
    atomic_int counter;
    atomic_init(&counter, 0);

    struct nested_arg na = { pool, &counter };

    tp_group_t *g = tp_group_create(pool);
    for (int i = 0; i < 10; i++)
        tp_group_submit(g, outer_task, &na);
    tp_group_wait(g);
    tp_group_destroy(g);

    /* 10 outer * 5 inner = 50 */
    CHECK(atomic_load(&counter) == 50);
    tp_destroy(pool);
    return 0;
}

/* ── 8. group wait with no tasks ──────────────────────────────── */

static int test_empty_group_wait(void)
{
    threadpool_t *pool = tp_create(2);
    tp_group_t *g = tp_group_create(pool);
    tp_group_wait(g);   /* should return immediately */
    tp_group_destroy(g);
    tp_destroy(pool);
    return 0;
}

/* ── 9. large parallel for ────────────────────────────────────── */

static int test_large_parallel_for(void)
{
    /* N must be ≤ 46340 to avoid signed int overflow in i*i
     * (sqrt(INT_MAX) ≈ 46340).  The compiler exploits overflow UB. */
    int N = 46340;
    int *arr = (int *)calloc((size_t)N, sizeof(int));
    CHECK(arr != NULL);

    threadpool_t *pool = tp_create(8);
    tp_parallel_for(pool, 0, N, fill_squares, arr);

    for (int i = 0; i < N; i++)
        CHECK(arr[i] == i * i);

    tp_destroy(pool);
    free(arr);
    return 0;
}

/* ── 10. deep recursive tree (V2 stack-safety proof) ──────────── */

/* Depth 16 = 65536 leaves.  V1 (global FIFO queue) would need
 * O(65536) stack frames in the worst case, overflowing 8 MB.
 * V2 (Chase-Lev LIFO) needs O(16) frames per thread. */
static int test_deep_recursive_tree(void)
{
    int id = 0;
    struct tree_node *root = build_tree(16, &id);  /* 65536 leaves */

    /* Serial reference */
    float serial_result = serial_tree(root);
    free_tree(root);
    id = 0;
    root = build_tree(16, &id);

    threadpool_t *pool = tp_create(4);
    struct tree_task_arg targ = { pool, root };
    tree_task(&targ);

    float diff = root->result - serial_result;
    if (diff < 0) diff = -diff;
    CHECK(diff < 0.01f);

    tp_destroy(pool);
    free_tree(root);
    return 0;
}

/* ── 11. deep recursive fibonacci ─────────────────────────────── */

static int test_deep_recursive_fib(void)
{
    threadpool_t *pool = tp_create(4);

    /* fib(20) = 6765 — creates ~13K groups, ~26K tasks. */
    struct fib_arg f = { .pool = pool, .n = 20, .result = -1 };
    fib_task(&f);
    CHECK(f.result == 6765);

    tp_destroy(pool);
    return 0;
}

/* ── 12. multiple sequential groups ───────────────────────────── */

static int test_sequential_groups(void)
{
    threadpool_t *pool = tp_create(4);

    for (int round = 0; round < 20; round++) {
        atomic_int counter;
        atomic_init(&counter, 0);

        tp_group_t *g = tp_group_create(pool);
        for (int i = 0; i < 100; i++)
            tp_group_submit(g, increment, &counter);
        tp_group_wait(g);
        tp_group_destroy(g);

        CHECK(atomic_load(&counter) == 100);
    }

    tp_destroy(pool);
    return 0;
}

/* ── main ─────────────────────────────────────────────────────── */

int main(void)
{
    printf("threadpool tests\n");
    printf("================\n\n");

    RUN(create_destroy);
    RUN(create_auto);
    RUN(parallel_for_correctness);
    RUN(parallel_for_single_thread);
    RUN(parallel_for_empty_range);
    RUN(parallel_for_distance_matrix);
    RUN(fork_join_two);
    RUN(fork_join_four);
    RUN(recursive_fib);
    RUN(recursive_tree);
    RUN(many_tasks);
    RUN(nested_groups);
    RUN(empty_group_wait);
    RUN(large_parallel_for);
    RUN(deep_recursive_tree);
    RUN(deep_recursive_fib);
    RUN(sequential_groups);

    printf("\n%d/%d passed", g_pass, g_total);
    if (g_fail > 0) printf(", %d FAILED", g_fail);
    printf("\n");

    return g_fail > 0 ? 1 : 0;
}
