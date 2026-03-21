/*
 * stress_threadpool.c — Stress tests for production readiness.
 *
 * Covers scenarios the unit tests don't:
 *   1. High contention (many submitters, few workers)
 *   2. Rapid pool create/destroy cycles
 *   3. Deque overflow → ext queue fallback
 *   4. Oversubscription (more workers than cores)
 *   5. Sustained mixed-pattern load
 *   6. Shutdown during active work (tp_request_shutdown)
 *   7. Correctness under sustained parallel-for pressure
 */

#include "threadpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdatomic.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* ── TSan detection ──────────────────────────────────────────────
 * TSan can't trace happens-before through the threadpool's fence-
 * based synchronization to user data.  We add explicit annotations
 * for tree-structured tests where the main thread writes data that
 * workers later read through a chain of task submissions. */
#if defined(__SANITIZE_THREAD__)
#define STRESS_TSAN 1
#elif defined(__has_feature)
#if __has_feature(thread_sanitizer)
#define STRESS_TSAN 1
#endif
#endif

#ifdef STRESS_TSAN
void __tsan_acquire(void *addr);
void __tsan_release(void *addr);
#define TSAN_RELEASE(addr) __tsan_release(addr)
#define TSAN_ACQUIRE(addr) __tsan_acquire(addr)
#else
#define TSAN_RELEASE(addr) ((void)0)
#define TSAN_ACQUIRE(addr) ((void)0)
#endif

/* ── Test harness ─────────────────────────────────────────────── */

static int g_pass, g_fail, g_total;

#define CHECK(cond) do {                                              \
    if (!(cond)) {                                                    \
        fprintf(stderr, "    FAIL: %s  (line %d)\n",                  \
                #cond, __LINE__);                                     \
        return -1;                                                    \
    }                                                                 \
} while (0)

#define RUN(name) do {                                                \
    printf("  %-50s ", #name);                                        \
    fflush(stdout);                                                   \
    g_total++;                                                        \
    if (test_##name() == 0) { printf("PASS\n"); g_pass++; }          \
    else                    { printf("FAIL\n"); g_fail++; }           \
} while (0)

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ── 1. High contention: many tasks from main thread ─────────── */
/* All tasks go through ext queue since the main thread is not a
 * worker.  Stresses ext queue locking under contention. */

static void atomic_inc(void *arg)
{
    atomic_int *c = (atomic_int *)arg;
    atomic_fetch_add(c, 1);
}

static int test_high_contention_ext_queue(void)
{
    int N = 100000;
    threadpool_t *pool = tp_create(8);
    CHECK(pool != NULL);

    atomic_int counter;
    atomic_init(&counter, 0);

    tp_group_t *g = tp_group_create(pool);
    CHECK(g != NULL);
    for (int i = 0; i < N; i++) {
        int rc = tp_group_submit(g, atomic_inc, &counter);
        CHECK(rc == 0);
    }
    tp_group_wait(g);
    tp_group_destroy(g);

    CHECK(atomic_load(&counter) == N);
    tp_destroy(pool);
    return 0;
}

/* ── 2. Rapid pool create/destroy cycles ─────────────────────── */
/* Catches resource leaks (threads, mutexes, memory). */

static int test_rapid_create_destroy(void)
{
    for (int i = 0; i < 200; i++) {
        threadpool_t *pool = tp_create(4);
        CHECK(pool != NULL);

        /* Do a tiny bit of work each time to exercise the full path. */
        atomic_int counter;
        atomic_init(&counter, 0);
        tp_group_t *g = tp_group_create(pool);
        CHECK(g != NULL);
        tp_group_submit(g, atomic_inc, &counter);
        tp_group_wait(g);
        tp_group_destroy(g);
        CHECK(atomic_load(&counter) == 1);

        tp_destroy(pool);
    }
    return 0;
}

/* ── 3. Deque overflow → ext queue fallback ──────────────────── */
/* A worker task submits > DEQUE_DEFAULT_CAP (4096) children.
 * Once the deque fills, remaining tasks must go to the ext queue.
 * All tasks must still execute correctly. */

struct burst_ctx {
    threadpool_t *pool;
    atomic_int   *counter;
    int           burst_size;
};

static void burst_parent(void *arg)
{
    struct burst_ctx *ctx = (struct burst_ctx *)arg;
    tp_group_t *g = tp_group_create(ctx->pool);
    if (!g) return;

    for (int i = 0; i < ctx->burst_size; i++) {
        tp_group_submit(g, atomic_inc, ctx->counter);
    }
    tp_group_wait(g);
    tp_group_destroy(g);
}

static int test_deque_overflow_fallback(void)
{
    threadpool_t *pool = tp_create(2);
    CHECK(pool != NULL);

    atomic_int counter;
    atomic_init(&counter, 0);

    /* Submit 8000 tasks from a worker context — exceeds deque cap of 4096. */
    struct burst_ctx ctx = { pool, &counter, 8000 };

    tp_group_t *g = tp_group_create(pool);
    CHECK(g != NULL);
    tp_group_submit(g, burst_parent, &ctx);
    tp_group_wait(g);
    tp_group_destroy(g);

    CHECK(atomic_load(&counter) == 8000);
    tp_destroy(pool);
    return 0;
}

/* ── 4. Oversubscription ─────────────────────────────────────── */
/* More workers than CPU cores.  Must not deadlock or starve. */

static void fill_val(int start, int end, void *arg)
{
    int *a = (int *)arg;
    for (int i = start; i < end; i++)
        a[i] = i + 1;
}

static int test_oversubscription(void)
{
    int ncpus = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpus <= 0) ncpus = 4;
    int nthreads = ncpus * 4;  /* 4x oversubscription */

    threadpool_t *pool = tp_create(nthreads);
    CHECK(pool != NULL);
    CHECK(tp_get_nthreads(pool) == nthreads);

    /* Parallel for with many chunks. */
    int N = 10000;
    int *arr = (int *)calloc((size_t)N, sizeof(int));
    CHECK(arr != NULL);

    tp_parallel_for(pool, 0, N, fill_val, arr);

    for (int i = 0; i < N; i++)
        CHECK(arr[i] == i + 1);

    free(arr);
    tp_destroy(pool);
    return 0;
}

/* ── 5. Sustained mixed-pattern load ─────────────────────────── */
/* Run parallel-for + recursive tasks + fork-join concurrently
 * on the same pool over many iterations. */

struct tree_node_s {
    struct tree_node_s *left;
    struct tree_node_s *right;
    int    leaf_id;
    atomic_int result;
};

static struct tree_node_s *build_small_tree(int depth, int *next_id)
{
    struct tree_node_s *n = (struct tree_node_s *)calloc(1, sizeof(*n));
    if (!n) return NULL;
    atomic_init(&n->result, 0);
    if (depth == 0) {
        n->leaf_id = (*next_id)++;
        return n;
    }
    n->leaf_id = -1;
    n->left  = build_small_tree(depth - 1, next_id);
    n->right = build_small_tree(depth - 1, next_id);
    return n;
}

static void free_small_tree(struct tree_node_s *n)
{
    if (!n) return;
    free_small_tree(n->left);
    free_small_tree(n->right);
    free(n);
}

struct stree_arg {
    threadpool_t     *pool;
    struct tree_node_s *node;
    struct tree_node_s *root; /* sync point for TSan */
};

static void stree_task(void *arg)
{
    struct stree_arg *ta = (struct stree_arg *)arg;
    struct tree_node_s *n = ta->node;
    TSAN_ACQUIRE(ta->root); /* pairs with TSAN_RELEASE after tree build */

    if (n->leaf_id >= 0) {
        atomic_store(&n->result, n->leaf_id + 1);
        return;
    }

    struct stree_arg left  = { ta->pool, n->left,  ta->root };
    struct stree_arg right = { ta->pool, n->right, ta->root };

    tp_group_t *g = tp_group_create(ta->pool);
    if (!g) return;
    tp_group_submit(g, stree_task, &left);
    tp_group_submit(g, stree_task, &right);
    tp_group_wait(g);
    tp_group_destroy(g);

    atomic_store(&n->result,
                 atomic_load(&n->left->result) +
                 atomic_load(&n->right->result));
}

static int serial_tree_sum(struct tree_node_s *n)
{
    if (n->leaf_id >= 0) return n->leaf_id + 1;
    return serial_tree_sum(n->left) + serial_tree_sum(n->right);
}

static void pfor_fill(int start, int end, void *arg)
{
    int *a = (int *)arg;
    for (int i = start; i < end; i++)
        a[i] = i * i;
}

static int test_sustained_mixed_patterns(void)
{
    threadpool_t *pool = tp_create(8);
    CHECK(pool != NULL);

    int N = 5000;
    int *arr = (int *)calloc((size_t)N, sizeof(int));
    CHECK(arr != NULL);

    for (int iter = 0; iter < 50; iter++) {
        /* Parallel for */
        memset(arr, 0, (size_t)N * sizeof(int));
        tp_parallel_for(pool, 0, N, pfor_fill, arr);
        for (int i = 0; i < N; i++)
            CHECK(arr[i] == i * i);

        /* Recursive tree */
        int id = 0;
        struct tree_node_s *root = build_small_tree(8, &id);  /* 256 leaves */
        CHECK(root != NULL);
        int expected = serial_tree_sum(root);

        TSAN_RELEASE(root); /* tree is fully built; pairs with ACQUIRE in stree_task */
        struct stree_arg targ = { pool, root, root };
        tp_group_t *tg = tp_group_create(pool);
        CHECK(tg != NULL);
        tp_group_submit(tg, stree_task, &targ);
        tp_group_wait(tg);
        tp_group_destroy(tg);
        CHECK(atomic_load(&root->result) == expected);
        free_small_tree(root);

        /* Fork-join */
        atomic_int counter;
        atomic_init(&counter, 0);
        tp_group_t *g = tp_group_create(pool);
        CHECK(g != NULL);
        for (int i = 0; i < 100; i++)
            tp_group_submit(g, atomic_inc, &counter);
        tp_group_wait(g);
        tp_group_destroy(g);
        CHECK(atomic_load(&counter) == 100);
    }

    free(arr);
    tp_destroy(pool);
    return 0;
}

/* ── 6. tp_request_shutdown during active work ───────────────── */
/* Verify that tp_group_wait returns after shutdown is requested,
 * and that tp_destroy completes without hanging. */

static void slow_task(void *arg)
{
    atomic_int *counter = (atomic_int *)arg;
    /* Simulate some work. */
    volatile int sink = 0;
    for (int i = 0; i < 10000; i++)
        sink += i;
    (void)sink;
    atomic_fetch_add(counter, 1);
}

static int test_shutdown_during_work(void)
{
    threadpool_t *pool = tp_create(4);
    CHECK(pool != NULL);

    atomic_int counter;
    atomic_init(&counter, 0);

    /* Submit a batch of tasks. */
    tp_group_t *g = tp_group_create(pool);
    CHECK(g != NULL);
    for (int i = 0; i < 1000; i++)
        tp_group_submit(g, slow_task, &counter);

    /* Request shutdown immediately — some tasks may not finish. */
    tp_request_shutdown(pool);

    /* tp_group_wait should return (possibly early). */
    tp_group_wait(g);
    tp_group_destroy(g);

    /* tp_destroy should not hang. */
    tp_destroy(pool);

    /* At least some tasks should have run. */
    int completed = atomic_load(&counter);
    CHECK(completed > 0);

    return 0;
}

/* ── 7. Parallel-for correctness under sustained pressure ────── */
/* Many iterations of parallel-for on the same pool, verifying
 * every element every time. */

static void fill_with_check(int start, int end, void *arg)
{
    int *a = (int *)arg;
    for (int i = start; i < end; i++)
        a[i] = i * 3 + 7;
}

static int test_parallel_for_sustained(void)
{
    threadpool_t *pool = tp_create(8);
    CHECK(pool != NULL);

    int N = 50000;
    int *arr = (int *)calloc((size_t)N, sizeof(int));
    CHECK(arr != NULL);

    for (int iter = 0; iter < 100; iter++) {
        memset(arr, 0, (size_t)N * sizeof(int));
        tp_parallel_for(pool, 0, N, fill_with_check, arr);
        for (int i = 0; i < N; i++)
            CHECK(arr[i] == i * 3 + 7);
    }

    free(arr);
    tp_destroy(pool);
    return 0;
}

/* ── 8. Single-element ranges and edge cases ─────────────────── */

static void count_calls(int start, int end, void *arg)
{
    atomic_int *c = (atomic_int *)arg;
    for (int i = start; i < end; i++)
        atomic_fetch_add(c, 1);
}

static int test_edge_cases(void)
{
    threadpool_t *pool = tp_create(4);
    CHECK(pool != NULL);

    /* Single element */
    atomic_int counter;
    atomic_init(&counter, 0);
    tp_parallel_for(pool, 0, 1, count_calls, &counter);
    CHECK(atomic_load(&counter) == 1);

    /* Two elements */
    atomic_init(&counter, 0);
    tp_parallel_for(pool, 0, 2, count_calls, &counter);
    CHECK(atomic_load(&counter) == 2);

    /* Range smaller than thread count */
    atomic_init(&counter, 0);
    tp_parallel_for(pool, 0, 3, count_calls, &counter);
    CHECK(atomic_load(&counter) == 3);

    /* Range equal to thread count */
    atomic_init(&counter, 0);
    tp_parallel_for(pool, 0, 4, count_calls, &counter);
    CHECK(atomic_load(&counter) == 4);

    /* Negative/empty ranges */
    atomic_init(&counter, 0);
    tp_parallel_for(pool, 5, 5, count_calls, &counter);
    CHECK(atomic_load(&counter) == 0);
    tp_parallel_for(pool, 10, 3, count_calls, &counter);
    CHECK(atomic_load(&counter) == 0);

    /* NULL pool is UB, so we don't test it. */

    tp_destroy(pool);
    return 0;
}

/* ── 9. Many sequential groups (lifecycle churn) ─────────────── */

static int test_group_lifecycle_churn(void)
{
    threadpool_t *pool = tp_create(4);
    CHECK(pool != NULL);

    for (int i = 0; i < 10000; i++) {
        atomic_int counter;
        atomic_init(&counter, 0);

        tp_group_t *g = tp_group_create(pool);
        CHECK(g != NULL);
        tp_group_submit(g, atomic_inc, &counter);
        tp_group_wait(g);
        tp_group_destroy(g);

        CHECK(atomic_load(&counter) == 1);
    }

    tp_destroy(pool);
    return 0;
}

/* ── 10. Deeply nested recursive tasks ───────────────────────── */
/* Depth 20 = 1M leaves.  LIFO pop bounds stack to O(20) per
 * thread.  Would overflow with a naive FIFO approach. */

struct deep_node {
    struct deep_node *left;
    struct deep_node *right;
    int    leaf_id;
    int    result;
};

static struct deep_node *build_deep_tree(int depth, int *next_id)
{
    struct deep_node *n = (struct deep_node *)calloc(1, sizeof(*n));
    if (!n) return NULL;
    if (depth == 0) {
        n->leaf_id = (*next_id)++;
        return n;
    }
    n->leaf_id = -1;
    n->left  = build_deep_tree(depth - 1, next_id);
    n->right = build_deep_tree(depth - 1, next_id);
    return n;
}

static void free_deep_tree(struct deep_node *n)
{
    if (!n) return;
    free_deep_tree(n->left);
    free_deep_tree(n->right);
    free(n);
}

struct deep_arg {
    threadpool_t   *pool;
    struct deep_node *node;
    int             cutoff;
};

static int serial_deep_sum(struct deep_node *n)
{
    if (n->leaf_id >= 0) return 1;
    return serial_deep_sum(n->left) + serial_deep_sum(n->right);
}

static void deep_serial_task(struct deep_node *n)
{
    if (n->leaf_id >= 0) { n->result = 1; return; }
    deep_serial_task(n->left);
    deep_serial_task(n->right);
    n->result = n->left->result + n->right->result;
}

static void deep_task(void *arg)
{
    struct deep_arg *da = (struct deep_arg *)arg;
    struct deep_node *n = da->node;

    if (n->leaf_id >= 0) { n->result = 1; return; }

    if (da->cutoff <= 0) {
        /* Below cutoff: run serially to avoid task-creation overhead. */
        deep_serial_task(n);
        return;
    }

    struct deep_arg left  = { da->pool, n->left,  da->cutoff - 1 };
    struct deep_arg right = { da->pool, n->right, da->cutoff - 1 };

    tp_group_t *g = tp_group_create(da->pool);
    if (!g) { deep_serial_task(n); return; }
    tp_group_submit(g, deep_task, &left);
    tp_group_submit(g, deep_task, &right);
    tp_group_wait(g);
    tp_group_destroy(g);

    n->result = n->left->result + n->right->result;
}

static int test_deep_tree_stress(void)
{
#ifdef STRESS_TSAN
    /* TSan inflates stack frames ~5x; reduce depth to avoid overflow
     * during recursive tp_group_wait → execute_task → deep_task chains. */
    int depth = 14;
    int cutoff = 5;
#else
    int depth = 20;  /* 1,048,576 leaves */
    int cutoff = 10;
#endif
    int id = 0;
    struct deep_node *root = build_deep_tree(depth, &id);
    CHECK(root != NULL);
    CHECK(id == (1 << depth));

    int expected = serial_deep_sum(root);
    CHECK(expected == (1 << depth));

    threadpool_t *pool = tp_create(8);
    CHECK(pool != NULL);

    struct deep_arg da = { pool, root, cutoff };
    deep_task(&da);

    CHECK(root->result == expected);

    tp_destroy(pool);
    free_deep_tree(root);
    return 0;
}

/* ── main ─────────────────────────────────────────────────────── */

int main(void)
{
    printf("threadpool stress tests\n");
    printf("=======================\n\n");

    double t0 = now_sec();

    RUN(high_contention_ext_queue);
    RUN(rapid_create_destroy);
    RUN(deque_overflow_fallback);
    RUN(oversubscription);
    RUN(sustained_mixed_patterns);
    RUN(shutdown_during_work);
    RUN(parallel_for_sustained);
    RUN(edge_cases);
    RUN(group_lifecycle_churn);
    RUN(deep_tree_stress);

    double elapsed = now_sec() - t0;
    printf("\n%d/%d passed (%.1fs)", g_pass, g_total, elapsed);
    if (g_fail > 0) printf(", %d FAILED", g_fail);
    printf("\n");

    return g_fail > 0 ? 1 : 0;
}
