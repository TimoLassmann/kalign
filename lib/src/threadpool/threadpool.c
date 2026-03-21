/*
 * threadpool.c — Chase-Lev work-stealing thread pool.
 *
 * See threadpool.h for the public API and usage examples.
 *
 * Internal architecture:
 *
 *   Queues:
 *   - Per-worker Chase-Lev deques (lock-free, LIFO pop, FIFO steal).
 *   - Global "external" queue (mutex-protected) for non-worker submissions.
 *
 *   Work-finding priority (per worker):
 *     1. Own deque (LIFO) — preserves DFS order, bounds stack depth.
 *     2. External queue   — picks up tasks submitted from non-workers.
 *     3. Steal from peer  — random victim, FIFO steals shallowest work.
 *
 *   Sleeping:
 *     Event-count protocol (wake_gen + condvar).  Workers spin briefly,
 *     then park.  Submitters bump wake_gen and signal if anyone sleeps.
 *
 *   Group recycling:
 *     Per-worker free lists for tp_group_t objects.  Groups are 16 bytes
 *     and follow strict LIFO create/destroy ordering, so thread-local
 *     recycling eliminates virtually all malloc/free after warmup.
 *
 *   Stack depth:
 *     LIFO pop means a recursive tree traversal nests O(tree_depth) wait
 *     frames per thread, not O(tree_size).  For a balanced tree of 1M
 *     leaves (depth ~20), each frame is ~500 bytes → ~10 KB total.
 *
 * Ref: Chase & Lev, "Dynamic Circular Work-Stealing Deque", SPAA 2005.
 */

#include "threadpool.h"

#include <pthread.h>
#include <stdatomic.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sched.h>

/* ── TSan annotations for Chase-Lev deque ────────────────────────
 * TSan can't track happens-before through atomic_thread_fence,
 * so we annotate the deque push/pop/steal with explicit acquire/
 * release on the buffer slot.  Zero cost when TSan is off. */
#if defined(__SANITIZE_THREAD__)
#define TP_TSAN 1
#elif defined(__has_feature)
#if __has_feature(thread_sanitizer)
#define TP_TSAN 1
#endif
#endif

#ifdef TP_TSAN
void __tsan_acquire(void *addr);
void __tsan_release(void *addr);
#define TSAN_RELEASE(addr) __tsan_release(addr)
#define TSAN_ACQUIRE(addr) __tsan_acquire(addr)
#else
#define TSAN_RELEASE(addr) ((void)0)
#define TSAN_ACQUIRE(addr) ((void)0)
#endif

/* Abort on pthread errors that indicate programming bugs (EINVAL,
 * EDEADLK, etc.).  These never fail in correct code, but catching
 * them immediately beats silent corruption. */
#define PTHREAD_CHECK(call) do {                                      \
    int _rc = (call);                                                 \
    if (_rc != 0) {                                                   \
        fprintf(stderr, "threadpool: %s failed (%d)\n", #call, _rc); \
        abort();                                                      \
    }                                                                 \
} while (0)

/* ── Internal task ────────────────────────────────────────────── */

typedef struct {
    void       (*fn)(void *);
    void        *arg;
    atomic_int  *pending;   /* &group->pending, or NULL */
} tp_task_t;

/* ══════════════════════════════════════════════════════════════════
 * Chase-Lev work-stealing deque
 *
 * Owner pushes/pops at bottom (LIFO).  Thieves steal from top (FIFO).
 * Lock-free: owner ops use fences, steal uses CAS on top.
 * Fixed capacity (power of 2).
 * ══════════════════════════════════════════════════════════════════ */

#define DEQUE_DEFAULT_CAP 4096   /* per worker; holds O(tree_depth) tasks */

typedef struct {
    tp_task_t  *buf;
    long        cap;            /* power of 2 */
    atomic_long bottom;         /* modified by owner              */
    char        _pad[64];       /* keep bottom and top on separate cache lines */
    atomic_long top;            /* CAS'd by thieves              */
} ws_deque_t;

static int deque_init(ws_deque_t *dq, long cap)
{
    dq->buf = (tp_task_t *)calloc((size_t)cap, sizeof(tp_task_t));
    if (!dq->buf) return -1;
    dq->cap = cap;
    atomic_store_explicit(&dq->bottom, 0, memory_order_relaxed);
    atomic_store_explicit(&dq->top, 0, memory_order_relaxed);
    return 0;
}

static void deque_destroy(ws_deque_t *dq)
{
    free(dq->buf);
}

/* Owner push.  Returns 0 on success, -1 if full. */
static int deque_push(ws_deque_t *dq, tp_task_t task)
{
    long b = atomic_load_explicit(&dq->bottom, memory_order_relaxed);
    long t = atomic_load_explicit(&dq->top,    memory_order_acquire);
    if (b - t >= dq->cap)
        return -1;                          /* full */
    long slot = b & (dq->cap - 1);
    dq->buf[slot] = task;
    TSAN_RELEASE(&dq->buf[slot]);
    atomic_thread_fence(memory_order_release); /* task visible before bottom bumps */
    atomic_store_explicit(&dq->bottom, b + 1, memory_order_relaxed);
    return 0;
}

/* Owner pop (LIFO).  Returns 0 on success, -1 if empty. */
static int deque_pop(ws_deque_t *dq, tp_task_t *out)
{
    long b = atomic_load_explicit(&dq->bottom, memory_order_relaxed) - 1;
    atomic_store_explicit(&dq->bottom, b, memory_order_relaxed);
    atomic_thread_fence(memory_order_seq_cst);
    long t = atomic_load_explicit(&dq->top, memory_order_relaxed);

    if (t <= b) {
        long slot = b & (dq->cap - 1);
        *out = dq->buf[slot];
        TSAN_ACQUIRE(&dq->buf[slot]);
        if (t == b) {
            /* Last element — race with a stealer. */
            if (!atomic_compare_exchange_strong_explicit(
                    &dq->top, &t, t + 1,
                    memory_order_seq_cst, memory_order_relaxed)) {
                /* Stealer won. */
                atomic_store_explicit(&dq->bottom, b + 1, memory_order_relaxed);
                return -1;
            }
            atomic_store_explicit(&dq->bottom, b + 1, memory_order_relaxed);
        }
        return 0;
    }
    /* Was already empty. */
    atomic_store_explicit(&dq->bottom, b + 1, memory_order_relaxed);
    return -1;
}

/* Thief steal (FIFO).  Returns 0 on success, -1 if empty or contention. */
static int deque_steal(ws_deque_t *dq, tp_task_t *out)
{
    long t = atomic_load_explicit(&dq->top,    memory_order_acquire);
    atomic_thread_fence(memory_order_seq_cst);
    long b = atomic_load_explicit(&dq->bottom, memory_order_acquire);

    if (t < b) {
        long slot = t & (dq->cap - 1);
        *out = dq->buf[slot];
        if (!atomic_compare_exchange_strong_explicit(
                &dq->top, &t, t + 1,
                memory_order_seq_cst, memory_order_relaxed))
            return -1;                      /* another stealer won */
        TSAN_ACQUIRE(&dq->buf[slot]);
        return 0;
    }
    return -1;                              /* empty */
}

/* ══════════════════════════════════════════════════════════════════
 * External queue — mutex-protected, for non-worker submissions
 * ══════════════════════════════════════════════════════════════════ */

typedef struct {
    tp_task_t      *buf;
    int             cap;
    int             head;
    int             tail;
    int             count;
    pthread_mutex_t lock;
} ext_queue_t;

static int ext_init(ext_queue_t *q, int cap)
{
    q->buf = (tp_task_t *)calloc((size_t)cap, sizeof(tp_task_t));
    if (!q->buf) return -1;
    q->cap = cap;
    q->head = q->tail = q->count = 0;
    if (pthread_mutex_init(&q->lock, NULL) != 0) {
        free(q->buf);
        q->buf = NULL;
        return -1;
    }
    return 0;
}

static void ext_destroy(ext_queue_t *q)
{
    free(q->buf);
    PTHREAD_CHECK(pthread_mutex_destroy(&q->lock));
}

/* Caller must hold q->lock. */
static int ext_grow_locked(ext_queue_t *q)
{
    if (q->cap > INT_MAX / 2) return -1;  /* overflow guard */
    int new_cap = q->cap * 2;
    tp_task_t *nb = (tp_task_t *)calloc((size_t)new_cap, sizeof(tp_task_t));
    if (!nb) return -1;
    for (int i = 0; i < q->count; i++)
        nb[i] = q->buf[(q->head + i) % q->cap];
    free(q->buf);
    q->buf  = nb;
    q->head = 0;
    q->tail = q->count;
    q->cap  = new_cap;
    return 0;
}

static int ext_push(ext_queue_t *q, tp_task_t task)
{
    PTHREAD_CHECK(pthread_mutex_lock(&q->lock));
    if (q->count == q->cap && ext_grow_locked(q) != 0) {
        PTHREAD_CHECK(pthread_mutex_unlock(&q->lock));
        return -1;
    }
    q->buf[q->tail] = task;
    q->tail = (q->tail + 1) % q->cap;
    q->count++;
    PTHREAD_CHECK(pthread_mutex_unlock(&q->lock));
    return 0;
}

static int ext_try_pop(ext_queue_t *q, tp_task_t *out)
{
    PTHREAD_CHECK(pthread_mutex_lock(&q->lock));
    if (q->count == 0) {
        PTHREAD_CHECK(pthread_mutex_unlock(&q->lock));
        return -1;
    }
    *out = q->buf[q->head];
    q->head = (q->head + 1) % q->cap;
    q->count--;
    PTHREAD_CHECK(pthread_mutex_unlock(&q->lock));
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * Pool, worker, group structures
 * ══════════════════════════════════════════════════════════════════ */

#define MAX_FREE_GROUPS 64

typedef struct worker {
    ws_deque_t         deque;
    struct threadpool *pool;
    pthread_t          thread;
    int                id;
    unsigned           rng;     /* xorshift state for random steal */
    struct tp_group   *free_groups; /* per-worker group free list */
    int                n_free;
} worker_t;

struct threadpool {
    worker_t       *workers;
    int             nworkers;
    ext_queue_t     ext;

    atomic_int      wake_gen;       /* event count for wakeup       */
    atomic_int      n_sleeping;
    pthread_mutex_t wake_lock;
    pthread_cond_t  wake_cond;

    atomic_int      shutdown;
};

struct tp_group {
    union {
        threadpool_t   *pool;  /* valid when active */
        struct tp_group *next; /* valid when on free list */
    };
    atomic_int pending;
};

/* ── Thread-local: current worker (NULL on non-worker threads) ── */
static _Thread_local worker_t *tl_worker = NULL;

/* ── RNG for steal-target selection ───────────────────────────── */

static unsigned xorshift32(unsigned *s)
{
    unsigned x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return *s = x;
}

/* ── Task execution ───────────────────────────────────────────── */

static void execute_task(tp_task_t *t)
{
    t->fn(t->arg);
    if (t->pending)
        atomic_fetch_sub(t->pending, 1);
}

/* ── Worker wakeup ────────────────────────────────────────────── */

static void notify_workers(threadpool_t *pool)
{
    atomic_fetch_add_explicit(&pool->wake_gen, 1, memory_order_release);
    if (atomic_load_explicit(&pool->n_sleeping, memory_order_relaxed) > 0) {
        PTHREAD_CHECK(pthread_mutex_lock(&pool->wake_lock));
        PTHREAD_CHECK(pthread_cond_signal(&pool->wake_cond));
        PTHREAD_CHECK(pthread_mutex_unlock(&pool->wake_lock));
    }
}

/* ── Work-finding protocols ───────────────────────────────────── */

/* Called by worker threads: own deque (LIFO) → ext queue → steal. */
static int try_find_work(worker_t *w, tp_task_t *out)
{
    /* 1. Own deque — LIFO preserves DFS order, bounds stack depth. */
    if (deque_pop(&w->deque, out) == 0) return 0;

    /* 2. External queue — picks up tasks from non-worker submitters. */
    if (ext_try_pop(&w->pool->ext, out) == 0) return 0;

    /* 3. Steal from a random peer — FIFO takes shallowest work. */
    int n = w->pool->nworkers;
    if (n > 1) {
        int start = (int)(xorshift32(&w->rng) % (unsigned)n);
        for (int i = 0; i < n; i++) {
            int v = (start + i) % n;
            if (v == w->id) continue;
            if (deque_steal(&w->pool->workers[v].deque, out) == 0)
                return 0;
        }
    }
    return -1;
}

/* Called by non-worker threads (main thread in tp_group_wait). */
static int try_find_work_ext(threadpool_t *pool, tp_task_t *out,
                             unsigned *rng)
{
    if (ext_try_pop(&pool->ext, out) == 0) return 0;

    int n = pool->nworkers;
    int start = (int)(xorshift32(rng) % (unsigned)n);
    for (int i = 0; i < n; i++) {
        int v = (start + i) % n;
        if (deque_steal(&pool->workers[v].deque, out) == 0)
            return 0;
    }
    return -1;
}

/* ── Worker thread ────────────────────────────────────────────── */

static void *worker_main(void *arg)
{
    worker_t *w = (worker_t *)arg;
    tl_worker = w;
    threadpool_t *pool = w->pool;
    tp_task_t task;
    unsigned spins = 0;

    for (;;) {
        if (try_find_work(w, &task) == 0) {
            execute_task(&task);
            spins = 0;
            continue;
        }

        if (atomic_load(&pool->shutdown))
            break;

        if (++spins < 64)
            continue;           /* brief spin before parking */

        /* Park until new work arrives (event-count protocol). */
        int gen = atomic_load_explicit(&pool->wake_gen,
                                       memory_order_acquire);
        atomic_fetch_add(&pool->n_sleeping, 1);
        PTHREAD_CHECK(pthread_mutex_lock(&pool->wake_lock));

        while (atomic_load_explicit(&pool->wake_gen,
                                    memory_order_acquire) == gen
               && !atomic_load(&pool->shutdown))
            PTHREAD_CHECK(pthread_cond_wait(&pool->wake_cond,
                                            &pool->wake_lock));

        PTHREAD_CHECK(pthread_mutex_unlock(&pool->wake_lock));
        atomic_fetch_sub(&pool->n_sleeping, 1);
        spins = 0;
    }
    return NULL;
}

/* ══════════════════════════════════════════════════════════════════
 * Public API: pool lifecycle
 * ══════════════════════════════════════════════════════════════════ */

threadpool_t *tp_create(int nthreads)
{
    if (nthreads <= 0) {
        nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
        if (nthreads <= 0) nthreads = 4;
    }

    threadpool_t *pool = (threadpool_t *)calloc(1, sizeof(*pool));
    if (!pool) return NULL;

    pool->nworkers = nthreads;
    atomic_store_explicit(&pool->shutdown,   0, memory_order_relaxed);
    atomic_store_explicit(&pool->wake_gen,   0, memory_order_relaxed);
    atomic_store_explicit(&pool->n_sleeping, 0, memory_order_relaxed);

    if (pthread_mutex_init(&pool->wake_lock, NULL) != 0) {
        free(pool);
        return NULL;
    }
    if (pthread_cond_init(&pool->wake_cond, NULL) != 0) {
        PTHREAD_CHECK(pthread_mutex_destroy(&pool->wake_lock));
        free(pool);
        return NULL;
    }

    if (ext_init(&pool->ext, 1024) != 0) {
        PTHREAD_CHECK(pthread_cond_destroy(&pool->wake_cond));
        PTHREAD_CHECK(pthread_mutex_destroy(&pool->wake_lock));
        free(pool);
        return NULL;
    }

    pool->workers = (worker_t *)calloc((size_t)nthreads, sizeof(worker_t));
    if (!pool->workers) {
        ext_destroy(&pool->ext);
        PTHREAD_CHECK(pthread_cond_destroy(&pool->wake_cond));
        PTHREAD_CHECK(pthread_mutex_destroy(&pool->wake_lock));
        free(pool);
        return NULL;
    }

    for (int i = 0; i < nthreads; i++) {
        worker_t *w = &pool->workers[i];
        if (deque_init(&w->deque, DEQUE_DEFAULT_CAP) != 0) {
            for (int j = 0; j < i; j++)
                deque_destroy(&pool->workers[j].deque);
            free(pool->workers);
            ext_destroy(&pool->ext);
            PTHREAD_CHECK(pthread_cond_destroy(&pool->wake_cond));
            PTHREAD_CHECK(pthread_mutex_destroy(&pool->wake_lock));
            free(pool);
            return NULL;
        }
        w->pool = pool;
        w->id   = i;
        w->rng  = (unsigned)(i + 1) * 2654435761u;
    }

    /* Default stacks are fine: LIFO pop bounds stack depth to
     * O(tree_depth), which is a few KB at most. */
    for (int i = 0; i < nthreads; i++) {
        if (pthread_create(&pool->workers[i].thread, NULL,
                           worker_main, &pool->workers[i]) != 0) {
            atomic_store(&pool->shutdown, 1);
            PTHREAD_CHECK(pthread_mutex_lock(&pool->wake_lock));
            PTHREAD_CHECK(pthread_cond_broadcast(&pool->wake_cond));
            PTHREAD_CHECK(pthread_mutex_unlock(&pool->wake_lock));
            for (int j = 0; j < i; j++)
                PTHREAD_CHECK(pthread_join(pool->workers[j].thread, NULL));
            for (int j = 0; j < nthreads; j++)
                deque_destroy(&pool->workers[j].deque);
            free(pool->workers);
            ext_destroy(&pool->ext);
            PTHREAD_CHECK(pthread_mutex_destroy(&pool->wake_lock));
            PTHREAD_CHECK(pthread_cond_destroy(&pool->wake_cond));
            free(pool);
            return NULL;
        }
    }

    return pool;
}

void tp_destroy(threadpool_t *pool)
{
    if (!pool) return;

    atomic_store(&pool->shutdown, 1);
    PTHREAD_CHECK(pthread_mutex_lock(&pool->wake_lock));
    PTHREAD_CHECK(pthread_cond_broadcast(&pool->wake_cond));
    PTHREAD_CHECK(pthread_mutex_unlock(&pool->wake_lock));

    for (int i = 0; i < pool->nworkers; i++)
        PTHREAD_CHECK(pthread_join(pool->workers[i].thread, NULL));

    for (int i = 0; i < pool->nworkers; i++) {
        /* Drain per-worker group free lists. */
        tp_group_t *g = pool->workers[i].free_groups;
        while (g) {
            tp_group_t *next = g->next;
            free(g);
            g = next;
        }
        deque_destroy(&pool->workers[i].deque);
    }
    free(pool->workers);
    ext_destroy(&pool->ext);
    PTHREAD_CHECK(pthread_mutex_destroy(&pool->wake_lock));
    PTHREAD_CHECK(pthread_cond_destroy(&pool->wake_cond));
    free(pool);
}

int tp_get_nthreads(threadpool_t *pool)
{
    return pool ? pool->nworkers : 0;
}

void tp_request_shutdown(threadpool_t *pool)
{
    if (!pool) return;
    /* Only atomic store — safe to call from a signal handler. */
    atomic_store(&pool->shutdown, 1);
}

/* ══════════════════════════════════════════════════════════════════
 * Public API: task groups
 * ══════════════════════════════════════════════════════════════════ */

tp_group_t *tp_group_create(threadpool_t *pool)
{
    tp_group_t *g = NULL;
    worker_t *w = tl_worker;

    /* Try per-worker free list — zero contention. */
    if (w && w->pool == pool && w->free_groups) {
        g = w->free_groups;
        w->free_groups = g->next;
        w->n_free--;
    }

    if (!g) {
        g = (tp_group_t *)calloc(1, sizeof(*g));
        if (!g) return NULL;
    }

    g->pool = pool;
    atomic_store_explicit(&g->pending, 0, memory_order_relaxed);
    return g;
}

int tp_group_submit(tp_group_t *group, void (*fn)(void *), void *arg)
{
    atomic_fetch_add(&group->pending, 1);
    tp_task_t task = { .fn = fn, .arg = arg, .pending = &group->pending };

    worker_t *w = tl_worker;
    int pushed = -1;

    if (w && w->pool == group->pool)
        pushed = deque_push(&w->deque, task);

    if (pushed != 0) {
        /* Non-worker or deque full — fall back to global queue. */
        if (ext_push(&group->pool->ext, task) != 0) {
            atomic_fetch_sub(&group->pending, 1);
            return -1;
        }
    }

    notify_workers(group->pool);
    return 0;
}

void tp_group_wait(tp_group_t *group)
{
    threadpool_t *pool = group->pool;
    worker_t *w = tl_worker;
    unsigned ext_rng = (unsigned)(uintptr_t)group * 2654435761u;
    if (ext_rng == 0) ext_rng = 1;
    unsigned spins = 0;

    while (atomic_load(&group->pending) > 0) {
        if (atomic_load_explicit(&pool->shutdown, memory_order_relaxed))
            break;

        tp_task_t task;
        int found;

        if (w && w->pool == pool)
            found = (try_find_work(w, &task) == 0);
        else
            found = (try_find_work_ext(pool, &task, &ext_rng) == 0);

        if (found) {
            execute_task(&task);
            spins = 0;
        } else if (++spins < 128) {
            /* spin */
        } else {
            sched_yield();
            spins = 0;
        }
    }
}

void tp_group_destroy(tp_group_t *group)
{
    worker_t *w = tl_worker;

    if (w && w->pool == group->pool && w->n_free < MAX_FREE_GROUPS) {
        group->next = w->free_groups;
        w->free_groups = group;
        w->n_free++;
    } else {
        free(group);
    }
}

/* ══════════════════════════════════════════════════════════════════
 * Public API: parallel for
 * ══════════════════════════════════════════════════════════════════ */

struct pfor_chunk {
    void (*fn)(int, int, void *);
    void *arg;
    int   start;
    int   end;
};

static void pfor_worker(void *arg)
{
    struct pfor_chunk *c = (struct pfor_chunk *)arg;
    c->fn(c->start, c->end, c->arg);
}

#define MAX_STACK_CHUNKS 64

void tp_parallel_for(threadpool_t *pool, int start, int end,
                     void (*fn)(int, int, void *), void *arg)
{
    if (start >= end) return;

    int n       = end - start;
    int nchunks = pool->nworkers + 1;
    if (nchunks > n) nchunks = n;

    struct pfor_chunk stack_chunks[MAX_STACK_CHUNKS];
    struct pfor_chunk *chunks = stack_chunks;
    int heap = 0;

    if (nchunks > MAX_STACK_CHUNKS) {
        chunks = (struct pfor_chunk *)calloc((size_t)nchunks,
                                             sizeof(struct pfor_chunk));
        if (!chunks) { fn(start, end, arg); return; }
        heap = 1;
    }

    tp_group_t *group = tp_group_create(pool);
    if (!group) {
        if (heap) free(chunks);
        fn(start, end, arg);
        return;
    }

    int chunk_size = (n + nchunks - 1) / nchunks;
    int serial_from = end;  /* if submit fails, run remaining chunks serially */
    for (int i = 0; i < nchunks; i++) {
        int cs = start + i * chunk_size;
        int ce = cs + chunk_size;
        if (ce > end) ce = end;
        if (cs >= end) break;
        chunks[i].fn    = fn;
        chunks[i].arg   = arg;
        chunks[i].start = cs;
        chunks[i].end   = ce;
        if (tp_group_submit(group, pfor_worker, &chunks[i]) != 0) {
            serial_from = cs;
            break;
        }
    }

    tp_group_wait(group);
    tp_group_destroy(group);

    /* Execute any chunks that failed to submit. */
    if (serial_from < end)
        fn(serial_from, end, arg);

    if (heap) free(chunks);
}
