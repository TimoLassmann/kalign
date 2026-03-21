/*
 * threadpool.h — Lightweight work-stealing thread pool.
 *
 * A POSIX-threads-based thread pool using Chase-Lev work-stealing deques.
 * Designed for recursive divide-and-conquer workloads (guide-tree traversal,
 * Hirschberg alignment) and data-parallel loops (distance matrices).
 *
 * Supports three parallelism patterns:
 *
 *   1. PARALLEL FOR — tp_parallel_for()
 *      Splits [start, end) into chunks, one per worker.  The calling thread
 *      participates.  Best for embarrassingly parallel loops.
 *
 *   2. FORK-JOIN — tp_group_submit() + tp_group_wait()
 *      Submit N independent tasks, then block until all complete.  The
 *      waiting thread executes queued work while it waits.
 *
 *   3. RECURSIVE TASKS — groups created inside tasks (nested wait)
 *      A task can create a new group, submit children, and wait.  LIFO
 *      deque ordering ensures the C stack depth is O(tree_depth), not
 *      O(tree_size).  Safe for trees with millions of nodes.
 *
 * LIFECYCLE:
 *   threadpool_t *pool = tp_create(0);          // 0 = auto-detect CPUs
 *   // ... use pool for parallel work ...
 *   tp_destroy(pool);                           // joins all workers, frees all memory
 *
 * THREAD SAFETY:
 *   - tp_create / tp_destroy: call from one thread only.
 *   - tp_group_*: a group must be used by one owner thread at a time.
 *     The owner creates it, submits tasks, waits, and destroys it.
 *     Tasks submitted to the group may run on any worker thread.
 *   - tp_parallel_for: safe to call from any thread, including from
 *     inside a task (nested parallelism).
 *   - tp_request_shutdown: safe to call from any thread or signal handler.
 *
 * ERROR HANDLING:
 *   - tp_create returns NULL on failure (out of memory, pthread errors).
 *   - tp_group_create returns NULL on out of memory.
 *   - tp_group_submit returns -1 on failure (out of memory); the task
 *     is NOT executed.  Check the return value.
 *   - tp_parallel_for never fails: it falls back to serial execution
 *     if allocation or submission fails.
 *
 * GRACEFUL SHUTDOWN (signal handling):
 *   Call tp_request_shutdown(pool) to make tp_group_wait() return early.
 *   This function is async-signal-safe (it only writes an atomic flag).
 *   After shutdown, call tp_destroy(pool) from the main thread.
 *
 * USAGE EXAMPLES:
 *
 *   // Pattern 1: Parallel for (distance matrix)
 *   void compute_row(int start, int end, void *arg) {
 *       float **dm = (float **)arg;
 *       for (int i = start; i < end; i++)
 *           dm[i] = compute_distances(i);
 *   }
 *   tp_parallel_for(pool, 0, num_seqs, compute_row, dm);
 *
 *   // Pattern 2: Fork-join (two independent tasks)
 *   tp_group_t *g = tp_group_create(pool);
 *   tp_group_submit(g, forward_pass,  &fwd_arg);
 *   tp_group_submit(g, backward_pass, &bwd_arg);
 *   tp_group_wait(g);
 *   tp_group_destroy(g);
 *
 *   // Pattern 3: Recursive tasks (guide-tree traversal)
 *   void align_node(void *arg) {
 *       struct node *n = (struct node *)arg;
 *       if (is_leaf(n)) { align_leaf(n); return; }
 *       tp_group_t *g = tp_group_create(pool);
 *       tp_group_submit(g, align_node, n->left);
 *       tp_group_submit(g, align_node, n->right);
 *       tp_group_wait(g);
 *       tp_group_destroy(g);
 *       merge_alignments(n);
 *   }
 */

#ifndef THREADPOOL_H
#define THREADPOOL_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct threadpool threadpool_t;
typedef struct tp_group   tp_group_t;

/* ── Pool lifecycle ─────────────────────────────────────────────── */

/* Create a pool with nthreads workers (0 = auto-detect CPU count).
 * Returns NULL on failure. */
threadpool_t *tp_create(int nthreads);

/* Destroy the pool.  Joins all worker threads and frees all memory.
 * All task groups MUST have been waited on and destroyed first. */
void tp_destroy(threadpool_t *pool);

/* Request shutdown without blocking.  Causes tp_group_wait() to
 * return early and workers to exit.  Follow with tp_destroy().
 * Async-signal-safe (writes an atomic flag only). */
void tp_request_shutdown(threadpool_t *pool);

/* Number of worker threads in the pool. */
int tp_get_nthreads(threadpool_t *pool);

/* ── Task groups ────────────────────────────────────────────────── */

/* Create a group for coordinating related tasks.
 * Returns NULL on out of memory. */
tp_group_t *tp_group_create(threadpool_t *pool);

/* Submit fn(arg) to the pool under this group.
 * Returns 0 on success, -1 on failure (task NOT executed).
 * May be called from any thread, including from inside a task. */
int tp_group_submit(tp_group_t *group, void (*fn)(void *), void *arg);

/* Block until every task in the group has completed.
 * The calling thread participates in executing queued work while waiting.
 * Returns early if tp_request_shutdown() has been called. */
void tp_group_wait(tp_group_t *group);

/* Free a group.  Must be called after tp_group_wait(). */
void tp_group_destroy(tp_group_t *group);

/* ── Parallel for ───────────────────────────────────────────────── */

/* Partition [start, end) into chunks and call fn(chunk_start, chunk_end, arg)
 * for each chunk.  The calling thread participates in the work.
 * Falls back to serial execution on allocation failure (never fails). */
void tp_parallel_for(threadpool_t *pool, int start, int end,
                     void (*fn)(int start, int end, void *arg), void *arg);

#ifdef __cplusplus
}
#endif

#endif /* THREADPOOL_H */
