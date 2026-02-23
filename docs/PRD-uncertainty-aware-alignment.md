# PRD: Uncertainty-Aware Multiple Sequence Alignment

## Overview

Two novel features for kalign that exploit its unique architecture to improve
alignment accuracy while maintaining speed:

1. **Confidence-guided targeted refinement** -- Extract alignment confidence
   from Hirschberg forward-backward margins; refine only uncertain edges.
2. **Stochastic anchor ensembles** -- Use multiple random anchor sets to build
   more robust guide trees with built-in distance confidence estimates.

Together these form a coherent contribution: "uncertainty-aware MSA" where
both the guide tree and the alignment itself are informed by confidence
estimates, and compute is focused where uncertainty is highest.

---

## Feature 1: Confidence-Guided Targeted Refinement

### Motivation

Kalign currently performs a single progressive alignment pass with no
refinement. Competitors (MUSCLE, MAFFT) do blind iterative refinement --
re-aligning at every tree edge or random partition. This is wasteful because
most edges produce confident alignments that refinement does not improve.

Kalign's Hirschberg algorithm already computes forward and backward DP scores
meeting at a midpoint. The gap between the best and second-best meetup scores
is a free confidence signal that nobody currently exploits.

### Design

#### 1.1 Confidence Extraction

**Where:** Inside the three `*_meetup()` functions in `aln_seqseq.c`,
`aln_seqprofile.c`, and `aln_profileprofile.c`.

**What:** Track both the best score (`max`) and the second-best score
(`max2`) across all candidate meetup columns and transitions. The confidence
margin is `max - max2`.

**How:** Each meetup function already iterates over all columns `i` in
`startb..endb` testing 7 transition types. Currently only `max` is tracked.
Add a `max2` variable that tracks the second-highest score. At the end:

```
margin = max - max2;
```

The margin should be normalised by the alignment region size to produce a
per-column confidence. This normalisation accounts for the fact that longer
alignments have larger absolute scores.

**Accumulation:** The Hirschberg algorithm recurses, calling meetup at every
split level. We accumulate a per-edge confidence score by summing the
normalised margins across all meetup calls for one progressive alignment step.
This is stored in the task's `score` field (currently unused after alignment).

**Data flow:**
- `aln_mem` gets two new fields: `float margin_sum` and `int margin_count`
- Each meetup call adds its margin to `margin_sum` and increments
  `margin_count`
- After alignment, `task->score = margin_sum / margin_count` (average
  per-split confidence for this progressive alignment step)

#### 1.2 Refinement Engine

**New file:** `lib/src/aln_refine.c` (with `aln_refine.h`)

**Core function:**

```c
int refine_alignment(struct msa* msa, struct aln_param* ap,
                     struct aln_tasks* tasks);
```

**Algorithm (tree-dependent targeted refinement):**

1. **Rank edges by confidence.** Sort `tasks->list[]` by ascending `score`
   (lowest confidence = most uncertain). The first N edges (configurable,
   default: `min(numseq/2, 50)`) are candidates for refinement.

2. **For each candidate edge** (task `t` with children `a`, `b`):

   a. **Save current state.** Deep-copy `gaps[]` arrays for all sequences in
      `sip[a]` and `sip[b]`. Record current alignment length `plen[c]`.

   b. **Strip gaps between groups.** For each sequence in `sip[a]` union
      `sip[b]`, remove the gap columns that were introduced at this edge and
      all ancestor edges. This is done by calling a new function
      `strip_group_gaps()` that zeros only the inter-group gaps while
      preserving intra-group gaps (gaps introduced by merges below nodes
      `a` and `b` in the tree).

   c. **Rebuild profiles.** Build a fresh profile for group `a` by calling
      `build_group_profile()` which walks the aligned sequences in `sip[a]`
      column-by-column, summing residue counts and substitution scores. Same
      for group `b`.

   d. **Re-align.** Call the existing alignment machinery
      (`aln_runner` / `aln_runner_serial`) with `mode = ALN_MODE_FULL` on the
      two group profiles. This produces a new path.

   e. **Score comparison.** Compare old vs new alignment quality. Two options:
      - Use the alignment score from `aln_runner` (DP score)
      - Use a fast SP score sample (compare a random sample of sequence pairs)

      If the new alignment is not better, restore the saved `gaps[]` and
      continue to the next edge. If better, commit the new gaps.

   f. **Propagate gaps.** Call `make_seq()` with the new path to update
      `gaps[]` for all member sequences of `a` and `b`. Update `plen[c]`.

3. **Convergence.** Repeat steps 1-2 for up to `max_iterations` rounds
   (default: 2). Stop early if no edge improved in the last round.

#### 1.3 Helper Functions

**`strip_group_gaps(struct msa* msa, int node_a, int node_b, int** saved_gaps)`**

Removes inter-group gap columns from the gap arrays of sequences under
`node_a` and `node_b`. Works by:
1. Saving current `gaps[]` arrays for all sequences in both groups.
2. Reconstructing the alignment of group A alone (from the subtree below A)
   and group B alone.
3. The intra-group gaps (from merges within A's subtree and B's subtree) are
   preserved; only the gap columns added when A and B were merged (and at
   all ancestor merges) are stripped.

Implementation approach: Walk the gap arrays and the path stored during the
original merge. Identify which gaps were added at this edge vs inherited from
below.

**Alternative simpler approach:** Instead of surgical gap stripping, fully
strip all gaps from both groups (call `clean_aln()` logic on the subset) and
rebuild from scratch using the subtree below each node. This is simpler and
more robust, though slightly slower. For the first implementation, prefer this
approach.

**`build_group_profile(struct msa* msa, struct aln_param* ap, int node, struct aln_tasks* tasks, float** profile_out, int* len_out)`**

Builds a merged profile for an already-aligned group of sequences. Two
strategies:

*Strategy A (simple, recommended for v1):* Replay the progressive alignment
for the subtree below `node`. Walk the tree bottom-up, calling
`make_profile_n()` for leaves and `update_n()` at each internal merge, using
the existing gap arrays (which encode the intra-group alignment). This reuses
all existing code and guarantees consistency.

*Strategy B (faster, for v2):* Directly construct a profile from the aligned
columns. For each alignment column, sum residue counts and substitution
scores across all sequences in the group. This avoids the tree replay but
requires a new column-iteration function.

#### 1.4 Integration with `aln_wrap.c`

The refinement step is called after the initial progressive alignment and
before `finalise_alignment()`:

```c
int kalign_run(struct msa* msa, int n_threads, int type,
               float gpo, float gpe, float tgpe)
{
        /* ... existing code up to create_msa_tree() ... */
        RUN(create_msa_tree(msa, ap, tasks));
        msa->aligned = ALN_STATUS_ALIGNED;

        /* NEW: Targeted refinement */
        RUN(refine_alignment(msa, ap, tasks));

        RUN(finalise_alignment(msa));
        /* ... rest unchanged ... */
}
```

A command-line flag `--norefine` (or `--refine N` to set max iterations)
should be added to `run_kalign.c` to control this.

---

## Feature 2: Stochastic Anchor Ensembles

### Motivation

Kalign represents each sequence as a vector of BPM edit distances to 32
anchor sequences, then clusters in this 32-dimensional space. The anchor
choice is deterministic and length-stratified. A bad anchor set produces a
bad embedding, a bad tree, and a bad alignment. There is no way to detect
this.

### Design

#### 2.1 Multiple Anchor Sets

**Modified file:** `lib/src/pick_anchor.c`

**New function:**

```c
int* pick_anchor_diverse(struct msa* msa, int* n, struct rng_state* rng);
```

In addition to the existing length-stratified `pick_anchor()`, add a
diversity-aware variant that:

1. Picks the first anchor randomly (using `rng`).
2. Picks subsequent anchors greedily by maximum minimum BPM distance to all
   already-chosen anchors (farthest-first traversal).

This requires computing O(K * N) BPM distances during anchor selection
(where K=32, N=numseq), which is cheap relative to the full distance matrix.

**New function:**

```c
int pick_anchor_random(struct msa* msa, int* n, struct rng_state* rng);
```

Pure random selection of K anchors. Fast, provides maximum variance between
ensemble members.

#### 2.2 Ensemble Distance Matrix

**New file:** `lib/src/anchor_ensemble.c` (with `anchor_ensemble.h`)

**Core function:**

```c
int build_ensemble_distance(struct msa* msa, int n_ensemble,
                            struct rng_state* rng,
                            float*** dm_out, int* num_anchors_out,
                            float** variance_out);
```

**Algorithm:**

1. Generate `n_ensemble` (default: 5) different anchor sets. Use a mix of
   strategies:
   - 1x length-stratified (current method, for reproducibility)
   - 2x diversity-aware (farthest-first with different random seeds)
   - 2x pure random

2. For each anchor set `k`, compute the N x K_k distance matrix `D_k` using
   the existing `d_estimation()` with `pair=0`.

3. For each pair of sequences `(i, j)`, compute the Euclidean distance in
   each embedding:
   ```
   d_k(i,j) = euclidean_dist(D_k[i], D_k[j])
   ```

4. The consensus distance is the **median** across ensemble members:
   ```
   d(i,j) = median(d_1(i,j), ..., d_E(i,j))
   ```
   Median is preferred over mean because it is robust to outlier anchor sets.

5. The distance **variance** across ensemble members:
   ```
   var(i,j) = variance(d_1(i,j), ..., d_E(i,j))
   ```
   This is a free confidence signal for each pairwise distance.

**Note:** We do NOT compute full N x N distance matrices. The Euclidean
distances are only computed during bisecting k-means (inside `split2()`) and
UPGMA (for small subclusters). The ensemble produces `n_ensemble` separate
N x K distance matrices. The k-means splitting then uses the median of the
per-ensemble Euclidean distances.

**Efficient implementation:** Instead of computing full pairwise Euclidean
distances, modify `split2()` in `bisectingKmeans.c` to compute distances to
centroids using all ensemble embeddings and take the median. This keeps the
O(N * K * E) complexity per split (where E is the ensemble size).

#### 2.3 Variance-Weighted Tree Construction

The distance variance from the ensemble can be used in two ways:

**Option A (simple):** Feed the median distances into the existing bisecting
k-means unchanged. The median is more robust than a single anchor set, so
tree quality improves without algorithmic changes to the clustering.

**Option B (advanced):** Weight sequences during k-means assignment by the
inverse of their distance variance. Sequences with high variance (unreliable
distances) get less influence on centroid computation, preventing them from
distorting cluster assignments.

Recommend Option A for v1, Option B as a follow-up experiment.

#### 2.4 Integration into Guide Tree Construction

**Modified file:** `lib/src/bisectingKmeans.c`

The `build_tree_kmeans()` function currently:
1. Calls `pick_anchor()` once
2. Calls `d_estimation()` once
3. Calls `bisecting_kmeans()` with the single distance matrix

With the ensemble, this becomes:
1. For `k = 0..n_ensemble-1`:
   a. Call `pick_anchor_variant(msa, &n, rng, strategy[k])` to get anchor set
   b. Call `d_estimation(msa, dm[k], anchors[k], n, 0)` to compute distances
2. Store all `n_ensemble` distance matrices in the task structure
3. Modify `split2()` to compute centroid distances using median across
   ensembles

The UPGMA fallback for small clusters uses full pairwise BPM distances (via
`d_estimation(..., pair=1)`), which is independent of anchors. No change
needed there.

#### 2.5 Variance as Confidence Signal for Refinement

The per-edge distance variance from the ensemble provides a second
confidence signal (complementary to the Hirschberg margin). Edges in the
guide tree that were constructed from high-variance distances are more likely
to be suboptimal and therefore benefit more from refinement.

In the refinement step (Feature 1), use a combined confidence score:

```
edge_confidence = w1 * hirschberg_margin + w2 * (1 / tree_edge_variance)
```

Where `tree_edge_variance` is the average distance variance between the two
groups of sequences at that tree edge. Low confidence edges are refined
first. The weights `w1`, `w2` are tunable (default: equal weight after
normalisation to [0,1] range).

---

## New Files

| File | Purpose |
|------|---------|
| `lib/src/aln_refine.c` | Targeted refinement engine |
| `lib/src/aln_refine.h` | Header for refinement |
| `lib/src/anchor_ensemble.c` | Ensemble anchor distance computation |
| `lib/src/anchor_ensemble.h` | Header for anchor ensemble |

## Modified Files

| File | Change |
|------|--------|
| `lib/src/aln_struct.h` | Add `margin_sum`, `margin_count` to `aln_mem` |
| `lib/src/aln_seqseq.c` | Track `max2` in meetup, write margin to `aln_mem` |
| `lib/src/aln_seqprofile.c` | Same |
| `lib/src/aln_profileprofile.c` | Same |
| `lib/src/aln_wrap.c` | Call `refine_alignment()` after `create_msa_tree()` |
| `lib/src/bisectingKmeans.c` | Accept ensemble distance matrices, median in `split2()` |
| `lib/src/pick_anchor.c` | Add `pick_anchor_diverse()`, `pick_anchor_random()` |
| `lib/src/pick_anchor.h` | Declare new functions |
| `lib/src/task.h` | Add `float confidence` to `struct task` |
| `lib/CMakeLists.txt` | Add `aln_refine.c`, `anchor_ensemble.c` to `source_files` |
| `src/run_kalign.c` | Add `--norefine` / `--refine N` CLI flags |
| `CMakeLists.txt` | Add `KALIGN_REFINE_MAX_ITER` tunable (default: 2) |

---

## Code Style Requirements

All new code must follow existing kalign conventions:

### Header Pattern

```c
#ifndef HEADER_NAME_H
#define HEADER_NAME_H

#ifdef HEADER_NAME_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct forward_decl;

EXTERN int function_name(struct type* arg);

#undef HEADER_NAME_IMPORT
#undef EXTERN

#endif
```

### Source Pattern

```c
#include "tldevel.h"
/* other includes */

#define HEADER_NAME_IMPORT
#include "header_name.h"

/* static function declarations */
static int helper_func(int arg);

/* public functions first */
int function_name(struct type* arg)
{
        int local_var = 0;
        /* ... */
        MMALLOC(ptr, sizeof(type) * count);
        RUN(other_func(arg));
        RUNP(ptr = alloc_func());
        /* ... */
        return OK;
ERROR:
        return FAIL;
}

/* static functions after */
static int helper_func(int arg)
{
        /* ... */
        return OK;
ERROR:
        return FAIL;
}
```

### Conventions

- **Indentation:** 8-space tabs (actual tab characters)
- **Braces:** Opening brace on same line as `if`/`for`/`while`; function
  opening brace on next line (K&R variant)
- **Memory:** Use `MMALLOC()`, `MREALLOC()`, `MFREE()` macros (never raw
  `malloc`/`free`)
- **Errors:** Every function that can fail returns `int` (`OK`/`FAIL`).
  Use `RUN()` to call and check. Use `RUNP()` for functions returning
  pointers. Use `ASSERT()` for preconditions. Every function with `MMALLOC`
  or `RUN` must have an `ERROR:` label at the end that returns `FAIL`.
- **Logging:** `LOG_MSG()`, `WARNING_MSG()`, `ERROR_MSG()` macros
- **Naming:** `snake_case` for functions and variables. Uppercase for macros
  and constants. Struct names are lowercase.
- **OpenMP:** Guard all OpenMP pragmas with `#ifdef HAVE_OPENMP`. Use
  `task` parallelism (not `parallel for`) for tree operations.
- **Types:** Use `int` for most integers, `float` for DP scores and profile
  data, `uint8_t` for encoded sequences, `register` keyword for DP loop
  variables (matches existing style even though modern compilers ignore it).
- **Profile stride:** 64 floats per position (`<< 6` for pointer
  arithmetic).
- **Float sentinel:** `-FLT_MAX` for impossible DP states.
- **GPL header:** Include the GPL v3 license header in every new file.
  Copyright: `Copyright 2006, 2019, 2024, 2025 Timo Lassmann`.

---

## Implementation Plan

### Phase 1: Confidence Extraction (small, testable)

1. Add `margin_sum` and `margin_count` to `struct aln_mem`
2. Modify the three `*_meetup()` functions to track `max2` and accumulate
   margins
3. Store per-task confidence in `task->score` after alignment
4. Add a test that verifies confidence values are sensible (high for easy
   alignments, low for difficult ones)

### Phase 2: Refinement Engine (core feature)

1. Implement `strip_group_gaps()` -- zero gaps for a subset of sequences
2. Implement `build_group_profile()` -- replay subtree to build profile
3. Implement `refine_alignment()` -- the main refinement loop
4. Integrate into `aln_wrap.c`
5. Add `--norefine` CLI flag
6. Benchmark on BAliBASE: compare baseline vs 1-round vs 2-round refinement

### Phase 3: Anchor Ensemble (independent of Phase 2)

1. Implement `pick_anchor_diverse()` using farthest-first traversal
2. Implement `pick_anchor_random()`
3. Implement `build_ensemble_distance()` -- multi-anchor distance matrices
4. Modify `split2()` in `bisectingKmeans.c` to use median distances
5. Benchmark: single anchor set vs ensemble of 5

### Phase 4: Combined System

1. Connect ensemble variance to refinement edge selection
2. Tune weights for combined confidence score
3. Full benchmark: baseline vs ensemble-only vs refine-only vs both
4. Timing comparison: measure overhead of each feature
5. Test on BAliBASE, BaliFam100, BRAliBASE

### Phase 5: Polish

1. Ensure deterministic output when `n_ensemble=1` and `--norefine` (exact
   backward compatibility)
2. Thread safety for all new code
3. Memory leak testing (valgrind/ASAN)
4. Update Python bindings to expose refinement control
5. Documentation

---

## Testing Strategy

### Unit Tests

- **Confidence extraction:** Align two identical sequences; margin should be
  very high. Align two random sequences; margin should be low.
- **Refinement:** Align a BAliBASE test case. Verify that refinement does
  not decrease the SP score. Verify that `--norefine` produces identical
  output to current kalign.
- **Anchor ensemble:** Verify that median of 5 distance matrices is
  element-wise correct. Verify that variance is non-negative.
- **Determinism:** With fixed RNG seed, ensemble produces identical results
  across runs. Without ensemble (`n_ensemble=1`), output matches current
  kalign exactly.

### Integration Tests

- **BAliBASE benchmark:** Run full suite. Report SP scores for:
  - Baseline (current kalign)
  - +Refinement only (1 round, 2 rounds)
  - +Ensemble only (3, 5, 7 ensemble members)
  - +Both
- **Timing:** Wall-clock times for each configuration
- **Memory:** Peak RSS for each configuration

### Regression Tests

- All existing tests must pass unchanged when refinement is disabled
- The DSSIM determinism test (different thread counts, shuffled input)
  must still pass

---

## Tunable Parameters

| Parameter | Default | CMake Variable | Description |
|-----------|---------|----------------|-------------|
| `max_refine_iter` | 2 | `KALIGN_REFINE_MAX_ITER` | Max refinement rounds |
| `max_refine_edges` | numseq/2 | -- | Max edges to refine per round |
| `n_ensemble` | 5 | `KALIGN_ENSEMBLE_SIZE` | Number of anchor sets |
| `num_anchors` | 32 | -- | Anchors per set (unchanged) |

These should also be controllable via the C API (`kalign_run` parameters)
and exposed through the Python bindings.

---

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| Refinement increases runtime significantly | Targeted refinement limits work to uncertain edges; `--norefine` for speed |
| Ensemble overhead for large datasets | O(N * K * E) is still subquadratic; E=5 is small |
| Profile reconstruction from subtree is complex | Phase 2 uses simple "replay subtree" approach, not novel code |
| Numerical issues with float accumulation | Monitor score magnitudes; consider double for refinement scoring |
| Breaking determinism | Fixed RNG seed for ensemble; refinement is deterministic given input |
