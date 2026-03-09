# PRD: Unified C Entry Point — `kalign_align_full`

## Goal

Replace the 7+ separate alignment entry points in kalign's C library with a single comprehensive function `kalign_align_full()` that accepts an array of per-run configs. All callers (CLI, Python bindings, benchmark optimizer) route through this one function, eliminating duplicated routing logic, silent parameter dropping, and inconsistent sentinel handling.

## Current State: Detailed Audit

### C Library Entry Points (`lib/src/aln_wrap.c` + `lib/src/ensemble.c`)

#### 1. `kalign()` — aln_wrap.c line 110

```c
int kalign(char **seq, int *len, int numseq, int n_threads, int type,
           float gpo, float gpe, float tgpe, char ***aligned, int *out_aln_len)
```

- **What it does**: Wraps `kalign_run()` with arr→msa→arr conversion.
- **Hardcoded**: `refine=NONE`, `adaptive_budget=0`, `quiet=1`
- **Missing**: Everything else (vsm, seq_weights, consistency, realign, ensemble, dist_scale)
- **Used by**: Legacy C API consumers only

#### 2. `kalign_run()` — aln_wrap.c line 263

```c
int kalign_run(struct msa *msa, int n_threads, int type,
               float gpo, float gpe, float tgpe, int refine, int adaptive_budget)
```

- **What it does**: Thin wrapper around `kalign_run_seeded()`.
- **Hardcoded**: `tree_seed=0, tree_noise=0.0, dist_scale=0.0, vsm_amax=-1.0, use_seq_weights=-1.0, consistency_anchors=0, consistency_weight=2.0`
- **Missing**: All advanced features
- **Used by**: `kalign()` wrapper above, and the pybind11 router as the "most basic" fallback

#### 3. `kalign_run_seeded()` — aln_wrap.c line 133

```c
int kalign_run_seeded(struct msa *msa, int n_threads, int type,
                      float gpo, float gpe, float tgpe,
                      int refine, int adaptive_budget,
                      uint64_t tree_seed, float tree_noise,
                      float dist_scale, float vsm_amax,
                      float use_seq_weights,
                      int consistency_anchors, float consistency_weight)
```

- **What it does**: Full single-run alignment with all knobs.
- **Sentinel handling**:
  - `use_seq_weights >= 0.0` → override (default is 0.0 from `aln_param_init`)
  - `dist_scale > 0.0` → override (GUARDED: 0.0 means "keep default", which is also 0.0)
  - `vsm_amax >= 0.0` → override (default is 2.0 for protein, 0.0 for DNA)
- **Missing**: realign iterations (separate function), ensemble
- **Consistency**: YES — builds anchor table if `consistency_anchors > 0`, frees it after refinement
- **Used by**: pybind11 router (when consistency > 0), CLI fallback, `kalign_run()`, ensemble per-run alignments

#### 4. `kalign_run_dist_scale()` — aln_wrap.c line 268

```c
int kalign_run_dist_scale(struct msa *msa, int n_threads, int type,
                          float gpo, float gpe, float tgpe,
                          int refine, int adaptive_budget,
                          float dist_scale, float vsm_amax,
                          float use_seq_weights)
```

- **What it does**: Like `kalign_run_seeded` but WITHOUT tree noise and WITHOUT consistency.
- **Sentinel handling**:
  - `dist_scale` → UNCONDITIONAL assignment (differs from seeded!)
  - `vsm_amax >= 0.0` → override
  - `use_seq_weights >= 0.0` → override
- **Missing**: `tree_seed`, `tree_noise`, `consistency_anchors`, `consistency_weight`
- **BUG/GOTCHA**: No consistency support. If pybind11 router picks this path, consistency params are silently dropped.
- **Used by**: pybind11 router (when `dist_scale > 0 || vsm_amax >= 0 || seq_weights >= 0` but no consistency and no realign)

#### 5. `kalign_run_realign()` — aln_wrap.c line 361

```c
int kalign_run_realign(struct msa *msa, int n_threads, int type,
                       float gpo, float gpe, float tgpe,
                       int refine, int adaptive_budget,
                       float dist_scale, float vsm_amax,
                       int realign_iterations,
                       float use_seq_weights,
                       int consistency_anchors, float consistency_weight)
```

- **What it does**: Full single-run + iterative tree rebuild from alignment distances.
- **Sentinel handling**:
  - `dist_scale` → UNCONDITIONAL assignment
  - `vsm_amax >= 0.0` → override
  - `use_seq_weights >= 0.0` → override
- **Key behavior**: Builds initial tree from BPM, then after first alignment: finalise → compute pairwise distances → strip gaps → rebuild UPGMA tree → re-align. Consistency table built ONCE before first alignment (not rebuilt per iteration).
- **Missing**: `tree_seed`, `tree_noise` (always deterministic initial tree)
- **Used by**: pybind11 router (when `realign > 0`), CLI (when `realign > 0`), ensemble per-run alignments (when `realign > 0`)

#### 6. `kalign_ensemble()` — ensemble.c line 223

```c
int kalign_ensemble(struct msa* msa, int n_threads, int type,
                    int n_runs, float gpo, float gpe, float tgpe,
                    uint64_t seed, int min_support,
                    const char* save_poar_path,
                    int refine, float dist_scale, float vsm_amax,
                    int realign, float use_seq_weights,
                    int consistency_anchors, float consistency_weight)
```

- **What it does**: N runs with hardcoded scale-factor diversity table, POAR consensus/selection.
- **Hardcoded values**:
  - `use_seq_weights`: forced to 0.0 when `-1.0` (sentinel). Comment says "hurts ensemble performance".
  - Per-run diversity from `run_params[]` table: 12 entries with `gpo_scale`, `gpe_scale`, `tgpe_scale`, `noise` multipliers.
  - Run 0 always uses base params, no noise.
  - Post-selection refinement: HARDCODES `KALIGN_REFINE_CONFIDENT` for the refinement re-run (line 421), regardless of the `refine` param passed for per-run alignments.
  - Post-selection refinement: HARDCODES `dist_scale=0.0f` (line 423).
  - Auto min_support when `min_support=0`: `min_sup = (n_runs + 2) / 3`, clamped to >= 2.
- **Missing**: Per-run matrix types (always uses same `type` for all runs)
- **Used by**: pybind11 router (when `ensemble > 0`), CLI (when `ensemble > 0`)

#### 7. `kalign_ensemble_custom()` — ensemble.c line 515

```c
int kalign_ensemble_custom(struct msa* msa, int n_threads, int type,
                           int n_runs,
                           const float* run_gpo, const float* run_gpe,
                           const float* run_tgpe, const int* run_types,
                           const float* run_noise,
                           uint64_t seed, int min_support,
                           int refine, float vsm_amax,
                           int realign, float use_seq_weights,
                           int consistency_anchors, float consistency_weight)
```

- **What it does**: Like `kalign_ensemble` but with per-run gap penalties, types, and noise.
- **Hardcoded values**:
  - `use_seq_weights`: forced to 0.0 when `-1.0` (sentinel), same as `kalign_ensemble`.
  - Post-selection refinement: HARDCODES `KALIGN_REFINE_CONFIDENT` (line 676).
  - Post-selection refinement: HARDCODES `dist_scale=0.0f` (line 678).
  - Post-selection refinement: HARDCODES `adaptive_budget=0` (line 675).
  - Auto min_support: same `(n_runs + 2) / 3` formula.
- **Missing**: Per-run `dist_scale`, `vsm_amax`, `seq_weights`, `consistency` (all shared across runs). `save_poar` (removed from this function).
- **Used by**: pybind11 `ensemble_custom_file_to_file()` only (not accessible from CLI or Python wrapper)

#### 8. `kalign_post_realign()` — aln_wrap.c line 539

```c
int kalign_post_realign(struct msa *msa, int n_threads, int type,
                        float gpo, float gpe, float tgpe,
                        int refine, int adaptive_budget,
                        float dist_scale, float vsm_amax,
                        int realign_iterations,
                        float use_seq_weights)
```

- **What it does**: Takes an already-finalized MSA (e.g., from ensemble) and does realign iterations.
- **Missing**: `consistency_anchors`, `consistency_weight` — no consistency support in post-realign!
- **Used by**: Was intended for post-ensemble realign, but this was found to HURT performance (0.758→0.708) so it's currently unused.

### Routing Logic Comparison

#### pybind11 `run_alignment()` (`_core.cpp:73-104`)

```
if load_poar         → kalign_consensus_from_poar()
elif ensemble > 0    → kalign_ensemble()              ← old hardcoded diversity, NOT ensemble_custom
elif realign > 0     → kalign_run_realign()            ← OK, has consistency
elif consistency > 0 → kalign_run_seeded()             ← OK
elif dist_scale > 0  → kalign_run_dist_scale()         ← BUG: drops consistency!
  OR vsm_amax >= 0
  OR seq_weights >= 0
else                 → kalign_run()                    ← drops everything advanced
```

**Known bugs in the router**:
1. If you pass `vsm_amax=2.0` and `consistency=5`, the router picks `kalign_run_seeded` (correct, because `consistency > 0` is checked first). But if you pass `vsm_amax=2.0` and `consistency=0`, the router picks `kalign_run_dist_scale` which doesn't support consistency — OK in this case since consistency=0, but fragile logic.
2. `kalign_ensemble()` is the old hardcoded version. `ensemble_custom` is only accessible via a separate function.

#### CLI `run_kalign()` (`src/run_kalign.c:409-464`)

```
if load_poar         → kalign_consensus_from_poar()
elif ensemble > 0    → kalign_ensemble()              ← old hardcoded diversity only
elif realign > 0     → kalign_run_realign()
else                 → kalign_run_seeded()             ← always this for single-run
```

**Key differences from pybind11 router**:
1. CLI always uses `kalign_run_seeded` for non-ensemble/non-realign — never falls through to `kalign_run_dist_scale` or `kalign_run`. This means the CLI consistently handles vsm_amax and seq_weights correctly.
2. CLI has `--fast` (consistency=0) and `--precise` (ensemble=3, realign=1) modes.
3. CLI has no access to `ensemble_custom` at all.

### Sentinel Value Summary

| Parameter | Sentinel | Meaning | Where checked |
|-----------|----------|---------|---------------|
| `gpo` | -1.0 | Use matrix default | `aln_param_init`: `if(gpo >= 0.0)` |
| `gpe` | -1.0 | Use matrix default | `aln_param_init`: `if(gpe >= 0.0)` |
| `tgpe` | -1.0 | Use matrix default | `aln_param_init`: `if(tgpe >= 0.0)` |
| `vsm_amax` | -1.0 | Use biotype default (2.0 protein, 0.0 DNA) | `if(vsm_amax >= 0.0)` |
| `use_seq_weights` | -1.0 | Use biotype default (0.0 for all) | `if(use_seq_weights >= 0.0)` |
| `dist_scale` | 0.0 | Off | Inconsistent: guarded in seeded, unconditional in others |
| `tree_seed` | 0 | Deterministic tree | `if(tree_seed != 0 && tree_noise > 0.0f)` |
| `tree_noise` | 0.0 | No perturbation | `if(tree_seed != 0 && tree_noise > 0.0f)` |
| `consistency_anchors` | 0 | Off | `if(consistency_anchors > 0)` |
| `consistency_weight` | 2.0 | Default weight | Only used when consistency > 0 |
| `adaptive_budget` | 0 | Off | Direct flag check |
| `refine` | 0 (NONE) | No refinement | Switch/if chain |
| `realign` | 0 | No tree rebuild | `if(realign > 0)` triggers `kalign_run_realign` |

### Inconsistencies Found

1. **`dist_scale` sentinel handling**: In `kalign_run_seeded`, `dist_scale` is GUARDED (`if(dist_scale > 0.0f)`), so 0.0 means "keep default" (which is also 0.0 — benign). In `kalign_run_dist_scale`, `kalign_run_realign`, and `kalign_post_realign`, it's UNCONDITIONAL (`ap->dist_scale = dist_scale`). Effect: none in practice (default is 0.0), but inconsistent code.

2. **`use_seq_weights` in ensemble**: Both `kalign_ensemble` and `kalign_ensemble_custom` force `use_seq_weights = 0.0` when the sentinel `-1.0` is passed (lines 249-250, 545-546). But if an explicit positive value is passed, it flows through. Net effect: seq_weights is always 0 in ensemble mode through all current code paths, but it's not technically blocked at the C level for explicit positive values.

3. **Post-selection refinement mode**: In both ensemble functions, the post-selection refinement HARDCODES `KALIGN_REFINE_CONFIDENT` regardless of the `refine` parameter. The `refine` parameter is used for the per-run alignments, but the post-selection re-run always uses CONFIDENT. This is intentional but undocumented.

4. **`kalign_run_dist_scale` drops consistency**: This function doesn't accept consistency parameters at all. The pybind11 router can fall into this path when `vsm_amax` or `seq_weights` is set but `consistency=0` and `realign=0`. Currently benign (consistency=0 means no consistency anyway), but fragile.

5. **`ensemble_custom` missing per-run params**: `vsm_amax`, `use_seq_weights`, `consistency_anchors`, `consistency_weight`, `dist_scale` are all shared across runs in `kalign_ensemble_custom`. This is an arbitrary limitation — the optimizer would benefit from per-run control of all these.

## Proposed Design

### Core principle: every run is fully described by one config

Instead of "shared params + per-run overrides", each run gets its own complete config. For ensemble, you pass an array. No ambiguity about what's shared vs per-run.

### Per-run config struct

```c
/* kalign_config.h */

/* Describes everything needed for a single alignment run. */
struct kalign_run_config {
    /* Sequence type / substitution matrix */
    int type;                    /* KALIGN_TYPE_* constant (UNDEFINED = auto-detect) */

    /* Gap penalties (-1.0 = use matrix defaults) */
    float gpo;                   /* gap open penalty */
    float gpe;                   /* gap extend penalty */
    float tgpe;                  /* terminal gap extend penalty */

    /* Scoring modifiers */
    float vsm_amax;              /* variable scoring matrix amplitude (-1.0 = biotype default) */
    float dist_scale;            /* distance-dependent gap scaling (0.0 = off) */
    float use_seq_weights;       /* profile rebalancing pseudocount (-1.0 = biotype default) */

    /* Consistency transform */
    int consistency_anchors;     /* number of anchor sequences K (0 = off) */
    float consistency_weight;    /* bonus scale for consistency (default: 2.0) */

    /* Refinement */
    int refine;                  /* KALIGN_REFINE_* constant (default: NONE) */
    int adaptive_budget;         /* scale refinement trials by uncertainty (0 = off) */

    /* Realign (iterative tree rebuild) */
    int realign;                 /* number of realign iterations (0 = off) */

    /* Guide tree perturbation */
    uint64_t tree_seed;          /* random seed for noisy tree (0 = deterministic) */
    float tree_noise;            /* tree perturbation sigma (0.0 = none) */
};
```

### Ensemble config struct

```c
/* Controls orchestration when n_runs > 1. */
struct kalign_ensemble_config {
    uint64_t seed;               /* base RNG seed for diversity generation */
    int min_support;             /* POAR consensus threshold (0 = auto) */
    const char* save_poar;       /* path to save POAR table (NULL = don't save) */
};
```

### Initialization functions

```c
/* Returns a run config with all sentinel/default values. */
struct kalign_run_config kalign_run_config_defaults(void);

/* Returns an ensemble config with defaults. */
struct kalign_ensemble_config kalign_ensemble_config_defaults(void);
```

Default values for `kalign_run_config_defaults()`:
```c
struct kalign_run_config kalign_run_config_defaults(void) {
    struct kalign_run_config cfg;
    cfg.type = KALIGN_TYPE_UNDEFINED;
    cfg.gpo = -1.0f;              /* sentinel: use matrix default */
    cfg.gpe = -1.0f;
    cfg.tgpe = -1.0f;
    cfg.vsm_amax = -1.0f;         /* sentinel: use biotype default */
    cfg.dist_scale = 0.0f;        /* off */
    cfg.use_seq_weights = -1.0f;   /* sentinel: use biotype default */
    cfg.consistency_anchors = 0;   /* off */
    cfg.consistency_weight = 2.0f;
    cfg.refine = KALIGN_REFINE_NONE;
    cfg.adaptive_budget = 0;
    cfg.realign = 0;
    cfg.tree_seed = 0;             /* deterministic */
    cfg.tree_noise = 0.0f;         /* no perturbation */
    return cfg;
}
```

Default values for `kalign_ensemble_config_defaults()`:
```c
struct kalign_ensemble_config kalign_ensemble_config_defaults(void) {
    struct kalign_ensemble_config ens;
    ens.seed = 42;
    ens.min_support = 0;       /* auto: (n_runs + 2) / 3 */
    ens.save_poar = NULL;
    return ens;
}
```

### The unified function

```c
int kalign_align_full(
    struct msa* msa,
    const struct kalign_run_config* runs,     /* array of run configs */
    int n_runs,                               /* 1 = single-run, >1 = ensemble */
    const struct kalign_ensemble_config* ens,  /* NULL when n_runs == 1 */
    int n_threads
);
```

### Diversity table helper

The old `kalign_ensemble()` auto-generates per-run diversity (scaled gap penalties + tree noise) from a base config. This becomes a standalone helper that **generates** the run config array:

```c
/* Expand one base config into n_runs configs using the built-in diversity table.
   Copies base into each slot, then applies gap penalty scale factors and tree noise.
   Run 0 always gets the base config unchanged.
   Caller allocates out[n_runs]. */
int kalign_generate_ensemble_runs(
    const struct kalign_run_config* base,
    int n_runs,
    uint64_t seed,
    struct kalign_run_config* out
);
```

### Internal logic (pseudocode)

```
kalign_align_full(msa, runs, n_runs, ens, n_threads):
    essential_input_check(msa)
    detect_alphabet(msa)

    if n_runs > 1:
        // ENSEMBLE PATH
        for k = 0..n_runs-1:
            copy = deep_copy(msa)
            resolve sentinels in runs[k] (vsm_amax, seq_weights, gpo/gpe/tgpe)
            if runs[k].realign > 0:
                kalign_run_realign(copy, runs[k] params...)
            else:
                kalign_run_seeded(copy, runs[k] params...)
            extract POAR from alignment
            keep alignment for scoring

        score all alignments against POAR
        decide: consensus vs selection (based on ens->min_support)
        if selection: post-selection refinement (REFINE_CONFIDENT)
        copy winner back to msa
        compute confidence from POAR

    else:
        // SINGLE-RUN PATH
        resolve sentinels in runs[0]
        if runs[0].realign > 0:
            kalign_run_realign(msa, runs[0] params...)
        else:
            kalign_run_seeded(msa, runs[0] params...)

    sort by original rank
    return OK
```

## Usage Examples

### Single run with defaults

```c
struct kalign_run_config run = kalign_run_config_defaults();
kalign_align_full(msa, &run, 1, NULL, 4);
```

### Single run with custom params

```c
struct kalign_run_config run = kalign_run_config_defaults();
run.gpo = 8.0f;
run.vsm_amax = 2.0f;
run.consistency_anchors = 5;
run.realign = 1;
kalign_align_full(msa, &run, 1, NULL, 4);
```

### Ensemble with auto-diversity (replaces `kalign_ensemble`)

```c
struct kalign_run_config base = kalign_run_config_defaults();
base.vsm_amax = 2.0f;
base.realign = 1;

struct kalign_run_config runs[3];
kalign_generate_ensemble_runs(&base, 3, 42, runs);

struct kalign_ensemble_config ens = kalign_ensemble_config_defaults();
kalign_align_full(msa, runs, 3, &ens, 4);
```

### Ensemble with fully custom per-run params (optimizer use case)

```c
struct kalign_run_config runs[3];

runs[0] = kalign_run_config_defaults();
runs[0].gpo = 8.0f;  runs[0].vsm_amax = 2.0f;  runs[0].tree_noise = 0.0f;

runs[1] = kalign_run_config_defaults();
runs[1].gpo = 6.0f;  runs[1].vsm_amax = 1.5f;  runs[1].tree_noise = 0.3f;
runs[1].consistency_anchors = 5;  /* different consistency per run! */

runs[2] = kalign_run_config_defaults();
runs[2].gpo = 10.0f; runs[2].type = KALIGN_TYPE_CORBLOSUM;
runs[2].tree_noise = 0.5f;  runs[2].tree_seed = 12345;

struct kalign_ensemble_config ens = { .seed = 42, .min_support = 0, .save_poar = NULL };
kalign_align_full(msa, runs, 3, &ens, 4);
```

### CLI mode presets (using the config structs)

```c
/* --fast */
struct kalign_run_config run = kalign_run_config_defaults();
run.consistency_anchors = 0;
kalign_align_full(msa, &run, 1, NULL, n_threads);

/* default */
struct kalign_run_config run = kalign_run_config_defaults();
run.consistency_anchors = 5;
kalign_align_full(msa, &run, 1, NULL, n_threads);

/* --precise */
struct kalign_run_config base = kalign_run_config_defaults();
base.consistency_anchors = 5;
base.realign = 1;
struct kalign_run_config runs[3];
kalign_generate_ensemble_runs(&base, 3, 42, runs);
struct kalign_ensemble_config ens = kalign_ensemble_config_defaults();
kalign_align_full(msa, runs, 3, &ens, n_threads);
```

## Python API (Option A: single function, list of dicts)

### pybind11 binding

```cpp
// Single Python function that handles both single-run and ensemble
m.def("align_full", [](const std::string& input_path,
                        const std::string& output_path,
                        py::list run_configs,       // list of dicts
                        py::dict ensemble_config,   // {} for single-run
                        int n_threads,
                        const std::string& format) {
    int n_runs = py::len(run_configs);
    std::vector<kalign_run_config> runs(n_runs);
    for (int i = 0; i < n_runs; i++) {
        runs[i] = kalign_run_config_defaults();
        py::dict d = run_configs[i];
        if (d.contains("gpo"))       runs[i].gpo = d["gpo"].cast<float>();
        if (d.contains("vsm_amax"))  runs[i].vsm_amax = d["vsm_amax"].cast<float>();
        // ... etc for all fields, only override what's present
    }
    // ... build ensemble config if n_runs > 1, call kalign_align_full
});
```

### Python wrapper

```python
def align_full(input_path, output_path, *,
               runs,                      # list of dicts, one per run
               min_support=0,
               ensemble_seed=42,
               save_poar=None,
               expand_ensemble=None,      # int: generate N runs from runs[0] using diversity table
               n_threads=1,
               format="fasta"):
    """
    Unified alignment function.

    Args:
        runs: List of run config dicts. Each dict can contain:
            gpo, gpe, tgpe, vsm_amax, dist_scale, use_seq_weights,
            consistency_anchors, consistency_weight, refine, adaptive_budget,
            realign, tree_seed, tree_noise, type
            Missing keys use defaults.
        expand_ensemble: If set, takes runs[0] as base and generates
            this many runs using the built-in diversity table.
            Ignores runs[1:] if present.
        min_support: POAR consensus threshold (0 = auto). Only for ensemble.
        ensemble_seed: Base RNG seed for ensemble. Only for ensemble.
    """
```

### Python usage examples

**Single run:**
```python
kalign.align_full("in.fa", "out.fa",
    runs=[{"vsm_amax": 2.0, "consistency_anchors": 5}])
```

**Ensemble with diversity table (easy mode):**
```python
kalign.align_full("in.fa", "out.fa",
    runs=[{"vsm_amax": 2.0, "realign": 1}],
    expand_ensemble=3)
```

**Ensemble with fully custom per-run params (optimizer):**
```python
kalign.align_full("in.fa", "out.fa",
    runs=[
        {"gpo": 8.0, "vsm_amax": 2.0, "tree_noise": 0.0},
        {"gpo": 6.0, "vsm_amax": 1.5, "tree_noise": 0.3, "consistency_anchors": 5},
        {"gpo": 10.0, "type": "corblosum", "tree_noise": 0.5},
    ],
    min_support=0)
```

### Compatibility with pymoo optimizer

The optimizer's decision vector maps directly to the array-of-runs:

```python
def decode_to_runs(x, n_runs):
    """Map pymoo decision vector → list of run config dicts."""
    DIMS_PER_RUN = 12  # gpo, gpe, tgpe, vsm_amax, dist_scale, seq_weights,
                       # consistency_anchors, consistency_weight, refine,
                       # adaptive_budget, realign, tree_noise
    runs = []
    for i in range(n_runs):
        off = i * DIMS_PER_RUN
        runs.append({
            "gpo":                  x[off + 0],
            "gpe":                  x[off + 1],
            "tgpe":                 x[off + 2],
            "vsm_amax":             x[off + 3],
            "dist_scale":           x[off + 4],
            "use_seq_weights":      x[off + 5],
            "consistency_anchors":  int(x[off + 6]),
            "consistency_weight":   x[off + 7],
            "refine":               int(x[off + 8]),
            "adaptive_budget":      int(x[off + 9]),
            "realign":              int(x[off + 10]),
            "tree_noise":           x[off + 11],
        })
    return runs

def evaluate(x, input_path, n_threads):
    n_runs = int(x[60])  # from ensemble block at end of vector
    runs = decode_to_runs(x, n_runs)

    # Same call regardless of n_runs
    kalign.align_full(input_path, output_path,
                      runs=runs,
                      min_support=int(x[61]),
                      n_threads=n_threads)

    return f1, tc, wall_time  # three NSGA-II objectives
```

**Decision vector layout** (fixed size, `max_runs=5`):
```
┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐
│ runs[0]  │ runs[1]  │ runs[2]  │ runs[3]  │ runs[4]  │ ensemble │
│ 12 dims  │ 12 dims  │ 12 dims  │ 12 dims  │ 12 dims  │ 3 dims  │
└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
                                               ↑ masked if n_runs < 5
```

Benefits for NSGA-II:
- **No shared-vs-per-run ambiguity** — each run block is self-contained
- **Run-level crossover** — swap entire 12-dim run blocks between individuals
- **Clean masking** — unused run blocks are simply ignored
- **One code path** — evaluation function doesn't branch on mode

## What gets removed

| Function | Action |
|----------|--------|
| `kalign_run()` | Deprecated wrapper: `cfg = defaults(); kalign_align_full(msa, &cfg, 1, NULL, n_threads)` |
| `kalign_run_seeded()` | **Keep as internal helper** — called by `kalign_align_full` for each single run |
| `kalign_run_dist_scale()` | **Remove entirely** — redundant subset of `kalign_run_seeded` |
| `kalign_run_realign()` | **Keep as internal helper** — called by `kalign_align_full` when `realign > 0` |
| `kalign_ensemble()` | Deprecated wrapper: calls `kalign_generate_ensemble_runs` + `kalign_align_full` |
| `kalign_ensemble_custom()` | Deprecated wrapper: builds run config array + `kalign_align_full` |
| `kalign_post_realign()` | **Remove** — unused, was found to hurt performance |
| `kalign()` | Deprecated wrapper: arr→msa conversion + `kalign_align_full` |

**Important**: `kalign_run_seeded` and `kalign_run_realign` remain as internal implementation — they do the actual alignment work. `kalign_align_full` orchestrates WHEN to call them and with WHAT parameters.

## Sentinel handling (unified)

All sentinel resolution happens once in `kalign_align_full()`, before passing params to internal helpers:

```c
/* For each run config: */
for (int k = 0; k < n_runs; k++) {
    resolved[k] = runs[k];  /* copy */

    /* Resolve vsm_amax sentinel */
    if (resolved[k].vsm_amax < 0.0f)
        resolved[k].vsm_amax = (biotype == PROTEIN) ? 2.0f : 0.0f;

    /* Resolve use_seq_weights sentinel */
    if (resolved[k].use_seq_weights < 0.0f)
        resolved[k].use_seq_weights = 0.0f;  /* biotype default is 0.0 for all */

    /* gpo/gpe/tgpe: -1.0 sentinel passed through to aln_param_init which handles it */
    /* dist_scale: 0.0 = off, always passed through (no sentinel) */
}
```

## Ensemble post-selection refinement

Currently hardcoded to `KALIGN_REFINE_CONFIDENT` in both ensemble functions. In the unified design, this stays hardcoded with a clear comment — it's always beneficial and changing it has not been tested. If we ever want to make it configurable, it would go into `kalign_ensemble_config`.

## Backward compatibility

Old function signatures remain in `kalign.h` as **deprecated wrappers** that construct run configs and call `kalign_align_full()`:

1. Existing C code that calls `kalign_run()` or `kalign_ensemble()` continues to work unchanged.
2. The wrappers can be removed in a future major version.
3. The old `kalign()` function (char** API) remains for external consumers.

## Testing Strategy

### Phase 1: Verify behavioral equivalence

For each old function, create a test that:
1. Runs alignment with the old function
2. Runs the same alignment with `kalign_align_full` using equivalent config
3. Asserts the output alignments are IDENTICAL (byte-for-byte)

Test cases from `tests/data/`:
- BB11001 (small, protein)
- BB12006 (medium, protein)
- BB30014 (protein with insertions)

Test matrix:
```
kalign_run(gpo=-1, gpe=-1, tgpe=-1, refine=NONE, adaptive=0)
  == kalign_align_full(&defaults, 1, NULL, n_threads)

kalign_run_seeded(gpo=8.0, gpe=1.0, tgpe=0.5, vsm_amax=2.0, seq_weights=1.0, consistency=5)
  == kalign_align_full(&{same params}, 1, NULL, n_threads)

kalign_run_dist_scale(dist_scale=0.5, vsm_amax=2.0, seq_weights=1.0)
  == kalign_align_full(&{same params, consistency=0}, 1, NULL, n_threads)

kalign_run_realign(realign=2, consistency=5, vsm_amax=2.0)
  == kalign_align_full(&{same params}, 1, NULL, n_threads)

kalign_ensemble(n_runs=3, seed=42, refine=CONFIDENT, vsm_amax=2.0, realign=1)
  == kalign_generate_ensemble_runs(&base, 3, 42, runs) + kalign_align_full(runs, 3, &ens, n_threads)

kalign_ensemble_custom(n_runs=3, run_gpo=[...], ...)
  == kalign_align_full(custom_runs, 3, &ens, n_threads)
```

### Phase 2: Benchmark regression check

Run BAliBASE 218 cases with:
- Old CLI path → scores
- New `kalign_align_full` path → scores
- Assert identical F1/TC (to 6 decimal places)

### Phase 3: Python API regression

Run existing Python test suite (`tests/python/`). All C tests + Python tests must pass.

## Implementation Order

1. Add `kalign_run_config`, `kalign_ensemble_config` structs and `kalign_*_defaults()` to new header `kalign_config.h`
2. Add `kalign_generate_ensemble_runs()` — extracts the existing diversity table into a config generator
3. Implement `kalign_align_full()` in `aln_wrap.c`, calling existing `kalign_run_seeded` and `kalign_run_realign` internally
4. Move ensemble orchestration from `ensemble.c` into `kalign_align_full()` (or keep as helper called by it)
5. Write equivalence tests (Phase 1)
6. Convert CLI to use `kalign_align_full()` — verify all existing tests pass
7. Convert pybind11 to single `align_full` function — verify Python tests pass
8. Run BAliBASE regression (Phase 2)
9. Add deprecated wrappers for old functions
10. Remove `kalign_run_dist_scale` and `kalign_post_realign`
11. Expose `align_full` in Python wrapper (`__init__.py`) with `expand_ensemble` convenience param

## Files Modified

| File | Change |
|------|--------|
| `lib/include/kalign/kalign_config.h` | **NEW**: `kalign_run_config`, `kalign_ensemble_config`, defaults functions |
| `lib/include/kalign/kalign.h` | Add `kalign_align_full()`, `kalign_generate_ensemble_runs()`. Mark old functions deprecated. |
| `lib/src/aln_wrap.c` | Implement `kalign_align_full()`. Convert old functions to wrappers. Remove `kalign_run_dist_scale`. |
| `lib/src/ensemble.c` | Extract diversity table into `kalign_generate_ensemble_runs()`. Move orchestration to `kalign_align_full()`. |
| `src/run_kalign.c` | Replace if/else chain with config struct + `kalign_align_full()` |
| `python-kalign/_core.cpp` | Replace router + `ensemble_custom_file_to_file` with single `align_full` binding |
| `python-kalign/__init__.py` | Add `align_full()` wrapper with `expand_ensemble` convenience |
| `tests/` | Add equivalence tests |

## Risks

1. **Ensemble diversity table**: The old `kalign_ensemble` uses a hardcoded scale-factor table. `kalign_generate_ensemble_runs` must reproduce this exactly. We keep the table as-is and just expose it as a config generator.

2. **Post-selection refinement**: Currently hardcodes `REFINE_CONFIDENT`. Keeping this hardcoded is safe; making it configurable would need testing.

3. **Thread safety**: `kalign_run_config` is a plain value struct with no heap pointers. `kalign_ensemble_config` has `save_poar` (read-only string) and `seed`. Multiple threads can safely use different configs.

4. **Binary compatibility**: Adding new functions and structs is fine. Old functions become wrappers → no ABI break. Can be removed in a future major version.

5. **Struct size stability**: If we add fields to `kalign_run_config` later, code that uses `kalign_run_config_defaults()` gets the new defaults automatically. Code that initializes with `= {0}` or partial initialization would get zeros for new fields, which should be safe (sentinels and "off" values are 0 or -1.0).
