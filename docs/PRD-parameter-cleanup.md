# PRD: Kalign Parameter Architecture Cleanup

> **Status note:** this document captures the original three-mode design
> (`fast` / `default` / `accurate`). A fourth mode, `recall`, was added
> during implementation. The architecture below is otherwise current;
> mentally substitute "four modes × three biotypes" for the preset grid.

## Goal

Replace the current parameter specification system — which conflates sequence
types with matrix choices, uses sentinel values resolved at multiple depths, and
exposes internals to users — with a clean two-layer design:

1. **User layer**: one parameter (`mode`) that selects a fully optimised preset.
2. **Engine layer**: one function (`kalign_align_full`) that takes an array of
   fully concrete run configurations.

The optimizer populates the engine layer directly. Users never see it.

---

## 1. User-Facing API

### Python

```python
# The only parameter most users ever need:
result = kalign.align(sequences, mode="default")
result = kalign.align(sequences, mode="fast")
result = kalign.align(sequences, mode="accurate")

# Force biotype when auto-detection is wrong (rare):
result = kalign.align(sequences, mode="default", seq_type="protein")

# Threading:
result = kalign.align(sequences, mode="default", n_threads=4)
```

### Parameters exposed to users

| Parameter              | Type          | Values                              | Default     |
|------------------------|---------------|-------------------------------------|-------------|
| `mode`                 | str           | `"fast"`, `"default"`, `"accurate"` | `"default"` |
| `seq_type`             | str           | `"auto"`, `"dna"`, `"rna"`, `"protein"` | `"auto"` |
| `n_threads`            | int           | >= 1                                | 1           |
| `gap_open`             | float or None | penalty value                       | None        |
| `gap_extend`           | float or None | penalty value                       | None        |
| `terminal_gap_extend`  | float or None | penalty value                       | None        |

Mode handles everything for the common case. The gap penalty parameters
exist for backward compatibility and expert use.

### Gap penalty override rule

When the user provides **any** gap penalty (`gap_open`, `gap_extend`, or
`terminal_gap_extend`), the following happens:

1. The `mode` parameter is **ignored** (regardless of what was passed).
2. The `fast` preset for the detected biotype is loaded (single run).
3. The user's gap penalties **replace** the corresponding preset values.
4. Unspecified penalties keep the `fast` preset defaults.

This gives expert users direct control over gap scoring while keeping
everything else (matrix, VSM, seq_weights, etc.) at sensible optimised
values.

```python
# Uses default mode (5-run ensemble):
kalign.align(sequences)

# Uses fast preset, overrides gap open only:
kalign.align(sequences, gap_open=5.0)

# Uses fast preset, overrides all three:
kalign.align(sequences, gap_open=5.0, gap_extend=0.5,
             terminal_gap_extend=0.3)

# mode is ignored when gap penalties are set:
kalign.align(sequences, mode="accurate", gap_open=5.0)  # → fast + gpo=5.0
```

### File-based variants

```python
# Read from file, return AlignedSequences (names + sequences)
result = kalign.align_from_file("input.fasta", mode="default")

# Read from file, write to file
kalign.align_file_to_file("input.fasta", "output.fasta",
                          mode="default", format="fasta")
```

`align_file_to_file` adds one extra parameter:

| Parameter | Type | Values                          | Default   |
|-----------|------|---------------------------------|-----------|
| `format`  | str  | `"fasta"`, `"msf"`, `"clu"`    | `"fasta"` |

### Output formats

- `align()` returns `List[str]` (aligned sequences) or ecosystem objects
  via `fmt="biopython"` / `fmt="skbio"`.
- `align_from_file()` returns `AlignedSequences` (names + sequences +
  optional confidence).
- `align_file_to_file()` writes to disk, returns nothing.

### CLI

```bash
kalign -i input.fasta -o output.fasta                    # mode=default
kalign -i input.fasta -o output.fasta --mode fast
kalign -i input.fasta -o output.fasta --mode accurate
kalign -i input.fasta -o output.fasta --type protein     # force biotype
kalign -i input.fasta -o output.fasta --gpo 5.0          # fast + gpo override
```

---

## 2. Mode Presets

Each mode is a static lookup table of fully concrete run configurations. The
number of runs, the matrix per run, the gap penalties per run, the tree seed
per run — everything is baked in. No computation, no scaling, no indirection.

Presets are **biotype-specific**: protein, DNA, and RNA each have their own
optimised configurations.

### Mode semantics

Three modes trade speed against accuracy. The user's only decision is where
on that tradeoff they want to be.

| Mode       | Intent                                              |
|------------|-----------------------------------------------------|
| `fast`     | Fastest possible alignment, acceptable quality       |
| `default`  | Good balance of speed and accuracy for routine use   |
| `accurate` | Best achievable accuracy, speed secondary            |

### Preset grid

There are **9 preset slots**: 3 modes × 3 biotypes.

|            | Protein | DNA | RNA |
|------------|---------|-----|-----|
| `fast`     | slot    | slot| slot|
| `default`  | slot    | slot| slot|
| `accurate` | slot    | slot| slot|

Each slot is filled by the optimizer independently. The number of runs,
the choice of matrices, the gap penalties, whether realignment is used —
all of this is determined per slot by multi-objective optimization, not
prescribed by the PRD.

Constraints on the optimizer:
- `n_runs` must be 1, 3, or 5 (no other ensemble sizes).
- `fast` should be Pareto-optimal for speed; `accurate` for quality;
  `default` for balanced.
- Within each biotype, `fast` must be strictly faster than `default`,
  and `default` strictly faster than `accurate`.

### Preset structure (per mode × biotype)

Each preset specifies:

**Ensemble-level parameters** (structural — not per-run):

| Field          | Description                                     |
|----------------|-------------------------------------------------|
| `n_runs`       | Number of alignment runs (1, 3, or 5)           |
| `min_support`  | POAR consensus column threshold (ensemble only) |

These are inherently about the ensemble as a whole, not about any
individual alignment run.

**Per-run parameters** (each run in the ensemble has its own values):

| Field          | Description                                     |
|----------------|-------------------------------------------------|
| `matrix`       | Substitution matrix for this run                |
| `gpo`          | Gap open penalty                                |
| `gpe`          | Gap extend penalty                              |
| `tgpe`         | Terminal gap extend penalty                     |
| `vsm_amax`     | Variable scoring matrix amplitude (0 = off)     |
| `seq_weights`  | Profile rebalancing pseudo-count (0 = off)      |
| `dist_scale`   | Distance-dependent gap scaling (0 = off)        |
| `realign`      | Alignment-guided tree-rebuild iterations (0+)   |
| `refine`       | Post-alignment refinement strategy              |
| `adaptive_budget` | Scale refinement trials by uncertainty (0=off)|
| `tree_seed`    | RNG seed for guide tree construction            |
| `tree_noise`   | Guide tree perturbation sigma                   |

The optimizer may choose identical values across runs for some
parameters (e.g. the same `vsm_amax` for all 5 runs) or vary them
(e.g. `realign=2` on run 1 but `realign=0` on run 3). That is the
optimizer's decision, not an architectural constraint.

This is the **complete parameter space** exposed to the optimizer.

---

## 3. Engine Layer (C API)

### Substitution matrix constants

```c
#define KALIGN_MATRIX_AUTO         0  /* auto-select for biotype          */
#define KALIGN_MATRIX_PFASUM43     1  /* 1/3 bit, divergent protein       */
#define KALIGN_MATRIX_PFASUM60     2  /* 1/3 bit, moderate protein        */
#define KALIGN_MATRIX_CORBLOSUM66  3  /* 1/3 bit, close protein           */
#define KALIGN_MATRIX_DNA          4  /* DNA match/mismatch (+5/-4)       */
#define KALIGN_MATRIX_DNA_INTERNAL 5  /* DNA internal (tgpe=8)            */
#define KALIGN_MATRIX_RNA          6  /* RNA RIBOSUM-like (~160-383)      */
```

Every value maps to exactly one scoring table. No duplicates.

GONNET (gon250) remains in the source as dead code but is not assigned a
constant and is unreachable through any API.

### Matrix default penalties

Each matrix has intrinsic default gap penalties. These are used as starting
points when constructing configs, never resolved at alignment time.

| Matrix        | gpo    | gpe   | tgpe   | Score range   |
|---------------|--------|-------|--------|---------------|
| PFASUM43      | 7.0    | 1.25  | 1.0    | -6 to 13      |
| PFASUM60      | 7.0    | 1.25  | 1.0    | -6 to 14      |
| CorBLOSUM66   | 5.5    | 2.0   | 1.0    | -4 to 13      |
| DNA           | 8.0    | 6.0   | 0.0    | -4 to 5       |
| DNA_INTERNAL  | 8.0    | 6.0   | 8.0    | -4 to 5       |
| RNA           | 217.0  | 39.4  | 292.6  | ~160 to 383   |

PFASUM43, PFASUM60, and CorBLOSUM66 are all in 1/3-bit units. Their gap
penalties are directly comparable. A penalty value like `gpo=7.0` means the
same thing across all three matrices.

### Refinement constants

```c
#define KALIGN_REFINE_NONE      0  /* no post-alignment refinement       */
#define KALIGN_REFINE_ALL       1  /* refine all columns (two-pass)      */
#define KALIGN_REFINE_CONFIDENT 2  /* refine high-confidence columns     */
#define KALIGN_REFINE_INLINE    3  /* per-node refinement during tree    */
```

### Run config struct

```c
struct kalign_run_config {
    int matrix;            /* KALIGN_MATRIX_*                            */
    float gpo;             /* gap open penalty (concrete, no sentinel)   */
    float gpe;             /* gap extend penalty                         */
    float tgpe;            /* terminal gap extend penalty                */
    float vsm_amax;        /* variable scoring matrix amplitude (0=off) */
    float seq_weights;     /* profile rebalancing pseudo-count (0=off)  */
    float dist_scale;      /* distance-dependent gap scaling (0=off)    */
    int refine;            /* KALIGN_REFINE_*                            */
    int adaptive_budget;   /* scale refinement by uncertainty (0=off)   */
    int realign;           /* tree-rebuild iterations (0=none)           */
    uint64_t tree_seed;    /* guide tree RNG seed                        */
    float tree_noise;      /* guide tree perturbation sigma (0=none)    */
};
```

No sentinel values. Every field is a concrete, usable value. A config
obtained from `kalign_run_config_defaults()` or `kalign_get_mode_preset()`
is directly usable without further resolution.

### Ensemble config struct

```c
struct kalign_ensemble_config {
    int min_support;       /* POAR consensus column threshold            */
};
```

Simplified from current struct. `seed` removed (each run has its own
`tree_seed`). `save_poar` removed (debug feature, not part of core API).

### Preset function

```c
int kalign_get_mode_preset(const char *mode,
                           int biotype,
                           struct kalign_run_config *runs,
                           int *n_runs,
                           struct kalign_ensemble_config *ens);
```

- `mode`: `"fast"`, `"default"`, `"accurate"` (case-insensitive).
- `biotype`: `ALN_BIOTYPE_PROTEIN`, `ALN_BIOTYPE_DNA`, or `ALN_BIOTYPE_RNA`.
  Determined by auto-detection before this call.
- `runs`: caller-allocated array (minimum 8 elements).
- `n_runs`: filled with the number of runs in the preset.
- `ens`: filled with ensemble config.
- Returns 0 on success, -1 on unknown mode.

### Alignment entry point

```c
int kalign_align_full(struct msa *msa,
                      const struct kalign_run_config *runs,
                      int n_runs,
                      const struct kalign_ensemble_config *ens,
                      int n_threads);
```

This is the single entry point for all alignment. It receives fully concrete
configs and executes them. No parameter resolution, no sentinel handling.

If `n_runs == 1`: single alignment using `runs[0]`.
If `n_runs > 1`: ensemble — run each config, build POAR, consensus.

When `runs[k].matrix == KALIGN_MATRIX_AUTO`:
- If `msa->biotype == DNA`: resolves to `KALIGN_MATRIX_DNA`, uses DNA
  default penalties (overriding the config's gpo/gpe/tgpe).
- If `msa->biotype == RNA`: resolves to `KALIGN_MATRIX_RNA`, uses RNA
  default penalties.
- If `msa->biotype == protein`: resolves to PFASUM43 or PFASUM60 using
  the length-ratio heuristic (ratio < 1.5 → PFASUM43, else PFASUM60).
  Keeps the config's gpo/gpe/tgpe.

This is the **only** place where AUTO is resolved, and it happens **once**
before the alignment loop.

---

## 4. Parameters Exposed to the Optimizer

The optimizer sees the engine layer directly. It constructs
`kalign_run_config` arrays and calls `kalign_align_full`.

### Optimizer search space

**Ensemble-level parameters** (one value per preset slot):

| Parameter      | Type   | Range              | Description                    |
|----------------|--------|--------------------|--------------------------------|
| `n_runs`       | int    | {1, 3, 5}          | Number of ensemble runs        |
| `min_support`  | int    | [0, n_runs]         | POAR consensus threshold       |

**Per-run parameters** (optimised independently for each of the n_runs):

| Parameter      | Type   | Range              | Description                    |
|----------------|--------|--------------------|--------------------------------|
| `matrix`          | int    | {1, 2, 3}          | PFASUM43, PFASUM60, CorBLOSUM66|
| `gpo`             | float  | [1.0, 20.0]        | Gap open penalty               |
| `gpe`             | float  | [0.1, 5.0]         | Gap extend penalty             |
| `tgpe`            | float  | [0.1, 5.0]         | Terminal gap extend penalty    |
| `vsm_amax`        | float  | [0.0, 3.0]         | VSM amplitude (0 = disabled)   |
| `seq_weights`     | float  | [0.0, 3.0]         | Pseudo-count (0 = disabled)    |
| `dist_scale`      | float  | [0.0, 2.0]         | Distance gap scaling (0 = off) |
| `realign`         | int    | {0, 1, 2, 3}       | Tree-rebuild iterations        |
| `refine`          | int    | {0, 1, 2, 3}       | Refinement strategy            |
| `adaptive_budget` | int    | {0, 1}             | Scale refinement by uncertainty|
| `tree_seed`       | uint64 | fixed per run index | Deterministic seed             |
| `tree_noise`      | float  | [0.0, 0.5]         | Tree perturbation              |

The matrix values {1, 2, 3} correspond to the three protein matrices on
the same 1/3-bit scale, so a single set of gap penalty ranges works for
all of them.

For DNA/RNA optimisation, the matrix field is fixed (DNA or RNA) and the
penalty ranges are adjusted to match those matrices' scales.

The optimizer may choose to use the same value for a parameter across all
runs (e.g. `vsm_amax=2.0` for every run) or vary it per run. That is a
search strategy decision, not an architectural constraint.

### What the optimizer produces

A JSON file with this structure:

```json
{
  "protein": {
    "fast": {
      "n_runs": 1,
      "min_support": 0,
      "runs": [
        {
          "matrix": "pfasum60",
          "gpo": 8.4087,
          "gpe": 0.5153,
          "tgpe": 0.4927,
          "vsm_amax": 1.448,
          "seq_weights": 1.063,
          "dist_scale": 0.0,
          "realign": 0,
          "refine": 0,
          "adaptive_budget": 0,
          "tree_seed": 42,
          "tree_noise": 0.1623
        }
      ]
    },
    "default": {
      "n_runs": 5,
      "min_support": 3,
      "runs": [ { ... }, { ... }, { ... }, { ... }, { ... } ]
    },
    "accurate": {
      "n_runs": 5,
      "min_support": 3,
      "runs": [ { ... }, { ... }, { ... }, { ... }, { ... } ]
    }
  },
  "dna": {
    "fast": { ... },
    "default": { ... },
    "accurate": { ... }
  },
  "rna": {
    "fast": { ... },
    "default": { ... },
    "accurate": { ... }
  }
}
```

Each run object contains every per-run parameter. The JSON maps 1:1 to the
C preset tables. No interpretation, no transformation. The values in the
JSON are the values used in alignment.

---

## 5. Removed Parameters

The following are **removed from the user-facing API**. They remain
accessible only through the engine layer (`kalign_align_full`) for
benchmarking and optimiser use.

| Removed parameter      | Reason                                        |
|------------------------|-----------------------------------------------|
| `ensemble`             | Set by mode preset (1, 3, or 5)               |
| `ensemble_seed`        | Replaced by per-run `tree_seed`               |
| `matrix` / `seq_type`  | Set by mode preset per run                    |
| `vsm_amax`             | Set by mode preset                            |
| `seq_weights`          | Set by mode preset                            |
| `realign`              | Set by mode preset                            |
| `refine`               | Set by mode preset                            |
| `min_support`          | Set by mode preset                            |
| `consistency`          | Set by mode preset (currently always 0)       |
| `consistency_weight`   | Set by mode preset (dead when consistency=0)  |
| `dist_scale`           | Set by mode preset                            |
| `adaptive_budget`      | Set by mode preset                            |
| `save_poar`            | Debug feature, not user-facing                |
| `load_poar`            | Debug feature, not user-facing                |

---

## 6. Backward Compatibility

### C header

Old `KALIGN_TYPE_*` constants become `#define` aliases:

```c
/* Deprecated — use KALIGN_MATRIX_* */
#define KALIGN_TYPE_DNA              KALIGN_MATRIX_DNA
#define KALIGN_TYPE_DNA_INTERNAL     KALIGN_MATRIX_DNA_INTERNAL
#define KALIGN_TYPE_RNA              KALIGN_MATRIX_RNA
#define KALIGN_TYPE_PROTEIN          KALIGN_MATRIX_PFASUM43
#define KALIGN_TYPE_PROTEIN_PFASUM43 KALIGN_MATRIX_PFASUM43
#define KALIGN_TYPE_PROTEIN_PFASUM60 KALIGN_MATRIX_PFASUM60
#define KALIGN_TYPE_PROTEIN_PFASUM_AUTO KALIGN_MATRIX_AUTO
#define KALIGN_TYPE_UNDEFINED        KALIGN_MATRIX_AUTO
#define KALIGN_TYPE_PROTEIN_CORBLOSUM66 KALIGN_MATRIX_CORBLOSUM66
```

`KALIGN_TYPE_PROTEIN_DIVERGENT` gets no alias — it mapped to GONNET
which is dead code. Any code using it gets a compile error, directing
the author to pick a specific matrix.

Old functions (`kalign_run`, `kalign_run_seeded`, `kalign_ensemble`,
etc.) remain as thin wrappers that build a config and call
`kalign_align_full`. They accept -1.0 sentinels for backward compat,
resolving them into the config immediately. New code should not use them.

### Python

`align()`, `align_from_file()`, and `align_file_to_file()` continue to
accept all current keyword arguments during a deprecation period. When
old-style parameters are used, a `DeprecationWarning` is raised. The
old parameters are resolved into a run config and passed to the engine.

After the deprecation period, the signatures shrink to:

```python
def align(sequences, *, mode="default", seq_type="auto",
          gap_open=None, gap_extend=None, terminal_gap_extend=None,
          n_threads=None, fmt="plain", ids=None)

def align_from_file(input_file, *, mode="default", seq_type="auto",
                    gap_open=None, gap_extend=None,
                    terminal_gap_extend=None, n_threads=None)

def align_file_to_file(input_file, output_file, *, mode="default",
                       seq_type="auto", gap_open=None, gap_extend=None,
                       terminal_gap_extend=None, format="fasta",
                       n_threads=None)
```

Gap penalty parameters are retained permanently (not deprecated). When
any is set, the gap penalty override rule applies: mode is ignored, the
`fast` preset is loaded, and user penalties replace preset defaults.

---

## 7. Implementation Order

1. **Add `KALIGN_MATRIX_*` constants** to `kalign.h`. Add backward-compat
   aliases for `KALIGN_TYPE_*`. Add `biotype` parameter to
   `kalign_get_mode_preset`.

2. **Rename `type` → `matrix`** in `kalign_run_config`. Update all code
   that reads the field.

3. **Remove sentinels from the new path.** `kalign_run_config_defaults()`
   returns concrete PFASUM43 protein values. `aln_param_init` uses values
   directly (no `if(gpo >= 0)` guard). Old functions keep sentinel
   resolution internally.

4. **Fix v2 presets.** Change all `KALIGN_TYPE_PROTEIN_DIVERGENT` →
   `KALIGN_MATRIX_PFASUM43` in `kalign_get_mode_preset()`. These runs
   were always PFASUM43; the "gonnet" label was a bug.

5. **Add DNA/RNA preset stubs** in `kalign_get_mode_preset()`. Initially
   use matrix defaults; optimise later.

6. **Simplify Python API.** New signatures with `mode` + `seq_type` +
   `n_threads`. Old parameters accepted with deprecation warnings.
   Remove `_MODE_PRESETS` dict and `-1.0` sentinel logic from Python.

7. **Update optimizer** to use `KALIGN_MATRIX_*` constants and produce
   the JSON format described in section 4.

8. **Benchmark** all three modes × all three biotypes. Verify protein
   results match optimizer predictions now that the matrix bug is fixed.

9. **Remove deprecated parameters** after one release cycle.
