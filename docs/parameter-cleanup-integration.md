# Kalign Python API: Integration Guide

This document describes kalign's Python API architecture and how to use
it from benchmark runners and the NSGA-III optimizer.

---

## 1. Architecture: Two Paths

Kalign's Python bindings have exactly two code paths:

### Path 1: Mode-based (public API)

For standard alignment using preset configurations:

```python
import kalign

# In-memory
result = kalign.align(sequences, mode="default")
result = kalign.align(sequences, mode="fast")
result = kalign.align(sequences, mode="accurate")

# File-based (returns AlignedSequences with names + sequences)
result = kalign.align_from_file("input.fasta", mode="default")

# File-to-file
kalign.align_file_to_file("input.fasta", "output.fasta", mode="default")
```

These call `_core.align_mode()` / `_core.align_from_file_mode()` /
`_core.align_file_to_file_mode()`, which delegate to the C function
`kalign_get_mode_preset()` → `kalign_align_full()`.

**Mode presets** are biotype-aware. Each mode × biotype (protein / DNA / RNA)
returns fully concrete per-run configurations optimized by NSGA-III (protein)
or using sensible defaults (DNA/RNA).

**Gap penalty override rule**: If any gap penalty (`gap_open`, `gap_extend`,
`terminal_gap_extend`) is set, mode is forced to `"fast"` and the specified
penalties replace the preset values. Value `-1.0` means "use preset default".

### Path 2: Optimizer (internal API)

For the NSGA-III optimizer that needs full per-run control:

```python
from kalign._core import ensemble_custom_file_to_file

ensemble_custom_file_to_file(
    input_file="input.fasta",
    output_file="output.fasta",
    run_gpo=[8.4, 6.2, 7.1],        # per-run gap open
    run_gpe=[0.5, 1.0, 0.8],        # per-run gap extend
    run_tgpe=[0.5, 0.9, 0.6],       # per-run terminal gap extend
    run_noise=[0.16, 0.25, 0.10],    # per-run tree noise sigma
    run_types=[1, 2, 3],             # per-run matrix constants
    format="fasta",
    seq_type=1,                      # KALIGN_MATRIX_PFASUM43
    seed=42,                         # base seed (run k gets seed+k)
    min_support=0,                   # POAR consensus threshold
    refine=2,                        # KALIGN_REFINE_CONFIDENT
    vsm_amax=2.0,                    # shared across runs
    realign=1,                       # shared across runs
    seq_weights=1.0,                 # shared across runs
    n_threads=4,
    consistency_anchors=0,           # shared across runs
    consistency_weight=2.0,          # shared across runs
)
```

This builds a `kalign_run_config[]` array and calls `kalign_align_full()`.

---

## 2. C Entry Point

```c
int kalign_align_full(struct msa *msa,
                      const struct kalign_run_config *runs,
                      int n_runs,
                      const struct kalign_ensemble_config *ens,
                      int n_threads);
```

Single entry point for all alignment. Receives fully concrete configs
(no sentinel values) and executes them.

### `kalign_run_config` (14 fields)

```c
struct kalign_run_config {
    int matrix;                  /* KALIGN_MATRIX_* constant             */
    float gpo;                   /* gap open penalty                     */
    float gpe;                   /* gap extend penalty                   */
    float tgpe;                  /* terminal gap extend penalty          */
    float vsm_amax;              /* variable scoring matrix amp (0=off)  */
    float seq_weights;           /* profile rebalancing (0=off)          */
    float dist_scale;            /* distance-dependent gap scaling       */
    int refine;                  /* KALIGN_REFINE_* constant             */
    int adaptive_budget;         /* scale refinement by uncertainty      */
    int realign;                 /* iterative tree-rebuild iters (0=off) */
    uint64_t tree_seed;          /* random seed for tree perturbation    */
    float tree_noise;            /* tree perturbation sigma (0.0=none)   */
    int consistency_anchors;     /* anchor consistency rounds (0=off)    */
    float consistency_weight;    /* anchor consistency bonus (default 2) */
};
```

### `kalign_ensemble_config`

```c
struct kalign_ensemble_config {
    int min_support;             /* POAR consensus threshold (0=auto)    */
};
```

### Mode presets

```c
int kalign_get_mode_preset(const char *mode,    /* "fast"/"default"/"accurate" */
                           int biotype,         /* ALN_BIOTYPE_PROTEIN or _DNA */
                           struct kalign_run_config *runs,  /* out: array[8] */
                           int *n_runs,         /* out */
                           struct kalign_ensemble_config *ens);  /* out */
```

---

## 3. Matrix Constants

| Constant                   | Value | Description         |
|---------------------------|-------|---------------------|
| `KALIGN_MATRIX_AUTO`       | 0     | Auto-detect         |
| `KALIGN_MATRIX_PFASUM43`   | 1     | PFASUM43 (protein)  |
| `KALIGN_MATRIX_PFASUM60`   | 2     | PFASUM60 (protein)  |
| `KALIGN_MATRIX_CORBLOSUM66`| 3     | CorBLOSUM66 (protein) |
| `KALIGN_MATRIX_DNA`        | 5     | DNA                 |
| `KALIGN_MATRIX_DNA_INTERNAL`| 6    | DNA internal        |
| `KALIGN_MATRIX_RNA`        | 7     | RNA                 |

Python access:
```python
from kalign._core import (
    MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66,
    MATRIX_DNA, MATRIX_DNA_INTERNAL, MATRIX_RNA, MATRIX_AUTO,
)
```

Legacy names (`PROTEIN`, `DNA`, etc.) still work as aliases.

### Critical bug fix: value 3 was GONNET, now CorBLOSUM66

Old code used `PROTEIN` (value 3) thinking it was a generic protein
matrix, but it actually triggered the GONNET scoring path (gon250:
scores ~0-200). This is now `MATRIX_CORBLOSUM66` (CorBLOSUM66,
1/3-bit units, scores -4 to 13).

The correct three-matrix mapping for protein optimization:

```python
from kalign._core import MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66

matrix_map_int = [MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66]
matrix_map_str = ["pfasum43", "pfasum60", "corblosum66"]
```

---

## 4. Python API Signatures

### `kalign.align()`

```python
kalign.align(
    sequences: list[str],
    seq_type: str | int = "auto",      # "auto", "dna", "rna", "protein"
    gap_open: float | None = None,     # overrides mode → forces "fast"
    gap_extend: float | None = None,
    terminal_gap_extend: float | None = None,
    n_threads: int | None = None,
    mode: str | None = None,           # "fast", "default", "accurate"
    fmt: str = "plain",                # "plain", "biopython", "skbio"
    ids: list[str] | None = None,
) -> list[str]
```

### `kalign.align_from_file()`

```python
kalign.align_from_file(
    input_file: str,
    seq_type: str | int = "auto",
    gap_open: float | None = None,
    gap_extend: float | None = None,
    terminal_gap_extend: float | None = None,
    n_threads: int | None = None,
    mode: str | None = None,
) -> AlignedSequences  # unpacks as (names, sequences)
```

### `kalign.align_file_to_file()`

```python
kalign.align_file_to_file(
    input_file: str,
    output_file: str,
    format: str = "fasta",
    seq_type: str | int = "auto",
    gap_open: float | None = None,
    gap_extend: float | None = None,
    terminal_gap_extend: float | None = None,
    n_threads: int | None = None,
    mode: str | None = None,
)
```

### `_core.ensemble_custom_file_to_file()` (optimizer only)

```python
from kalign._core import ensemble_custom_file_to_file

ensemble_custom_file_to_file(
    input_file: str,
    output_file: str,
    run_gpo: list[float],          # per-run
    run_gpe: list[float],          # per-run
    run_tgpe: list[float],         # per-run
    run_noise: list[float],        # per-run
    run_types: list[int] = [],     # per-run matrix constants
    format: str = "fasta",
    seq_type: int = 1,             # fallback if run_types empty
    seed: int = 42,                # base seed (run k → seed+k)
    min_support: int = 0,          # POAR threshold
    refine: int = 0,               # KALIGN_REFINE_* constant
    vsm_amax: float = -1.0,       # -1.0 = C default
    realign: int = 0,
    seq_weights: float = -1.0,     # -1.0 = C default
    n_threads: int = 1,
    consistency_anchors: int = 0,
    consistency_weight: float = 2.0,
)
```

---

## 5. Optimizer Search Space → `kalign_run_config` Mapping

```
Per-run parameter     | optimizer field      | run_config field
----------------------|----------------------|------------------
gap open              | run_gpo[k]           | .gpo
gap extend            | run_gpe[k]           | .gpe
terminal gap extend   | run_tgpe[k]          | .tgpe
tree noise            | run_noise[k]         | .tree_noise
matrix type           | run_types[k]         | .matrix
tree seed             | seed + k             | .tree_seed
VSM amplitude         | vsm_amax             | .vsm_amax
seq weights           | seq_weights          | .seq_weights
realign iters         | realign              | .realign
refinement mode       | refine               | .refine
consistency anchors   | consistency_anchors  | .consistency_anchors
consistency weight    | consistency_weight   | .consistency_weight
adaptive budget       | (not searched)       | .adaptive_budget = 0
dist scale            | (not searched)       | .dist_scale = 0.0
```

Ensemble-level:
```
min_support           | min_support          | ens.min_support
```

---

## 6. Refinement Constants

```python
from kalign._core import REFINE_NONE, REFINE_ALL, REFINE_CONFIDENT, REFINE_INLINE
```

| Constant           | Value | Description                            |
|-------------------|-------|----------------------------------------|
| `REFINE_NONE`      | 0     | No refinement                          |
| `REFINE_ALL`       | 1     | Refine all columns                     |
| `REFINE_CONFIDENT` | 2     | Refine only high-confidence columns    |
| `REFINE_INLINE`    | 3     | Inline refinement during alignment     |

---

## 7. Matrix Default Penalties

| Matrix        | gpo    | gpe   | tgpe   | Score range   |
|---------------|--------|-------|--------|---------------|
| PFASUM43      | 7.0    | 1.25  | 1.0    | -6 to 13      |
| PFASUM60      | 7.0    | 1.25  | 1.0    | -6 to 14      |
| CorBLOSUM66   | 5.5    | 2.0   | 1.0    | -4 to 13      |
| DNA           | 8.0    | 6.0   | 0.0    | -4 to 5       |
| DNA_INTERNAL  | 8.0    | 6.0   | 8.0    | -4 to 5       |
| RNA           | 217.0  | 39.4  | 292.6  | ~160 to 383   |

All three protein matrices are in 1/3-bit units with directly comparable
gap penalty ranges (optimizer searches `[2.0, 15.0]` for gpo across all three).

---

## 8. How Presets Get Into the C Library

Optimized presets are hardcoded in `lib/src/aln_wrap.c` in the functions
`preset_protein()`, `preset_dna()`, and `preset_rna()`.

Each function handles all three modes (fast/default/accurate) and fills
the caller-provided `runs[]` array and `*n_runs` count.

To update presets after re-optimization:
1. Run the optimizer to get the Pareto front
2. Select the fast/default/accurate configurations
3. Export to JSON:
   ```json
   {
     "protein": {
       "fast": {
         "n_runs": 1,
         "min_support": 0,
         "runs": [{
           "matrix": "pfasum60",
           "gpo": 8.4087, "gpe": 0.5153, "tgpe": 0.4927,
           "vsm_amax": 1.448, "seq_weights": 1.063,
           "dist_scale": 0.0, "realign": 0, "refine": 0,
           "adaptive_budget": 0, "tree_seed": 42, "tree_noise": 0.1623,
           "consistency_anchors": 0, "consistency_weight": 2.0
         }]
       }
     }
   }
   ```
4. Translate JSON → C `preset_run()` calls in `aln_wrap.c`
5. Rebuild and test
