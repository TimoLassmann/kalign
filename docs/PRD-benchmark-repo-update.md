# PRD: Benchmark Repo Update for Kalign v3.5

## Context

The kalign Python bindings have been simplified to a clean two-path
architecture. The old parameter-based API (`ensemble=`, `refine=`,
`vsm_amax=`, `consistency=`, etc.) has been removed from the public
`kalign.align*()` functions. These functions now only accept `mode`
as the configuration mechanism.

Full technical details: `docs/parameter-cleanup-integration.md` in the
kalign repo.

**Important**: The benchmark repo's `CLAUDE.md` says "Do not edit
external code". All changes here are ONLY to the benchmark repo. The
kalign library itself is already updated and installed.

---

## Prerequisites

Install the updated kalign from source:

```bash
cd /Users/timo/code/kalign
uv pip install -e .
```

Verify:
```bash
uv run python -c "import kalign; print(kalign.__version__); print(kalign.align(['ATCG','ATCGG']))"
```

---

## Kalign API Quick Reference

### Standard alignment (mode-based)

```python
import kalign

# Three modes: "fast", "default", "accurate"
result = kalign.align_from_file("input.fasta", mode="default")
names, sequences = result  # AlignedSequences unpacks as 2-tuple

# File-to-file
kalign.align_file_to_file("input.fasta", "output.fasta", mode="default")
```

Parameters accepted by all three `align*()` functions:
- `mode`: `"fast"` | `"default"` | `"accurate"` (default: `"default"`)
- `seq_type`: `"auto"` | `"dna"` | `"rna"` | `"protein"` (default: `"auto"`)
- `gap_open`, `gap_extend`, `terminal_gap_extend`: optional float overrides
  (when any is set, mode is forced to `"fast"`)
- `n_threads`: int (default: 1)

No other parameters exist. There is no `ensemble`, `refine`, `vsm_amax`,
`realign`, `consistency`, `seq_weights`, etc. in the public API. These
are handled internally by the mode presets.

### Optimizer path (for NSGA-III evaluation)

```python
from kalign._core import ensemble_custom_file_to_file

ensemble_custom_file_to_file(
    input_file, output_file,
    run_gpo=[...], run_gpe=[...], run_tgpe=[...], run_noise=[...],
    run_types=[...],          # per-run KALIGN_MATRIX_* constants
    format="fasta",
    seq_type=1,               # fallback matrix if run_types empty
    seed=42,                  # base seed (run k gets seed+k)
    min_support=0,            # POAR consensus threshold
    refine=0,                 # KALIGN_REFINE_* constant
    vsm_amax=-1.0,            # -1.0 = C default
    realign=0,
    seq_weights=-1.0,         # -1.0 = C default
    n_threads=1,
    consistency_anchors=0,
    consistency_weight=2.0,
)
```

This is the ONLY way to pass fine-grained per-run parameters. It is used
exclusively by the optimizer, not by benchmark runners.

### Scoring

```python
import kalign

# SP score (0-100)
sp = kalign.compare("reference.fasta", "test.fasta")

# Detailed: recall, precision, F1, TC
d = kalign.compare_detailed("reference.fasta", "test.fasta", max_gap_frac=-1.0)
# d = {"recall": ..., "precision": ..., "f1": ..., "tc": ..., ...}

# With BAliBASE XML column mask
d = kalign.compare_detailed("ref.fasta", "test.fasta", column_mask=[1,0,1,...])
```

---

## Step 1: Update `src/runners.py` — `run_kalign()`

### 1a. Simplify `run_kalign()`

Remove all deprecated parameter kwargs. The function should only accept
`mode` and `seq_type`:

```python
def run_kalign(
    input_fasta: Path,
    mode: str = "default",
    seq_type: str = "auto",
) -> AlignResult:
    """Run kalign via the Python API using mode presets."""
    import kalign as _kalign

    start = time.perf_counter()
    result = _kalign.align_from_file(
        str(input_fasta),
        seq_type=seq_type,
        mode=mode,
    )
    wall = time.perf_counter() - start
    ...
```

### 1b. Update `METHODS` registry

```python
METHODS = {
    "kalign": {
        "fn": run_kalign,
        "mode": "default",
    },
    "kalign_fast": {
        "fn": run_kalign,
        "mode": "fast",
    },
    "kalign_accurate": {
        "fn": run_kalign,
        "mode": "accurate",
    },
    # External tools unchanged
    "mafft": {"fn": run_mafft},
    "muscle": {"fn": run_muscle},
    "clustalo": {"fn": run_clustalo},
}
```

Remove `"kalign_precise"` — use `"kalign_accurate"` instead. (If you
need backward compatibility with existing result files, keep `"precise"`
as an alias that maps to `mode="accurate"`.)

### 1c. Update `METHOD_COLORS` and `METHOD_ORDER`

Replace `"kalign_precise"` with `"kalign_accurate"` everywhere.

### 1d. Update `run_method()` dispatcher

After simplifying `run_kalign()`, the dispatcher should pass only
`mode` and `seq_type`.

**Smoke test**:
```bash
uv run python -c "
from src.runners import run_method
from pathlib import Path
r = run_method('kalign_fast', Path('data/downloads/balibase/bb3_release/RV11/BB11001.tfa'), Path('/tmp/test'))
print(f'wall_time={r.wall_time:.2f}s')
"
```

---

## Step 2: Update `config/config.yaml`

Replace `kalign_precise` with `kalign_accurate`:

```yaml
methods:
  kalign:
    mode: default
  kalign_fast:
    mode: fast
  kalign_accurate:
    mode: accurate
  mafft: {}
  muscle: {}
  clustalo: {}
```

Find-and-replace `kalign_precise` → `kalign_accurate` in:
- `config/config.yaml`
- Any Snakemake rules that reference method names

---

## Step 3: Update optimizer imports (in kalign repo)

**Note**: This is in `/Users/timo/code/kalign/benchmarks/`, not the
benchmark repo. Listed here for completeness.

### 3a. Fix `matrix_map_int`

```python
from kalign._core import MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66

"matrix_map_int": [MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66],
"matrix_map_str": ["pfasum43", "pfasum60", "corblosum66"],
```

### 3b. Keep `consistency` and `consistency_weight` in search space

These are real per-run parameters passed through
`ensemble_custom_file_to_file()` → `kalign_run_config` → C engine.

### 3c. Export JSON format for optimized presets

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

---

## Step 4: Run benchmark to verify

```bash
cd /Users/timo/work/Documents/Manuscripts/2026_kalign_35

# Quick test
podman run --rm -v $(pwd):/work -w /work kalign-35-bench \
    snakemake --cores 4 results/summaries/alignment_accuracy.json
```

Or on host:
```bash
uv run python workflow/scripts/run_balibase_quick.py
```

Expected results for protein (BAliBASE, 218 families):
- `kalign` (default mode): F1 ~ 0.72-0.73, TC ~ 0.47
- `kalign_fast`: F1 ~ 0.70-0.72
- `kalign_accurate`: F1 ~ 0.76-0.77

---

## Step 5: Update figure scripts (cosmetic)

Replace "precise" → "accurate" in labels and legends:
- `figures/fig_balibase.py`
- `figures/fig_bralibase.py`
- `figures/fig_speed.py`
- `figures/fig_summary_heatmap.py`
- `figures/style.py`

---

## Implementation Order

1. **Step 1** (runners.py) — critical, unblocks everything
2. **Step 2** (config.yaml) — quick, do with step 1
3. **Step 3** (optimizer) — independent, in kalign repo
4. **Step 4** (verify) — after steps 1-2
5. **Step 5** (figures) — cosmetic, can wait

Steps 1-2 are in the benchmark repo.
Step 3 is in the kalign repo.
Steps 4-5 are in the benchmark repo.

---

## Risk: Container Rebuild

The benchmark container has kalign baked in. After the update, rebuild:

```bash
cd /Users/timo/work/Documents/Manuscripts/2026_kalign_35
bash containers/build.sh
```

Or for host-side testing:
```bash
cd /Users/timo/code/kalign && uv pip install -e .
```
