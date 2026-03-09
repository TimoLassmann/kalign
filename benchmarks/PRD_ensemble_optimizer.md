# PRD: Ensemble Parameter Optimizer

## Goal

Build `benchmarks/optimize_ensemble.py` — an NSGA-II optimizer that finds the best per-run parameters for kalign's ensemble alignment mode. This is the companion to `optimize_params.py` (single-run optimizer) but targets the ensemble pipeline where multiple diverse alignments are combined via POAR consensus.

## Background

The ensemble pipeline runs N independent alignments with different settings, then combines them using a Partial Order Alignment Representation (POAR) consensus. Currently, the per-run diversity is controlled by a hardcoded scale-factor table in `ensemble.c`. We've added `kalign_ensemble_custom()` which accepts fully independent per-run parameters, enabling proper optimization.

Key insight from prior experiments: the consistency transform (anchor-based alignment improvement) has never been tested inside ensemble runs. It should boost the quality of each individual input alignment without reducing diversity, since diversity comes from gap penalties/matrices/noise, not from anchoring.

## Search Space

### Per-run parameters (× N runs)

Each of the N runs gets its own:

| Parameter | Range | Type | Description |
|-----------|-------|------|-------------|
| `gpo` | [2.0, 15.0] | continuous | Gap open penalty |
| `gpe` | [0.5, 5.0] | continuous | Gap extend penalty |
| `tgpe` | [0.1, 3.0] | continuous | Terminal gap extend penalty |
| `matrix` | {PFASUM43, PFASUM60, CorBLOSUM66} | discrete | Substitution matrix |
| `noise` | [0.0, 0.5] | continuous | Tree perturbation noise sigma |

= 5 parameters per run

### Shared parameters (apply to all runs)

| Parameter | Range | Type | Description |
|-----------|-------|------|-------------|
| `vsm_amax` | [0.0, 5.0] | continuous | Variable scoring matrix strength |
| `consistency` | {0, 1, 2, 3, 5, 8} | discrete | Anchor consistency rounds per run |
| `consistency_weight` | [0.5, 5.0] | continuous | Consistency bonus weight |
| `realign` | {0, 1, 2} | discrete | Tree-rebuild iterations per run |
| `min_support` | {0, 1, 2, ..., N} | discrete | POAR consensus threshold (0 = auto) |

= 5 shared parameters

### Total parameter counts by N runs

| N runs | Per-run | Shared | Total | Recommended pop_size |
|--------|---------|--------|-------|---------------------|
| 3 | 15 | 5 | 20 | 100 |
| 5 | 25 | 5 | 30 | 150 |
| 8 | 40 | 5 | 45 | 200 |

## Architecture

### Script: `benchmarks/optimize_ensemble.py`

Follows the same architecture as `optimize_params.py`:

```
optimize_ensemble.py
├── Parameter space definition (encode/decode per-run arrays)
├── evaluate_ensemble_params()     — run one ensemble config on a case list
├── evaluate_ensemble_cv()         — stratified k-fold CV wrapper
├── _eval_one_ensemble()           — top-level function for ProcessPoolExecutor
├── Dashboard (rich live)          — same look & feel as optimize_params.py
├── EnsembleCVProblem (pymoo)      — NSGA-II problem with serial/parallel eval
├── GenerationCallback             — dashboard updates + checkpoint saving
├── load_checkpoint()              — resume from pickle
└── main()                         — CLI entry point
```

### Key differences from `optimize_params.py`

1. **Parameter encoding**: The decision vector has a variable-length per-run section. For N=3: `[gpo_0, gpe_0, tgpe_0, noise_0, ..., gpo_2, gpe_2, tgpe_2, noise_2, vsm_amax, consistency_weight, matrix_0, matrix_1, matrix_2, consistency, realign, min_support]`. The `decode_params()` function returns a dict with lists for per-run params.

2. **Alignment call**: Uses `kalign._core.ensemble_custom_file_to_file()` instead of `kalign.align_file_to_file()`.

3. **Fixed N runs**: The `--n-runs` CLI argument sets N (default: 3). Separate optimization runs for different N values, then compare Pareto fronts.

4. **Display**: The dashboard `format_params_short()` shows per-run params compactly, e.g. `R0: gpo=7.0/PFASUM43 R1: gpo=3.5/PFASUM60 R2: gpo=10.5/CorBLOSUM66 | vsm=2.0 c=3 re=1 ms=2`.

### CLI interface

```bash
# Quick test
uv run python -m benchmarks.optimize_ensemble --n-runs 3 --pop-size 20 --n-gen 5

# Production run on Threadripper (3 runs)
uv run python -m benchmarks.optimize_ensemble \
    --n-runs 3 --pop-size 100 --n-gen 50 \
    --n-workers 56 --n-threads 1 --n-folds 5

# Production run (5 runs)
uv run python -m benchmarks.optimize_ensemble \
    --n-runs 5 --pop-size 150 --n-gen 60 \
    --n-workers 56 --n-threads 1

# Production run (8 runs)
uv run python -m benchmarks.optimize_ensemble \
    --n-runs 8 --pop-size 200 --n-gen 80 \
    --n-workers 56 --n-threads 1

# Resume after interrupt
uv run python -m benchmarks.optimize_ensemble \
    --resume benchmarks/results/ensemble_optim/gen_checkpoint.pkl \
    --n-gen 80 --n-workers 56

# Add wall time as 3rd objective
uv run python -m benchmarks.optimize_ensemble --n-runs 3 --optimize-time
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--n-runs` | 3 | Number of ensemble runs (fixed per optimization) |
| `--pop-size` | 100 | NSGA-II population size |
| `--n-gen` | 50 | Total number of generations |
| `--n-folds` | 5 | Stratified CV folds |
| `--n-threads` | 1 | OpenMP threads per alignment |
| `--n-workers` | 1 | Parallel worker processes |
| `--seed` | 42 | Random seed |
| `--optimize-time` | false | Add wall time as 3rd objective |
| `--output-dir` | `benchmarks/results/ensemble_optim` | Output directory |
| `--no-dashboard` | false | Disable rich dashboard |
| `--resume` | None | Path to checkpoint file for resume |

## Objectives

Same as single-run optimizer:
1. **Maximize F1** (category-averaged, held-out CV folds)
2. **Maximize TC** (category-averaged, held-out CV folds)
3. *(Optional)* **Minimize wall time**

## Evaluation pipeline

For each individual in the population:

```
decode_params(x) → per-run arrays + shared params
    ↓
for each CV fold:
    for each test case:
        kalign._core.ensemble_custom_file_to_file(
            input, output,
            run_gpo=[gpo_0, ..., gpo_N],
            run_gpe=[gpe_0, ..., gpe_N],
            run_tgpe=[tgpe_0, ..., tgpe_N],
            run_noise=[noise_0, ..., noise_N],
            run_types=[matrix_0, ..., matrix_N],
            vsm_amax=vsm_amax,
            realign=realign,
            consistency_anchors=consistency,
            consistency_weight=consistency_weight,
            min_support=min_support,
            refine=REFINE_CONFIDENT,  # always on for ensemble
            seed=42,
        )
        score_alignment_detailed(reference, output)
    ↓
    category-averaged F1, TC for this fold
    ↓
mean across folds → CV F1, CV TC
```

## Baseline

The baseline for comparison is the current best ensemble result:
- `ens3+vsm+ref+ra1`: F1=0.768, TC=0.467
- Uses hardcoded scale-factor table, no consistency, PFASUM43 for all runs

This will be computed via `kalign.align_file_to_file(ensemble=3, vsm_amax=2.0, refine="confident", realign=1)` at startup.

## Dashboard

Same rich live dashboard as `optimize_params.py`:

```
┌─ Progress ──────────────────────────────────────────────────┐
│ Gen 12/50  Eval 480/5000  Elapsed 45.2m  ETA 142m          │
│ Current: R0:gpo=7.0 R1:gpo=3.5 R2:gpo=10.5 | vsm=2.0 c=3 │
├─ Best Solutions ────────────────────────────────────────────┤
│ Baseline:  F1=0.7680  TC=0.4670                             │
│                                                             │
│ Best F1:   0.7823 (+0.0143)                                 │
│   R0:7.0/P43 R1:3.5/P60 R2:10.5/CB66 | vsm=1.8 c=3 re=1  │
│                                                             │
│ Best TC:   0.4891 (+0.0221)                                 │
│   R0:6.5/P43 R1:4.0/P43 R2:9.0/P60 | vsm=2.1 c=5 re=1    │
├─ Pareto Front ──────────────────────────────────────────────┤
│ # │ CV F1  │ CV TC  │ Parameters                            │
│ 0 │ 0.7823 │ 0.4801 │ R0:7.0/P43 R1:3.5/P60 ...           │
│ 1 │ 0.7791 │ 0.4891 │ R0:6.5/P43 R1:4.0/P43 ...           │
│ ...                                                         │
├─ Trend ─────────────────────────────────────────────────────┤
│ Gen  1: F1=0.7512  Gen  5: F1=0.7634  Gen 10: F1=0.7789   │
├─ Recent ────────────────────────────────────────────────────┤
│ F1=0.7654 TC=0.4512 R0:8.1/P43 R1:5.2/CB66 R2:3.0/P60    │
│ F1=0.7823 TC=0.4801 R0:7.0/P43 R1:3.5/P60 R2:10.5/CB66   │
└─────────────────────────────────────────────────────────────┘
```

## Checkpoint / Resume

Identical to `optimize_params.py`:
- `GenerationCallback` saves `gen_checkpoint.pkl` after every completed generation
- Stores: population X and F arrays, generation count, evaluation history
- Atomic write (write to .tmp, rename)
- `--resume` loads checkpoint, reconstructs NSGA2 with saved population as initial sampling
- Remaining generations = `--n-gen` minus completed generations
- `--n-workers` and `--n-threads` can change between runs
- `--n-folds`, `--seed`, `--n-runs` must stay the same

## Clean interrupt

Identical to `optimize_params.py`:
- `_kill_pool()` sends SIGTERM to worker processes on KeyboardInterrupt
- `os._exit(1)` in main handler to skip atexit join-hangs
- Single Ctrl+C exits cleanly, prints checkpoint path and resume command

## Output files

Written to `--output-dir` (default: `benchmarks/results/ensemble_optim`):

1. **`gen_checkpoint.pkl`** — per-generation checkpoint for resume
2. **`pareto_front.txt`** — human-readable Pareto front with full per-run parameters
3. **`optim_checkpoint.pkl`** — final results pickle (Pareto configs, history, baselines)

## Post-optimization analysis

After optimization completes:
1. Print Pareto front as rich table
2. Re-evaluate best-F1 solution on full dataset (all 218 cases)
3. Show per-category breakdown (RV11–RV50) comparing optimized vs baseline
4. Overfit check: flag if full-dataset score exceeds CV score by >0.02

## Code reuse

The following can be imported from `optimize_params.py` or a shared module:
- `stratified_kfold()` — fold splitting logic
- `score_alignment_detailed()` — already in `scoring.py`
- `BenchmarkCase`, `balibase_cases`, etc. — already in `datasets.py`

The following must be new (ensemble-specific):
- `decode_ensemble_params()` / `encode_ensemble_params()`
- `evaluate_ensemble_params()` — calls `ensemble_custom_file_to_file`
- `evaluate_ensemble_cv()` — CV wrapper
- `format_ensemble_params_short()` / `format_ensemble_params_long()`
- `EnsembleCVProblem` — pymoo Problem subclass
- Dashboard and callback can be adapted from the existing ones

## Implementation order

1. Parameter space definition + encode/decode
2. `evaluate_ensemble_params()` + `evaluate_ensemble_cv()`
3. `EnsembleCVProblem` with serial and parallel evaluation
4. Dashboard (adapt from existing)
5. `GenerationCallback` with checkpoint saving
6. `main()` with CLI, baseline, optimization loop, results
7. Resume logic
8. Smoke test with `--pop-size 4 --n-gen 2 --n-runs 3`
