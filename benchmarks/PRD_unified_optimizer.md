# PRD: Unified Parameter Optimizer

## Goal

Build `benchmarks/optimize_unified.py` — a single NSGA-II optimizer that searches across kalign's entire operating range: from fast single-run alignment (no consistency, no ensemble) through consistency-enhanced single-run to multi-run ensemble with POAR consensus. The optimizer always uses three objectives: F1, TC, and wall time. The resulting 3D Pareto surface reveals the full speed/accuracy trade-off landscape in one run.

## Motivation

Previous optimization was split into separate scripts for single-run (`optimize_params.py`) and ensemble (`optimize_ensemble.py`), each producing isolated Pareto fronts that can't be directly compared. A unified optimizer solves this:

1. **Direct comparability**: All configurations live on the same Pareto surface.
2. **Discovery of hybrid regimes**: The optimizer might find unexpected sweet spots.
3. **One run, one answer**: Instead of running 4+ separate optimizations and manually comparing results, one run produces the complete picture.
4. **Mandatory wall time**: Time is always the 3rd objective, so the Pareto front naturally stratifies from fast/rough to slow/accurate.

## Complete Parameter Inventory

This section maps every lever that exists in kalign's C code, regardless of whether it's currently exposed through the Python API. The goal is to understand the full landscape before deciding what to optimize.

### C-level parameters (from `aln_param` struct + function signatures)

| Parameter | C field/arg | Type | Default (protein) | Where set | Description |
|-----------|------------|------|-------------------|-----------|-------------|
| `gpo` | `ap->gpo` | float | 7.0 (P43/P60), 5.5 (CB66), 55 (GON250) | `aln_param_init` | Gap open penalty |
| `gpe` | `ap->gpe` | float | 1.25 (P43/P60), 2.0 (CB66), 8 (GON250) | `aln_param_init` | Gap extend penalty |
| `tgpe` | `ap->tgpe` | float | 1.0 (P43/P60/CB66), 4 (GON250) | `aln_param_init` | Terminal gap extend penalty |
| `subm` | `ap->subm[23][23]` | float[][] | Matrix-dependent | `aln_param_init` | Substitution matrix (selected by `type`) |
| `type` | arg to `aln_param_init` | int | `PROTEIN_PFASUM43` | caller | Which substitution matrix: PFASUM43, PFASUM60, CorBLOSUM66, GON250 |
| `vsm_amax` | `ap->vsm_amax` | float | 2.0 protein, 0.0 DNA/RNA | `aln_param_init` | Variable scoring matrix: subtracts `max(0, amax-d)` from subst scores for close pairs |
| `use_seq_weights` | `ap->use_seq_weights` | float | 0.0 | `aln_param_init` | Profile rebalancing pseudocount (0=off) |
| `dist_scale` | `ap->dist_scale` | float | 0.0 | `aln_param_init` | Distance-dependent gap scaling (0=off) |
| `consistency_anchors` | `ap->consistency_anchors` | int | 0 | `kalign_run_seeded` | Number of anchor sequences K for consistency transform (0=off) |
| `consistency_weight` | `ap->consistency_weight` | float | 2.0 | `kalign_run_seeded` | Bonus scale for consistency anchors |
| `adaptive_budget` | `ap->adaptive_budget` | int | 0 | caller | Scale refinement trial count by uncertainty (0=off, 1=on) |
| `subm_offset` | `ap->subm_offset` | float | 0.0 | computed | VSM offset, computed per alignment step (not user-settable) |
| `refine` | arg | int | NONE=0 | caller | Post-alignment refinement: NONE(0), ALL(1), CONFIDENT(2), INLINE(3) |
| `realign` | arg to `kalign_run_realign` | int | 0 | caller | Alignment-guided UPGMA tree rebuild iterations (0=off) |
| `tree_seed` | arg to `kalign_run_seeded` | uint64 | 0 | caller | Random seed for noisy guide tree (0=deterministic) |
| `tree_noise` | arg to `kalign_run_seeded` | float | 0.0 | caller | Tree perturbation noise sigma (0=none) |
| `n_runs` | arg to `kalign_ensemble*` | int | 1 | caller | Number of ensemble runs (1=single-run) |
| `min_support` | arg to ensemble | int | 0 | caller | POAR consensus threshold (0=auto selection-vs-consensus) |

### What flows where

```
kalign_run_seeded():
    gpo, gpe, tgpe, type (→ subst matrix)
    vsm_amax, use_seq_weights, dist_scale
    consistency_anchors, consistency_weight
    refine, adaptive_budget
    tree_seed, tree_noise

kalign_run_realign():
    same as kalign_run_seeded EXCEPT:
    - tree_seed/tree_noise NOT supported (uses alignment-derived tree)
    - adds realign_iterations
    - consistency is rebuilt each iteration

kalign_ensemble_custom():
    PER-RUN: gpo[], gpe[], tgpe[], type[], noise[]
    SHARED:  vsm_amax, use_seq_weights, consistency_anchors,
             consistency_weight, realign, refine, min_support
    NOTE: use_seq_weights IS passed through to each per-run alignment.
          The -1.0 sentinel defaults to 0.0, but explicit positive values
          work fine. Prior finding that seq_weights "hurts ensemble" was
          with hardcoded scale-factor table — worth re-exploring.
```

### Parameters NOT currently in the ensemble_custom API

| Parameter | Status | Worth adding? |
|-----------|--------|--------------|
| `dist_scale` | Hardcoded to 0.0 in ensemble_custom | Low priority — VSM serves similar purpose |
| `adaptive_budget` | Hardcoded to 0 | Low priority — minor effect |
| `refine` per-run | Currently shared across all runs | Possible but complex — would need per-run refine[] array |

### Available substitution matrices

| `type` constant | Name | Default gaps | Origin | Notes |
|----------------|------|-------------|--------|-------|
| `KALIGN_TYPE_PROTEIN_PFASUM43` | PFASUM43 | gpo=7.0 gpe=1.25 tgpe=1.0 | Keul et al. 2017, 43% clustering | Current default |
| `KALIGN_TYPE_PROTEIN_PFASUM60` | PFASUM60 | gpo=7.0 gpe=1.25 tgpe=1.0 | Keul et al. 2017, 60% clustering | Optimization found this is best |
| `KALIGN_TYPE_PROTEIN` | CorBLOSUM66 | gpo=5.5 gpe=2.0 tgpe=1.0 | BLOSUM66 variant | Higher entropy |
| `KALIGN_TYPE_PROTEIN_DIVERGENT` | GON250 | gpo=55 gpe=8 tgpe=4 | Gonnet 1992 | Very different scale, for divergent seqs |

## Decisions for the Unified Optimizer

### What to optimize

| Parameter | Include? | Rationale |
|-----------|----------|-----------|
| `n_runs` | YES {1, 3, 5} | Core mode variable |
| `gpo_i` per-run | YES [2, 15] | Per-run diversity is key for ensemble |
| `gpe_i` per-run | YES [0.5, 5] | Per-run diversity |
| `tgpe_i` per-run | YES [0.1, 3] | Per-run diversity |
| `noise_i` per-run | YES [0, 0.5] | Tree perturbation per-run |
| `matrix_i` per-run | YES {P43, P60, CB66} | Matrix diversity per-run |
| `vsm_amax` shared | YES [0, 5] | Major effect on quality |
| `seq_weights` shared | YES [0, 5] | Works for single-run; *also* passed to per-run alignments in ensemble — re-explore |
| `consistency` shared | YES {0, 1, 2, 3, 5, 8, 10} | Huge effect, expensive |
| `consistency_weight` shared | YES [0.5, 5] | Tunes consistency strength |
| `realign` shared | YES {0, 1, 2} | Alignment-guided tree rebuild |
| `refine` shared | YES {NONE, CONFIDENT, ALL, INLINE} | Post-alignment refinement mode |
| `min_support` shared | YES {0, 1, ..., max_runs} | POAR consensus threshold |
| `dist_scale` | NO | Largely redundant with VSM, 0 extra dims |
| `adaptive_budget` | NO | Minor effect, adds noise to search |
| GON250 matrix | NO | Very different scale, would need separate gap ranges |

### Refine as a searchable parameter

Currently refine is hardcoded:
- Single-run: not used (NONE)
- Ensemble: REFINE_CONFIDENT for post-selection refinement

But refine={NONE, CONFIDENT, ALL, INLINE} could be optimized. INLINE is particularly interesting — it does refinement *during* progressive alignment rather than as a post-step. This has never been benchmarked in combination with consistency or ensemble.

Encode as discrete: {0=NONE, 1=ALL, 2=CONFIDENT, 3=INLINE}

### Masking rules (which params are inactive when)

The masking logic forces inactive parameters to neutral values during decode, so mutations in those dimensions are silent:

| Condition | Masked parameters | Forced value | Why |
|-----------|-------------------|-------------|-----|
| `n_runs == 1` | `noise_0` | 0.0 | No tree perturbation in single-run path (uses deterministic tree) |
| `n_runs == 1` | `min_support` | 0 | No POAR consensus |
| `n_runs == 1` | run_1..N params | Copy of run_0 | Dead dimensions |
| `n_runs > 1` | — | — | **All shared params remain active**, including `seq_weights` |
| `consistency == 0` | `consistency_weight` | 1.0 | No anchors → weight irrelevant |
| `realign > 0` | `noise_i` | 0.0 | `kalign_run_realign` doesn't use tree_seed/noise (uses alignment-derived tree) |

**Important**: `seq_weights` is NOT masked for ensemble mode. The C code passes `use_seq_weights` through to each per-run alignment via `kalign_run_seeded`/`kalign_run_realign`. The previous finding that "seq_weights hurts ensemble" was with the old hardcoded parameters. With optimized per-run params, seq_weights might help. The optimizer should be free to explore this.

**Important**: When `realign > 0`, the ensemble path calls `kalign_run_realign` which builds a deterministic alignment-derived tree. Tree noise has no effect in this path. So if the optimizer picks `realign > 0`, per-run noise values should be masked to 0. However, this creates a complex interaction: `realign=0` enables noise diversity, `realign>0` gives better trees but loses noise diversity. The optimizer should discover which trade-off wins.

### C API changes needed

The current `ensemble_custom_file_to_file` Python binding needs one addition:

1. **`refine` parameter**: Currently hardcoded to `REFINE_CONFIDENT` in the optimize_ensemble.py evaluation. The binding already accepts it as a parameter. We just need to make it a decision variable instead of hardcoding.

No C code changes needed. All the levers already exist.

## Parameter Space (final)

### Per-run parameters (allocated for max_runs slots)

| Parameter | Range | Type | Dims per run |
|-----------|-------|------|-------------|
| `gpo_i` | [2.0, 15.0] | continuous | 1 |
| `gpe_i` | [0.5, 5.0] | continuous | 1 |
| `tgpe_i` | [0.1, 3.0] | continuous | 1 |
| `noise_i` | [0.0, 0.5] | continuous | 1 |
| `matrix_i` | {0, 1, 2} | discrete | 1 |

= 5 per run × max_runs

### Shared parameters

| Parameter | Range | Type | Dims |
|-----------|-------|------|------|
| `vsm_amax` | [0.0, 5.0] | continuous | 1 |
| `seq_weights` | [0.0, 5.0] | continuous | 1 |
| `consistency_weight` | [0.5, 5.0] | continuous | 1 |
| `n_runs` | {0, 1, 2} for max5, {0, 1, 2, 3} for max8 | discrete | 1 |
| `consistency` | {0..6} → [0, 1, 2, 3, 5, 8, 10] | discrete | 1 |
| `realign` | {0, 1, 2} | discrete | 1 |
| `refine` | {0, 1, 2, 3} → [NONE, ALL, CONFIDENT, INLINE] | discrete | 1 |
| `min_support` | {0..max_runs} | discrete | 1 |

= 8 shared (3 continuous + 5 discrete)

### Total dimensions

| max_runs | Per-run | Shared | Total | Recommended pop_size |
|----------|---------|--------|-------|---------------------|
| 5 | 25 | 8 | 33 | 200 |
| 8 | 40 | 8 | 48 | 300 |

## Evaluation Pipeline

### Decision routing

```python
def evaluate_unified(params, cases, n_threads=1, quiet=True):
    n_runs = params["n_runs"]

    if n_runs == 1:
        # Single-run path — full parameter control
        kalign.align_file_to_file(
            input, output,
            gap_open=params["gpo"][0],
            gap_extend=params["gpe"][0],
            terminal_gap_extend=params["tgpe"][0],
            seq_type=matrix_api_name(params["matrix"][0]),
            vsm_amax=params["vsm_amax"],
            seq_weights=params["seq_weights"],
            consistency=params["consistency"],
            consistency_weight=params["consistency_weight"],
            realign=params["realign"],
            refine=refine_api_name(params["refine"]),
        )
    else:
        # Ensemble path — per-run arrays
        ensemble_custom_file_to_file(
            input, output,
            run_gpo=params["gpo"][:n_runs],
            run_gpe=params["gpe"][:n_runs],
            run_tgpe=params["tgpe"][:n_runs],
            run_noise=params["noise"][:n_runs],
            run_types=params["matrix"][:n_runs],
            vsm_amax=params["vsm_amax"],
            seq_weights=params["seq_weights"],  # NOT forced to 0!
            realign=params["realign"],
            consistency_anchors=params["consistency"],
            consistency_weight=params["consistency_weight"],
            refine=params["refine"],
            min_support=params["min_support"],
            seed=42,
        )
```

### CV evaluation

Same as existing: stratified k-fold (default 5) on BAliBASE 218 cases.

## Objectives (always 3)

1. **Maximize F1** — category-averaged, held-out CV folds → pymoo minimizes -F1
2. **Maximize TC** — category-averaged, held-out CV folds → pymoo minimizes -TC
3. **Minimize wall time** — total CV evaluation time in seconds

No `--optimize-time` flag. Wall time is always the 3rd objective.

## Baselines

Computed at startup (parallelized across folds + workers):

| Baseline | Description | Expected profile |
|----------|-------------|-----------------|
| **fast** | n_runs=1, consistency=0, realign=0, refine=NONE, PFASUM60, default gaps | F1~0.72, TC~0.47, ~10s |
| **accurate** | n_runs=1, consistency=8, realign=2, refine=NONE, optimized gaps from run 1 | F1~0.76, TC~0.47, ~180s |
| **ensemble** | n_runs=3, consistency=0, realign=1, refine=CONFIDENT, vsm=2.0 | F1~0.77, TC~0.47, ~193s |

## Parallelism

### Fine-grained (individual x fold) job distribution

Same as updated `optimize_params.py`:

- Each job = one (individual, fold) pair, running ~44 cases
- For pop=200, k=5 → 1000 jobs per generation
- Distributed across N workers

### Baseline parallelization

All baseline evaluations (3 baselines × (k folds + 1 full)) run in parallel at startup.

### Top-level function

```python
def _eval_one_unified_fold(args_tuple):
    ind_idx, fold_idx, x, test_cases, n_threads, max_runs = args_tuple
    params = decode_unified_params(x, max_runs)
    result = evaluate_unified(params, test_cases, n_threads, quiet=True)
    return ind_idx, fold_idx, params, result
```

## Dashboard

### Layout

```
┌─ Progress ──────────────────────────────────────────────────────┐
│ Gen 12/100  Eval 480/20000  Elapsed 2.1h  ETA 5.3h             │
│ Gen time 128s  Workers 56                                       │
├─ Baselines ─────────────────────────────────────────────────────┤
│ fast:      F1=0.7160  TC=0.4660  Time=10s    (single, c=0)     │
│ accurate:  F1=0.7557  TC=0.4713  Time=180s   (single, c=8)     │
│ ensemble:  F1=0.7680  TC=0.4670  Time=193s   (ens3, c=0)       │
├─ Best by Objective ─────────────────────────────────────────────┤
│ Best F1:   0.7823 (+0.0143 vs accurate)  ens3 c=5 re=1  210s   │
│ Best TC:   0.5102 (+0.0389 vs accurate)  ens5 c=3 re=2  340s   │
│ Fastest:   8.2s   (-1.8s vs fast)        single c=0 re=0       │
├─ Pareto Front (top 12 by F1) ──────────────────────────────────┤
│ # │ Mode  │ CV F1  │ CV TC  │ Time  │ c  re ref│ Key params    │
│ 0 │ ens5  │ 0.7823 │ 0.5001 │ 340s  │ 5  1  C │ vsm=1.4      │
│ 1 │ ens3  │ 0.7791 │ 0.4891 │ 210s  │ 3  1  C │ vsm=1.8      │
│ 2 │ ens3  │ 0.7756 │ 0.4812 │ 155s  │ 0  1  C │ vsm=2.0      │
│ 3 │ single│ 0.7557 │ 0.4713 │ 180s  │ 8  2  N │ vsm=1.4 P60  │
│ 4 │ single│ 0.7401 │ 0.4650 │  42s  │ 3  1  N │ vsm=2.0 P60  │
│ 5 │ single│ 0.7160 │ 0.4580 │   9s  │ 0  0  N │ vsm=1.8 P60  │
├─ Mode Distribution ────────────────────────────────────────────┤
│ Pareto: single=8  ens3=12  ens5=5                               │
│ Pop:    single=65  ens3=82  ens5=53                             │
├─ Trend ────────────────────────────────────────────────────────┤
│ Gen   1: F1=0.7012  Gen  5: F1=0.7434  Gen 10: F1=0.7689      │
├─ Recent ───────────────────────────────────────────────────────┤
│ F1=0.7654 TC=0.4512 t=145s ens3 c=3 re=1 ref=C vsm=1.8       │
│ F1=0.7201 TC=0.4380 t=12s  single c=0 re=0 ref=N vsm=2.1     │
└─────────────────────────────────────────────────────────────────┘
```

### Key dashboard features

1. **Mode column**: "single", "ens3", "ens5", "ens8" — most important column.
2. **Refine column**: N=None, A=All, C=Confident, I=Inline — compact 1-char code.
3. **Mode distribution panel**: Pareto and population counts per mode.
4. **Three baselines**: Context for all three speed regimes.
5. **Time always shown**: Every entry shows wall time.
6. **Compact params**: For ensemble, show shared params only (c, re, ref, vsm, sw). Full per-run details go to output file.

## CLI Interface

```bash
# Quick smoke test
uv run python -m benchmarks.optimize_unified --pop-size 20 --n-gen 5

# Production run on Threadripper (max 5 ensemble runs)
uv run python -m benchmarks.optimize_unified \
    --pop-size 200 --n-gen 100 \
    --n-workers 56 --n-threads 1

# Larger search (max 8 ensemble runs)
uv run python -m benchmarks.optimize_unified \
    --max-runs 8 --pop-size 300 --n-gen 120 \
    --n-workers 56

# Resume
uv run python -m benchmarks.optimize_unified \
    --resume benchmarks/results/unified_optim/gen_checkpoint.pkl \
    --n-gen 150 --n-workers 56
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--max-runs` | 5 | Max ensemble runs. Sets n_runs choices: {1,3,5} or {1,3,5,8} |
| `--pop-size` | 200 | NSGA-II population size |
| `--n-gen` | 100 | Total generations |
| `--n-folds` | 5 | Stratified CV folds |
| `--n-threads` | 1 | OpenMP threads per alignment |
| `--n-workers` | 1 | Parallel worker processes |
| `--seed` | 42 | Random seed |
| `--output-dir` | `benchmarks/results/unified_optim` | Output directory |
| `--run-name` | None | Subdirectory name |
| `--no-dashboard` | false | Disable rich dashboard |
| `--resume` | None | Checkpoint file path |

## Checkpoint / Resume

Same pattern as existing optimizers:
- Saves `pop_X`, `pop_F`, `n_gen_completed`, history, plus `max_runs` for validation
- Atomic write via temp file + rename
- On resume: validates `max_runs`, `n_folds`, `seed` match

## Clean Interrupt

- `_kill_pool()` sends SIGTERM to workers
- `os._exit(1)` to skip atexit join-hangs

## Output Files

1. **`gen_checkpoint.pkl`** — per-generation checkpoint
2. **`pareto_front.txt`** — human-readable with full per-run params for ensemble solutions
3. **`optim_results.pkl`** — full results pickle

### pareto_front.txt format

```
# Unified kalign optimization (NSGA-II, 3 objectives: F1, TC, time)
# pop_size=200 n_gen=100 max_runs=5 n_folds=5 seed=42
# Baselines:
#   fast:     F1=0.7160 TC=0.4660 Time=10s
#   accurate: F1=0.7557 TC=0.4713 Time=180s
#   ensemble: F1=0.7680 TC=0.4670 Time=193s

[0] mode=ens5 CV_F1=0.7823 CV_TC=0.5001 Time=340s
    n_runs=5
    vsm_amax=1.359 seq_weights=0.8
    consistency=5 consistency_weight=1.17
    realign=1 refine=CONFIDENT
    min_support=2
    run_0: gpo=7.0 gpe=0.55 tgpe=0.41 noise=0.10 PFASUM60
    run_1: gpo=3.5 gpe=1.20 tgpe=0.80 noise=0.25 PFASUM43
    ...

[1] mode=single CV_F1=0.7557 CV_TC=0.4713 Time=180s
    n_runs=1
    gap_open=8.472 gap_extend=0.554 terminal_gap_extend=0.409
    vsm_amax=1.359 seq_weights=3.407
    consistency=8 consistency_weight=1.167
    realign=2 refine=NONE
    matrix=PFASUM60
```

## Post-Optimization Analysis

1. Full Pareto front as rich table
2. Group by mode, show best F1/TC/time per mode
3. Re-evaluate top-3 on full dataset
4. Per-category breakdown (RV11-RV50)
5. Overfit check (full > CV by >0.02)
6. Print recommended settings for three tiers:
   - **Fast** (best F1 among solutions under 15s)
   - **Default** (best F1 among solutions under 60s)
   - **Accurate** (best F1 overall)

## Implementation Order

1. Parameter space definition + encode/decode with masking
2. `evaluate_unified()` with single-run / ensemble routing
3. `_eval_one_unified_fold()` for parallel eval
4. `UnifiedCVProblem` with 3 mandatory objectives
5. Dashboard with mode columns, multiple baselines, mode distribution
6. `GenerationCallback` with checkpoint saving
7. `main()` with CLI, parallel baselines, optimization loop, results
8. Resume logic with max_runs validation
9. Post-optimization analysis
10. Smoke test: `--pop-size 10 --n-gen 3 --max-runs 5 --n-workers 4`
