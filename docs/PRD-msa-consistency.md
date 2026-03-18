# PRD: MSA-Derived Consistency Merge for Ensemble Alignment

## Goal

Add an alternative ensemble merge strategy that extracts residue-residue
consistency scores from N completed ensemble MSAs and uses them as bonus
weights in a final progressive alignment, replacing the POAR consensus path.

## Background

The current ensemble pipeline:
1. Runs N progressive alignments with diverse parameters
2. Builds a POAR table, scores each run, picks the best
3. Either selects the best single run or builds a POAR consensus MSA

This works well for F1 (0.810) but TC lags (0.523 vs 0.58+ for MUSCLE/MAFFT).
TC measures column-level consistency — exactly what a consistency-transformed
re-alignment should improve.

The idea: instead of POAR consensus, extract pairwise residue correspondences
from the N ensemble MSAs (which residues land in the same column across
multiple runs?) and use those agreement counts as bonus scores during a fresh
progressive alignment. The best-scoring run's gap penalties and matrix are
reused for the final alignment — no extra parameters to optimize.

## Design

### Reuse `consistency_table` + `sparse_bonus`

The existing anchor consistency system provides:
- `consistency_table` — stores position maps indexed by [seq × K + slot]
- `sparse_bonus` — sparse bonus matrix queried during DP
- `anchor_consistency_get_bonus()` / `_get_bonus_profile()` — computes
  bonuses for seq-seq, seq-profile, and profile-profile merges
- DP integration — `sparse_bonus_lookup()` already wired into match scores

For MSA consistency, K = n_runs, and position maps come from MSA column
structure instead of pairwise alignments. The "anchor position" is replaced
by "column index" — the `get_bonus` functions work unchanged.

### Data flow (when `consistency_merge = 1`)

```
1. Run N alignments with diverse params        [unchanged]
2. Extract POARs from each alignment            [unchanged]
3. Score alignments, find best_k                [unchanged]
4. NEW: Build consistency_table from N aligned MSAs
5. NEW: Copy original MSA, attach consistency_table
6. NEW: Run progressive alignment with best_k's params + consistency bonus
7. Copy result back to original MSA
8. Compute per-residue confidence from POAR     [unchanged]
```

### Column map extraction

For each ensemble MSA k, for each sequence i, build:
```
col_map[pos] = column index in aligned MSA k
```
where `pos` is the ungapped residue position. Stored in
`consistency_table.pos_maps[seq_i * K + run_k]`.

When queried for pair (seq_a, seq_b): if both map to the same column in
run k, that's a vote. Bonus = weight × (votes / K).

## API Changes

### C: `kalign_ensemble_config`
- `int consistency_merge;` — 0 = POAR (default), 1 = MSA consistency
- `float consistency_merge_weight;` — bonus weight (default 2.0)

### Python: `ensemble_custom_file_to_file()`
- `consistency_merge` (int, default 0)
- `consistency_merge_weight` (float, default 2.0)

### Optimizer: `optimize_unified.py`
- `consistency_merge`: Choice({0, 1}), only when n_runs > 1
- `consistency_merge_weight`: Real([0.5, 10.0]), only when consistency_merge=1

## Weight Parameter

A weight parameter is needed. The bonus is added to the substitution score
in the DP. Too high → rigidly reproduces ensemble consensus. Too low → no
benefit. The optimizer finds the sweet spot. Default 2.0.

## Implementation

### New files
- `lib/src/msa_consistency.c` — `msa_consistency_build()`
- `lib/src/msa_consistency.h` — header

### Modified files
- `lib/include/kalign/kalign_config.h` — ensemble_config fields
- `lib/src/aln_wrap.c` — defaults function
- `lib/src/ensemble.c` — consistency merge path
- `lib/CMakeLists.txt` — new source file
- `python-kalign/_core.cpp` — expose params
- `benchmarks/optimize_unified.py` — optimizer variables
