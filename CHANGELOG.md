# Changelog

## [3.5.0] - Unreleased

### Headline: Three named modes

Kalign v3.5 introduces three named modes that package the best-performing
configurations into simple flags:

| Mode | CLI | Python | Description |
|------|-----|--------|-------------|
| **default** | `kalign` | `mode="default"` or omit | Consistency anchors + VSM. Best general-purpose. |
| **fast** | `kalign --fast` | `mode="fast"` | VSM only. Fastest, equivalent to kalign v3.4. |
| **precise** | `kalign --precise` | `mode="precise"` | Ensemble(3) + VSM + realign. Highest precision. |

Explicit parameters always override mode defaults.

**Breaking change:** The default behavior now includes consistency anchors
(K=5) and VSM. To get the old v3.4 behavior, use `--fast` or `mode="fast"`.

### New features

- **Three named modes** (`--fast`, `--precise`): Simple entry points inspired
  by minimap2's `-x preset` pattern. Most users only need to pick a mode.
  Expert parameters remain available and override mode defaults.
- **Ensemble alignment** (`--ensemble N`, or `--precise`): Runs N alignments
  with varied gap penalties and tree noise, combines results via POAR
  (Pairs of Aligned Residues) consensus. Improves F1 by ~5 points on
  BAliBASE at ~10x time cost.
- **POAR consensus save/load** (`--save-poar`, `--load-poar`): Save the
  consensus table after an ensemble run, then instantly re-threshold
  with different `--min-support` values without re-alignment.
- **Per-column and per-residue confidence scores** from ensemble mode.
  Accessible in Python via `result.column_confidence` and
  `result.residue_confidence`. Writable as Stockholm PP annotations.
- **Anchor consistency**: Computes K fast BPM anchor alignments and feeds
  consistency signals into the substitution matrix. Enabled by default
  (K=5) in the default mode.
- **Variable Scoring Matrix (VSM)**: Adapts substitution scores based on
  estimated pairwise distance. Enabled by default for protein; disabled
  for DNA/RNA.
- **Alignment-guided realignment**: Rebuilds the guide tree from pairwise
  identities in the aligned MSA, then re-aligns. Used by `--precise` mode.
- **Sequence weight rebalancing**: Pseudocount-based profile rebalancing
  to reduce bias from redundant sequences. Auto-disabled in ensemble
  mode where POAR consensus already handles imbalance.
- **Alignment refinement** (`--refine`): Post-alignment column-level
  refinement with modes "all", "confident", and "inline".
- **PFASUM substitution matrices**: New `--type pfasum`, `pfasum43`,
  `pfasum60` options for PFASUM-family scoring.
- **Python `align()` parity**: The in-memory `align()` function now
  supports `mode`, `vsm_amax`, `realign`, and `ensemble_seed` parameters
  (previously only available on file-based APIs).
- **New Python constants**: `MODE_DEFAULT`, `MODE_FAST`, `MODE_PRECISE`,
  `REFINE_INLINE`, `PROTEIN_PFASUM43`, `PROTEIN_PFASUM60`,
  `PROTEIN_PFASUM_AUTO`.
- **Container support**: Podman/Docker files for reproducible benchmarks.
- **Downstream benchmark suite** (`benchmarks/downstream/`): Measures
  alignment quality on downstream tasks: HMMER profile search,
  phylogenetic tree accuracy (RF distance), and positive selection
  detection (HyPhy BUSTED).

### Changed

- **License changed from GPL-3.0-or-later to Apache-2.0.** All
  dependencies are permissive-licensed, so no compatibility issues.
- **Default mode now uses consistency anchors (K=5) and VSM.** The
  previous default (no consistency, no VSM) is now `--fast`.
- Refine default is `none` (Python docstring previously said
  `confident` incorrectly; now fixed).
- `seq_weights` auto-disabled in ensemble mode.
- CLI `--ensemble` without a value defaults to 5 runs.
- C library API expanded with `kalign_run_seeded()`,
  `kalign_run_realign()`, `kalign_ensemble()`,
  `kalign_consensus_from_poar()`, `kalign_msa_compare_detailed()`,
  `kalign_msa_compare_with_mask()`.

### Fixed

- `kalign_arr_to_msa()` now initializes all MSA struct fields (biotype,
  alnlen, seq_distances, col_confidence, run_parallel, seq->confidence,
  seq->rank). Previously caused `free(): invalid pointer` on Linux
  where glibc does not zero malloc'd memory.
- Python `align()` (in-memory path) no longer crashes on Linux due to
  the uninitialized fields above.

### BAliBASE results (218 cases, protein)

| Mode | Recall | Precision | F1 | TC | Time |
|------|--------|-----------|-----|-----|------|
| fast (v3.4 equiv) | 0.804 | 0.656 | 0.716 | 0.466 | 10s |
| **default** | 0.810 | 0.665 | 0.724 | 0.472 | 29s |
| **precise** | 0.796 | 0.752 | 0.768 | 0.467 | 193s |
| ClustalO | 0.840 | 0.710 | 0.764 | 0.559 | -- |
| MAFFT | 0.867 | 0.715 | 0.778 | 0.590 | -- |
| MUSCLE5 | 0.870 | 0.721 | 0.783 | 0.581 | -- |
