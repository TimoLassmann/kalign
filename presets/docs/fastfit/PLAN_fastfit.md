# Kalign 4 fast mode: one setting for small and large protein alignments (the fast-fit plan)

Status: plan of the one pre-specified experiment that follows the explorations in
`k4/explore/EXPLORATION_LOG.md` (manuscript repository). It is frozen (commit, tag, sha256 of this file and of the
code in `presets/`) before any registered fit starts. Every number below that is not a rule was measured before
the freeze (section 11). Nothing in this plan was tuned on a fit result of this experiment.

Question: can a single-run (fast-mode) protein setting, found by nested cross-validation over small benchmarks and
large extHomFam families together, with the guide-tree construction searched alongside the alignment parameters,
match Kalign 4.0.0 fast on the small benchmarks and do better on large alignments, within a stated CPU limit?

## 1. Fixed parts of the engine

| Fixed | Value | Reason |
|---|---|---|
| Gap model | member gap model (build 99aecb2, branch exp/fast-tree-fit, on by default) | the engine's gap model; not a parameter |
| Runs | 1 (no ensemble) | fast mode; the ensemble is a later experiment that must beat the result of this one |
| Pairwise-alignment guide distance (`run_guide_dist = 1`) | off | aligns every pair; infeasible at large N |
| Anchors for k-means | 32 | other values gave non-deterministic output (PARKED) |
| Residues per sequence for tree distances | first 1,024 | known limitation: the Myers kernel caps its pattern at 1,024; kept for all distances so every option sees the same input; whole sequences changed 16 of 2,916 pool alignments |
| Threads per alignment | 1 | CPU objective (section 4) |

## 2. Search space

Alignment parameters: the stage-1 single-run variables of the gm1 refit (`experiments/redesign/tools/refit/space.py`,
same bounds and masking), without `guide_dist`.

| Variable | Range |
|---|---|
| gap open | [2, 15] |
| gap extend | [0.5, 5] |
| terminal gap extend | [0.1, 3] |
| matrix | PFASUM43, PFASUM60, CorBLOSUM66 |
| VSM | mode {legacy, per-merge}; amax [0, 5]; kappa [0, 1]; d0 [0.5, 10] |
| refinement | none, all, confident, inline |
| sequence weights | [0, 5] |
| consistency rounds | {0, 1, 2, 3, 5, 8, 10}; weight [0.5, 5] |
| realign | {0, 1, 2} |

Guide-tree variables (per-run options of build 99aecb2):

| Variable | Values |
|---|---|
| tree method | k-means (bisecting k-means, WPGMA below the handover size; 4.0.0's design) or hybrid (bisecting k-means, single linkage by Prim below the handover size) |
| handover size | 50 to 5,000, searched on a log scale (log10 in [1.699, 3.699], rounded to an integer) |
| tree distance | Myers raw, Myers normalised, LCS indel (FAMSA's), 1 - LCS / shorter |
| tree alphabet | reduced (13 letters) or full |
| k-means distance | 4.0.0's BPM distance to the anchors, or the tree distance |
| realign-pass tree | the 4.0.0 tree or the chosen tree |

An input smaller than the handover size gets a tree built entirely by WPGMA or Prim (no k-means), so pure WPGMA and
pure single-linkage trees are in the space for inputs up to 5,000 sequences. Kalign's "UPGMA" is WPGMA (plain mean
of the two merged distances, plus 0.001) and is reported as WPGMA.

Masking (one cache key per distinct output): k-means method with the BPM k-means distance never reads the tree
distance or alphabet (both set to canonical values); the realign-pass tree is canonical when realign = 0.

## 3. Data

**Small benchmarks.** The gm1 refit's protein pool (2,916 cases: BAliBASE 3, PREFAB 4, OXBench, SABRE, BaliFam100),
its 5 grouped folds (`folds_protein.tsv`) and its per-fold training cores (800 + 40 guard cases) and screening
subsets (300 cases), unchanged.

**Large alignments.** extHomFam v37 (390 families; reference alignments from HOMSTRAD, 3-4 reference sequences in
most families). Families with at least 1,000 sequences: 372; at least 10,000: 310.

- Subsets: all reference sequences plus randomly drawn other sequences of the family (seed 1), in file order, to
  1,000 and 10,000 sequences (the subset rule of the explorations).
- Scoring: SP and TC of the reference sequences only (all-gap columns of those rows dropped), against the
  reference alignment (`kalign.compare_detailed`, all columns). q = (SP + TC) / 2.
- **Lockbox (held out, scored once, at the end):** 60 families drawn first (seed 1234, stratified by extHomFam
  size category) from the families not used by any exploration, plus the four 200,000-sequence inputs of the paper
  benchmark (cox, helicase_NC, oat, TyrKc; those families are kept out of the fit). Each lockbox family is scored at
  1,000 sequences, at 10,000 where it has them, and 10 of them (drawn the same way from those with at least 50,000)
  at 50,000.
- **Fit families:** all other families with at least 1,000 sequences, including the 23 used by explorations (listed in `recipe/largefam.py` EXPLORED). Each
  is scored at 1,000 and, where it has them, at 10,000. A family is one group: both its subsets are in the same
  fold.
- Folds: fit families are assigned to the 5 folds by seeded hash order (seed 1234), stratified by size category;
  the small pool keeps its folds. Fold k's training data = the small training core of fold k + the fit families of
  the other folds. Fold "all" (final fit) = the small core of fold "all" + all fit families.
- Per-fold screening data: the small screening subset (300) + 40 training families at 1,000 + 12 at 10,000
  (seeded, stratified by size category). Core data (re-scoring and selection): the small core + every training
  family at 1,000 and at 10,000.

## 4. Objectives and constraints

Computed on the same set of cases for every candidate:

1. Small accuracy Q_small: per-case q = (SP + TC) / 2, averaged within benchmark, benchmarks weighted equally.
2. Large accuracy Q_large: per-family q, averaged over families, at 1,000 and at 10,000; Q_large = mean of the two.
3. CPU: R = sqrt( (CPU_small / CPU_small,ref) x (CPU_large / CPU_large,ref) ), total single-thread user + system CPU
   seconds of the aligner process, divided by Kalign 4.0.0 fast on the same cases.

Reference ("4.0.0 fast"): the shipped 4.0.0 fast preset run by the 4.0.0 build (daff955), single-threaded, through
the same worker (identity with the CLI preset checked, section 11). Constraints (infeasible otherwise): no failed
case (error or timeout, 600 s per case); CPU_small <= 1.5 x reference; CPU_large <= 2 x reference.

CPU seconds are measured at full machine load (64 hardware threads, up to 60 aligners at once). Noise is measured,
not removed: a calibration process re-runs 4.0.0 fast (member gap model build) on a fixed set (20 small cases, 3
families at 1,000, 1 at 10,000) every 30 minutes without the cache, and its coefficient of variation is reported.
Speed claims in the paper come from the controlled timing block, not from the fit.

## 5. Optimiser

NSGA-III (pymoo 0.6.1.6), population 56, Das-Dennis reference directions for 3 objectives (9 partitions, 55
directions), seed 1234, the mixed-variable operators of the gm1 refit. Initial population: 6 seeds (4.0.0 fast
parameters and the gm1 fast preset, each with the 4.0.0 tree, with the E21 tree (hybrid, handover 2,000, LCS indel,
full alphabet, k-means on the tree distance), and with the E21 tree and realign 0) + random individuals. Each
generation is scored on the screening data. Stop: at least 20 generations and archive hypervolume gain < 0.1 % over
the last 5, or 100 generations. Objectives are normalised for NSGA-III and the hypervolume by a monotone transform:
accuracy axes (1 - Q) / (1 - (Q_ref - 0.10)); CPU log(R) mapped from R = 0.25 to R = 3; reference point 1.1 on
every axis. A candidate is stopped early (infeasible) at its first failed case or once its CPU on the evaluated
cases exceeds 3 x the reference.

Re-scoring: the non-dominated feasible screening points with Q_small >= Q_small,ref - 0.01 (screening data) are
re-scored on the core data, at most 40 of them: if more qualify, the 20 with the highest Q_large, then the 20 with
the lowest R among the rest, then by Q_large. The two references (4.0.0 fast; gm1 fast preset) are always scored on
the core.

## 6. Selection (identical in every fold and in the final fit; core data)

- **Fast pick:** feasible, Q_small >= Q_small,ref - 0.002; the highest Q_large; ties: lower R, then key. No eligible
  point: the gm1 fast preset.
- **Ultra-fast pick:** feasible, Q_small >= Q_small,ref - 0.01 and Q_large >= Q_large,ref; the lowest R; ties:
  higher Q_large, then key. No eligible point: none.

## 7. Evaluation and gate

Out of fold: for each fold k, the fold's picks, 4.0.0 fast and the gm1 fast preset on all test cases of fold k
(small cases of the fold; its fit families at 1,000 and 10,000). Confidence intervals: 95 %, bootstrap over groups
(small) and families (large), 10,000 resamples, seed 1234.

The fast pick of the final fit ships as the new fast preset only if all hold:
1. out-of-fold small SP minus 4.0.0 fast: lower confidence bound >= -0.005;
2. out-of-fold Q_large minus 4.0.0 fast: lower confidence bound > 0;
3. lockbox: Q at 1,000 and Q at 10,000 both above 4.0.0 fast;
4. every 50,000 and 200,000 lockbox input completes (8 threads, 4 h per input).
Otherwise the gm1 fast preset ships.

The ultra-fast pick of the final fit is reported as a candidate extra mode only if: out-of-fold small SP minus 4.0.0
fast, lower bound >= -0.015; out-of-fold Q_large minus 4.0.0 fast >= 0; out-of-fold R <= 0.8; lockbox Q at 1,000 and
10,000 >= 4.0.0 fast; the 50,000 and 200,000 inputs complete.

Reported regardless of the gate: every evaluated point (fold, generation, parameters, three objectives), the
fronts of every fold and of the final fit, where the references and both picks sit, lockbox results of both picks,
FAMSA on the lockbox inputs, the calibration noise, and a repeat of the final fit with seed 4321 and a fully random
initial population (reported, not a gate).

## 8. Known limitations (stated in the paper)

- tree distances use the first 1,024 residues of each sequence;
- k-means uses 32 anchors;
- most extHomFam families have 3-4 reference sequences, so each family's score is noisy; the result rests on
  many families;
- CPU seconds measured under load carry noise (section 4).

## 9. Leak guards

A fit receives only its training case ids and training families; the evaluator refuses anything else. No lockbox
family or input is reachable by any fit or by the out-of-fold step; lockbox scoring runs once, after the final pick
is written and hashed.

## 10. Outputs

Per fit: run.json, refs.json, evaluated.jsonl (every candidate), hv.json, front_core.json, selection.json.
Pooled: out-of-fold rows, gate table, lockbox rows. Report: HTML (fronts per fold and final, front over
generations, references and picks, parameter table of the front).

## 11. Measured before the freeze

- single-thread CPU of 4.0.0 fast (member gap model build), 1,000 sequences: 0.3-1.6 s (cyt5, rnr, WW); 10,000:
  2.7-31 s (kex, lacI, ghf1, mdd, EGF_Lam); hybrid tree adds about 40 %.
- identity checks run before the freeze: (I1) the 4.0.0 build's worker with the 4.0.0 fast parameters gives the
  same alignments as its CLI preset; (I2) explicit 4.0.0 tree values give the same alignments as the defaults on
  build 99aecb2; (I3) every subset contains all reference sequences and each reference alignment scores SP = TC = 1
  against itself.
