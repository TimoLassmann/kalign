# Kalign 4 accurate mode (protein): an ensemble that aligns only what it can support (the ensemble plan)

Status: FROZEN 2026-10-15 (commit and tag `ens-plan-freeze`, sha256 of this file in the tag message), BEFORE the
fast-fit lockbox is scored, so that no decision below can depend on lockbox results (section 9). Measurements made
before the freeze are in section 11; nothing is tuned on a fit result of this experiment.

Question: can an ensemble protein setting, found by nested cross-validation and built from five diverse members of
which one is the new fast setting, align residues more precisely (a larger share of its aligned residue pairs are
correct) than the fast setting, while recovering as many correct pairs (SP) as the fast setting, within a stated CPU
limit, and with per-column confidence that tells which columns are correct?

Kalign 4 has four settings: fast and accurate, for protein and for nucleotides. Protein fast is fixed (fast-fit pick
0976900c, tag k4-fast-0976900c). This experiment finds protein accurate. The pre-release accurate preset (parameter
key 2852d1c2, called "pre-release accurate" below) is a starting point for the search and a reported row, not a
reference to be beaten or kept.

## 1. Fixed parts of the engine

| Fixed | Value | Reason |
|---|---|---|
| Build | `release/fast-preset` at 5b93950 (per-member settings, below), tagged at the freeze | the shipped engine |
| Members | K = 5 | as the pre-release accurate preset; fixed |
| Member 1 | the fast setting 0976900c, every parameter frozen (gaps, matrix, VSM, refinement, sequence weights 2.86, no realign, no consistency, VSM mode legacy, its guide tree) | accurate mode contains the fast alignment; the ensemble must beat it |
| Gap model | member gap model in every member | the engine's gap model |
| Consensus | POAR consensus with per-column confidence, its options searched (section 2) | the source of precision and of confidence |
| Size | no size cut-off: accurate mode runs the ensemble at every input size; the user chooses the mode | |
| Pairwise-alignment guide distance | off | aligns every pair; infeasible at large N |
| Anchors for k-means and for realign distances | 32 | as the fast fit |
| Residues per sequence for tree distances | first 1,024 | as the fast fit |
| Threads per alignment | 1 | CPU objective (section 4) |

Engine behaviour that is part of the setting and stays as it is (ensemble.c, consensus_fill.c): the consistency-merge
rerun uses the winning member's gaps and matrix with sequence weights 0, no consistency and no realign; fill-in segment
runs use no realign and no consistency. Without a consensus (min_support 0) the engine returns the best-scoring member.

## 2. Search space

Diversity is the point of the ensemble: members 2..5 are searched per member, so they can differ from member 1 and
from each other.

Ensemble-level:

| Variable | Range |
|---|---|
| min_support | 0..5 (0: best member, no consensus) |
| consistency_merge | {0, 1}; weight [0.5, 10] |
| fill-in | {off, realign, project, attach}; anchor fraction [0, 1]; support 1..5 |

Per member k = 2..5 (member 1 fixed):

| Variable | Range |
|---|---|
| gap open, extend, terminal extend | [2, 15], [0.5, 5], [0.1, 3] |
| matrix | PFASUM43, PFASUM60, CorBLOSUM66 |
| VSM amax / kappa / d0 | [0, 5] / [0, 1] / [0.5, 10] |
| VSM mode | legacy, per-merge |
| sequence weights | [0, 5] |
| realign | {0, 1, 2} |
| consistency rounds; weight | {0, 1, 2, 3, 5, 8, 10}; [0.5, 5] |
| refinement | none, all, confident, inline |
| tree noise | [0, 0.5] |
| tree method | k-means or hybrid (k-means + Prim) |
| handover size | 50 to 5,000 (log scale) |
| tree distance | Myers raw, Myers normalised, LCS indel, 1 - LCS / shorter |
| tree alphabet | reduced or full |
| k-means distance | BPM to anchors, or the tree distance |

Tree seeds are 42 + k. Masking as in the fast fit (one cache key per distinct output).

## 3. Data

Identical to PLAN_fastfit.md section 3, so fast and accurate are reported on the same cases: the 2,916-case small
protein pool with its 5 grouped folds, cores and screening subsets; the 308 extHomFam fit families in the same folds;
the same lockbox (64 families + the four 200,000-sequence inputs).

| Size | Role |
|---|---|
| small pool | search (screening subset of 300), re-scoring and selection (core), out-of-fold evaluation |
| 1,000 | search (40 training families per fold, as the fast fit), re-scoring and selection (every training family), out-of-fold evaluation (every held-out family) |
| 10,000 | in the search only as one timing guard input per fold (section 4); re-scoring and selection on 24 training families per fold (seeded 1234, stratified by size category, drawn from the families of the fold's training data that have 10,000 sequences); out-of-fold evaluation on every held-out 10,000 input |
| 50,000 | measured before the freeze on 3 fit families (section 11); lockbox (10 inputs): completion, wall time, peak memory reported |
| 200,000 | lockbox (4 inputs): completion, wall time, peak memory reported |

The 10,000 inputs are not searched for accuracy because one ensemble costs hundreds of CPU-s per input there
(section 11: the pre-release design took 296-1,496 CPU-s and timed out on 2 of 3 families at 1,800 s); scoring 12 of
them for every candidate of six fits would take weeks.

## 4. Objectives and constraints

On the same cases for every candidate (screening data during the search):

1. Small precision P_small: per-case precision (correct aligned pairs / aligned pairs), averaged within benchmark,
   benchmarks weighted equally.
2. Large precision P_1k: per-family precision at 1,000 sequences.
3. CPU: R = sqrt((CPU_small / ref) x (CPU_1k / ref)), single-thread CPU, reference = the fast setting on the same
   cases.

Constraints (infeasible otherwise):
- recall floor: small SP >= small SP of the fast setting - 0.005, and SP at 1,000 >= SP at 1,000 of the fast setting
  - 0.005 (same cases). Precision can be raised by aligning fewer residues; the floor rules that out. The share of
  residues left unaligned (projected width ratio, test width / reference width) is reported for every pick;
- no failed case (error or timeout, 1,800 s per case);
- 10,000 timing guard: every candidate also aligns one 10,000-sequence training family of its fold (single thread),
  fixed per fold: of the fold's 24 selection families at 10,000 (section 3), the one with the median CPU of the fast
  setting (ties: family name); more than 600 s CPU, or an error, makes the candidate infeasible. Its accuracy is not
  used;
- CPU_small <= 7.1 x fast setting and CPU_1k <= 53.6 x fast setting: twice the ratio measured for the pre-release
  accurate design with the fast setting as member 1 (section 11: 3.55 and 26.8).

References scored on every data set: the fast setting (the bar), the pre-release accurate preset (reported), FAMSA
(reported).

## 5. Optimiser

NSGA-III (pymoo 0.6.1.6), population 56, Das-Dennis directions for 3 objectives (9 partitions), seed 1234, the
mixed-variable operators of the fast fit. Initial population: the pre-release accurate members 2..5 and ensemble options
with the fast setting as member 1; the same with the fast setting's guide tree in every member; random individuals.
Stop: at least 20 generations and archive hypervolume gain < 0.1 % over 5, or 60 generations. Screening, re-scoring
(at most 40 points per fit on the core) and objective normalisation as in PLAN_fastfit.md section 5, with precision
in place of Q and the fast setting as the reference.

## 6. Selection (identical in every fold and in the final fit; core data)

Knee rule of DECISION_fastfit_selection.md, with precision:
1. eligible: feasible on the core; small precision > fast setting; small SP >= fast setting - 0.005; SP at 1,000 >=
   fast setting - 0.005; on the 24 families at 10,000: every input completes (1,800 s, 1 thread), precision >= fast
   setting, SP >= fast setting - 0.005;
2. b = the eligible point with the highest P_1k;
3. keep every eligible c with P_1k(b) - P_1k(c) <= 1.96 x SE_c (paired family bootstrap on the core, 2,000
   resamples, seed 1234);
4. pick the lowest R; ties: higher P_1k, then key.
No eligible point: no accurate setting comes out of this experiment, and that is reported.

## 7. Evaluation and gate

Out of fold: each fold's pick, the fast setting, the pre-release accurate preset and FAMSA on all test cases of the
fold (small, 1,000, 10,000). 95 % bootstrap CIs over groups (small) and families (large), 10,000 resamples, seed 1234.

The final pick becomes protein accurate only if all hold (differences are pick minus fast setting):
1. out-of-fold small precision: lower bound > 0;
2. out-of-fold small SP: estimate >= -0.005 and lower bound >= -0.015;
3. out-of-fold precision at 1,000: lower bound > 0; SP at 1,000: estimate >= -0.005;
4. out-of-fold at 10,000: precision estimate > 0, SP estimate >= -0.005;
5. every out-of-fold input at 1,000 and 10,000 completes (8 threads, 1 h per input);
6. lockbox (second use, declared in section 9): precision at 1,000 and at 10,000 above the fast setting, SP not below
   the fast setting - 0.005; the 50,000 and 200,000 inputs: completion, wall time and peak memory reported.
Otherwise no accurate setting comes out of this experiment, and that is reported.

## 8. Confidence and the precision dial (secondary outcomes, reported, not a gate)

On the out-of-fold cases, for the pick and the pre-release accurate preset:
- AUROC of per-column confidence for column correctness, per benchmark and pooled; precision of the columns with
  confidence >= 0.5, 0.7, 0.9 and the share of columns kept (the `--confidence-threshold` use); the calibration
  analysis of the paper (fig_calibration) on the pick;
- the pick with min_support set to 1, 2, 3, 4, 5 and fill-in off (members unchanged): precision, SP and width ratio
  at each, so the documentation can state what a stricter support threshold buys.

## 9. Leak guards and the second use of the lockbox

As PLAN_fastfit.md section 9. The lockbox was created for the fast fit and is scored for it once. It is used a second
time, once, for this experiment's gate (criterion 6). This is declared here, before the fast fit's lockbox results
exist, so no choice in this plan can depend on them. After the fast-fit lockbox is opened, this plan is not changed;
any change is a new plan that cannot use the lockbox.

## 10. Outputs

As the fast fit: per fit run.json, refs.json, evaluated.jsonl, hv.json, front_core.json, selection.json; pooled
out-of-fold rows, gate table, lockbox rows; HTML report.

## 11. Measured before the freeze

Done on 2026-10-15 (manuscript k4/explore/ensemble_premeasure; workstation ~/tmp/kalign-fastfit/runs/ens_premeasure;
build 5b93950), at moderate machine load, single thread through the fit's Python path:

| CPU-s (total) | fast setting | pre-release design, fast as member 1 | pre-release accurate |
|---|---|---|---|
| small screening subset (300 cases) | 40.4 | 143.2 | 164.5 |
| 3 families at 1,000 | 3.1 | 83.2 | 77.6 |
| 3 families at 10,000 | 25.3 | GAL4 320; Extradiol_dioxy timed out at 1,800 (1,496 without a limit); serpin timed out | GAL4 296; 2 timed out |
| GAL4 at 50,000, 8 threads | | 15,712 CPU-s, 4,609 s wall, 2.1 GB | |

Ratios to the fast setting (pre-release design, fast as member 1): small 3.55, 1,000 26.8 -> the CPU limits of
section 4. Python path = C CLI byte for byte for the pre-release accurate preset at 1,000 (3/3).
Still running at the freeze and reported later (they change no rule above): the remaining 10,000 and 50,000 timings,
and the per-member identity checks below.

Planned checks (fit data only, never a lockbox input):
- identity: a 1-member ensemble with the fast setting equals `--mode fast` byte for byte; the pre-release accurate
  preset through the Python path equals `--mode accurate`; a per-member value reaches the output (a 5-member run whose
  members differ only in that value, with the winning member logged), for sequence weights, realign, consistency and
  VSM mode;
- single-thread CPU, wall time and peak memory of the fast setting and of the pre-release accurate design with the
  fast setting as member 1: small-pool screening subset, 3 families at 1,000, 3 at 10,000, and (8 threads) 3 at 50,000
  -> C, the 10,000 cost estimate, and the projected wall time of the fit.
