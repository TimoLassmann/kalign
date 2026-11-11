# Kalign 4 protein accurate: confidence re-fit (the ensconf plan)

Status: written 2026-10-28, before any result of this experiment exists (the identity checks of section 6 ran
first; they score nothing that this plan selects on). Code: `recipe/ensfit_conf.py`, `recipe/scoring_conf.py`,
`recipe/_align_worker_conf.py`. Decision record: manuscript `k4/DECISIONS.md` D5; evidence that led here:
`k4/explore/ensemble_width/FINDINGS.md` (collapse at 1,000 sequences from the members' k-means guide tree; Tests A and
B of 2026-10-28). Those tests used the fit's own training/test data out of fold and in-sample; they motivated the
design, so this experiment's out-of-fold result is not a blind confirmation. That is stated wherever it is reported.

Question: does an ensemble whose members share the fast setting's guide-tree algorithm, perturbed per member,
produce a complete alignment that is at least as accurate as the fast setting (SP, TC) and whose column confidence
lowers the pair Brier score below that of the fast alignment, out of fold?

## 1. Fixed

| Fixed | Value |
|---|---|
| Engine | kalign `exp/sl-tree-noise` ea02ca9 (adds `tree_sl_noise`; byte-identical to release 9a49d37 at noise 0, kat_slnoise.json K1-K3), registered as `WS/kalign/ea02ca9` |
| Members | K = 5; per fit, the members of that fit's ensemble pick (`runs/ensfit/<fit>/selection_a4.json`, PLAN_ensemble.md + amendments 1-4) with their scoring, refinement, consistency, realign, sequence-weight, VSM and k-means tree-noise values unchanged |
| Guide tree | every member uses member 1's tree keys (the fast setting's: hybrid k-means + single linkage, handover 1,163, LCS indel, reduced alphabet) |
| Tree seeds | 42 + k (engine); member 1 has no noise, so member 1 is the fast alignment |
| Ensemble options | consistency merge, fill-in as in the fit's pick |
| Threads | 1 per alignment; 1,800 s per case |

## 2. Grid (per fit)

| Variable | Values |
|---|---|
| single-linkage tree noise sigma of members 2-5 (`run_tree_sl_noise` = [0, s, s, s, s]) | 0, 0.05, 0.1, 0.2, 0.3, 0.5 |
| consensus support `min_support` | 1, 2 |

12 points per fit. Six fits: folds 0-4 (held-out estimate) and the final fit on all fit data (the shipped setting).

## 3. Data

The fit data of PLAN_ensemble.md / PLAN_fastfit.md: per fit, the training core (small core of 840 cases + every
training family at 1,000); out of fold, each fold's test cases (small; families at 1,000). Not 10,000: accurate mode
refuses more than 5,000 sequences (release 6e691c6). No lockbox input is read (leak guards of `load_ens_fold` and the
evaluator's allowed set).

## 4. Objective and selection

Pair Brier score B of a case (`recipe/scoring_conf.py`). Population: every pair of scored residues (core columns of
the reference for the small pool, every column for extHomFam) from two different reference sequences. y = 1 if the
reference aligns the pair; p = the column confidence of the test column holding both residues, or 0 if the test does
not align them. B = sum (p - y)^2 / number of reference pairs. A pair neither alignment aligns contributes 0, so all
alignments of a case are scored on the same population and leaving residues unaligned cannot lower B (it costs 1 per
missed reference pair). A single alignment without confidence (the fast setting) is scored with p = 1:
B = (missed + wrong pairs) / reference pairs.

B_small = mean over benchmarks of the mean B of their cases (benchmarks weighted equally); B_1k = mean over families.
J = (B_small + B_1k) / 2. SP (recall) and TC as the pool scorers.

Selection on the training core: eligible = no failed case; SP_small, TC_small, SP_1k, TC_1k each >= the fast
setting's on the same cases; no family at 1,000 with SP more than 0.2 below the fast setting's. Pick = lowest J; ties:
smaller sigma, then smaller support. No eligible point: the fit has no pick and counts as the fast setting out of fold.
CPU is recorded and reported, not used (members are fixed, so the grid barely changes CPU).

## 5. Out of fold and gate

Arms on each fold's test cases: the fold's pick, the fast setting, the fold's ensemble pick as fitted ("asfit"), and
asfit with member 1's tree in every member ("tree1"). 95 % bootstrap intervals over groups (small) and families
(1,000), 10,000 resamples, seed 1234; a failed case scores SP = TC = 0, B = 1.

Gate (pooled folds 0-4, pick minus fast setting): the final fit's pick becomes protein accurate only if all hold:
1. small: SP and TC difference >= 0, lower bound >= -0.005;
2. 1,000: SP and TC difference >= 0, lower bound >= -0.005;
3. at most 3 of the out-of-fold families at 1,000 with SP more than 0.2 below the fast setting's (asfit: 9 of 308,
   FAMSA 2: 8);
4. B below the fast setting's: upper bound of the difference < 0, small and 1,000;
5. every out-of-fold input completes.
Otherwise no accurate setting comes out of this experiment, and that is reported. Reported beside the gate: pick
minus asfit for every measure, per-family drop counts, CPU. Confirmation beyond out of fold (the ensemble lockbox,
already aligned once with the shipped accurate on 2026-10-20 and reported descriptively, and/or fresh held-out
families) is decided separately (decision D7) and is not part of this gate.

## 6. Identity checks (run before the fit; `python -m recipe.ensfit_conf check`, runs/ensconf/identity.json)

On 20 small cases and 2 families at 1,000 of the final fit's core:
- I1: the confidence worker with build ea02ca9 reproduces the shipped evaluator's cached rows (build 5b93950,
  recipe/_align_worker.py) for the fast setting and the final ensemble pick: same alignment sha256, SP, precision, TC;
- I2: sigma 0 with `run_tree_sl_noise` gives the same alignment as the same dict without it;
- I3: correct pairs / reference pairs of the confidence scorer equals the pool scorer's recall (1e-12);
- I4: the fast setting is scored as a hard alignment and the ensemble with its confidence; B of the hard alignment
  equals (missed + wrong pairs) / reference pairs.
Result 2026-10-28: I1 44/44 (old rows found 44/44), I2 22/22, I3 66/66, I4 22/22.

## 7. Cost

Measured in the identity check: about 55 CPU-s per family at 1,000 and 0.36 CPU-s per small case for a shared-tree
ensemble. 12 points x 6 fits x (840 small + 246-308 families) is about 1.0 M CPU-s; out of fold about 0.05 M. About
8-10 h on the shared 60-slot pool. Chain: fit (fits 0-4, all) -> oof -> gate; log runs/ensconf/chain.log.
