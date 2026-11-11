# Implementation notes: the ensemble fit (PLAN_ensemble.md)

The plan (`PLAN_ensemble.md`, tag `ens-plan-freeze`, sha256 77d2f157...9a23) is not changed. This file records every
detail the plan leaves open and the reading chosen for it (the most literal one consistent with `PLAN_fastfit.md`),
plus two places where the plan's wording and the engine disagree. Code: `presets/recipe/ensfit.py` (the fit),
`ensfit_prep.py` (10,000 selection and guard families), `ensfit_identity.py` (identity checks), `ensfit_run.py`
(launcher); `_align_worker.py` forwards the per-member lists when a dict has them (older dicts and keys unchanged).

## Engine, settings, data

1. **Engine.** Build 5b93950 on the workstation (`~/tmp/kalign-rel/5b93950`: CLI and the venv with the Python
   module), registered as `WS/kalign/5b93950/build.json` (build id = sha256 of the CLI). Every alignment goes through
   the fast fit's evaluator path (`_align_worker.py` -> `_core.ensemble_custom_file_to_file`, 1 thread, seed 42).
2. **Cache.** `WS/cache/ensfit.sqlite`, separate from the fast fit's cache; keys carry the build id as well.
3. **Member 1 / reference.** `WS/runs/fastfit/final_pick.json` (key 0976900c, re-checked) plus `run_gap_model [1]`.
   Identity check a: equal to `kalign --mode fast -n 1` of the f739560 release CLI.
4. **Pre-release accurate.** `presets/results/v2/presets_v2.json` of the 5b93950 tree (sha256 040f3af2...), protein
   accurate, key 2852d1c2 (equal to the param key of its dict), used exactly as the release header carries it (gap model
   0, shipped tree) for the reported row. Identity check b: equal to `kalign --mode accurate -n 1` of the 5b93950 CLI.
5. **Seeds (plan section 5).** Seed 1 = the pre-release accurate members 2..5 (their gaps, matrices, VSM, refinement,
   sequence weights 3.27, realign 2, consistency 8 / 3.558, VSM legacy, noise 0) and ensemble options (min_support 2,
   fill-in realign, anchor fraction 0.0796, support 1), member 1 = the fast setting, member gap model in every member;
   the shipped tree of members 2..5 is written as explicit values (k-means, handover 50, Myers raw, reduced alphabet,
   BPM k-means distance), which is byte-identical to leaving it at the engine default (identity check d). Seed 2 = the
   same with the fast setting's tree in members 2..5 (k-means + Prim, handover 1,163, LCS indel, reduced alphabet,
   k-means on the tree distance). The other 54 individuals: pymoo `MixedVariableSampling`, `default_rng(1234)`, as the
   fast fit.
6. **Variable encodings** as the fast fit: handover size searched as log10 in [1.699, 3.699], rounded; realign
   Integer 0..2; consistency rounds Choice {0,1,2,3,5,8,10}; min_support Integer 0..5; fill support Integer 1..5;
   matrix Choice {PFASUM43, PFASUM60, CorBLOSUM66}; tree noise Real [0, 0.5]. 82 variables (19 per member x 4 + 6).
7. **Not searched, fixed.** The realign-pass tree (a fast-fit variable) is not in the plan's member list: fixed at 1
   (the chosen tree also in the realign pass), the engine default and the fast setting's value. Guide distance BPM,
   1,024 residues, member gap model, 32 anchors (engine constant).
8. **Tree seeds.** `_core` sets tree seed = 42 + k for member index k = 0..4 (seed 42); member j (1-based) gets
   42 + j - 1. Tree noise only acts when it is > 0.
9. **Masking (one cache key per distinct output).** Applied: consistency rounds 0 -> consistency weight 1.0; VSM
   legacy -> kappa = d0 = 1.0; VSM per-merge -> amax 0 (ignored by the engine); k-means tree with the BPM k-means
   distance -> tree distance and alphabet canonical; consistency merge 0 -> its weight 2.0; consistency merge 1 ->
   min_support 0 and fill-in off (the merge path in ensemble.c never reads them); anchor fraction only for fill-in
   realign/project; fill support only for attach. Not applied: the refit space's "tree noise 0 when realign > 0" rule,
   because it does not preserve the output (the realign pass starts from the alignment built on the noisy tree).
10. **Data per fold.** Screening = the small screening subset (300) + the fast fit's 40 screening families at 1,000
    (the `k1000` entries of `exthomfam_fold<f>_screen.txt`) + the fold's guard input at 10,000; the fast fit's 12
    screening inputs at 10,000 are not used (10,000 is not searched for accuracy). Core = the small core (800 + 40) +
    every training family at 1,000 + the fold's 24 selection families at 10,000; the fold's other training inputs at
    10,000 are not used. Fold "all" = all fit families. Leak guards: `load_fold`, the core/test/lockbox check and the
    evaluator's allowed set (the core) as in the fast fit; no lockbox input is reachable.
11. **24 selection families at 10,000** (`ensfit_prep select10k`): the fast fit's draw rule
    (`largefam.strat_pick`: largest-remainder allocation over extHomFam size categories, at least 1 per category,
    hash order sha256("1234|select10k<fold>|<family>")) over the fold's training families that have 10,000
    sequences.
12. **Guard family.** The fast setting aligns the 24 inputs of each fold through the fit's evaluator (single thread,
    the ensfit cache, so the fits reuse these rows); the guard is the lower median (12th smallest CPU of 24; ties by
    family name). Written to `WS/runs/ensfit/select10k.json` before launch; each fit records its sha256.

## Objectives, constraints, search

13. **Objectives.** P_small = per-case precision of the protein scorer (core columns / PREFAB pair core), mean within
    benchmark, benchmarks weighted equally; P_1k = per-family precision at 1,000, mean over families; R from CPU_small
    and CPU_1k only (the guard and the 10,000 inputs are not in R). CPU = in-call RUSAGE delta of the worker, as the
    fast fit.
14. **Constraints** (plan section 4) are computed the same way on the screening data (search) and on the core
    (re-scoring, selection), each time against the fast setting on the same cases: small SP and SP at 1,000 >= fast
    - 0.005; CPU_small <= 7.1 x and CPU_1k <= 53.6 x; guard CPU <= 600 s and no guard error (on the core the guard
    input is one of the 24 and its row is the same alignment); no failed case (error or timeout) on any evaluated
    input, the 24 at 10,000 included on the core. Timeout 1,800 s wall per case for every case (the fast fit's wall
    timeout standing in for CPU).
15. **Guard enforcement.** A search candidate's guard alignment is stopped once its CPU (from /proc) exceeds 600 s
    (row error "cpu_limit"): the candidate is infeasible exactly as if it had finished above 600 s. References run
    without this limit. A cached "cpu_limit" row is re-run in full whenever the same input is requested without the
    limit.
16. **Early stop.** PLAN_fastfit's stop "once the CPU exceeds 3 x the reference" would stop feasible ensemble
    candidates (limits 7.1 x and 53.6 x), so its principle is kept instead (stop only a candidate that is already
    certainly infeasible): at the first failed case, or once its CPU over the evaluated cases exceeds the sum of all
    limits (7.1 x fast small CPU + 53.6 x fast CPU at 1,000 + 600 s for the guard on the screening data; + 1,800 s per
    10,000 input on the core).
17. **Normalisation (literal).** Precision axes (1 - P) / (1 - (P_fast - 0.10)); CPU axis log R mapped from R = 0.25
    to R = 3, unchanged from PLAN_fastfit; hypervolume reference point 1.1. Consequence: a point with R above about
    3.85 lies outside the hypervolume box. If no feasible screening point has R below that, the hypervolume stays 0
    and the fast fit's stop rule (which needs a positive hypervolume) never fires: the search then runs to the cap of
    60 generations. The NSGA-III search itself is unaffected (the transform is monotone). See the pilot below.
18. **NSGA-III** as the fast fit (pymoo 0.6.1.6, population 56, Das-Dennis 3 objectives / 9 partitions, seed 1234,
    mixed-variable mating and duplicate elimination); an invalid candidate gets F = (2, 2, 2); 6 inequality
    constraints (the two recall floors, the two CPU limits, the guard, failed cases). Stop: >= 20 generations and HV
    gain < 0.1 % over 5, or 60.
19. **Re-scoring** as PLAN_fastfit section 5 with precision: the non-dominated feasible screening points on
    (P_small, P_1k, R) with P_small >= P_small,fast - 0.01; at most 40 (the 20 with the highest P_1k, then the 20
    with the lowest R among the rest, then by P_1k).

## References and selection

20. **References.** The fast setting on screen and core at the start (every objective and constraint is relative to
    it); the pre-release accurate preset on screen and core after re-scoring (reported row; not a selection candidate,
    since its member 1 is not the fast setting). FAMSA is not run inside the fits (the fast fit also scored FAMSA with
    a separate comparison harness); it is a reported row for the out-of-fold step.
21. **Selection** (plan section 6) runs at the end of each fit on the re-scored candidates: eligibility as listed in
    the plan, the 10,000 conditions on the means over the 24 inputs, all 24 complete; b = the highest P_1k (ties: lower
    R, then key); SE_c = paired family bootstrap of P_1k(b) - P_1k(c) over the core families at 1,000 (only P_1k is
    involved), 2,000 resamples, `default_rng(1234)`, the same resample for every c; pick = the lowest R of the points
    within 1.96 SE (ties: higher P_1k, then key). As in DECISION_fastfit_selection.md, the picks at fixed delta 0.005
    and 0.01 are reported alongside (not used). `selection.json` carries the pick's sha256 (canonical parameter JSON).
22. **Not in the ensemble plan, not done:** a seed-robustness repeat and the fast fit's 30-minute CPU calibration
    loop. Six fits: folds 0-4 and the final fit.

## Plan text vs engine

23. **min_support 0.** The plan's table says "0: best member, no consensus". The engine (ensemble.c, unchanged by plan
    section 1) does this with min_support 0: it builds a POAR consensus at max(2, (n + 2) / 3) = 2 and uses it when its
    POAR score beats the best member's; otherwise it re-runs the best member with confident refinement and keeps
    whichever scores higher. The fit searches min_support 0..5 with the engine's behaviour (section 1: "engine
    behaviour ... stays as it is"); no masking assumes the plan's description. This affects only the reading of
    min_support 0 (and section 8's secondary analysis, not part of the fits).

## Before launch (2026-10-15, workstation)

24. **Guard families** (`results/ensfit/select10k.json`; fast-setting CPU, single thread): fold 0 cyclin (4.0 s), 1
    ATP-synt (5.7 s), 2 PTB (4.1 s), 3 slectin (4.7 s), 4 gpdh (4.6 s), final cyclin (4.0 s).
25. **Identity** (`results/ensfit/identity_summary.json`): a 22/22 byte-identical (20 small + 2 at 1,000); b 20/20;
    d 22/22; aligner constants as expected. c (5 members = the fast setting, member k alone changed, k = 2..5, 22
    inputs each = 88 runs per setting; winner logged by ens_probe, whose output equals the evaluator's 88/88 in every
    row): with min_support 5 / fill-in off the output differs from all-equal in 64/88 (sequence weights), 64/88
    (realign), 72/88 (consistency), 32/88 (VSM mode) runs, exactly the runs where member k's POAR score differs and
    the other four are equal; with the pre-release options (min_support 2, fill-in realign) the output never differs
    (0/88: the four identical members support every pair of their alignment); member 1 wins every run (tie order).
    VSM mode reaches member k but changes it on fewer inputs (32/88).
26. **Pilot** (fold 0, generation 1 only, `runs/ensfit/pilot_fold0`, unregistered; its alignments are cached and are
    reused by the registered fold-0 fit): 56 candidates, 0 feasible: every guard alignment (seeds 1 and 2 included)
    passed 600 s CPU and was stopped there; the early stop then left most other screening cases unaligned. 48,006 CPU-s,
    1,977 s wall (one fit alone). The hypervolume is 0 (no feasible point).
27. **Projection** (6 fits sharing 60 slots, 60 generations, since a zero hypervolume never triggers the stop rule):
    search at the pilot's cost 6 x 60 x 48,006 CPU-s = 17.3 M CPU-s / 60 slots = 288,000 s = 3.3 days; a feasible candidate
    costs at most the sum of its limits (about 33 + 1,730 + 600 CPU-s on fold-0 screening), so the search stays within
    about 3-5 days if feasible points appear; re-scoring worst case 40 points x (97 + 10,300 + about 24,000) CPU-s x 6 /
    60 = 1.6 days; pre-release accurate on the cores about 0.1 day. Total under 7 days (limit 10): launched.

## Amendment 1 (2026-10-15)

28. `AMENDMENT_ensemble_1.md` changes the plan after the pilot: no 10,000 guard in the search (screening = 300 small +
    40 families at 1,000), the 10,000 check only at selection (section 6, unchanged), generation cap 30, gate out of
    fold only (no lockbox). Code: `ensfit.py` drops the guard input, its constraint (5 inequality constraints) and the
    guard CPU limit; the early-stop cap on the screening data no longer has a 10,000 term; the normaliser never used
    the guard. Notes 12, 15, 17 (60-generation consequence) and 27 are superseded where they concern the guard.
29. Cache keys: the registered fold-0 fit of the aborted launch reused the pilot's 56 generation-1 keys from the cache
    (initial-population keys reproduce exactly across processes; a pilot key's cached screening rows are found by the
    evaluator); the new keys seen then were generation-2 offspring. Core re-scoring reads the same (build, key, case)
    rows, so screening rows are reused.
30. The aborted launch (15:21-15:23) and the pilot are kept in `runs/ensfit/_aborted_20261015/`.

## Amendment 2 (2026-10-15)

31. `AMENDMENT_ensemble_2.md`: no CPU limits (constraint and early stop removed; 3 inequality constraints: the two SP
    floors and failed cases), two fits (fold 0 and final; `ensfit_run.FITS`), CPU axis log R from 30 to 50. The
    amendment-1 runs are kept in `runs/ensfit/_aborted_amend1/`; their alignments are in the cache and are reused.
