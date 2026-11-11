# Kalign presets by pooled nested cross-validation (version 1): the Freeze-1 plan

Status: final plan, to be committed and tagged as Freeze 1 before any fitting starts. It supersedes the draft
`experiments/redesign/docs/PLAN_pooled_refit_v1.md` (2026-10-05) and describes what is implemented in `presets/`,
not what was intended. Everything below that carries a number was measured or counted in step 0 and the build
tasks (`presets/docs/step0/`, `presets/docs/build/`). Nothing in this plan has been tuned on a fit result: no fit has
been run with registered settings.

Goal: new Kalign presets obtained in a way that is reproducible, defensible and honest, and a paper that can say
exactly how they were obtained. End state: the presets and one recipe (code) that regenerates them.

How to run (one command, workstation): `presets/run_all.sh start --freeze FREEZE1.json`. The command refuses to
start if `presets/` differs from the frozen commit (section 12).

## 1. The procedure in one picture

```
 FREEZE 1  (this document + code, committed and tagged before anything runs)
 ----------------------------------------------------------------------------------------------
 PROTEIN POOL  2,916 cases                         NUCLEOTIDE POOL  3,225 cases
 BAliBASE 3 (BB+BBS) 385, PREFAB 4 1,666,          BRAliBASE II 574, MDSA 299,
 OXBench 395, SABRE 411, BaliFam100 59             Rfam seed (all eligible families) 2,352
        |                                                   |
 [1] features from the REFERENCE alignment only (identity, #seqs, length, gap structure)
        |
 [2] groups of related cases (all benchmarks of a pool together)  ->  5 grouped folds
        |            (stratified by benchmark x difficulty; seed 1234)
        |
        +------ for each fold k = 0..4 (all five folds at once, sharing 30 aligner slots) -----+
        |  TRAINING = the other 4 folds             TEST = fold k, ALL its cases, never reduced |
        |    [3] training core (cap per group, difficulty bins, guard set, CPU cap)             |
        |    [4] NSGA-III, two stages, two speeds                                               |
        |    fixed rule selects fast / default / recall / accurate  ---> aligned on TEST        |
        +---------------------------------------------------------------------------------------+
        |   every alignment goes through ONE cache: (build, parameters, case) -> scores + CPU
        v
 POOLED OUT-OF-FOLD RESULTS on all cases (BCa CIs; per benchmark / homolog_in_train / difficulty)
   new presets, 3.7.0 presets, MAFFT L-INS-i, MUSCLE 5, Clustal Omega on the same cases
        |
 FREEZE 2  sha256 of every out-of-fold row file, recorded by run_all.sh BEFORE the final fit starts
        v
 FINAL FIT: the identical procedure on the whole pool (fold "all")  ->  candidate presets
        |  re-timing at 4 workers (finalists and 3.7.0), optimism estimate
        v
 DECISION per pool x mode (rule of this document): new preset ships only if it beats 3.7.0 out of fold by
 the stated margin within budget; otherwise 3.7.0 stays
```

## 2. Decisions that define the experiment

1. Real reference alignments only; simulations are version 2.
2. Two pools, protein and nucleotide (Kalign's two sequence types; its nucleotide type covers DNA and RNA). A
   separate DNA-only / RNA-only preset is a version-2 question.
3. Kalign never chooses its own benchmark: training cases are chosen from properties of the reference alignments and
   of the grouping only. No sensitivity screen on Kalign output. The one exception is a runtime cap, not a quality
   screen: protein cases whose 3.7.0 `accurate` run took more than 10 s CPU are left out of the training cores
   (section 6); they are still tested.
4. Parameters are never averaged across folds. The shipped presets come from the same procedure on the whole pool.
5. No wall-clock budget: each NSGA-III stage stops when its hypervolume stops improving (section 7).
6. The out-of-fold numbers are the paper's numbers. The final fit is evaluated in sample only to estimate optimism.

## 3. Data

Every benchmark is downloaded from its original source and checked against a SHA-256 in
`presets/resources/REGISTRY.yaml` (source, version, licence, size). The recipe never redistributes data.

| Pool | Benchmark | Cases | Source |
|---|---|---|---|
| protein | BAliBASE 3 (BB + BBS) | 385 | drive5.com/bench (Edgar), common FASTA format |
| protein | PREFAB 4 | 1,666 | same archive (reference = the structural pair) |
| protein | OXBench | 395 | same archive |
| protein | SABRE | 411 | same archive |
| protein | BaliFam100 | 59 | github.com/rcedgar/balifam (CC0) |
| nucleotide | BRAliBASE II | 574 | registry (RNA; data set 1 families tRNA, 5S rRNA, U5, group II intron, SRP; data set 2 pairs) |
| nucleotide | MDSA (balibase, oxbench, smart `_all`) | 299 | dna.cs.byu.edu/mdsas (coding DNA aligned by protein structure) |
| nucleotide | Rfam 15.1 seed alignments | 2,352 | ftp.ebi.ac.uk/pub/databases/Rfam (CC0); ALL eligible families, gaps as given |

Eligibility rules for Rfam are those of step 0 (`recipe/cases_nucleotide.py`); the draft's 600-family sample is
replaced by all 2,352 eligible families because the sample was high-identity (median 0.80) and the cost is small.
Left out: HomFam, BaliFam1000/10000, BRAliBASE 2.1, all simulated data.

Duplicates. Exact duplicates would count one test twice. Found by content (identical multiset of input sequences and
of reference rows): BB40017 (= BB40016), the `twi_` copy of each of 12 identical SABRE `sup_/twi_` pairs, and 15 PREFAB
cases whose reference pair is the same pair in reverse order (e.g. 1aoiD_1aoiB / 1aoiB_1aoiD). 28 dropped in total
(`cases_protein.final.dropped.tsv`): the protein pool is 2,916 cases, not the 2,931 of the earlier decision text.
One pair of Rfam families with identical references (RF00929, RF04288) is kept; both sit in one group, one fold.
Two protein pairs share their input sequences but have different reference alignments and are kept as distinct tests:
BB30019 / BB40034 (identical input file) and PREFAB 1aoa_1bkrA / 1bkrA_1aoa (same two sequences, different reference
files); each pair sits in one group and one fold.

Data differences from the old pipeline that the paper states: the drive5 copy of BAliBASE differs from the original
in two cases (BB40002 has an extra input sequence; BB20002 has a different ABL_DROME sequence); on the other 216 of 218
cases the new scorers agree with the old pipeline to 1e-6.

## 4. Scoring

Per case, from the alignment Kalign produces and the reference: SP (the scorer's recall: fraction of reference residue
pairs recovered, on the reference's core/scored columns), TC (fraction of core columns reproduced), precision
(fraction of aligned pairs that are in the reference), F1 (reported, not an objective). Notes the paper must contain:
PREFAB references are pairs, so TC equals SP on all PREFAB cases; precision is computed on all columns, SP and TC on
core columns (`compare_detailed`); BRAliBASE is scored on all columns (no annotation); both scorers verify that every
input sequence is present and unchanged and raise otherwise (the nucleotide scorer raises on `integrity_bad > 0`).
Scorer versions (`protein-v1`, `nucleotide-v1`) are part of the cache key. A failed alignment (error or timeout) scores
0 on SP, TC, precision and F1 out of fold and keeps its CPU, so failures cannot flatter a preset.

## 5. Grouping and folds (leakage control)

Two cases are in one group when they are linked; groups go whole into one fold.

* Protein: linked when a REFERENCE sequence of one case is identical to, or contained in, a reference sequence of the
  other (semi-global identity 1.0 over the shorter sequence, plus an exact substring test so that an ambiguous residue
  X / B / Z, which the identity kernel never counts as identical, cannot hide a containment; sequences under 30 aa by
  exact hash), union-find over all benchmarks. (The substring test was added in review: without it 4 PREFAB pairs
  sharing a PDB chain, e.g. 1ubpC, sat in different folds; it merges 6 group pairs.) 543 groups; largest 438 cases
  (15.0%, below the 20% limit); 267 singletons; 189 groups span more than one benchmark.
  The draft's rule (50% identity over input and reference sequences) put 1,511 of 2,944 cases (51.3%) into one group
  through real shared PDB chains, so grouped 5-fold CV was impossible; the sensitivity of the rule was measured (0.9
  identity: 26.4%; 0.99: 18.2%; identical-or-contained: 14.8%; exact identity only: 3.5%). Because the rule leaves
  homologs across folds, every out-of-fold result is reported split by the covariate `homolog_in_train` (the case has a
  link of at least 50% identity, step-0 rule, to a case in another fold): 921 of 2,916 cases (31.6%; by benchmark
  bali3 143/385, balifam 22/59, ox 140/395, prefab4 544/1,666, sabre 72/411).
* Nucleotide: step-0 rule (input and reference sequences; length-dependent identity: shorter sequence under 50 nt needs
  0.95, under 100 nt 0.80, otherwise 0.5; k-mer prefilter k = 7, 4 shared) plus Rfam clan links (`#=GF CL`; 241 links,
  174 of which merged two components). 1,976 groups; largest 552 cases (17.1%). The covariate `homolog_in_train` is 0
  for every nucleotide case by construction (the groups use the rule that defines it; 74.9 M candidate pairs verified,
  no cross-fold link), so it carries no information. The paper claims leakage control at the level of sequence links
  and Rfam clans only: the prefilter misses about 12% of true 60-80% identity links and relatedness beyond clan level is
  not measured.
* Folds: 5, stratified by benchmark and difficulty, seed 1234; greedy assignment of whole groups, largest first, to the
  fold that reduces the squared deviation of benchmark x difficulty cell counts the most. Difficulty = identity bin of
  the reference core: protein <20 / 20-30 / 30-40 / >=40%; nucleotide <50 / 50-70 / 70-85 / >=85%. BRAliBASE: the four
  family groups with at least 20 cases are in four different folds.
* Fold sizes. Protein 585 / 584 / 583 / 582 / 582. Nucleotide 814 / 636 / 622 / 613 / 540: the 552-case group is in fold 0,
  and fold 0 holds most of the hard nucleotide cases (112 under 50% identity against 26-31 elsewhere). These are
  properties of the data; per-fold numbers are descriptive and the pooled numbers are the result.

## 6. Training core, screening subset, cap

Training uses the four other folds only; test folds are never reduced. All draws rank cases by
`sha256(1234 | purpose | case_id)`: deterministic, stable across folds, asserted for every core (no test case in a
training core, cap respected, screen inside core, guard disjoint from screen, no protein case over 10 s CPU).

* Cap: at most 5 cases per group (protein) / 25 (nucleotide) in the training core. When a group exceeds the cap, the
  cases of the scarcer benchmarks survive first (balifam, bali3, ox, sabre, prefab4), then hash order. A purely random
  cap thinned BaliFam to 17 and bali3 to 103 cases in fold 0 and was replaced.
* Protein core per fold: 800 cases + a fixed guard set of 40 cases from the >=40% identity bin (840 in every fold). All
  BaliFam that survive the cap; the rest split over bali3, prefab4, ox and sabre by water-filling; inside each benchmark
  the identity bins (<20, 20-30, 30-40, >=40%) are weighted 1 : 1 : 1 : 0.5. Cases with more than 10 s CPU for 3.7.0
  accurate (3-5 per fold, 30-worker timing) are excluded from cores (still tested). The guard set is drawn first and
  kept out of the screening subset, so it is unseen during the search.
  **The cap binds**: only 956-1,106 cases per fold survive it (1,255 for the final fit), so the equal split over four
  benchmarks is not reachable. Realised cores (folds 0-4): bali3 124-173, prefab4 257-344, ox 103-131, sabre 221-246,
  balifam 37-45 cases; final fit: bali3 175, prefab4 242, ox 141, sabre 231, balifam 51. "All BaliFam" therefore means all that survive the cap. The paper says so.
* Nucleotide core: all training cases after the cap (fold 0: 2,201; folds 1-4: 1,931-1,948; final fit 2,488).
* Screening subset (two-speed search): 300 cases per fold and pool, a stratified subset of the core (proportional over
  benchmark x difficulty cells, at least 3 per non-empty cell; nucleotide: split equally over the three benchmarks
  where possible; BRAliBASE has only 71-96 core cases after the cap and gives the rest to MDSA and Rfam).
  300, not the draft's 150: the standard error of a paired SP difference at 150 cases is 0.008-0.009, above the 0.005
  decision margin.

## 7. Search

* Two stages per fit, NSGA-III (pymoo 0.6.1.6), population 56, Das-Dennis p = 5 (56 reference directions for 4
  objectives), seed 1234. Stage 1: single-run presets; stage 2: ensembles of 3 or 5 runs seeded from the stage-1 front.
  ("Up to 5" includes 3: the 3-run nucleotide default of 3.7.0 must be representable. With 5 runs the protein decode
  equals `space.py`.)
* Search spaces unchanged: protein `tools/refit/space.py` (62 variables), nucleotide `tools/nested_nuc/nuc_space.py`
  (42 variables). Their sha256 are pinned in `recipe/presets_370.py` (`SPACE_SHA256`, so they are under Freeze 1) and a
  fit refuses to start with any other file. Limit: `nuc_space.to_worker` fixes vsm_mode 0, BPM guide tree, no fill-in, matching the 3.7.0
  nucleotide presets, so the nucleotide search cannot represent anything else.
* Seeds: every 3.7.0 single run seeds stage 1 (protein fast, default and the five accurate runs; nucleotide fast, recall
  and the runs of default and accurate). The stage-2 member pool is the CPU-filtered stage-1 front plus these runs,
  so the 3.7.0 ensembles are exact points; stage 2 starts from the 3.7.0 ensembles plus one ensemble of the five front
  members with the best mean z. Each 3.7.0 preset decodes to exactly its own parameter key (tested).
* Objectives, both pools: maximise SP, TC, precision (per-case means, averaged within each benchmark, benchmarks weighted
  equally); minimise total single-thread CPU. F1 is reported. NSGA-III sees the same monotone transform as the
  hypervolume (log CPU, normalised quality), so Pareto dominance is unchanged. Failed individuals are infeasible
  (constraint), not given fake worst values.
* Generations: stage 1 min 20, max 100; stage 2 min 10, max 140 (v1: 40; see Amendment v1.1); stop (from the minimum on) when the hypervolume of the
  stage archive improves by less than 0.1% over 5 generations. Hypervolume: objectives normalised against the 3.7.0
  presets on the screening subset; quality axes `(1 - q) / (1 - (q_fast - 0.10))` (from fast minus 0.10 up to 1); CPU on a
  log scale from 0.5x fast to 2x default (stage 1) or 2x accurate (stage 2); reference point (1.1, 1.1, 1.1, 1.1).
* Two speeds: every generation is scored on the screening subset; the final front of each stage is re-scored on the
  full core before selection. Only front points with screening CPU at most 1.5 x 1.2 x the 3.7.0 accurate CPU are
  re-scored (the only ones that can fit any budget); the same cut defines the stage-2 member pool.
* Runaway guard (search candidates only): an individual's evaluation stops at its first failed case, or when its CPU
  exceeds 3x the 3.7.0 accurate CPU on the same cases (the individual is then infeasible; it would lie outside every
  budget and outside the hypervolume box anyway). Never applied to the 3.7.0 references or out-of-fold tests.
  Timeouts: 600 s per case during fitting, 1,800 s in out-of-fold tests and re-timing. A timed-out case stores
  the error "timeout" and its process group is killed. A "timeout" cached under the shorter fitting limit (possibly
  by another fold's search) is never reused by an out-of-fold test or the re-timing: it is re-run under 1,800 s and
  replaced.
* Leak guard: the fit accepts training ids only and stops (`LeakError`) if any is a test case; the evaluator refuses any
  case outside its allowed set (the core during a fit, the test fold during testing); `nested test` verifies that the
  test cases are exactly fold k and disjoint from the fit's `cases_core.txt` (hash in `run.json`).

## 8. Selection (identical in every fold and in the final fit)

Candidates: the core-scored stage-1 and stage-2 fronts plus the 3.7.0 presets. Take the non-dominated set on (SP, TC,
precision, CPU). Within budget = 1.2 x the CPU of the 3.7.0 preset of that mode on the same core, measured through
the SAME parameter-set worker path (protein recall uses accurate; the C presets protein recall and accurate are
identical), pick the maximum mean z-score of SP, TC and precision (z over that set); recall additionally needs SP >= the
accurate choice. Ties go to lower CPU, then key. If nothing is within budget the mode falls back to 3.7.0 (flagged).

Known property, stated in the paper: equal-weight z of SP, TC and precision can pick a point that gains precision at the
cost of SP/TC; the decision rule then rejects it out of fold.

## 9. Evaluation of the 3.7.0 presets and of the external aligners

* 3.7.0 presets are decoded parameter dicts taken from `lib/src/aln_wrap.c`; `tests/test_presets_370.py` re-parses the C
  source so a later edit fails the test. The worker path reproduces the CLI byte for byte on 290 of 290 alignments.
  Only individuals with tree noise > 0 are seed-sensitive; every 3.7.0 preset has noise 0.
* External aligners, for comparison only (never tuned), same cases, same scorers, single thread: MAFFT 7.526 L-INS-i
  (`--localpair --maxiterate 1000`), MUSCLE 5.3 (`-align` up to 200 sequences, `-super5` above; the pools have 3 cases
  over 200), Clustal Omega 1.2.4 (Debian package binary: clustal.org refuses scripted downloads; 1.2.4 is the latest
  release). Timeout 3 h per case.

## 10. CPU protocol

CPU is single-thread user + system seconds of the aligner call. CPU timed at 30 concurrent workers is inflated by 5-30%
relative to 4 workers (step 0), and the cache keeps the first measurement. Therefore:
* all fitting evaluations, the budgets and the 3.7.0 references inside a fit are measured at the same machine-wide
  concurrency (30 aligner slots shared by all fit and test processes);
* the finalists and the 3.7.0 references are re-timed at 4 workers (`recipe.decide retime` holds the other 26 slots) in
  a separate cache; reported CPU and the final budget check use this re-timing. The 1.2x budget check at 4 workers is
  binding;
* the batched worker's CPU differed by up to 8% from the per-case worker in one 8-worker check; this is accepted noise
  at the scale of the budget and is why the final check is repeated at 4 workers.
Alignments and quality scores are deterministic; CPU is not, so a fit is reproducible from the cache but not
bit-identical from an empty cache.

## 11. Reporting, decision, honesty safeguards

Pooled out-of-fold tables (`recipe.pool_results`, `recipe.report`): n, errors, mean and 95% BCa interval (10,000
resamples, seed 1234; benchmark-weighted statistics resample within benchmark, jackknife acceleration over cases) of SP,
TC, precision, F1, CPU; overall, per benchmark, per `homolog_in_train`, per difficulty bin; paired differences new -
3.7.0; optimism (final-fit in-sample minus out-of-fold); selected preset per fold; search effort per fit.

Decision per pool and mode: the final-fit preset ships iff (1) out of fold, pooled, benchmark-weighted (the fit's own
estimand), dSP >= -0.005 and dTC >= -0.005 and at least one >= +0.005, d = new - 3.7.0; (2) its CPU on the final-fit core
re-timed at 4 workers is <= 1.2 x that of the 3.7.0 preset of the mode; (3) it differs from 3.7.0. Otherwise 3.7.0
stays. The result is `decision.json` and the generated `presets.json`.

Safeguards: this plan, the code and the seeds are tagged before any result (Freeze 1); the out-of-fold rows are hashed
before the final fit starts and rechecked after it (Freeze 2); `run_all.sh` refuses to run if the code or the data
tables differ from what was frozen; non-registered (smoke) settings are recorded in `run.json` and
`nested`, `pool_results` and `decide` refuse them without `--allow-smoke` (a reduced core, `--limit-core`, also makes a
fit non-registered); test folds are complete in every
fold and every failure counts as 0; every setting run, every deviation and every stopped run is listed in the
supplement. Statistical limit: the bootstrap resamples cases, not groups; cases of one group are correlated, so the
intervals can be somewhat too narrow (the decision rule uses point estimates and has no interval criterion). Limits
stated up front: the 5-per-group cap leaves protein cores unequal across benchmarks; the protein
homolog covariate shows that 31.6% of protein test cases have a >=50% identity homolog in training; nucleotide
leakage control is at sequence-link and clan level; Rfam seed alignments are curated alignments, not
structure-derived references (scored as given, gaps included); the nucleotide folds are unequal in size.

## 12. Running it and keeping it honest

`presets/run_all.sh` runs, detached and resumable, with one line per finished step in `WS/runs/full/progress.log`:

| step | content |
|---|---|
| cv_protein, cv_nucleotide | five fits per pool (all at once, sharing 30 slots), then the out-of-fold test of each fold on ALL its cases |
| pool_oof | pooled out-of-fold tables; FREEZE 2 record (`freeze2.json`); `--stop-after pool_oof` stops here for tagging |
| final_protein, final_nucleotide | the identical fit on the whole pool, in-sample scores |
| externals | MAFFT L-INS-i, MUSCLE 5, Clustal Omega on all cases of both pools |
| retime | 4-worker re-timing of finalists and 3.7.0 presets |
| pool_full | pooled tables with 4-worker CPU, externals, optimism; verifies Freeze 2 |
| decide | decision per pool x mode, `presets.json` |
| report | `report/report.md` and `report/tables/*.tsv` |

Freeze check: `python -m recipe.freeze record --rev <tagged commit> --out FREEZE1.json` (laptop, inside the repository)
writes the sha256 of every frozen file (`presets/recipe/*.py`, `tests/*.py`, `resources/*`, `run_all.sh`,
`docs/PLAN.md`) at that commit; `run_all.sh --freeze FREEZE1.json` checks the workstation copy against it
(different, missing or extra file: refused) at every start and resume, copies the manifest into the run directory, and
refuses to resume under a different manifest. The data tables the fit reads (case manifests, folds, cores, features,
groups) are hashed at the first start and verified at every resume.

## 13. Step-0 findings that shaped this plan

1. Protein grouping by 50% identity put 51.3% of the pool in one group (real shared PDB chains): replaced by identical-or-
   contained reference sequences plus the `homolog_in_train` covariate.
2. 13 exact-duplicate protein cases found in step 0, 15 more reversed PREFAB pairs found when the manifests were rebuilt.
3. The old Rfam sample (600 families) was high identity (median 0.80): all 2,352 eligible families, plus clan links in the
   grouping.
4. The paired SP difference has a standard error of 0.008-0.009 at 150 cases, above the 0.005 margin: screening subset 300.
5. CPU timed at 30 workers is 5-30% above 4-worker CPU, comparable to the 1.2x budget: all references and candidates at
   one concurrency, finalists re-timed at 4 workers.
6. The cache key lacked input/reference hashes and a scorer version: key is now (build, parameters, case, sha256 of
   input, sha256 of reference, scorer version), schema 2; 17,668 step-0 rows migrated after re-hashing every file.
7. The nucleotide scorer did not raise on integrity failures: it does now. Evaluate calls had no timeout: now every path has
   one and kills the process group.
8. A process per case costs about 0.2 s: the batched long-lived worker is used (identical alignments on 200 of 200
   cases; CPU within about 8%).
9. Only `fast` was shown byte-identical between CLI and wheel in step 0: all seven distinct 3.7.0 presets are now decoded
   from `aln_wrap.c` and verified (290 of 290).
10. 3.7.0 protein default is no better than fast on the pool (SP 0.681 vs 0.681, TC 0.589 vs 0.594) at 5.5x the CPU;
    accurate makes alignments more than 2x the reference width on 117 nucleotide and 64 protein cases (up to 8.7x). These
    are aligner behaviours, not pipeline errors; they are the reason precision is an objective.
11. Step-0 baselines (case means; CPU fast / default / accurate): protein SP 0.681 / 0.681 / 0.698, TC 0.594 / 0.589 / 0.610,
    CPU 75 / 417 / 2,419 s; nucleotide SP 0.836 / 0.830 / 0.824, TC 0.689 / 0.688 / 0.683, CPU 15 / 51 / 107 s (on the
    step-0 pools, before Rfam was completed).
12. The mean-z selection rule can pick a precision-for-SP trade; the out-of-fold decision rule is the safeguard.

## 14. Differences from the draft plan

| Draft | Implemented | Reason |
|---|---|---|
| protein link at 50% identity over input + reference | identical-or-contained reference sequences + covariate | 51% giant group |
| 2,931 protein cases (13 duplicates) | 2,916 (28 duplicates) | 15 reversed PREFAB pairs |
| Rfam sample of 600 / "a sample" | all 2,352 eligible families | high-identity sample; negligible cost |
| screening subset ~150 | 300 | standard error vs margin |
| equal split over benchmarks | cap-limited, scarce benchmarks first | the cap binds |
| plain NSGA-III on CPU | same Pareto relation on log CPU / normalised quality, constraint for failures | niching over 5-30x CPU spacing |
| ensembles of 5 | ensembles of 3 or 5 | the 3-run nucleotide default |
| out-of-fold results only | + final-fit in-sample (optimism), Freeze 2 hash | the plan's Freeze 2 made mechanical |
| snakemake workflow | `run_all.sh` + resumable Python steps | the engine needs a machine-wide slot pool and checkpoints |

## 15. Expected cost

See `presets/docs/build/integrate.md` (measured smoke, per-generation timing and the wall-time estimate).

## 16. After the decision (unchanged from the draft)

Phase R: tag the last commit that contains `experiments/redesign/` and `benchmarks/`
(`archive/redesign-experiment-<sha>`), inventory dependants, port what is needed, remove both in one commit; the
presets compiled into Kalign come from the generated `presets.json`. The 3.7.0 release is decided after the results.


## Amendment v1.1 (2026-10-07; written after the v1 results, before any v1.1 computation)

**Reason.** In v1 (tag `exp/refit-v1-results`) stage 2 (the ensemble search that produces default, recall and
accurate) stopped at its 40-generation cap in every fit (5 folds x 2 pools and both final fits) while the
hypervolume was still improving; stage 1 converged by the stop rule in 10 of 12 fits (protein fold 3 and the
protein final fit hit the 100-generation cap). The v1 ensemble presets are therefore search-limited. v1 found no
meaningful overfitting (optimism within +0.006 on SP and TC), so more search is not expected to hurt
generalisation, but this is measured, not assumed.

**Change.** The stage-2 generation cap is raised from 40 to 140. Nothing else changes: same pools, folds, cores,
screening subsets, seeds, NSGA-III settings, stop rule (hypervolume < 0.1% over 5 generations, from the minimum on),
selection rule, budgets, CPU protocol, decision rule and reporting.

**How.** Every v1 fit (10 folds + 2 final fits) is continued from its saved stage-2 state (`algorithm.pkl`) in a
copy of the v1 run directory (`WS/runs/v1_1`; v1 in `WS/runs/full` is left untouched). Stage 1 is not re-run (its
checkpoint is complete and its stop rule is unchanged). The fit's resume check accepts exactly this one settings
change and records it in each fit's `run.json` (`amendments`). Then, exactly as in v1: stage-2 front re-scoring on
the full core, selection, out-of-fold tests on all test cases, pooling, Freeze 2 (before the final fits), final-fit
continuation, external aligners (cached), re-timing at 4 workers, decision, report. Derived v1 outputs
(selection.json, oof_rows.tsv, stage-2 front_core.json, pooled results, re-timing) are removed from the copy so
that every v1.1 number is recomputed. All alignments go through the same cache, so v1 evaluations are reused, not
repeated.

**Reporting.** v1 and v1.1 are both reported, side by side: out-of-fold quality vs 3.7.0, optimism, convergence,
decision. v1.1 supersedes v1 for the shipped presets only if v1.1 is run to completion under this amendment.


## Addendum A1 (2026-10-07; after the v1.1 results): Kalign 3.6.1 and 3.7.0 command-line comparison

**Purpose.** Descriptive only: no fit, selection or decision uses it. It adds the previous release, Kalign 3.6.1
(commit c007ad0, the 3.6.1 version bump), to the comparison, and measures the single-thread speed of 3.7.0
against 3.6.1 under one harness.

**How.** `recipe.externals` gained two tool families, `kalign361_<mode>` and `kalign370_<mode>` (mode = fast, default,
recall, accurate): the command-line program of each build (`kalign -i IN -o OUT --mode <mode> -n 1`), run on every
case of both pools with the same scorer, cache-key scheme and CPU measure (user + sys from `wait4`) as the external
aligners. Builds: 3.6.1 built from `git archive c007ad0` on the workstation (Release); 3.7.0 = the fitting build
`WS/kalign/960a4f6`. Both at 8 concurrent workers, the two builds alternating per mode and pool, separate cache
(`WS/runs/addendum_a1/cache_a1.sqlite`). Pooled tables come from `recipe.pool_results` with these rows as external
arms (`WS/runs/addendum_a1/refit/<pool>/results/`), reading the v1.1 out-of-fold rows unchanged.

**Reading.** 3.6.1 and 3.7.0 presets were not re-fitted per fold, and both were partly tuned earlier on BAliBASE,
BRAliBASE and MDSA, which are in the pools; their scores are therefore in-sample for part of the pool. Differences
between 3.6.1 and 3.7.0 combine the 3.7.0 bug fixes, the exact speed-ups and the 3.7.0 preset changes.


## Amendment v2 (2026-10-07; written after the v1.1 results and the Addendum A1 comparison, before any v2 number)

**Reason (decision of the investigator, 2026-10-07).** The paper can only motivate presets that come out of this
procedure. v1 and v1.1 kept a 3.7.0 preset unless a new one beat it by a margin, so four of seven shipped presets
would have come from the earlier September procedure. In addition, the v1 selection rule (equal-weight z of SP, TC
and precision within 1.2 x the 3.7.0 CPU) tied the modes to the CPU of the old presets and chose precision-heavy
nucleotide presets. Inspecting the final-fit fronts (`presets/results/v1_1/figures/fig_mode_choice_*`) showed that
modes should sit at the knees of the quality curve, not in CPU slices, and that the curves support two modes per
pool (the middle protein step gains evenly, the nucleotide ensemble gains < 0.002 per CPU doubling) and no precision
mode (precision above the accurate mode is reached only by aligning fewer residues; SP falls to 0.53-0.65).

**Change.** Selection and decision only; fits, fronts, folds, cores and the cache are unchanged (v1.1 fronts are
re-used as they are). Code: `recipe/select_v2.py`, run by `presets/run_v2.sh`.
1. Quality q = (SP + TC) / 2 on the fit's training core (benchmark-weighted as before). Precision is reported, not
   selected on.
2. A front point p is acceptable within a CPU limit L when no front point o with cpu(p) < cpu(o) <= L gains
   q(o) - q(p) >= 0.002 x log2(cpu(o) / cpu(p)) (more CPU is paid only for at least 0.002 per doubling).
3. accurate = the cheapest acceptable point with L = infinity.
4. fast = the cheapest acceptable point with L = (total CPU on the same training core of the fastest of MAFFT
   L-INS-i, MUSCLE 5 and Clustal Omega) / 5 (fast stays at least 5 times faster than every other aligner). If no
   point lies within L, the cheapest point is used and flagged.
5. default = recall = accurate. The command line keeps all four names.
6. No decision rule: the selected presets ship. Out of fold they are compared with 3.7.0, 3.6.1 and the external
   aligners and every difference is reported, including losses.
7. Same procedure in every fold and in the final fit; out-of-fold rows hashed (Freeze 2) before the final
   selection is scored in sample; optimism = in-sample minus out of fold, as before. CPU is reported at the fitting
   concurrency (no separate 4-worker re-timing, which only served the dropped decision rule).
The constants 0.002 per doubling and 5 x are fixed here before any v2 out-of-fold number exists. They were chosen
after looking at the final-fit fronts (in-sample), which the paper states.

**Reporting.** v1, v1.1 and v2 are all reported. v2 defines the shipped presets. Added: out-of-fold SP and TC per
BAliBASE 3 reference category (RV11, RV12, RV20, RV30, RV40, RV50) for the v2 modes and the external aligners.

**Deviation v2-1 (2026-10-07 16:39).** `run_v2.sh` stopped at the first pooling step with `KeyError: 'budget_cpu'`:
`recipe.pool_results` read the v1 selection field `budget_cpu`, which v2 selections do not have (v2 has no CPU
budget). Fix: the field is read with `.get` (reported empty for v2). No selection, out-of-fold row or Freeze 2 hash is
affected: selection, the out-of-fold tests and Freeze 2 had completed before the error; the remaining steps
(in-sample scores, pooling, Freeze 2 check) are run from tag `exp/refit-v2-fix1`.


## Addendum A2 (2026-10-07; after the v2 results): Kalign 3.1 and 3.3.5

**Purpose.** Descriptive only, like Addendum A1: show where the shipped presets stand against Kalign 3 as
published (3.1, the release of the Kalign 3 paper, October 2019) and as widely installed (3.3.5, November 2022).
Prompted by a laptop check (`presets/results/v2/repro_laptop/`) in which a native arm64 build of 3.1 without AVX2
gave implausibly low scores; this addendum runs both releases as x86_64 builds with their SIMD paths.

**How.** `recipe.externals` tools `kalign31` (`kalign -i IN -o OUT`, defaults; the release has no thread option)
and `kalign335` (`kalign -i IN -o OUT -n 1`, defaults), every case of both pools, 8 workers, separate cache;
pooled with the v2 out-of-fold rows by `recipe.pool_results` (`presets/run_a2.sh`). Builds: 3.1 compiled directly
from its sources with `-O2 -mavx2` and `HAVE_AVX2` (its autotools configuration is not used); 3.3.5 with its CMake
defaults (SSE4.1/AVX/AVX2 on, OpenMP off). Neither release was tuned on this pool's folds; both default
parameters were set long before, partly on BAliBASE.

Deviation A2-1 (2026-10-07 18:15): the 3.1 build failed to link with the workstation's gcc ("multiple definition
of logsum_lookup": gcc >= 10 defaults to -fno-common, which the 2019 sources predate). Fix: `-fcommon` added to the
3.1 compile line. No other change.
