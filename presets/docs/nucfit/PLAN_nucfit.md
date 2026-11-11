# Kalign 4 nucleotide settings (DNA/RNA fast and accurate): the nucfit plan

Status: DRAFT 2026-11-02; decisions N1-N3 taken (k4/DECISIONS.md D12). Nothing in it has been run. It is frozen (commit, tag `nucfit-freeze`, sha256 of
this file and of `presets/` in the tag message) before any fit starts; every number marked "measured before the
freeze" is filled in by the pre-freeze checks of section 10, which score no case the selection uses. Changes after
the freeze are dated amendments that say what had been seen when they were made.

This plan mirrors the protein procedure (`presets/docs/fastfit/PLAN_fastfit.md` + `DECISION_fastfit_selection.md`
for fast; `fastfit/PLAN_ensemble.md` + `ensconf/PLAN_ensconf.md` for accurate). Where it differs, the difference is
listed in section 9 with the reason.

Questions:
- **Fast.** Can a single-run nucleotide setting, found by nested cross-validation with the guide tree searched
  alongside the alignment parameters, beat the pre-release nucleotide fast setting (key 4ff7e4d1) on held-out
  small and large nucleotide alignments within a CPU limit?
- **Accurate.** Can a five-member ensemble whose member 1 is the new nucleotide fast setting give a complete
  alignment at least as accurate as fast (SP, TC) whose column confidence lowers the pair Brier score below fast's,
  on held-out alignments?

## 1. The overfitting test, stated first

Every fit is run six times: on the training part of folds 0-4 (five fold fits) and on all fit data (the final fit,
which ships). Three estimates are reported side by side for every measure, each with a 95 % bootstrap interval:

| Estimate | What it is | Biased? |
|---|---|---|
| In-sample (IS) | the final fit's pick scored on all fit data (the data it was chosen on) | optimistic |
| Out of fold (OOF) | each fold fit's pick scored on that fold's held-out cases, pooled over folds 0-4 | unbiased for the procedure |
| Lockbox (LB) | the final fit's pick scored once on data no fit or check has read | unbiased for the shipped setting |

Overfitting criteria (pre-set; a failure is reported as overfitting and blocks shipping, section 7):
- **O1, optimism.** For each primary measure, IS minus OOF (the same cases: OOF uses the fold fits, IS the final fit,
  both on all fit data) must be at most 0.01, or within the paired bootstrap interval of zero. The difference is
  reported for every measure whether or not it passes.
- **O2, lockbox agreement.** The lockbox difference to the reference (pick minus reference) must lie inside the OOF
  95 % interval of the same difference. A lockbox result below that interval means the final fit did worse on new
  data than the procedure predicts.
- **O3, stability.** The six picks are reported side by side (all parameters). Reported, not a gate: the spread of
  each parameter across fold fits and the OOF score of every fold pick on every other fold's test cases.

Note: the final fit cannot be scored out of fold, because it was trained on everything. That is why both the OOF
estimate (of the procedure) and a lockbox (of the shipped setting) are needed.

## 2. Data

All inputs are real reference alignments unless marked as simulated. Manifest code: `recipe/cases_nucleotide.py`
(step 0, 1,473 small cases), extended for the large tier and the lockbox.

**Small pool:** the full gm1 nucleotide pool of the pre-release fit (3,225 cases: BRAliBASE II 574, MDSA 299,
Rfam 15.1 seed 2,352), unchanged (decision D12; protein also used its full gm1 pool). The pre-release nucleotide
presets were fitted on this pool, so the reference has seen every small case, including lockbox ones; this favours
the reference and makes the gate conservative.

**Large tier** (decision N1 = D12; code `recipe/nucfit_data.py`):
- Rfam seed families with more than 200 sequences (the families the small pool's size filter left out; 201-1,613
  sequences), at full size, small-pool filters otherwise (ragged, > 2,000 columns, empty sequence, gapless);
  benchmark `rfam_large`;
- simulated DNA (benchmark `sim_dna`): INDELible 1.03, HKY (kappa 2, base frequencies 0.3 0.2 0.2 0.3), gamma 4
  categories alpha 0.5, root length 500, indel model POW 1.7 max 10; birth-death trees as downstream_v2; grid
  1,000 and 10,000 sequences x depth 0.25 / 0.5 / 1.0 x indel rate 0.005 / 0.01 / 0.02 (the same rates at both
  sizes, so every root-to-tip lineage sees the same expected number of indels and only the number of sequences
  changes); 5 replicates at 1,000 and 3 at 10,000 per cell, the last replicate of every cell is lockbox.
  Amendment 2026-11-03, before any alignment was run on simulated data: the first grid (rates 0.02 / 0.05 / 0.10)
  gave true alignments up to 95,700 columns for 500-base sequences (190x) at 10,000 sequences, beyond the
  workstation's memory in INDELible and far from any real nucleotide family (Rfam large seed alignments are
  1.2-2x the sequence length). The two higher rates were replaced by two lower ones; the 0.02 sets are unchanged
  (same seeds, reproduced byte for byte). Peak INDELible memory at the widest remaining set: 23 GB;
- Q_large = mean of three category means (Rfam large, simulations at 1,000, simulations at 10,000), as the protein
  per-size means.
SILVA and Rfam full alignments are not used: they are made by aligners or covariance models rather than curated by
hand, so scoring against them partly rewards agreeing with another program.

**Grouping and leak guards.** BRAliBASE was built from Rfam seed alignments, and Rfam families of one clan are
related. Groups are built before folds:
1. sequences are matched across all cases (exact, or contained, T = U, case ignored), the rule of the protein
   grouping (`recipe/grouping.py`);
2. cases sharing a sequence, the same Rfam family or the same Rfam clan form one group;
3. MDSA cases are grouped by their source alignment.
Folds: 5, by seeded hash order of groups (seed 1234), stratified by benchmark and size. Largest group share is
reported; if one group exceeds 15 % of the pool it is split by family and the split is reported (as protein step 0).

**Lockbox (drawn first, before anything else is run on the new data):** 15 % of the groups (seed 1234, stratified
by benchmark, size and identity), plus 20 % of the large-tier families and 20 % of the simulation grid cells. The
lockbox file is hashed and committed at the freeze; the evaluator refuses lockbox ids until the final pick is
written and hashed. It is scored once, for both settings (fast and accurate), declared here.

**Observed at prepare (2026-11-03, before any fit; fit data only):** with the large Rfam families added, the
identity and clan links join one group holding 39 % of all cases (17 % in the gm1 refit), which contains every
BRAliBASE case and most large Rfam families. Folds therefore read: fold 0 = all 574 BRAliBASE cases and 72 of 92
non-lockbox large Rfam families; folds 1-4 = Rfam, MDSA and simulations only; the lockbox holds no BRAliBASE case
(the group is above the lockbox cap). Kept as is: the grouping is leak-free and the rule was fixed before the data
were seen; the cost is that fold 0's fit never trains on BRAliBASE (its BRAliBASE scores are a held-out-benchmark
test) and that the lockbox cannot test BRAliBASE. Reported with every OOF number.

**Amendment A4 (2026-11-04, user decision after the first fast evaluation; k4/DECISIONS.md D14): folds without an
oversized group.** The first complete fast fit failed its gate (OOF SP_small -0.0036 [-0.0067, -0.0008], TC_small
-0.0011; O1 failed) and the diagnosis traced most of the small-input loss to the group above, not to the search:
BRAliBASE's whole OOF score came from fold 0's pick, which had trained on no BRAliBASE; the cap of 25 cases per
group left 15 BRAliBASE cases in the training cores while Q_small weights BRAliBASE one third; and 72 of 92 Rfam
large families sat in fold 0. The split rule written above ("split by family") was not applied at prepare; it is
applied now, in a form that cannot leak (`recipe/nucfit_refold.py`):
1. every case-level link inside a group over the limits is computed (the grouping keeps only the links that merged
   components): shared sequence, the length-dependent identity rule on input + reference sequences with every
   prefilter candidate verified, every pair of Rfam cases (seed or large) in one clan;
2. while a component holds more than 10 % of the fit cases (small + Rfam large) or more than 35 % of one
   benchmark's fit cases, one case of it is removed: Rfam seed and MDSA first, then Rfam large, then BRAliBASE;
   inside a tier the case with the most links in the component (ties: case id);
3. removed cases get role `excluded` (no core, no test fold, no lockbox); the lockbox, all other groups and the
   simulations are unchanged; folds are re-assigned with `recipe/folds.assign_folds` (unchanged rule); cores by the
   rules of `recipe/nucfit_data.py` (unchanged). 5 folds (D14).
Result: 271 cases excluded (Rfam seed 232 of 2,063, MDSA 26 of 259, Rfam large 13 of 92, BRAliBASE 0); the group
splits into the five BRAliBASE families (191, 111, 93, 87, 85 cases) and one cluster of 30 Rfam large families;
largest group 7 % of the fit cases; every fold tests 82-191 BRAliBASE cases and trains on 106-113 (final fit 138).
Rfam large per test fold: 11 / 12 / 13 / 30 / 13. Checks (`nucfit_refold check`): the data checks of
`nucfit_data check`, no excluded case in any core, lockbox file unchanged, BRAliBASE in every test fold, and an
exhaustive search finding no exact, identity or clan link between folds or between fit data and the lockbox.
All fits are re-run from scratch on the new folds (cached alignments reused); the folds-v1 runs are archived
(runs/archive/folds_v1_2026-11-04) and reported as such, never pooled with the new ones.

## 3. Scoring

`recipe/scoring_nucleotide.py`: every reference column scored (`max_gap_frac = -1`, no core mask; validated
574/574 against the old BRAliBASE pipeline). SP (= recall), precision, TC, width ratio. q = (SP + TC) / 2.
Amendment 2026-11-03 (before any fit; scorer version nucleotide-v2): simulations at 10,000 sequences are scored on
1,000 of their sequences drawn with a seed from the simulation id (`ref_subset.fasta`, the true alignment of those
rows with all-gap columns removed); the test alignment of all 10,000 is projected onto those rows (score_rule
`projected`), exactly as the extHomFam reference rows inside a protein family at 10,000. Reason: scoring the full
10,000-row truth took more than 11 minutes per alignment (pairs grow with N squared) against 15 s to align it.
Known-answer tests: the full true alignment scores SP = TC = 1 against its subset; on a 1,000-sequence simulation
the projected rule with every row kept gives the all-columns scores exactly.
Small accuracy: per-case q averaged within benchmark, benchmarks weighted equally (decision N2: equal weight gives
Rfam one third). Large accuracy: per-family q, mean over families and sizes.

## 4. Fast fit (mirrors PLAN_fastfit.md)

Fixed: 1 run, no pairwise guide distance, 32 k-means anchors, 1 thread, first 1,024 residues for tree distances
(as protein). One build for candidates and reference: 9137658 (the release build; its nucleotide presets are the
reference). Code `recipe/nucfast.py` (recipe.fastfit unchanged except as listed in its docstring).

Search space: the pinned nucleotide stage-1 space of the pre-release fit (`nested_nuc/nuc_space.py`, sha256
869c22f5; adapter `recipe.fit.NucleotideSpace`), plus the protein fast fit's guide-tree variables and the gap model:

| Variable | Range |
|---|---|
| gap open, extend, terminal extend | [1, 20], [0.2, 5], [0.05, 3] (the pre-release values lie inside) |
| scoring | Kalign's nucleotide matrices: 1 PAM, 20 PAM, 200 PAM |
| VSM amax, refinement, sequence weights, consistency (+ weight), realign | the pinned space |
| gap model | 0 (Kalign gap rule, the pre-release nucleotide setting) or 1 (member gap model, protein fast) |
| guide tree | method (k-means or hybrid), handover size 50-5,000 (log), tree distance (4 options), alpha, k-means distance, realign tree |

Seeds: the pre-release nucleotide fast and accurate settings x guide trees {4.0.0, E21, protein fast} x gap model
{0, 1} (12), the rest of the population of 56 random.

Objectives on the same cases for every candidate: Q_small, Q_large, CPU ratio R to the pre-release nucleotide
fast setting (geometric mean of small and large). Constraints: no failed case (600 s); CPU_small <= 1.5x and
CPU_large <= 2x the reference.

Data per fit (recipe.core rules, as the pre-release nucleotide fit): small core = every training case, at most 25 per
group; small screen = 300 cases, equal share per benchmark, proportional over identity bins; large core = every
large training case; large screen = 24 cases, Rfam large and each simulation size proportional to their training
counts (at least 2 each). The search runs on the screen, re-scoring and selection on the core.

Amendment A1, 2026-11-03 (after generation 1 of the first launch, before any selection; only costs were looked at,
no accuracy): the simulations at 10,000 sequences are removed from the search screen (as the protein ensemble fit,
AMENDMENT_ensemble_1 (a)). Measured over the first 50 minutes of generation 2: they took 98 % of the search CPU
(154,281 of 157,874 s; median 209 s per alignment; all 63 timeouts), Rfam large 0.3 s and simulations at 1,000
8 s median; one generation would have taken about 5 h, the fast fit at least 4 days. They stay in the large core,
so the re-scoring of the front, Q_large, the CPU constraint and the selection still include them. The search
screen is now small screen + Rfam large + simulations at 1,000 (Q_large on the screen = mean of those two category
means). Everything else in this section is unchanged. The first launch is archived under
runs/archive/nucfast_preA1_2026-11-03; cached alignments are reused.

Amendment A3, 2026-11-03 (user decision: "relax the CPU limits ... better than being quicker with patchy
results"; after generation 2 of the A1 launch, before any selection; only feasibility counts were looked at): at
generation 2 only 0-1 of 56 candidates per fit met the CPU constraints, so speed, not accuracy, decided the search.
New values for the nucleotide fast fit: feasible = no failed case and CPU <= 10x the reference on small and on
large cases (was 1.5x, 2x); a candidate is stopped beyond 15x the reference CPU (was 3x); per-alignment limit
1,800 s (was 600 s; as the accurate fit), cached 600 s timeouts are re-run at 1,800 s; CPU axis of the hypervolume
log R over [0.25, 10] (was [0.25, 3]). The CPU objective R and the knee rule (lowest R within 1.96 SE of the best
Q_large) are unchanged, so the pick still prefers the cheapest of the most accurate points. The A1 launch is
archived under runs/archive/nucfast_A1_2026-11-03.

Optimiser: NSGA-III (pymoo 0.6.1.6), population 56, 3-objective Das-Dennis directions, seed 1234, initial
population seeded with the pre-release nucleotide fast and accurate parameters and the protein fast guide tree;
stop at >= 20 generations and hypervolume gain < 0.1 % over 5, or 100 generations. Re-scoring of up to 40 front
points on the full training core.

Selection (the knee rule of DECISION_fastfit_selection.md, written in from the start, identical in every fit):
eligible = feasible, small SP and TC >= the reference; b = highest Q_large; keep points within 1.96 paired-bootstrap
SE of b; pick the lowest R.

Gate (pick minus reference; OOF pooled over folds 0-4, 95 % bootstrap over groups and families, 10,000 resamples):
1. OOF small SP: lower bound >= -0.005; OOF small TC: estimate >= 0;
2. OOF Q_large: lower bound > 0;
3. overfitting criteria O1 and O2 (section 1);
4. lockbox Q_small and Q_large above the reference; every large lockbox input completes.
Otherwise the pre-release nucleotide fast setting ships unchanged, and that is reported.

## 5. Accurate fit (mirrors the shipped protein design: ensemble fit, then confidence step; both pre-specified)

Fixed: K = 5 members; member 1 = the new nucleotide fast setting, frozen; every member uses member 1's guide-tree
algorithm with single-linkage tree noise (the protein accurate design, 9137658); POAR consensus with confidence;
1 thread; 1,800 s per case; same 5,000-sequence limit as protein accurate.

Step 1, ensemble fit (`recipe/nucens.py`, as protein `ensfit.py`): NSGA-III over members 2-5 (alignment
parameters of section 4 per member, member noise [0, 0.5], guide-tree variables) and the ensemble options;
objectives precision P_small, P_large and CPU ratio R to the fast pick (as protein), constraint SP >= fast pick -
0.005 (small, large); knee rule on P_large (bootstrap stratified by large category). Member 1 = the same fit's
fast pick (fold k uses the fold-k pick; protein used the final pick in every fold). Large tier of the accurate fit:
Rfam large + simulations at 1,000 (accurate refuses inputs above 5,000 sequences, as protein). Seeds: members 2-5 =
the pre-release nucleotide accurate setting (a single run) x 4, and accurate / fast / accurate / fast.
Step 2, confidence step (`recipe/nucconf.py`, as protein `ensfit_conf.py`): on the step-1 pick of each fit, every
member takes member 1's guide tree; a grid over tree noise sigma {0, 0.05, 0.1, 0.2, 0.3, 0.5} x min_support
{1, 2, 3}, scored as below on the fit's core.

Step 2 score: pair Brier score B_small, B_large (definition of PLAN_ensconf.md section 4: unaligned reference
pairs cost 1, so leaving residues out cannot lower B). Eligible: no failed case; SP and TC >= the fast pick on small
and large (same cases); no large case with SP more than 0.2 below the fast pick. No CPU limit (as protein ensconf;
CPU is reported).

Selection: lowest J = (B_small + B_large) / 2 among eligible points, ties smaller sigma then smaller min_support
(the ensconf rule), on the fit's core. No eligible point: that fit has no accurate pick (counts as fast out of fold).

Gate (pick minus the new fast setting, OOF): SP and TC estimate >= 0 and lower bound >= -0.005 (small, large);
B upper bound < 0 (small, large); O1 and O2; lockbox SP and TC not below fast - 0.005 and B below fast; every input
completes. Otherwise the pre-release nucleotide accurate setting (key bd03cdb7) ships unchanged, and that is
reported.

Secondary outcomes (reported, not a gate): AUROC and calibration of column confidence; the pick at min_support 1-5.

## 6. Comparators (reported on OOF and lockbox cases; never tuned)

MAFFT L-INS-i and X-INS-i (X-INS-i uses predicted RNA structure; RNA cases only), MUSCLE 5, Clustal Omega,
Kalign 3.1, the pre-release nucleotide presets, and FAMSA if its nucleotide mode runs on the pool (checked before
the freeze). Downstream comparisons in section 8.

## 7. Shipping rule

A setting ships only if its gate passes, including the overfitting criteria. Both fits report IS, OOF and lockbox
numbers for every measure, the six picks' parameters, every evaluated point, and a repeat of the final fit with seed
4321 (reported, not a gate). The C presets are regenerated by `presets/release/emit_c_presets.py`; identity checks
(release CLI = fit worker, byte for byte, on all fit and lockbox inputs) run before any benchmark.

**Reporting rule A2: real and simulated data reported separately** (added 2026-11-03, before any selection, any
OOF evaluation or any lockbox look; user decision). The large tier mixes real alignments (Rfam large) and
simulations (INDELible at 1,000 and 10,000), and simulated data can reward settings that fit the simulator's
assumptions rather than real sequences. Therefore, for both settings (fast, accurate), in IS, OOF and lockbox:

1. Every accuracy measure (SP, TC, Q, precision where defined) and its pick-minus-reference difference with the 95 %
   bootstrap interval is reported in three columns: **real** = the small pool (BRAliBASE, MDSA, Rfam seed;
   benchmark-weighted, as Q_small) and Rfam large, each shown on its own; **simulated** = mean of the two
   simulation size categories, with 1,000 and 10,000 also shown on their own; and the combined Q_large of the gate.
   Real and simulated are never pooled into one number without the separate columns beside it.
2. Pre-specified check (reported, not a gate): the pick's OOF difference to the reference on **Rfam large alone**.
   If its 95 % interval lies entirely below 0 while Q_large passes the gate, the report says plainly that the
   large-tier gain comes from simulated data and does not hold on real families. The same statement is made for
   the lockbox.
3. Figures that show large-tier accuracy show real and simulated in separate panels or colours.

The gate (section 4, section 5) is unchanged.

## 8. Downstream applications (benchmark after the settings are fixed; the downstream_v2 design)

| Application | Test | Truth | Mirrors protein |
|---|---|---|---|
| RNA consensus secondary structure | RNAalifold on each alignment; base-pair F1 / MCC | Rfam SS_cons, BRAliBASE reference structure | new (RNA-specific) |
| Structural homology search | Infernal cmbuild on a training subset of each family, cmsearch on held-out members + shuffled and unrelated decoys; sensitivity, sensitivity at first false positive | family membership | HMMER block |
| Covariation | R-scape on each alignment: number of significantly covarying base pairs that are true pairs, and covariation on non-pairs | reference structure | new |
| Phylogeny | IQ-TREE (GTR+G) on simulated DNA; normalised RF to the true tree | simulated tree | trees block |
| Positive selection on coding DNA aligned as DNA | FUBAR / HyPhy on simulated codon data aligned directly as nucleotides (and as translate-align-backtranslate, for contrast) | simulated selected sites | codon block |
| Confidence | calibration and AUROC of column confidence on simulated DNA and reference sets | true / reference alignment | calibration block |

Not proposed: read mapping or genome alignment (not multiple sequence alignment tasks).

## 9. Differences from the protein procedure, and why

1. **The overfitting test (section 1) is a written gate.** For protein, the OOF gate and (fast only) the lockbox
   were the protection against overfitting; IS minus OOF was not a registered criterion. Here it is.
2. **Both accurate steps are pre-specified.** Protein accurate went through an ensemble fit and then a confidence
   re-fit (with amendments) that was motivated by out-of-fold results (stated in PLAN_ensconf.md as not blind).
   Here the same two steps are written down before either runs, and no OOF result exists before the fit.
3. **The knee rule is in the plan from the start** (for protein it replaced the registered rule after launch, as a
   dated deviation).
4. **One lockbox, declared for both settings.** For protein the lockbox was declared for fast and then for accurate,
   and the shipped accurate setting was not gated on it.
5. **Large tier.** No nucleotide counterpart of extHomFam with curated references exists; the large tier relies on
   large Rfam seed families and simulations (decision N1).

## 10. Before the freeze (fit data only, never the lockbox)

- Identity: the fit worker reproduces the release CLI for the pre-release nucleotide fast and accurate presets byte
  for byte; every reference scores SP = TC = 1 against itself; groups contain no shared sequence across folds.
- Cost: CPU of the reference settings on the small pool, 3 large families and 2 simulated 10,000 inputs ->
  projected fit wall time.
- Comparator availability (FAMSA nucleotide mode, MAFFT X-INS-i).

## 11. Decisions (taken 2026-11-02, k4/DECISIONS.md D12)

- **N1, large tier:** Rfam seed families of more than 200 sequences plus simulated DNA at 1,000 and 10,000; no SILVA
  or Rfam full alignments.
- **N2, weighting:** benchmarks weighted equally, as protein.
- **N3:** one setting for DNA and RNA (Kalign has one nucleotide biotype); results reported separately.
