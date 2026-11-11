# Step 0 / case features and grouping

Status 2026-10-05. Owner files: `recipe/features.py`, `recipe/grouping.py`, `tests/test_grouping.py` (it also holds
two small feature tests), this note. No fitting, no folds, nothing committed. Outputs on the workstation
(WS = `~/tmp/kalign-refit`, not in git): `data/features_{protein,nucleotide}.tsv`, `data/groups_{protein,nucleotide}.tsv`,
plus evidence files `data/groups_*.edges.tsv` (the links that merged two components) and `data/groups_*.stats.json`.
Sensitivity runs and scratch scripts: `scratch/grouping/`.

Run (workstation; numba needed for grouping, it is in `WS/venv-grouping`, a new venv with numpy, numba, pytest; the
kalign build venv has no numba):

    cd ~/tmp/kalign-refit/presets
    ../venv-grouping/bin/python -m recipe.features --pool protein          # 0.4 s for 2,944 cases (24 procs)
    ../venv-grouping/bin/python -m recipe.grouping --pool protein --validate 3000   # about 75 s, 24 threads
    ../venv-grouping/bin/python -m recipe.grouping --pool nucleotide --validate 3000 # about 3 min

Without numba the code falls back to pure python (tests pass on the laptop in 2 s; real runs would take hours).
Tests: 11 pass on laptop and workstation.

## 1. Features (from the reference alignment only)

One row per case, columns listed in the module docstring: `n_seq_input` (manifest), `n_seq`, `ref_len`,
`len_mean/sd/min/max` (ungapped length of the reference sequences), `mean_pid` (mean over sequence pairs of
identical / both-present columns), `mean_pid_core` (same on the scoring-mask columns: `core_columns`/`pair_core` for
protein, all columns for nucleotides), `core_cols`, `core_frac`, `gap_frac`, internal gap runs (`n_gap_runs`,
`gap_runs_per_seq`, `gap_run_len_mean/max`, `internal_gap_frac`) and terminal extensions (`lead_mean/max`,
`trail_mean/max`, `term_ext_frac`). Every non-letter is a gap. A gap run is a maximal stretch of gap cells in one
row strictly inside the row's first and last residue (so a column gap that is "all other sequences have a
residue" and one that is a block insertion are both counted; a distinction by column pattern is not made).
Case marks (lower case = not assessed) only enter through the mask of `mean_pid_core`/`core_*`.

Decisions: PREFAB references are pairs, so `n_seq`, gap structure and identity describe the pair, while
`n_seq_input` (about 50) is kept as the cost driver. No feature reads an aligner output.

Validation: `mean_pid_core` against the drive5 `info/pctids.txt` for the three benchmarks that have it (prefab4
1,681, ox 395, sabre 423 cases): correlation 1.000, mean absolute difference 0.02-0.03 percentage points. So the
identity definition equals drive5's. No other feature was checked against an external value (only hand-built
unit tests).

Medians [min-max] per benchmark (`mean_pid` = all columns, fraction):

| benchmark | cases | n_seq | ref_len | mean length | mean_pid | gap_frac | gap runs/seq | terminal frac |
|---|---|---|---|---|---|---|---|---|
| bali3 | 386 | 18 [4-142] | 504 [52-8481] | 309 | 0.29 [0.12-0.60] | 0.30 | 11.9 | 0.07 |
| balifam | 59 | 11 [4-142] | 134 [30-545] | 118 | 0.32 [0.13-0.50] | 0.16 | 4.6 | 0.00 |
| ox | 395 | 4 [3-122] | 130 [43-598] | 119 | 0.52 [0.07-1.00] | 0.07 | 1.8 | 0.03 |
| prefab4 | 1681 | 2 | 242 [54-1016] | 210 | 0.20 [0.04-0.97] | 0.10 | 3.5 | 0.03 |
| sabre | 423 | 4 [3-25] | 202 [38-807] | 144 | 0.16 [0.05-0.50] | 0.24 | 5.0 | 0.11 |
| bralibase | 574 | 5 [2-5] | 118 [63-377] | 110 | 0.66 [0.10-0.96] | 0.05 | 2.6 | 0.00 |
| mdsa | 299 | 5 [3-15] | 300 [66-399] | 266 | 0.54 [0.28-1.00] | 0.07 | 1.3 | 0.01 |
| rfam | 600 | 12 [4-200] | 108 [32-961] | 91 | 0.79 [0.42-1.00] | 0.08 | 3.3 | 0.00 |

No case has an undefined `mean_pid`. Note for the stratification: the BAliBASE `mean_pid` (all columns) is low
partly because the lower-case, non-core regions are included; `mean_pid_core` is the one comparable with drive5.

## 2. Grouping

Rule as in the plan: two cases are linked if a sequence of one is identical to, contained in, or >= 50 % identical
(over the shorter length) to a sequence of the other; groups = connected components, all benchmarks of a pool
together. Sequences compared = every sequence in the case's unaligned input plus its reference (upper case, T=U).

Implementation (`recipe/grouping.py`): unique sequences (hash) -> cases sharing a sequence are merged (union-find);
then a k-mer inverted index (protein k=5, >= 3 shared distinct k-mers; nucleotide k=7, >= 4) proposes pairs; pairs
whose cases are already in one component are dropped; the rest are verified best-first in rounds of 50,000 by a
numba semi-global alignment (free end gaps, linear gap cost; identity = identical residues on the best path /
length of the shorter sequence; X/N never match), and components are refreshed between rounds. Containment is the
identity = 1.0 case of the same machinery (no separate substring search was needed). Sequences shorter than 30
take part only in the exact-hash rule.
Timing: protein 63,310 unique sequences, 1.8 M candidates, 0.8 M verified, 75 s; nucleotide 15,109 unique
sequences, 9.4 M candidates, 9.4 M verified, about 3 min.

### 2.1 Things I had to decide (and measured)

1. **Nucleotide identity needs a length-dependent threshold.** First attempt (plain 0.5 for both pools) produced a
   nucleotide giant component with 1,183 cases mixing BRAliBASE, Rfam and MDSA coding DNA, including 144 links
   between Rfam families and MDSA coding sequences, which cannot be homology. Cause: with a 4-letter alphabet a
   short sequence reaches 50 % identity at some offset of a long one by chance. A shuffle null (2,000 random
   sequences of the pool, 2 M pairs): protein never exceeds 0.38 (also not at 30-40 residues, 0 pairs >= 0.5), so
   protein keeps 0.5; nucleotides: max identity 0.90 at shorter length 30-40, 0.83 at 40-60, 0.75 at 60-100, 0.50 at
   >= 100 (no pair >= 0.5 at >= 100). Rule now (nucleotides only, `IDENT_SCHEDULE`): shorter length < 50 needs
   0.95, < 100 needs 0.80, otherwise 0.50. This is my decision, not in the plan; it is conservative for short RNAs
   such as tRNA (about 75 nt) which now link only at >= 0.80 or by exact sequence.
2. **Residual false links.** Negative control: MDSA (coding DNA) and Rfam/BRAliBASE (RNA) should share nothing.
   After the fix 4 identity links between MDSA and Rfam families remain (all in groups of size 2-3, they do not
   matter for fold building but they show the rule is not clean). I did not inspect these 4 after the fix. In the
   earlier (pre-fix) run I looked at three similar pairs: one was a repetitive 40-nt sequence (low complexity), the
   other two (75 % over 76 nt, 51 % over 78 nt) I could not classify as chance or homology.
3. **Prefilter sensitivity** (all-pairs alignment of 3,000 random usable sequences, fraction of true links the
   k-mer prefilter proposed): protein 95.1 % overall (91.7 % for identity 0.5-0.6, 99.5 % for 0.6-0.8, 100 % above);
   nucleotide k=7/T=4: 88.2 % (153 true links in the sample, so imprecise). Alternatives tried for nucleotides
   (k=8,T=4: 71 %, k=7,T=3: 91 % at 16 M candidates); the number of components changed only between 769 and 754, so
   a few missed pairs rarely matter (bridges between components are redundant in big families). Missed links can
   only make groups smaller, never larger.
4. Identical results on re-runs (md5 of both `groups_*.tsv` equal, protein and nucleotide).

### 2.2 Result: protein pool (2,944 cases): ONE GROUP HOLDS HALF THE POOL

354 groups: 153 singletons, 66 pairs, 43 triples, 70 of size 4-12, a tail up to 55 and **G0001 with 1,511 cases (51 %)**:
bali3 199 of 386, prefab4 925 of 1681, ox 230 of 395, sabre 127 of 423, balifam 30 of 59. The next groups have
55, 52, 43, 40, 36 ... cases. Size buckets: 1: 153, 2: 66, 3-5: 75, 6-20: 45, 21-100: 14, > 100: 1.

Benchmarks spanning groups: 148 groups contain more than one benchmark. **384 of 386 BAliBASE-3 cases and all 59
BaliFam100 cases and 389 of 395 OXBench cases are in a multi-benchmark group.** BaliFam joins BAliBASE as expected
(58 bali3-balifam merging edges; only 5 of the 59 BaliFam cases share an identical sequence with another benchmark,
the rest are linked by >= 50 % identity, e.g. `balifam.PF00625.100` is identical to `bali3.BB30021`).
Merging edges by benchmark pair: bali3-prefab4 296, prefab4-sabre 262, ox-prefab4 120, bali3-sabre 69, bali3-ox 67,
bali3-balifam 58, ox-sabre 34, balifam-prefab4 5.

Sensitivity of the giant group (largest group / groups / singletons; runs in `scratch/grouping/protein_*`):

| compared sequences | identity | largest group | groups | singletons |
|---|---|---|---|---|
| input + reference (official) | 0.5 | 1511 | 354 | 153 |
| input + reference | 0.7 | 1476 | 366 | 158 |
| reference sequences only | 0.5 | 1346 | 414 | 186 |
| reference sequences only | 0.7 | 1120 | 430 | 191 |
| reference sequences only | 0.9 | 778 | 447 | 197 |
| input + reference | exact sequence only | 408 | 832 | 491 |

So the giant group is not an artefact of the PREFAB full-set inputs (50 homologs per case) or of the 50 % threshold:
even near-identical links (0.9, reference sequences only) leave 778 cases in one group, and the benchmarks
(BAliBASE, OXBench, SABRE, PREFAB, BaliFam) are drawn from overlapping sets of families and structures. **A 5-fold
grouped split cannot be built with this rule as written**: one fold would hold >= 51 % of the pool (>= 1,511 of
2,944), the others are tiny, and test and training distributions would differ completely per fold. Options for
Freeze 1 (a decision for the owner, I did not apply any): (a) split G0001 by a finer rule (higher identity, or
sub-clustering of the group by e.g. average-linkage on the same alignments) and accept near-neighbours across
folds, (b) leave G0001 in a single always-train block and test on the other groups, (c) keep the rule but fold
by benchmark-defined families (SCOP family for PREFAB/ox/sabre/balifam, BAliBASE reference) if such metadata is
added. The `edges.tsv` files list which links hold the group together if one wants to inspect bridges
(not done).

### 2.3 Result: nucleotide pool (1,473 cases)

756 groups: 692 singletons, 25 pairs, 16 triples, 4 of size 4, 7 of size 5, ..., and five large ones
which are the five BRAliBASE-I families: 191, 109, 95 (93 BRAliBASE + 2 Rfam), 93 (92 + 1 Rfam), 81 cases. Size
buckets: 1: 692, 2: 25, 3-5: 27, 6-20: 7, 21-100: 3, > 100: 2. MDSA cases are mostly isolated (165 of 169 groups
touched are MDSA-only; largest MDSA groups 13, 11, 9; MDSA-balibase/oxbench/smart overlaps are small).

BRAliBASE joins Rfam only in two groups (5 exact-sequence edges, 3 Rfam cases): the 600-case Rfam sample is a random
draw over families, so it contains only a few of the BRAliBASE families (SRP/tRNA/rRNA...) , not a design fault of the
rule. The expectation "BRAliBASE joins Rfam families" holds for exactly those few. Practical consequence: the
BRAliBASE cases form five groups of 81-191 cases each, so with 5 folds a fold will hold one whole family
(rRNA, SRP, tRNA, U5 or group-II intron): fold = RNA family, i.e. a leave-one-family-out evaluation, which is much
harsher than a random split and uneven in size (81 to 191). `data-set2` pairwise cases (116 of 574) are inside these
groups. Alternatively the Rfam sample could include the BRAliBASE families on purpose.

## 3. Open issues / honesty

* The protein giant group (section 2.2) is the main outcome and blocks fold assignment as planned.
* The nucleotide length schedule is my choice (calibrated on a shuffle null, not on real unrelated sequences); 4
  suspicious MDSA-Rfam links remain.
* Prefilter recall is 88 % (nucleotide) and 95 % (protein) on a random sample; this was not measured on the full
  data, and the nucleotide sample has only 153 true links.
* Single linkage: a few chance or multi-domain links can merge families (probably part of the protein giant group;
  not investigated: which edges are bridges).
* Features: `mean_pid` for bali3 uses all columns including non-core ones; only the core-restricted version was
  validated (against drive5).
* Fold assignment, stratification bins and the final choice of `IDENT_SCHEDULE` are not done (Freeze 1).
* Do not `rsync --delete` `presets/` to the workstation while other agents work there. I used plain `rsync -a` of
  my own files.
