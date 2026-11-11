# Step 0 / nucleotide case manifest

Status 2026-10-05. Owner files: `recipe/cases_nucleotide.py`, `recipe/scoring_nucleotide.py`,
`tests/test_scoring_nucleotide.py`, this note. No fitting. Nothing committed.

Output (on the workstation, not in git): `WS/data/cases_nucleotide.tsv` (1,473 cases, sha256 `ebba56df...`,
rebuilding gives the same hash), `WS/data/cases_nucleotide.stats.json` (all counts below), derived files in
`WS/data/derived/nucleotide/{mdsa,rfam}/`. Rebuild: `cd WS/presets && WS/kalign/960a4f6/venv/bin/python -m recipe.cases_nucleotide`
(2.7 s; `--rfam-sample N`, `--seed S`). Tests: 14 pass on laptop and workstation.

## Manifest format

The 11 shared columns (`case_id, pool, benchmark, subset, unaligned, reference, score_rule, n_seq, ref_len,
sha256_unaligned, sha256_reference`) plus extras `seq_kind` (rna / dna), `family`, `ref_identity`, `stratum`, `note`.
`recipe.scoring.read_manifest` accepts it. Paths are absolute workstation paths. `pool` = `nucleotide`,
`score_rule` = `all_columns` for every row (see scoring). `benchmark` is `bralibase`, `mdsa` or `rfam`;
`subset` is the BRAliBASE family (or `data-set2`), the MDSA database (balibase / oxbench / smart) or `seed`.
`ref_identity` = column-pooled mean pairwise identity of the reference (identical residue pairs / residue pairs
over all columns; T = U, case ignored; exact, no pair loop). It is a property of the reference only.
`seq_kind` is for the DNA-versus-RNA preset question of the plan (BRAliBASE, Rfam = rna; MDSA = dna). Note that
MDSA is coding DNA of real genes and Rfam/BRAliBASE contain both T and U spellings.

## Counts

| Benchmark | Found | Dropped | Kept |
|---|---|---|---|
| BRAliBASE II | 599 (481 set1 + 118 set2) | 25 gapless references | **574** (set1 g2intron 92, rRNA 81, SRP 93, tRNA 83, U5 109; set2 116) |
| MDSA | 1,869 (498 + 672 + 699) | 1,129 by sequence count (not 3..15), 415 by length (> 400), 26 gapless | **299** (balibase 29, oxbench 178, smart 92) |
| Rfam 15.1 seed | 4,227 families | 1,572 n_seq < 4; 98 n_seq > 200 (skipped, not subsampled); 5 aln len > 2000; 200 gapless | 2,352 eligible, **600 sampled** |
| **Total** | | | **1,473** |

BRAliBASE (574) and MDSA (299) reproduce the old DEV counts exactly (574 and 325 - 26 = 299). The 25 BRAliBASE
cases dropped are the ones listed in `stats.json` (`rRNA` 8, `tRNA` 15, `data-set2/aln66`, `aln118`). The old rule for gapless
is "a `-` in any sequence line"; the same test is used (`.` also counts, which changes nothing here).

### MDSA filter, reproduced
Old rule (`nested_nuc.discover`): count `>` lines, 3 <= n <= 15; longest non-header line (stripped) <= 400
(the files hold one line per sequence, so this is the alignment length). Gapless dropped afterwards.
Unaligned input: the reference with `-` and `.` removed, headers and case verbatim, one line per sequence
(port of `datasets._strip_gaps_fasta`), written to `derived/nucleotide/mdsa/<db>/<case>.afa`.
Counts agree with the data note: 29 + 192 + 104 = 325 pass, 0 + 14 + 12 gapless dropped.

### Rfam sample
Parse Stockholm (handles interleaved blocks; latin-1), `.` and `-` -> `-`, upper-case, U/T as given. Filters in
this order: 4..200 sequences, equal row lengths, alignment length <= 2000, no empty sequence, at least one gap.
Drops are counted by the first rule that fires (so 200 gapless is among families that passed the others).
Identity of the eligible set: min 0.41, quartiles 0.72 / 0.80 / 0.87, max 1.0 (median 0.80 is high: seed
alignments of related sequences).
Stratification: 4 n_seq bins (4-5, 6-9, 10-24, 25-200) x 4 identity bins (<0.60, 0.60-0.75, 0.75-0.90, >=0.90);
bin edges were fixed after looking at the distribution (cell sizes 2 .. 392). Allocation: proportional to cell
size with at least min(cell size, 3) per non-empty cell, largest remainder to reach exactly 600; within a cell
`random.Random(1234).sample`, cells visited in sorted label order, so the sample is a pure function of the file,
the seed and the bins. Sample vs eligible identity quartiles: 0.72 / 0.80 / 0.88 vs 0.72 / 0.80 / 0.87 (the
minimum per cell oversamples the sparse low-identity cells: 24 of 90 for n25-200 x id<0.60). The 16 cells and
their eligible / sampled counts are in `stats.json`. The final size is a later decision: change `--rfam-sample`.
References written to `derived/nucleotide/rfam/RFxxxxx.ref.fa`, inputs `RFxxxxx.fa` (one line per sequence,
names `name/start-end` as in Rfam).

### Size distributions (min / median / max)

| Benchmark | n_seq | ref_len | ref_identity |
|---|---|---|---|
| bralibase 574 | 2 / 5 / 5 | 63 / 118 / 377 | 0.10 / 0.66 / 0.96 |
| mdsa 299 | 3 / 5 / 15 | 66 / 300 / 399 | 0.28 / 0.54 / 1.00 |
| rfam 600 | 4 / 12 / 200 | 32 / 108 / 961 | 0.43 / 0.80 / 1.00 |

Rfam is the only source of larger families (25+ sequences: 8 + 24 + 58 + 61 = 151 of the sample). BRAliBASE
and MDSA are small (5 sequences median); BRAliBASE is nearly all five-sequence cases. Nucleotide pool weighting
"benchmarks equally" would give Rfam only one third despite being the broadest set; the weighting is not decided here.

## Scoring (`recipe/scoring_nucleotide.py`)

`score_case(case_row, aln_path)` -> `kalign.compare_detailed(reference, aln, max_gap_frac=-1.0)`: every reference
column is scored, no gap-fraction filter, no column mask (the old BRAliBASE rule; applied to MDSA as the old
nested pipeline did, i.e. MDSA's upper/lower-case core marking is **not** used, and to Rfam). Returns recall (= SP),
precision, f1, tc, `width_ratio` (alignment width / reference width, as in `bench_bralibase.py`), `n_seq`, and extras
`integrity_bad` (sequences whose residues differ from the input, missing or extra; old `krcommon.integrity`),
`ref_pairs`, `test_pairs`, `common_pairs`. Any `score_rule` other than `all_columns` raises.

### Validation against the old BRAliBASE pipeline
Ran the **old** `experiments/redesign/tools/bench_bralibase.py` (unchanged, copied to `WS/scratch_nuc/oldtools`) with
`--mode fast` on the full data (574 cases; cases listed by its own glob and gapless filter), then scored the same 574
alignment files with `score_case` through the manifest. Result: **574 of 574 bit-identical** in recall, precision,
f1, tc and width_ratio (max abs difference 0.0), n_seq and integrity_bad equal; the old script's case set equals the
manifest's BRAliBASE set. Old aggregate: R 0.807, P 0.800, F1 0.803, TC 0.711, width/ref 0.98; the harness
(`recipe.evaluate`, CLI `--mode fast -n 1`) reproduces the same means (SP 0.807, P 0.800, TC 0.711, width 0.98).
Script: a temporary validation script (`nuc_validate_bra.py`, not kept) (not in the repo; it is ~25 lines: read old TSV, call `score_case`, compare).
This validates BRAliBASE only; MDSA and Rfam have no old-pipeline reference beyond the old case selection.

## Cost measurement (step-0 smoke, scratch cache `WS/scratch_nuc/cache.sqlite`, not the project cache)

`recipe.evaluate` on all 1,473 cases, 14 workers, build 960a4f6 (CLI, `-n 1`): 0 errors, `integrity_bad` = 0 on all
alignments, `max` below is a single case.

| Preset | Benchmark | SP | Prec | TC | CPU total | CPU median / max per case | width ratio |
|---|---|---|---|---|---|---|---|
| fast | bralibase | 0.807 | 0.800 | 0.711 | 2.0 s | 0.003 / 0.01 s | 0.98 |
| fast | mdsa | 0.811 | 0.812 | 0.668 | 2.0 s | 0.006 / 0.02 s | 0.99 |
| fast | rfam (600) | 0.877 | 0.877 | 0.678 | 6.0 s | 0.005 / 0.15 s | 0.98 |
| default | bralibase | 0.800 | 0.844 | 0.708 | 7.4 s | 0.012 / 0.04 s | 1.10 |
| default | mdsa | 0.775 | 0.877 | 0.635 | 7.4 s | 0.021 / 0.08 s | 1.40 |
| default | rfam (600) | 0.886 | 0.907 | 0.696 | 30.7 s | 0.024 / 1.05 s | 1.16 |

Whole pool: fast 10.1 s, default 45.5 s of single-thread CPU (wall about 1 minute each with 14 workers; the
per-process start cost dominates wall). **The nucleotide pool is very cheap**: one full evaluation of a parameter
set with the 5-run accurate ensemble should cost on the order of 5 x default = about 4 minutes of CPU (not measured:
`accurate` / `recall` were not run), against hours for the protein pool. CPU values for sub-10 ms cases are
dominated by process start and are noisy (see infra.md issue 3); the CPU cap rule of the plan is irrelevant for this pool.

## Observations and open issues

1. `default` over-gaps on MDSA (width ratio 1.40, SP 0.775 below fast 0.811; precision higher). Not investigated
   (outside step 0), but it is exactly the kind of signal a nucleotide fit should be given.
2. Rfam eligibility bias: families with fewer than 4 sequences (37%) and gapless seed alignments (200) are excluded
   by the plan's filters; the sample is high-identity (median 0.80) compared with BRAliBASE (0.66) and MDSA (0.54),
   so Rfam adds little hard-case signal. Consider oversampling the < 0.75 identity cells in the final sample (families with
   n >= 10 and identity < 0.75 are 567 of the 2,352 eligible, 24%, from the `stats.json` cells).
3. Overlap for the grouping rule: BRAliBASE was derived from Rfam seeds, so some Rfam sample families (tRNA, SRP,
   U5, group II intron, 5S rRNA) share sequences with BRAliBASE cases; `family` holds the Rfam accession and the BRAliBASE
   family name, but **sequence-level overlap is not computed here**. MDSA balibase/oxbench cases are DNA versions
   of cases in the protein pool (different pool, no leak across pools but kept in mind). Rfam families in the
   same clan (e.g. many miRNA/snoRNA families) are not independent either; clan information is in the seed file
   (`#=GF CL`) and was not used.
3a. Rfam 15.1 also has duplicate sequences across families; not checked.
4. The Rfam identity stratum uses pooled identity (over the whole alignment), not the mean of per-pair identities
   over each pair's own aligned columns; they differ slightly for gappy alignments.
5. Mixed T/U inside BRAliBASE files and IUPAC letters are passed to the aligner as given; kalign handled all 1,473
   cases without error, and the integrity check (`A-Z` residues, case-insensitive) passed on all. Scoring is a
   residue-pair (position) measure and is not expected to depend on the T/U spelling; this was not tested separately.
6. The workstation copy of `tests/` was also rsynced (no `--delete`).
