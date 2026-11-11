# Step 0 / protein case manifest and scorer

Status 2026-10-05. Files: `recipe/cases_protein.py`, `recipe/scoring_protein.py`,
`tests/test_scoring_protein.py`, this note, plus one validation helper `recipe/_validate_bali3.py`
(not part of the pipeline). No fitting, nothing committed. Manifest and validation outputs live on the
workstation (WS = `~/tmp/kalign-refit`): `data/cases_protein.tsv`, `data/cases_protein.dropped.tsv`,
`data/cases_protein.stats.json`, normalised inputs in `data/protein/in/<benchmark>/<id>.fa`,
validation material in `scratch/protein/` (alignments, per-case tables, the extracted BAliBASE release).

Rebuild: `cd ~/tmp/kalign-refit/presets && <venv>/bin/python -m recipe.cases_protein` (deterministic: a second
build gives byte-identical TSV, tested on a synthetic tree).

## 1. Manifest: 2,944 cases, 0 dropped

Columns are the 11 shared ones (`case_id, pool, benchmark, subset, unaligned, reference, score_rule, n_seq,
ref_len, sha256_unaligned, sha256_reference`). `case_id = <benchmark>.<id>` (e.g. `bali3.BB11001`,
`balifam.PF00009.100`). `n_seq` = sequences in the unaligned input (= the cost driver), `ref_len` = reference
alignment length. `reference` points to the untouched raw file in `data/raw` (sha256 recorded); `unaligned` points
to a normalised copy (upper-case letters only, name = first header token), because the drive5 `in/` files are
extension-less and the BaliFam/MDSA-style inputs are not guaranteed clean. Gap stripping was not needed for any
protein input (drive5 and BaliFam `in/` files are already gap-free); the copy only normalises.

| benchmark | source | cases | subsets (cases) | n_seq input min/med/max | ref seqs min/med/max | ref_len min/med/max | gapless refs |
|---|---|---|---|---|---|---|---|
| bali3 | drive5 `bali3` (BB 218 + BBS 168) | 386 | BB-RV11 38, BB-RV12 44, BB-RV20 41, BB-RV30 30, BB-RV40 49, BB-RV50 16; BBS-RV11 38, BBS-RV12 44, BBS-RV20 41, BBS-RV30 30, BBS-RV50 15 | 4 / 18 / 142 | same as input | 52 / 504 / 8481 | 0 |
| prefab4 | drive5 `prefab4` (input = full set) | 1681 | pid<20 805, pid20-30 534, pid30-40 124, pid>=40 218 | 2 / 50 / 50 | 2 / 2 / 2 (pair) | 54 / 242 / 1016 | 7 |
| ox | drive5 `ox` | 395 | pid>=40 291, pid30-40 56, pid20-30 32, pid<20 16 | 3 / 4 / 122 | same | 43 / 130 / 598 | 16 |
| sabre | drive5 `sabre` | 423 | pid<20 177, pid20-30 153, pid30-40 65, pid>=40 28 | 3 / 4 / 25 | same | 38 / 202 / 807 | 0 |
| balifam | BaliFam100 (CC0) | 59 | balifam100 59 | 104 / 111 / 242 | 4 / 11 / 142 | 30 / 134 / 545 | 0 |

(Per-subset distributions: `data/cases_protein.stats.json`.) Subsets for prefab4/ox/sabre are percent-identity
bins of the core columns taken from the drive5 `info/pctids.txt`; BBS-RV40 does not exist (BBS has no RV40 cases).
Note the `ox` identity bins are very unbalanced (291 of 395 in `pid>=40`), and `ox` has 16 gapless references.

### Decisions
* **Which variants.** bali3 (BB and BBS; they overlap in sequences, the grouping rule handles that; BB/BBS are
  told apart in `subset`, so the old 218-case DEV set is exactly the `BB-*` subsets), prefab4, ox, sabre,
  balifam100. Not used: `oxx` (106 references have no upper-case residue), the MUSTANG variants, `bali3pdb`,
  `prefab4ref` (identical references to prefab4 with a 2-sequence input; the full-set input is the standard
  PREFAB protocol and what `score_prefab_case` did), `bali2dna*`.
* **PREFAB input = the full set (median 50 sequences).** This is the expensive choice (the 1,681 cases carry most of
  the pool's alignment cost); using `prefab4ref` would make them pairwise alignments. The scorer works for either,
  because it projects onto the reference names. The pool owner may switch by changing `DRIVE5_SETS`.
* **Gapless references are kept** (ox 16, prefab4 7): they have scorable residue pairs and test over-gapping;
  `validate_case` drops a case only if there is no core column, no assessed residue pair, a reference with <2
  sequences, duplicate names, ragged rows, or reference residues that differ from the input. None of the 2,944
  triggered any of these (`cases_protein.dropped.tsv` is empty). The old BRAliBASE rule (skip gapless) does not apply
  to proteins here.
* **BaliFam / BAliBASE overlap** is not resolved here (grouping step).

## 2. Scoring rules (`scoring_protein.score_case(case_row, aln_path)`)

All cases: project the test alignment onto the reference sequences (others removed, all-gap columns dropped),
rewrite the reference with `-` gaps (`.`/`-`/`!` all non-letters become `-`), call
`kalign.compare_detailed(ref, projected, column_mask=M)`.

| score_rule | benchmarks | mask M (on reference columns) | old code it ports |
|---|---|---|---|
| `core_columns` | bali3, ox, sabre, balifam | >= 1 upper-case residue and no lower-case residue (mixed-case columns are not core, count reported as `mixed_cols`) | `score_balifam_case` in experiments/redesign/tools/score_confirmation.py |
| `pair_core` | prefab4 | both residues of the reference pair upper case | `score_prefab_case` |

Returned: `recall` (SP), `precision`, `f1`, `tc`, `width_ratio` (= width of the projected test alignment /
reference width; for full-reference cases identical to the old `bench_balibase.py` ratio), `n_seq`, plus extras
`core_cols, mixed_cols, ref_pairs, test_pairs`. Integrity is enforced: every input sequence must be in the test
alignment with unchanged residues (case-insensitive), no extra sequences, no duplicate names, reference residues
reproduced by the projection, non-empty mask; otherwise `ValueError` (the evaluate harness stores it as an error row,
a case is never silently scored on a damaged alignment).

Differences from the old `score_*_case` functions: they only checked residues for the reference sequences and
reported `integrity_bad`; this version raises on any input/output mismatch (stricter), and `core_columns` is
also applied to BAliBASE 3 and OXBench and SABRE (the older code used it only for BaliFam).

## 3. Validation

**Scorer self-test (all 2,944 cases).** Each reference, padded with the extra input sequences in private
insertion columns, scored against itself gives recall = precision = F1 = TC = 1.0 for every case
(`test_every_reference_scores_one_against_itself`, 40 s on the workstation; skipped on the laptop). This checks
parsing, name handling, projection and mask handling on the whole pool, not accuracy.

**BAliBASE 3, 218 BB cases vs known numbers.** `kalign --mode M -n 1` run by `recipe.evaluate` (CLI of build
960a4f6, sha256 d0a20557..., 28 workers) on the drive5-format BB cases, scored with `scoring_protein`:

| mode | SP | precision | F1 | TC | known (old pipeline, 3.7.0, XML core blocks) |
|---|---|---|---|---|---|
| default | 0.8147 | 0.6739 | 0.7321 | 0.4761 | SP 0.814, F1 0.732, TC 0.476 |
| fast | 0.8221 | 0.6755 | 0.7351 | 0.4953 | 0.822 / 0.735 / 0.495 |

All six quoted numbers are reproduced to the stated precision (the known SP 0.814 is 0.8147 truncated rather than
rounded, or a difference below 0.001; I do not know which). Whole 218-case runs take seconds on 28 workers.

**Why drive5 core marking vs XML core blocks differ so little (explanation, measured by
`recipe/_validate_bali3.py` against the original `BAliBASE_R1-5.tar.gz` release at `~/tmp/BAliBASE_R1-5.tar.gz`,
sha256 5378cf4e...):**
* Column masks: for 217 of 218 BB cases the drive5 upper-case mask equals the XML `colsco-data` mask exactly
  (same reference columns, same core count). Upper-case marking in the drive5 files therefore is the XML core
  blocks. The one other case, BB20002, has a drive5 reference of 2036 columns versus 2228 in the MSF (192 fewer
  columns; I did not find out which, but both masks contain 16 core columns and the scores of the same alignment
  are identical to 1e-9).
* Sequence content: drive5 BB40002 has 56 sequences, the release 55 (extra `APP1_SCHPO` in drive5). Same mask, but
  recall of the same alignment differs by 0.0049 for this case (0.0031 precision), which moves the 218-case means
  by < 0.0001 (scored on the same alignments, drive5-mask vs XML-mask: SP 0.8147 vs 0.8147, F1 0.7321 vs 0.7321,
  TC 0.4761 vs 0.4761 for default; fast 0.8221 vs 0.8222, 0.7351 vs 0.7352, 0.4953 vs 0.4953). The other
  217 cases agree to 1e-9. A first version of the helper compared a 56-sequence alignment against the 55-sequence
  MSF directly and produced a spurious 0.05 for BB40002 (I did not investigate how the C comparison pairs sequences when the sets differ); the helper now projects onto the release names. This also means anyone scoring an alignment against
  an MSF must give it exactly the reference's sequence set.
* The 218 drive5 BB inputs are the same 218 case ids as the old DEV set; the alignments are not byte-identical to
  the old runs (different input file, sequence order) but scores agree.

**Smoke run on other benchmarks.** A random sample (seed 1) of 25 prefab4, 25 ox, 25 sabre, 10 BBS and 8 balifam
cases was aligned with `--mode default` and scored without error (sample means are not meaningful, n is tiny:
prefab4 SP 0.668, ox 0.828, sabre 0.518, BBS 0.895, balifam 0.859). No error rows. This was not compared with any
published number: **PREFAB, OXBench, SABRE and BaliFam scoring is validated only by the self-test and by being
a direct port of the score_prefab_case/score_balifam_case code, not against known results** (the earlier BaliFam/
PREFAB results were produced with a different kalign version and set of cases).

## 4. Open issues

1. The PREFAB input choice (full set vs pair) decides most of the CPU cost of the pool; Step 0's cost measurement
   should show it. If too expensive, `prefab4ref` is a drop-in (same references; a one-line change).
2. `pid` subsets use drive5's `pctids.txt`; the OX bins are skewed, SABRE `pid>=40` has only 28 cases. The
   grouping/stratification step may want its own difficulty bins from the references.
3. The metric on references with no gap at all (ox 16, prefab4 7) is valid but trivially easy for gap-free
   alignments; keep or drop is a Freeze-1 decision. They are tagged only through the manifest stats, not a column.
4. BAliBASE BB/BBS and BaliFam overlap: left to the grouping step.
5. `recipe/_validate_bali3.py` needs the original BAliBASE release (not registered in REGISTRY.yaml: it was found at
   `~/tmp/BAliBASE_R1-5.tar.gz` and used only for this check).
6. The validation caches (`scratch/protein/val_cache.sqlite`) are separate from the project cache on purpose;
   the project cache `WS/cache/cache.sqlite` was not written to.
