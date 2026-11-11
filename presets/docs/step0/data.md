# Step 0 / data: downloaded benchmarks

Status 2026-10-05. Everything below was downloaded on the workstation into `~/tmp/kalign-refit/data/raw/<benchmark>/`
(390 MB) by `recipe/download.py`; every file and the balifam clone are registered in `resources/REGISTRY.yaml`.
No fitting, no aligner runs. Counts below were measured by script on the unpacked files
(`recipe/_layout_stats.py` for the drive5 archive; the others by ad-hoc scripts of the same kind).

Re-run: `cd ~/tmp/kalign-refit/presets && ~/.local/bin/uv run --with pyyaml python -m recipe.download`
(`--only drive5_bench rfam bralibase mdsa balifam`, `--verify-only`). A file that exists is hashed and compared
with the registry (or the sha256 pinned in the source list for MDSA); a mismatch aborts, nothing is overwritten.
The first run's registry was overwritten once by another agent's `rsync --delete` of `presets/` to the workstation;
it was regenerated from the files on disk (retrieval times = file mtimes, i.e. the original download time).
**Do not `rsync --delete` presets/ to the workstation while agents work in parallel**: it removes files the others
created there.

## Registry summary

| Entry | Source | Version pinned | Licence statement as found |
|---|---|---|---|
| drive5_bench_tar | https://www.drive5.com/bench/bench.tar.gz, 21,414,017 B | sha256 0728653f...; server Last-Modified 2025-09-07 | none found: no licence file in the archive (only BENCH_README.txt), none on the drive5.com/bench page |
| balifam_repo | git clone https://github.com/rcedgar/balifam | commit f2a46ac425834b4acecce7e0d7cd1439090ffca2 (2021-06-29); tree sha256 of `balifam100/` a09128c1... | repo `LICENSE`: "Creative Commons Legal Code, CC0 1.0 Universal" |
| rfam_seed (+README, COPYING) | https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.seed.gz | **Rfam 15.1** (latest numbered release dir on 2026-10-05; CURRENT not used) | COPYING: "Rfam is freely available under the Creative Commons Zero ("CC0") licence." README: "freely available and in the public domain under the Creative Commons Zero licence." |
| bralibase_data_set1 / _set2 | Wayback copies (2016-04-08) of projects.binf.ku.dk/pgardner/bralibase/data-set{1,2}.tar.gz (the original host is dead); 131,201 B and 11,509 B | sha256 e93f1206..., f1bdede4... | none found in the archives; the paper is CC BY, data derived from Rfam (as in the plan) |
| mdsa_balibase / _oxbench / _smart | https://dna.cs.byu.edu/mdsas/tarred/{name}_mdsa_all.tar.gz | sha256 **identical** to the values in experiments/redesign/resources/REGISTRY.yaml (499cee4a..., be8b064b..., e4e1fc4a...); the script pins them and aborts on a difference | none found (index page gone, files live) |

BRAliBASE check: all 1,198 files of the two extracted archives match the per-file SHA-256 list of
`experiments/redesign/resources/manifests/bralibase.sha256` (1198 OK, 0 failures), i.e. the Wayback copies are
byte-identical to what the earlier experiments used.

## Common facts

- All alignments: FASTA, one sequence per record, **sequence lines may be wrapped** (balifam wrapped at 80, drive5 bali3 ref at 96 columns;  the MDSA files are one line per sequence). Parse records, do not assume one line.
- Gap characters differ: drive5 / balifam references use `.` (and a handful of `-`: bali3 4, bali3pdb 2,
  bali2dnaf also `!` x1,682 = the induced frame-shift markers); the **MUSTANG-realigned sets (`*m`: bali3pdbm, oxm,
  prefab4refm, sabrem) use `-`**; BRAliBASE, MDSA, Rfam use `-`. Normalise `.` and `-` before comparing.
- Core columns are marked by case in drive5, balifam and MDSA refs: **upper case = core/assess, lower case = not
  assessed** (drive5 README: "Upper-case regions indicate conservative regions that are intended for use in
  assessment. Lower-case regions should not be used."). A gap column inside a lower-case region is written `.`;
  inside an upper-case region it is `-` or `.` (use the column's case from the residues in it). No column
  without any upper-case letter was found in the multi-sequence refs except oxx (below).
  BRAliBASE and Rfam have no core marking: all columns are scored (`max_gap_frac=-1.0` as in the 3.5 manuscript).
- Input (unaligned) files hold the gap-free sequences in upper case; they may contain **more sequences than the ref**
  (prefab4, oxx, balifam: the ref is a subset of the input by name).

## 1. drive5 archive (protein; one tarball, 13 benchmarks, 13,931 files)

Unpacked to `data/raw/drive5_bench/unpacked/bench/bench1.0/` (note the extra `bench/` level and the `bench1.0/`
directory). Files have **no extension**. Layout per benchmark `<b>`:

```
bench1.0/<b>/in/<id>          input FASTA (unaligned)
bench1.0/<b>/ref/<id>         reference alignment, FASTA, case marks core
bench1.0/<b>/info/ids.txt     whitespace-separated case ids (= file names in in/ and ref/)
bench1.0/<b>/info/nseqs.txt   "<id>;<n>" per case (n = number of sequences in in/)
bench1.0/<b>/info/pctids.txt  "<id>;<percent identity in core columns>" per case
bench1.0/<b>/qscore/          13 qscore outputs of other aligners (not needed; useful only as a cross-check)
bench1.0/BENCH_README.txt     description, references
```
(The README says nrseqs.txt; the file is called `nseqs.txt`.) Names in `in/` and `ref/` agree except as noted.
Measured (number of cases; ref sequences min/median/max; alignment length min/median/max):

| dir | what | cases | ref seqs | aln len | notes |
|---|---|---|---|---|---|
| bali3 | BAliBASE 3, **both BB and BBS** (ids `BB11001`..., `BBS11001`...) | 386 (168 BBS, 218 BB) | 4 / 18 / 142 | 52 / 504 / 8481 | no case without upper-case; 218 BB = the old DEV subset (RV11 38, RV12 44, RV20 41, RV30 30, RV40 49, RV50 16) |
| bali3pdb | BALIS, structural subset | 180 | 3 / 5 / 14 | 68 / 444 / 1196 | |
| bali3pdbm | MU-BALIS (MUSTANG) | 120 | 3 / 4 / 12 | 77 / 437 / 1547 | gap `-` |
| ox | OXBench | 395 | 3 / 4 / 122 | 43 / 130 / 598 | 16 refs gapless (conserved) |
| oxm | MU-OXBench | 336 | 3 / 4 / 49 | 46 / 142 / 601 | gap `-`, 3 gapless |
| oxx | OXBench-X (extended) | 395 | ref 3/4/122 (same ids as ox) | same as ox | **`in/` holds 4/93/668 sequences, a superset of the ref's; `ref/` differs from ox's in content (hash differs)**; 106 refs have no upper-case residue at all (nothing to assess); 16 gapless |
| prefab4 | PREFAB 4: **input = full set** (2..50 seqs, median 50), **ref = the 2-sequence pair** | 1681 | 2 / 2 / 2 | 54 / 242 / 1016 | ids like `12asA_1atiA`; ref names are a subset of in (1681/1681); 5 refs have no lower case, 7 gapless; score the pair columns of the full alignment (extract by name) |
| prefab4ref | PREFAB-R: same references, **input = only the 2 ref sequences** | 1681 | 2 | 54 / 242 / 1016 | pairwise alignment problem, identical refs to prefab4 |
| prefab4refm | MU-PREFAB-R | 692 | 2 | 54 / 203 / 954 | gap `-` |
| sabre | SABRE (from SABmark 1.65) | 423 | 3 / 4 / 25 | 38 / 202 / 807 | |
| sabrem | MU-SABRE | 303 | 3 / 4 / 24 | 42 / 224 / 1076 | gap `-` |
| bali2dna, bali2dnaf | BAliBASE 2 reverse-translated to DNA (`f` = with frame shifts) | 141 each | 3 / 5 / 28 | 180 / 915 / 3009 | **DNA, not part of the protein pool** (not named in plan 2a); listed for completeness |

Plan 2a names BAliBASE 3, PREFAB 4, OXBench and SABRE: these are `bali3`, `prefab4` (and/or `prefab4ref`), `ox`
(optionally `oxx`), `sabre`. Which variants enter the pool is for the manifest builder; observations relevant to
that decision: the MUSTANG variants and `bali3pdb` re-use the same sequences as `bali3`/`ox`/`sabre` (overlap
for the grouping rule); `prefab4` and `prefab4ref` are the same 1,681 cases (use one); `oxx` ids equal `ox` ids.

## 2. BaliFam100 (protein)

`data/raw/balifam/repo/` is a git clone (HEAD f2a46ac4..., a single commit state recorded in the registry); only
`balifam100/` is used (`balifam1000/`, `balifam10000/` exist in the clone, not used).

```
repo/balifam100/in/<PFid>.100   input FASTA, 104..242 sequences (median 111)
repo/balifam100/ref/<PFid>.100  BAliBASE-3 reference alignment of a subset (4..142 seqs, median 11), wrapped FASTA,
                                upper case = core (26% of residue+gap characters), lower case = ignore, gap `.`
repo/balifam100/info/ids.txt    59 ids (PF00009.100, PF00018.100, ...)
repo/LICENSE (CC0), repo/README.md
```
59 cases. Ref names are a subset of the input names (checked for all 59). Score only the ref sequences'
columns of the full alignment (as for prefab4). Source of the refs is BAliBASE 3 (overlap with `bali3`).

## 3. Rfam seed (nucleotide, RNA), release 15.1

`data/raw/rfam/Rfam.seed` (decompressed, 38 MB; also `Rfam.seed.gz`, `README`, `COPYING`). One concatenated
**Stockholm 1.0** file with **4,227 families** (README: 4227; the plan's "~4,100" is outdated). Each record starts
`# STOCKHOLM 1.0` and ends `//`. Metadata `#=GF AC RF0xxxx`, `ID`, `DE`, `TP` (type, e.g. `Gene; miRNA` 1598,
`Gene; snRNA` 813, `Gene; sRNA` 765, `Cis-reg;` 351, `Gene; lncRNA` 218, `Gene; CRISPR` 65, ...), `#=GF SQ n` (number of
sequences; equals the row count in all families). Alignment rows `name/start-end<spaces>SEQUENCE` on one line
(not interleaved); `#=GC SS_cons` structure and `#=GC RF` reference line (a few have `RNA_structural_elements`,
`CORE`). Sequences: upper-case RNA (A C G U, rare N R Y M W K S H V D B), gap `-` only. **No core marking.**
Sequences per family min/median/max = 2 / 5 / 1613; alignment length 19 / 93 / 8395. Latin-1-safe (read with
`encoding="latin-1"`: some GF lines are not UTF-8). The sample (fixed seed, stratified) is for the manifest builder.

## 4. BRAliBASE II (nucleotide, RNA)

```
data-set1/data-set1/<family>/structural/aln<N>.fa   reference alignment (FASTA, wrapped at 60 columns)
data-set1/data-set1/<family>/unaligned/aln<N>.fa    the same sequences without gaps (same file name)
data-set1/data-set1/shift_rename.sh
data-set2/data-set2/structural/aln<N>.fasta         ".fasta" extension, not ".fa"
data-set2/data-set2/unaligned/aln<N>.fasta
```
under `data/raw/bralibase/unpacked/`. Families and cases: g2intron 92, rRNA 89, SRP 93, tRNA 98, U5 109 (data-set1:
481 cases, 5 sequences each) and data-set2: 118 pairwise cases (2 sequences) = **599**. Alignment length
(min/median/max): SRP 300/318/377, U5 116/127/153, g2intron 63/88/167, rRNA 116/121/128, tRNA 69/74/97, set2 70/76/92.
Gap `-`; no lower case; no core marking (all columns scored). Alphabet mixes `T` and `U` between and within
files; rRNA also contains IUPAC `M`, `S`. **25 references are gapless** (rRNA 8, tRNA 15, set2 2) and were skipped
by the earlier scoring (gives the 574 DEV cases; `max_gap_frac=-1.0` metric is undefined or trivial for them).
Unaligned and structural files have the same number of sequences in every case (checked).

## 5. MDSA (nucleotide, coding DNA aligned by protein structure)

Under `data/raw/mdsa/unpacked/<set>_mdsa_all/<set>_mdsa_all/<case>.afa` (a doubled directory level):

| set | files | seqs min/median/max | aln len min/median/max | gapless | has lower case |
|---|---|---|---|---|---|
| balibase | 498 | 4 / 20.5 / 384 | 156 / 1711 / 25443 | 0 | 395 |
| oxbench | 672 | 2 / 3 / 122 | 126 / 387 / 1794 | 63 | 31 |
| smart | 699 | 2 / 32 / 1769 | 54 / 378 / 3801 | 22 | 109 |

Case names: balibase `RV11_BB11001.afa`, `RV20_BBS20030.afa`, `RV64_zf_2_ref6.afa` (RV6x are extra balibase
sets); oxbench `431s1.afa`, `468.afa`; smart `ANK.afa`. Each `.afa` is the **reference alignment only**
(there is no separate unaligned directory): the unaligned input must be made by stripping `-`
(the earlier code used `datasets._strip_gaps_fasta`). Headers look like `>1aab_.dna|ref|NM_010439.3|`. Upper case
is the core marking inherited from the protein alignment (as in drive5 refs); gap `-`; all sequences have a length
divisible by 3 in all 1,869 files; IUPAC codes N Y M S W R occur rarely. Filter that reproduces the 3.5 DEV set:
3..15 sequences and alignment length <= 400 gives 29 + 192 + 104 = **325** files, equal to the earlier registry
(then minus 26 gapless = 299). MDSA balibase/oxbench cases are DNA versions of cases of the protein sets (overlap
across pools; the plan keeps them in the nucleotide pool and the grouping rule handles overlap).

## Decisions and open issues

- Rfam release pinned to 15.1; the "latest numbered release dir" was read from the EBI listing (14.9, 14.10, 15.0,
  15.1 are the last four). Release 15.1 is also what CURRENT pointed to.
- `bali3` contains BBS cases (168) as well as BB; the old DEV used only BB with XML annotations (218). The drive5 format
  has no XML: core columns come from the case marking, which can differ slightly from BAliBASE's XML core blocks
  (bench_balibase.py scores with XML masks). Whether the pool uses drive5 case-marking or the original
  `BAliBASE_R1-5.tar.gz` is a decision for the manifest builder; this step did not download the BAliBASE tarball
  (plan says one drive5 archive).
- drive5 licence: none found. I looked for a licence file in the archive and read the bench page text; I did not
  contact the author.
- BRAliBASE licence: none found in the archives; the "article CC BY" statement of the plan is not verifiable from
  the files.
- The Wayback URLs for BRAliBASE are archive copies of a dead host; their content matches the earlier hashes.
- Not run: any scoring or alignment; case counts of final pools are for the manifest builder.
- I did not record SHA-256 of every unpacked file in the registry (only of each downloaded file plus the balifam100
  tree hash); a per-file manifest of `data/raw/` can be produced by the manifest agent if wanted.
