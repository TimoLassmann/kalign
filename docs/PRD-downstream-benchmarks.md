# PRD: Downstream Application Benchmarks

## Goal

Demonstrate that kalign's ensemble confidence scores improve results in real-world downstream analyses, not just alignment accuracy metrics. Three pipelines targeting different communities, plus a foundational confidence calibration analysis:

0. **Confidence calibration** — foundational validation (are confidence scores meaningful?)
1. **Positive selection analysis** — evolutionary biology
2. **Phylogenetic tree accuracy** — phylogenetics / systematics
3. **HMMER homology detection** — protein family annotation / Pfam

Each pipeline produces quantitative evidence that confidence-aware alignment (ensemble + masking/weighting) outperforms raw alignment, and that kalign's integrated confidence is competitive with or better than the established GUIDANCE2 + MAFFT pipeline — at a fraction of the computational cost.

---

## Key comparator: GUIDANCE2

GUIDANCE2 (Sela et al. 2015, ~2,500 citations) is the *de facto* standard for alignment confidence scoring and column masking. It wraps an external aligner (typically MAFFT), runs it N times with bootstrap guide trees, and computes per-column and per-residue confidence. It is the method every reviewer will ask about.

**Kalign's advantage over GUIDANCE2:**
- **Speed**: GUIDANCE2 runs the full aligner N times with Perl overhead + guide tree resampling. Kalign's ensemble is native C, ~100-300x faster for comparable N.
- **Integration**: Confidence is a native output, not a post-processing step.
- **Per-residue scores**: Available directly in the alignment result (Stockholm PP format, Biopython letter_annotations).

GUIDANCE2 must appear as a method in all three downstream pipelines. The speed comparison alone justifies the paper.

### GUIDANCE2 installation

```bash
# In container
RUN wget http://guidance.tau.ac.il/ver2/GUIDANCE2.tar.gz -O /tmp/guidance2.tar.gz \
    && tar xzf /tmp/guidance2.tar.gz -C /opt \
    && chmod +x /opt/guidance.v2.02/www/Guidance/guidance.pl \
    && ln -s /opt/guidance.v2.02/www/Guidance/guidance.pl /usr/local/bin/guidance2 \
    && rm /tmp/guidance2.tar.gz
```

### GUIDANCE2 runner (`utils.py`)

```python
def run_guidance2(
    input_fasta: Path,
    output_dir: Path,
    seq_type: str = "aa",  # "aa", "nuc", "codon"
    n_bootstrap: int = 100,
    msa_program: str = "MAFFT",
) -> GuidanceResult:
    """Run GUIDANCE2 + MAFFT. Returns alignment + confidence scores."""
```

---

## Provenance and result tracking

Every result JSON includes a `provenance` block recording all information needed to reproduce the result exactly:

```python
@dataclass
class Provenance:
    timestamp: str                    # ISO 8601
    kalign_version: str               # git describe --tags --dirty
    kalign_commit: str                # git rev-parse HEAD
    container_image: Optional[str]    # podman image digest (if containerised)
    hostname: str
    cpu_model: str                    # /proc/cpuinfo or sysctl
    cpu_cores: int
    ram_gb: float
    os_version: str                   # platform.platform()
    python_version: str
    tool_versions: dict               # {"mafft": "7.520", "muscle": "5.1", ...}
    parameters: dict                  # All CLI flags passed to this run
```

### Version collection

```python
def collect_tool_versions() -> dict:
    """Collect versions of all external tools."""
    tools = {
        "kalign":   ("kalign --version", r"(\d+\.\d+\.\d+)"),
        "mafft":    ("mafft --version", r"(\S+)"),
        "muscle":   ("muscle --version", r"(\S+)"),
        "clustalo": ("clustalo --version", r"(\S+)"),
        "hmmer":    ("hmmbuild -h", r"HMMER (\S+)"),
        "iqtree":   ("iqtree2 --version", r"IQ-TREE.*version (\S+)"),
        "hyphy":    ("hyphy --version", r"(\S+)"),
        "indelible": ("indelible", r"INDELible V(\S+)"),
        "guidance2": ("guidance2 --version", r"(\S+)"),
    }
    # Returns {"tool_name": "version_string"}
```

### Result file naming

Results are stored with timestamps and git short hashes for traceability:

```
benchmarks/results/<pipeline>/
├── run_20260215_143022_abc1234.json   # Full results
├── run_20260215_143022_abc1234.log    # Console output
└── latest.json -> run_20260215_143022_abc1234.json
```

The `latest.json` symlink allows scripts and the Dash app to find the most recent run without parsing timestamps.

### Comprehensive per-case recording

**Critical design principle**: Record everything per case. You can always aggregate later but you cannot go back and measure what you did not record. Each per-case result dict includes all raw metrics, not just the primary ones. Specific additions beyond the pipeline-specific metrics:

```python
# For all pipelines with simulated data (Pipelines 1 & 2)
per_case_base = {
    # Simulation identity
    "sim_id": str,
    "sim_params": dict,           # Full parameter set
    "replicate": int,

    # Alignment quality against true alignment
    "sp_score": float,            # Sum-of-pairs score vs true alignment
    "tc_score": float,            # Total-column score vs true alignment

    # Alignment geometry
    "aln_length": int,            # Columns in alignment
    "aln_length_after_mask": int, # Columns retained after masking (if applicable)
    "gap_fraction": float,        # Fraction of gap characters
    "mean_pairwise_identity": float,

    # Confidence (ensemble methods only)
    "conf_mean": float,           # Mean column confidence
    "conf_median": float,
    "conf_q25": float,            # 25th percentile
    "conf_q75": float,            # 75th percentile
    "conf_min": float,
    "n_cols_above_50": int,       # Columns with confidence >= 0.5
    "n_cols_above_70": int,
    "n_cols_above_90": int,

    # Timing
    "align_wall_time": float,     # Alignment only (seconds)
    "downstream_wall_time": float, # Downstream analysis (FUBAR/IQ-TREE/hmmsearch)
    "peak_memory_mb": float,      # Peak RSS during alignment
}
```

---

## Architecture

### Module layout

```
benchmarks/
├── (existing files unchanged)
├── downstream/
│   ├── __init__.py
│   ├── __main__.py              # Entry: python -m benchmarks.downstream
│   ├── provenance.py            # Version collection, hardware info, Provenance dataclass
│   ├── simulation.py            # Shared: INDELible wrapper + tree generation
│   ├── calibration.py           # Pipeline 0: confidence calibration
│   ├── positive_selection.py    # Pipeline 1
│   ├── phylo_accuracy.py        # Pipeline 2
│   ├── hmmer_detection.py       # Pipeline 3
│   ├── figures.py               # Publication figure generation (matplotlib)
│   └── utils.py                 # Shared helpers (tree comparison, masking, stats)
├── data/
│   └── downloads/
│       ├── (existing: bb3_release/, data-set1/, data-set2/)
│       ├── pfam_seed/           # Pfam-A.seed families (Stockholm)
│       ├── swissprot/           # Swiss-Prot FASTA for HMMER searches
│       └── selectome/           # Selectome curated positive selection cases
└── results/
    ├── (existing JSON files)
    ├── calibration/             # Per-run JSONs
    ├── positive_selection/      # Per-run JSONs
    ├── phylo_accuracy/          # Per-run JSONs
    └── hmmer_detection/         # Per-run JSONs
```

### Container

Separate `Containerfile.downstream` extending the existing image pattern. Installs additional tools (INDELible, HyPhy, IQ-TREE, HMMER3, GUIDANCE2) alongside the existing aligners (clustalo, mafft, muscle). Does not modify the base `Containerfile`.

All external tools run inside the container to ensure version-pinned reproducibility. Results are mounted out for analysis on the host.

### CLI pattern

All pipelines follow the existing convention:

```bash
# Individual pipelines
python -m benchmarks.downstream.calibration -j 4
python -m benchmarks.downstream.positive_selection -j 4 --max-sims 50
python -m benchmarks.downstream.phylo_accuracy -j 4 --max-sims 100
python -m benchmarks.downstream.hmmer_detection -j 4 --max-families 50

# Run all
python -m benchmarks.downstream --all -j 4

# Generate publication figures from existing results
python -m benchmarks.downstream --figures -o figures/

# Download data only
python -m benchmarks.downstream --download-only

# Quick smoke test (5 cases per pipeline)
python -m benchmarks.downstream --all -j 4 --quick
```

Standard flags: `-j/--parallel`, `--max-*` (limit cases for quick testing), `--output` (JSON path), `-v/--verbose`, `--figures` (generate figures from latest results).

### Shared design patterns

- **Worker function**: `_run_one(args) -> dict` returning results or `{"error": str}`
- **Parallel execution**: `ProcessPoolExecutor` with `as_completed()`
- **Tempfile cleanup**: `NamedTemporaryFile(delete=False)` + `finally: unlink()`
- **Output**: JSON to `benchmarks/results/<pipeline>/run_YYYYMMDD_HHMMSS_<commit>.json`
- **Console**: Summary tables with per-method aggregates
- **Lazy download**: Datasets fetched on first use if not cached
- **Provenance**: Every result JSON includes full `Provenance` block
- **Memory tracking**: `resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss` after each tool invocation

---

## Pipeline 0: Confidence Calibration

### Motivation

Before using confidence scores for masking or weighting, we must demonstrate that they are *calibrated*: a confidence of 0.8 should mean the residue pair is correct ~80% of the time. This is the foundational result — it validates the entire enterprise and should be Figure 1 of the paper.

### Method

Use the simulated datasets from Pipelines 1 and 2 (true alignments are known). For each ensemble alignment:

1. For every residue pair (i,j) at each alignment column, check whether the same pair exists in the true alignment.
2. Bin pairs by their confidence score (0.0-0.05, 0.05-0.15, ..., 0.85-0.95, 0.95-1.0) — matching the PP encoding bins.
3. For each bin: `observed_accuracy = n_correct_pairs / n_total_pairs`.
4. Plot `observed_accuracy` vs `predicted_confidence` (bin midpoint).

A perfectly calibrated method gives a diagonal line. Overconfident methods lie below the diagonal.

### Additional calibration metrics

- **Brier score**: Mean squared error between confidence and binary correctness. Lower is better.
- **Expected Calibration Error (ECE)**: Weighted average of |accuracy - confidence| per bin.
- **Resolution**: Variance of per-bin accuracies (higher = more informative).
- **Column-level calibration**: Same analysis but at column granularity (average confidence vs fraction of correctly aligned residue pairs in the column).

### Comparison with GUIDANCE2

Run GUIDANCE2+MAFFT on the same simulated datasets. Extract GUIDANCE2 per-residue scores. Perform the same calibration analysis. Overlay both calibration curves on the same plot.

### Key output

```python
calibration_result = {
    "method": "kalign_ens3",
    "bins": [
        {"bin_low": 0.0, "bin_high": 0.05, "n_pairs": 12340, "n_correct": 1850, "accuracy": 0.15},
        {"bin_low": 0.05, "bin_high": 0.15, "n_pairs": 8920, "n_correct": 1070, "accuracy": 0.12},
        # ...
        {"bin_low": 0.95, "bin_high": 1.0, "n_pairs": 45200, "n_correct": 44300, "accuracy": 0.98},
    ],
    "brier_score": 0.082,
    "ece": 0.034,
    "resolution": 0.21,
    "n_total_pairs": 1234567,
}
```

---

## Pipeline 1: Positive Selection Analysis

### Motivation

Alignment errors cause false positive signals in dN/dS site-model tests. This is well-documented (Markova-Raina & Petrov 2011, Jordan & Goldman 2012, Privman et al. 2012) and a real pain point for evolutionary biologists. If confidence-masking reduces false positives while retaining true positives, that is directly actionable.

### External tools

| Tool | Purpose | Install |
|------|---------|---------|
| INDELible v1.03 | Codon sequence simulation | Build from source |
| HyPhy >=2.5 | FUBAR site-model test | `apt install hyphy` or build from source |
| GUIDANCE2 v2.02 | Alignment confidence (comparator) | Download from guidance.tau.ac.il |

**Why HyPhy/FUBAR over PAML/codeml**: FUBAR is faster (minutes vs hours per dataset), handles large alignments, and produces posterior probabilities per site. codeml (M2a vs M1a) is the classic approach but prohibitively slow for a large simulation study. We can include a small codeml validation subset if reviewers demand it.

### Data generation

**A. Simulated data (primary)**

Use INDELible to simulate codon evolution under site models:

**Fixed parameters:**
- Codon model: M8 (beta + positive selection) for positive cases, M7 (beta, no selection) for negative controls
- Codon frequencies: F3x4 (estimated from empirical data)
- Kappa (ts/tv ratio): 2.0

**Varied parameters (full factorial):**

| Parameter | Values | Purpose |
|-----------|--------|---------|
| Number of taxa | 8, 16, 32 | Alignment difficulty scales with N |
| Tree depth | 0.3, 0.7, 1.5 (expected substitutions/site) | Low / medium / high divergence |
| Indel rate | 0.01, 0.05, 0.10 (relative to substitution rate) | Controls gap density |
| Fraction of sites under positive selection | 0%, 5%, 10% | 0% = negative control |
| omega for selected sites | 3.0 | dN/dS > 1 at selected sites |
| Replicates per condition | 10 | Statistical power |

Total simulations: 3 x 3 x 3 x 3 x 10 = 810 datasets. Each is small (8-32 codon sequences, ~300-500 codons), so simulation is fast.

**INDELible control file template** (`simulation.py`):

```
[TYPE] CODON 1
[MODEL] m8model
  [submodel] 2.0          // kappa
  [statefreq] ...         // F3x4
  [rates] 0 <p0> <p1>     // site classes
  [indelmodel] POW 1.7
  [indelrate] <indel_rate>
[TREE] t1 <newick_tree>
[PARTITIONS] p1
  [t1 m8model <n_codons>]
[EVOLVE] p1 1 output
```

Each simulation produces:
- `TRUE.phy` -- true alignment (with gaps as evolved)
- `output_TRUE_1.fas` -- true unaligned sequences (gaps stripped)
- Site-class assignments (which sites are under positive selection)

**B. Empirical data (validation)**

Selectome database (selectome.org) provides curated genes with experimentally-supported evidence of positive selection, mapped to specific branches and codons. Download ~50-100 gene families with:
- Known positive selection sites (from branch-site tests with manual curation)
- Unaligned coding sequences
- Reference phylogenetic trees

This validates that simulation findings transfer to real data.

### Alignment methods

For each dataset, align the unaligned sequences with:

| Method | Label | Description |
|--------|-------|-------------|
| kalign | `kalign` | Single run, default params |
| kalign ensemble 3 | `kalign_ens3` | Ensemble, auto threshold |
| kalign ensemble 3, mask >=0.5 | `kalign_ens3_m50` | Remove columns with confidence < 0.5 |
| kalign ensemble 3, mask >=0.7 | `kalign_ens3_m70` | Remove columns with confidence < 0.7 |
| kalign ensemble 3, mask >=0.9 | `kalign_ens3_m90` | Remove columns with confidence < 0.9 |
| GUIDANCE2 + MAFFT | `guidance2_mafft` | GUIDANCE2 (100 bootstrap replicates) + MAFFT |
| GUIDANCE2 + MAFFT, mask >=0.93 | `guidance2_m93` | GUIDANCE2 default column threshold |
| MAFFT | `mafft` | `mafft --auto` |
| MUSCLE | `muscle` | MUSCLE v5 |
| Clustal Omega | `clustalo` | Default |
| True alignment | `true` | Upper bound (oracle) |

Confidence masking: after ensemble alignment, remove all columns where `column_confidence[i] < threshold`. This is a simple `utils.mask_by_confidence(sequences, column_confidence, threshold)` function.

### Selection test

Run HyPhy FUBAR on each alignment:

```bash
hyphy fubar --alignment <aln.fa> --tree <tree.nwk> --output <results.json>
```

FUBAR outputs posterior probability of positive selection per site. A site is called "positively selected" if `P(dN/dS > 1) > 0.9` (standard threshold).

**Tree for FUBAR**: Use the true tree for all methods. This isolates the effect of alignment quality from tree inference quality. (The phylo_accuracy pipeline tests tree inference separately.)

### Metrics

For each alignment method and simulation condition:

| Metric | Definition |
|--------|------------|
| TPR (sensitivity) | True positively-selected sites detected / total true positives |
| FPR (false positive rate) | Sites falsely called positive / total true negatives |
| Precision | True positives detected / total sites called positive |
| F1 | Harmonic mean of TPR and precision |
| SP score | Sum-of-pairs alignment accuracy vs true alignment |
| TC score | Total-column alignment accuracy vs true alignment |
| Alignment length | Columns in alignment (before/after masking) |
| Column retention | Fraction of columns retained after masking |
| Align wall time | Seconds for alignment step |
| FUBAR wall time | Seconds for selection test |
| Peak memory (MB) | Alignment peak RSS |

Aggregate by simulation condition (divergence, indel rate, N taxa) to show where confidence masking helps most.

### Statistical tests

For paired comparisons (same dataset, different methods):
- **Wilcoxon signed-rank test** on per-case F1, FPR
- **Effect size**: Cliff's delta (non-parametric)
- **95% confidence intervals** on all reported means via bootstrap (10,000 resamples)
- Report exact p-values; apply Holm-Bonferroni correction for multiple comparisons

### Per-case JSON structure

```json
{
  "sim_id": "n16_d0.7_i0.05_p0.10_rep3",
  "n_taxa": 16,
  "tree_depth": 0.7,
  "indel_rate": 0.05,
  "psel_fraction": 0.10,
  "replicate": 3,
  "method": "kalign_ens3_m70",

  "n_true_positive_sites": 30,
  "n_detected_positive": 28,
  "n_false_positive": 3,
  "tpr": 0.80,
  "fpr": 0.012,
  "precision": 0.88,
  "f1": 0.84,

  "sp_score": 0.91,
  "tc_score": 0.78,
  "aln_length": 423,
  "aln_length_after_mask": 380,
  "column_retention": 0.90,
  "gap_fraction": 0.12,
  "mean_pairwise_identity": 0.65,

  "conf_mean": 0.82,
  "conf_median": 0.90,
  "conf_q25": 0.68,
  "conf_q75": 0.97,

  "align_wall_time": 1.1,
  "fubar_wall_time": 12.3,
  "peak_memory_mb": 45.2,

  "fubar_per_site": [
    {"site": 1, "alpha": 2.3, "beta": 0.5, "prob_positive": 0.001},
    {"site": 2, "alpha": 0.8, "beta": 3.2, "prob_positive": 0.95}
  ]
}
```

---

## Pipeline 2: Phylogenetic Tree Accuracy

### Motivation

The standard justification for alignment quality: better alignments produce more accurate phylogenies. Simulation gives us a known true tree to compare against.

### External tools

| Tool | Purpose | Install |
|------|---------|---------|
| INDELible v1.03 | Sequence simulation | (shared with Pipeline 1) |
| IQ-TREE 2 | Maximum-likelihood tree inference | `apt install iqtree2` or build from source |
| GUIDANCE2 v2.02 | Alignment confidence (comparator) | (shared) |

Tree comparison in Python: `dendropy` (pure Python, reliable RF distance implementation, no system dependencies).

### Data generation

**A. Simulated data (primary)**

Simulate protein evolution (amino acid model, not codon -- faster, complements Pipeline 1):

**Fixed parameters:**
- Substitution model: WAG + Gamma(4)
- Sequence length: 300 amino acids

**Varied parameters:**

| Parameter | Values | Purpose |
|-----------|--------|---------|
| Number of taxa | 16, 32, 64 | Scaling |
| Tree depth | 0.5, 1.0, 2.0, 4.0 (expected subs/site) | Low -> very high divergence |
| Indel rate | 0.02, 0.05, 0.10 | Gap density |
| Mean indel length | 2, 5 (geometric distribution) | Short vs long indels |
| Replicates | 20 | |

Total: 3 x 4 x 3 x 2 x 20 = 1,440 simulations. Each produces ~16-64 protein sequences of ~300 aa -- lightweight.

**Tree generation**: Random birth-death trees using `dendropy.simulate.treesim.birth_death_tree()`, scaled to target depth. This avoids hand-picking topologies and gives realistic tree shapes.

**B. Empirical datasets (validation)**

TreeBASE (treebase.org) archived empirical datasets:
- Select ~50 protein datasets with published ML trees
- Re-align with each method, re-infer tree, compare to published tree
- Not a definitive ground truth but shows consistency with prior analyses

### Alignment methods

Same as Pipeline 1 plus:

| Method | Label | Description |
|--------|-------|-------------|
| kalign ensemble 3, weighted | `kalign_ens3_wt` | Pass column confidence as site weights to IQ-TREE |
| GUIDANCE2 + MAFFT | `guidance2_mafft` | GUIDANCE2 (100 bootstrap replicates) + MAFFT |
| GUIDANCE2 + MAFFT, masked | `guidance2_m93` | GUIDANCE2 default column masking |
| GUIDANCE2 + MAFFT, weighted | `guidance2_wt` | GUIDANCE2 column scores as IQ-TREE site weights |
| True alignment | `true` | Upper bound (oracle, simulated data only) |

IQ-TREE supports site weights via `-a weights_file`, where each column gets a weight [0,1]. This uses the continuous confidence scores without binary masking -- a unique advantage of per-column confidence.

### Tree inference

```bash
iqtree2 -s <alignment.fa> -m WAG+G4 -nt 1 --prefix <output>
# Or with site weights:
iqtree2 -s <alignment.fa> -m WAG+G4 -a <weights.txt> -nt 1 --prefix <output>
```

Use `-nt 1` per job since we parallelize across simulations.

### Metrics

| Metric | Definition |
|--------|------------|
| nRF | Normalized Robinson-Foulds distance (0 = identical topology, 1 = maximally different) |
| Quartet distance | Fraction of quartets that disagree (more sensitive than RF for large trees) |
| Branch score distance | Sum of squared branch length differences (topology + lengths) |
| SP score | Alignment accuracy vs true alignment |
| TC score | Total-column accuracy vs true alignment |
| Alignment length | Columns before/after masking |
| Column retention | Fraction retained |
| Align wall time | Seconds |
| IQ-TREE wall time | Seconds |
| IQ-TREE log-likelihood | Final ML value (sanity check) |
| Peak memory (MB) | Alignment peak RSS |

Both tree distances computed with `dendropy.calculate.treecompare`.

### Statistical tests

Same framework as Pipeline 1:
- Wilcoxon signed-rank on paired nRF values
- Cliff's delta effect sizes
- Bootstrap 95% CIs on means
- Holm-Bonferroni correction

---

## Pipeline 3: HMMER Homology Detection

### Motivation

Profile HMMs built from MSAs are the foundation of protein family databases (Pfam, InterPro). Better seed alignments produce more sensitive HMMs. This directly measures practical bioinformatics utility.

### External tools

| Tool | Purpose | Install |
|------|---------|---------|
| HMMER 3.4 | `hmmbuild` + `hmmsearch` | `apt install hmmer` |
| GUIDANCE2 v2.02 | Alignment confidence (comparator) | (shared) |

### Data

**Pfam seed alignments:**
- Source: `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz`
- Format: Stockholm (multiple families concatenated, `//` delimited)
- Parse with Biopython `Bio.AlignIO.parse(..., "stockholm")`
- Each family gives: family accession (PF00001), seed sequences + their accessions

**Swiss-Prot (search database):**
- Source: `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`
- ~570K sequences, ~250MB compressed
- Downloaded once, cached in `data/downloads/swissprot/`

**Pfam domain annotations on Swiss-Prot:**
- Source: `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.regions.tsv.gz`
  - Tab-separated: family accession, UniProt accession, start, end
  - ~500MB compressed, ~2GB uncompressed
  - Parse into a dict: `family_accession -> set(uniprot_accessions)` = ground truth

**Family selection:**
- Pick families spanning a range of sizes and divergence levels
- Filter: seed size between 10 and 200 sequences (too small = unreliable, too large = slow)
- Filter: full family size between 50 and 10,000 members in Swiss-Prot
- Sample ~200 families across this range
- Store the selected family list so runs are reproducible

### Pipeline per family

```
1. Extract seed sequences from Pfam-A.seed (strip gaps -> unaligned FASTA)
2. Align with each method
3. For ensemble: also produce confidence-masked alignments
4. hmmbuild: build HMM from each alignment
5. hmmsearch: search HMM against Swiss-Prot
6. Parse hits above E-value threshold
7. Compare detected accessions to ground truth (Pfam-A.regions)
```

### Alignment methods

| Method | Label | Description |
|--------|-------|-------------|
| Pfam seed (original) | `pfam_original` | The curated Pfam alignment -- baseline to beat is hard |
| kalign | `kalign` | Single run, default |
| kalign ensemble 3 | `kalign_ens3` | Ensemble, auto threshold |
| kalign ensemble 3, mask >=0.7 | `kalign_ens3_m70` | Confidence-masked |
| GUIDANCE2 + MAFFT | `guidance2_mafft` | GUIDANCE2 unmasked |
| GUIDANCE2 + MAFFT, masked | `guidance2_m93` | GUIDANCE2 default masking |
| MAFFT | `mafft` | `mafft --auto` |
| MUSCLE | `muscle` | MUSCLE v5 |
| Clustal Omega | `clustalo` | Default |

Including `pfam_original` as a reference is important -- Pfam seed alignments are manually curated, so they set an upper bound. If any automated method approaches or matches the curated alignment's detection power, that's noteworthy.

### HMMER commands

```bash
# Build HMM from alignment
hmmbuild --amino <family.hmm> <alignment.sto>

# Search against Swiss-Prot
hmmsearch --tblout <hits.tbl> --noali -E 0.001 <family.hmm> <swissprot.fa>
```

Parse `--tblout` output: columns are target name, accession, E-value, score, etc.

### Metrics

At a fixed E-value threshold (1e-3, the Pfam gathering threshold ballpark):

| Metric | Definition |
|--------|------------|
| Sensitivity (recall) | True family members detected / total true family members in Swiss-Prot |
| Specificity (precision) | True family members detected / total hits returned |
| F1 | Harmonic mean |
| Mean log E-value (true positives) | Lower = more confident detection |
| HMM effective number of seqs | `hmmbuild` reported neff (alignment information content) |
| Alignment length | Columns before/after masking |
| Column retention | Fraction retained |

Also compute at multiple E-value thresholds (1e-1, 1e-3, 1e-5, 1e-10) to show the full sensitivity/specificity tradeoff.

### Per-case JSON structure

```json
{
  "family": "PF00001",
  "family_name": "7tm_1",
  "n_seed_seqs": 42,
  "mean_pairwise_identity": 0.31,
  "n_true_members_swissprot": 523,
  "method": "kalign_ens3_m70",

  "aln_length": 285,
  "aln_length_after_mask": 240,
  "column_retention": 0.84,
  "gap_fraction": 0.18,
  "hmm_neff": 8.2,

  "hits_at_1e1": {"tp": 510, "fp": 12, "sens": 0.975, "prec": 0.977, "f1": 0.976},
  "hits_at_1e3": {"tp": 498, "fp": 3, "sens": 0.952, "prec": 0.994, "f1": 0.972},
  "hits_at_1e5": {"tp": 485, "fp": 1, "sens": 0.927, "prec": 0.998, "f1": 0.961},
  "hits_at_1e10": {"tp": 450, "fp": 0, "sens": 0.860, "prec": 1.000, "f1": 0.925},

  "mean_log_evalue_tp": -42.3,
  "conf_mean": 0.85,
  "conf_median": 0.92,

  "align_wall_time": 1.8,
  "hmmbuild_wall_time": 0.3,
  "hmmsearch_wall_time": 18.5,
  "peak_memory_mb": 120.0
}
```

### Statistical tests

Same framework as Pipelines 1 & 2. Family is the unit of analysis for paired tests.

---

## Publication figures

The benchmark results should produce 6 publication-quality figures. All figures use consistent color coding across methods and consistent axis scaling. Generated by `figures.py` using matplotlib with `Nature`-style formatting (7-inch single-column or 3.5-inch half-column, 300 DPI, PDF output for vector graphics).

### Color scheme

```python
METHOD_COLORS = {
    "kalign":           "#1f77b4",  # Blue
    "kalign_ens3":      "#2ca02c",  # Green
    "kalign_ens3_m50":  "#98df8a",  # Light green
    "kalign_ens3_m70":  "#006400",  # Dark green
    "kalign_ens3_m90":  "#004d00",  # Darker green
    "kalign_ens3_wt":   "#ff7f0e",  # Orange
    "guidance2_mafft":  "#d62728",  # Red
    "guidance2_m93":    "#ff9896",  # Light red
    "guidance2_wt":     "#e377c2",  # Pink
    "mafft":            "#9467bd",  # Purple
    "muscle":           "#8c564b",  # Brown
    "clustalo":         "#7f7f7f",  # Gray
    "true":             "#000000",  # Black (oracle)
    "pfam_original":    "#bcbd22",  # Olive
}
```

### Figure 1: Confidence calibration (Pipeline 0)

**Layout**: 2 panels side by side (7 x 3.5 inches)

- **Panel A**: Calibration plot. X-axis: predicted confidence (0-1). Y-axis: observed accuracy (0-1). Lines for kalign_ens3 and guidance2_mafft. Perfect calibration (y=x) as dashed gray line. Error bars from bootstrap CIs. Inset text: Brier score and ECE for each method.
- **Panel B**: Confidence score distribution. Histogram of per-residue confidence scores for both methods. Shows how scores are distributed (bimodal? concentrated near 1?).

### Figure 2: Positive selection -- FPR reduction (Pipeline 1)

**Layout**: 3 panels (7 x 7 inches, 1 row x 3 columns)

- **Panel A**: FPR vs tree depth. X-axis: divergence (0.3, 0.7, 1.5). Y-axis: false positive rate. Lines for each method. Error bars: bootstrap 95% CI. Shows masking advantage grows with divergence.
- **Panel B**: ROC curve. TPR vs FPR across confidence thresholds (0.1 to 0.99) for kalign ensemble. Each point is a threshold. GUIDANCE2 ROC overlaid. Shows optimal operating point.
- **Panel C**: F1 heatmap. Rows: methods. Columns: conditions (divergence x indel rate). Color intensity = F1 score. Highlights which conditions benefit most from confidence masking.

### Figure 3: Positive selection -- speed vs accuracy (Pipeline 1)

**Layout**: 1 panel (3.5 x 3.5 inches)

- **Scatter plot**: X-axis: wall time (log scale). Y-axis: F1 score. Points for each method. Point size proportional to alignment accuracy (SP score). This is the "money plot" -- kalign ensemble should be in the upper-left quadrant (fast + accurate), GUIDANCE2 in the upper-right (accurate but slow).

### Figure 4: Phylogenetic tree accuracy (Pipeline 2)

**Layout**: 4 panels (7 x 7 inches, 2x2 grid)

- **Panel A**: nRF vs tree depth. Lines for each method. Error bars. Standard divergence curve.
- **Panel B**: nRF vs indel rate at fixed depth (1.0). Isolates indel effect.
- **Panel C**: Delta nRF (ensemble - single-run) by condition. Bar chart showing improvement magnitude. Grouped by divergence.
- **Panel D**: Weighted vs masked vs raw. Boxplot comparing `kalign_ens3`, `kalign_ens3_m70`, `kalign_ens3_wt`, `guidance2_m93`, `guidance2_wt` at high divergence (depth 2.0+). Shows that confidence weighting (continuous) beats binary masking.

### Figure 5: HMMER detection (Pipeline 3)

**Layout**: 4 panels (7 x 7 inches, 2x2 grid)

- **Panel A**: Sensitivity at E-value 1e-3 vs family divergence. Scatter with LOWESS trend lines per method. Shows where ensemble wins.
- **Panel B**: Per-family delta sensitivity. Paired scatter: kalign_ens3_m70 sensitivity vs kalign sensitivity. Points above diagonal = ensemble wins. Color by family divergence.
- **Panel C**: Sensitivity-specificity curves. Sensitivity vs E-value threshold per method. Standard ROC-like view.
- **Panel D**: HMM information content. hmmbuild neff vs sensitivity. Shows that confidence masking doesn't destroy alignment information.

### Figure 6: Speed comparison (cross-pipeline)

**Layout**: 2 panels (7 x 3.5 inches)

- **Panel A**: Bar chart of mean alignment wall time per method across all pipelines. Log scale. Group kalign variants together, GUIDANCE2 variants together, plain aligners together. Error bars: SD across datasets.
- **Panel B**: Speedup table as heatmap. Rows: dataset size (N taxa). Columns: methods. Cell value: speedup ratio vs GUIDANCE2+MAFFT. Expected: kalign ensemble is 50-300x faster.

### Figure generation code

```python
# figures.py
def generate_all_figures(results_dir: Path, output_dir: Path):
    """Generate all 6 publication figures from result JSONs."""
    cal = load_latest_results(results_dir / "calibration")
    psel = load_latest_results(results_dir / "positive_selection")
    phylo = load_latest_results(results_dir / "phylo_accuracy")
    hmmer = load_latest_results(results_dir / "hmmer_detection")

    fig1_calibration(cal, output_dir / "fig1_calibration.pdf")
    fig2_positive_selection(psel, output_dir / "fig2_positive_selection.pdf")
    fig3_speed_accuracy(psel, output_dir / "fig3_speed_accuracy.pdf")
    fig4_phylo_accuracy(phylo, output_dir / "fig4_phylo_accuracy.pdf")
    fig5_hmmer_detection(hmmer, output_dir / "fig5_hmmer_detection.pdf")
    fig6_speed_comparison(psel, phylo, hmmer, output_dir / "fig6_speed.pdf")
```

---

## Container: `Containerfile.downstream`

```dockerfile
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# -- System dependencies --
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake g++ git curl wget ca-certificates \
    python3 python3-pip python3-venv python3-dev \
    pkg-config zlib1g-dev \
    clustalo mafft hmmer iqtree \
    perl libwww-perl libbio-perl-perl \
    && rm -rf /var/lib/apt/lists/*

# -- MUSCLE v5 from source --
RUN git clone https://github.com/rcedgar/muscle.git /tmp/muscle \
    && cd /tmp/muscle/src \
    && sed -i 's/defined(__arm64__)/defined(__arm64__) || defined(__aarch64__)/' myutils.h \
    && bash build_linux.bash \
    && cp /tmp/muscle/bin/muscle /usr/local/bin/muscle \
    && rm -rf /tmp/muscle

# -- INDELible v1.03 from source --
RUN git clone https://github.com/timmassingham/INDELible.git /tmp/indelible \
    && cd /tmp/indelible/src \
    && make \
    && cp /tmp/indelible/src/indelible /usr/local/bin/indelible \
    && rm -rf /tmp/indelible

# -- HyPhy from source --
RUN git clone https://github.com/veg/hyphy.git /tmp/hyphy \
    && cd /tmp/hyphy \
    && cmake -DCMAKE_BUILD_TYPE=Release -DNOAVX=ON . \
    && make -j$(nproc) HYPHYMP \
    && cp /tmp/hyphy/HYPHYMP /usr/local/bin/hyphy \
    && cp -r /tmp/hyphy/res /usr/local/lib/hyphy \
    && rm -rf /tmp/hyphy
ENV HYPHY_LIB=/usr/local/lib/hyphy

# -- GUIDANCE2 --
RUN wget -q http://guidance.tau.ac.il/ver2/GUIDANCE2.tar.gz -O /tmp/guidance2.tar.gz \
    && tar xzf /tmp/guidance2.tar.gz -C /opt \
    && chmod +x /opt/guidance.v2.02/www/Guidance/guidance.pl \
    && ln -s /opt/guidance.v2.02/www/Guidance/guidance.pl /usr/local/bin/guidance2 \
    && rm /tmp/guidance2.tar.gz

# -- Kalign C build --
COPY . /kalign
WORKDIR /kalign
RUN mkdir build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make -j$(nproc)

# -- Record tool versions at build time --
RUN echo "kalign=$(build/src/kalign --version 2>&1 | head -1)" > /tool_versions.txt \
    && echo "mafft=$(mafft --version 2>&1 | head -1)" >> /tool_versions.txt \
    && echo "muscle=$(muscle --version 2>&1 | head -1)" >> /tool_versions.txt \
    && echo "clustalo=$(clustalo --version 2>&1 | head -1)" >> /tool_versions.txt \
    && echo "hmmer=$(hmmbuild -h 2>&1 | grep '^# HMMER' | head -1)" >> /tool_versions.txt \
    && echo "iqtree=$(iqtree2 --version 2>&1 | grep 'IQ-TREE' | head -1)" >> /tool_versions.txt \
    && echo "indelible=1.03" >> /tool_versions.txt

# -- Python environment --
RUN python3 -m venv /venv
ENV PATH="/venv/bin:/kalign/build/src:$PATH"
RUN pip install uv \
    && uv pip install -e ".[all]" \
    && uv pip install dendropy biopython pandas matplotlib scipy seaborn

# -- Data directories --
RUN mkdir -p benchmarks/data/downloads/pfam_seed \
             benchmarks/data/downloads/swissprot \
             benchmarks/data/downloads/selectome \
             benchmarks/results/calibration \
             benchmarks/results/positive_selection \
             benchmarks/results/phylo_accuracy \
             benchmarks/results/hmmer_detection

EXPOSE 8050

CMD ["python", "-m", "benchmarks.downstream", "--help"]
```

Build and run:

```bash
# Build
podman build -f Containerfile.downstream -t kalign-downstream .

# Run a specific pipeline
podman run -it \
  -v ./benchmarks/data:/kalign/benchmarks/data \
  -v ./benchmarks/results:/kalign/benchmarks/results \
  kalign-downstream \
  python -m benchmarks.downstream.positive_selection -j 4 --max-sims 20

# Run all pipelines
podman run -it \
  -v ./benchmarks/data:/kalign/benchmarks/data \
  -v ./benchmarks/results:/kalign/benchmarks/results \
  kalign-downstream \
  python -m benchmarks.downstream --all -j 4

# Generate figures from existing results (no re-computation)
podman run -it \
  -v ./benchmarks/results:/kalign/benchmarks/results \
  -v ./benchmarks/figures:/kalign/benchmarks/figures \
  kalign-downstream \
  python -m benchmarks.downstream --figures -o benchmarks/figures/

# Quick smoke test (5 cases per pipeline)
podman run -it \
  -v ./benchmarks/data:/kalign/benchmarks/data \
  -v ./benchmarks/results:/kalign/benchmarks/results \
  kalign-downstream \
  python -m benchmarks.downstream --all -j 4 --quick
```

---

## Shared module: `simulation.py`

### INDELible wrapper

```python
def generate_indelible_dataset(
    tree: str,                  # Newick string
    model: str,                 # "WAG", "M8", "M7"
    seq_length: int,            # Amino acids or codons
    indel_rate: float,
    indel_length_mean: float,
    output_dir: Path,
    seed: int = 42,
) -> SimulatedDataset
```

Returns:
```python
@dataclass
class SimulatedDataset:
    true_alignment: Path        # FASTA with gaps
    unaligned: Path             # FASTA without gaps
    true_tree: Path             # Newick
    site_classes: List[int]     # Per-site class (0=neutral, 1=positive)
    params: dict                # All simulation parameters
```

### Tree generation

```python
def random_birth_death_tree(
    n_taxa: int,
    target_depth: float,        # Expected subs/site on longest path
    seed: int = 42,
) -> str                        # Newick string
```

Uses `dendropy.simulate.treesim.birth_death_tree()`, then scales branch lengths to target depth.

---

## Shared module: `utils.py`

### Confidence masking

```python
def mask_alignment_by_confidence(
    sequences: List[str],
    column_confidence: List[float],
    threshold: float,
) -> Tuple[List[str], int]:
    """Remove columns where confidence < threshold.
    Returns (masked_sequences, n_columns_retained)."""
```

### Confidence weighting (for IQ-TREE)

```python
def write_site_weights(
    column_confidence: List[float],
    output_path: Path,
) -> None:
    """Write IQ-TREE -a format: one weight per column."""
```

### Alignment accuracy against true alignment

```python
def alignment_accuracy(
    test_sequences: List[str],
    true_sequences: List[str],
    test_names: List[str],
    true_names: List[str],
) -> dict:
    """Compute SP and TC scores vs true alignment.
    Returns {"sp_score": float, "tc_score": float}."""
```

Uses `kalign.compare()` internally.

### Tree comparison

```python
def compare_trees(
    true_tree: str,     # Newick
    inferred_tree: str, # Newick
) -> dict:
    """Returns {"nrf": float, "quartet_dist": float, "branch_score_dist": float}"""
```

### Statistical helpers

```python
def bootstrap_ci(values: List[float], n_bootstrap: int = 10000,
                 alpha: float = 0.05) -> Tuple[float, float]:
    """Bootstrap confidence interval."""

def wilcoxon_paired(a: List[float], b: List[float]) -> dict:
    """Returns {"statistic": float, "p_value": float, "cliffs_delta": float}."""

def holm_bonferroni(p_values: List[float]) -> List[float]:
    """Apply Holm-Bonferroni correction."""
```

### Alignment method runner

```python
METHODS = {
    "kalign":           {"fn": run_kalign, "ensemble": 0, "mask": None},
    "kalign_ens3":      {"fn": run_kalign, "ensemble": 3, "mask": None},
    "kalign_ens3_m50":  {"fn": run_kalign, "ensemble": 3, "mask": 0.5},
    "kalign_ens3_m70":  {"fn": run_kalign, "ensemble": 3, "mask": 0.7},
    "kalign_ens3_m90":  {"fn": run_kalign, "ensemble": 3, "mask": 0.9},
    "kalign_ens3_wt":   {"fn": run_kalign, "ensemble": 3, "mask": None, "weights": True},
    "guidance2_mafft":  {"fn": run_guidance2, "n_bootstrap": 100, "mask": None},
    "guidance2_m93":    {"fn": run_guidance2, "n_bootstrap": 100, "mask": 0.93},
    "guidance2_wt":     {"fn": run_guidance2, "n_bootstrap": 100, "mask": None, "weights": True},
    "mafft":            {"fn": run_mafft},
    "muscle":           {"fn": run_muscle},
    "clustalo":         {"fn": run_clustalo},
}
```

Centralises the alignment step so all pipelines share the same method definitions and produce comparable results.

---

## Data download strategy

### Pfam data (`hmmer_detection.py`)

```python
PFAM_SEED_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz"
PFAM_REGIONS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.regions.tsv.gz"
SWISSPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
```

Downloads are cached in `benchmarks/data/downloads/`. The `Pfam-A.seed.gz` file is ~250MB, `Pfam-A.regions.tsv.gz` is ~500MB, Swiss-Prot is ~250MB compressed. Total: ~1GB download, ~3GB on disk.

**Family selection**: On first run, parse `Pfam-A.seed.gz` and `Pfam-A.regions.tsv.gz` to build a family index. Select ~200 families matching the size criteria (10-200 seed seqs, 50-10,000 full members in Swiss-Prot). Save the selected list to `benchmarks/data/downloads/pfam_seed/family_selection.json` for reproducibility.

### Selectome data (`positive_selection.py`)

```python
SELECTOME_URL = "https://selectome.org/download"
```

Download curated gene family alignments with annotated positive selection sites. Cache in `benchmarks/data/downloads/selectome/`.

### INDELible data (Pipelines 0, 1 & 2)

No external data -- everything is simulated. INDELible control files are generated programmatically and executed in temp directories.

---

## Timing estimates

Per-pipeline wall-clock time (on an 8-core machine, `-j 4`):

| Pipeline | N cases | Bottleneck | Est. time |
|----------|---------|-----------|-----------|
| Calibration | Reuses Pipeline 1+2 sims | Analysis only | ~30 min |
| Positive selection | 810 sims x 11 methods | HyPhy FUBAR (~30s each) + GUIDANCE2 (~2 min each) | ~12 hours |
| Phylo accuracy | 1440 sims x 12 methods | IQ-TREE (~10s each) + GUIDANCE2 (~2 min each) | ~10 hours |
| HMMER detection | 200 families x 9 methods | hmmsearch (~20s each) + GUIDANCE2 (~3 min each) | ~3 hours |

GUIDANCE2 is the bottleneck -- it is ~50-300x slower than kalign ensemble per dataset. This is itself a key result.

Quick mode (`--quick`): 5 cases per pipeline, ~15 minutes total. Useful for CI and development.

Full run mode (`--all`): ~25 hours total. Run overnight in container.

---

## Success criteria

The benchmarks succeed if they demonstrate:

1. **Calibration**: Kalign confidence scores are reasonably calibrated (ECE < 0.05) and comparable to or better than GUIDANCE2 calibration
2. **Positive selection**: Confidence-masked kalign ensemble has >=10% lower FPR than raw kalign at matched TPR, competitive with GUIDANCE2+MAFFT
3. **Phylo accuracy**: Confidence-weighted kalign ensemble has >=5% lower nRF than raw kalign, especially at high divergence
4. **HMMER detection**: Kalign ensemble seed alignments match or approach Pfam curated alignment sensitivity, outperforming single-run aligners
5. **Speed**: Kalign ensemble is >=50x faster than GUIDANCE2+MAFFT for comparable downstream accuracy

The publication story requires at least 3 of these 5 criteria to be met, with the speed advantage being the most reliable.

---

## Implementation order

1. **`provenance.py` + `utils.py` + `simulation.py`** -- shared infrastructure, no external tool dependencies for unit testing
2. **`calibration.py`** -- simple analysis on simulated data, validates confidence scoring
3. **Pipeline 2 (phylo accuracy)** -- simplest downstream pipeline (simulate -> align -> tree -> compare), validates simulation and alignment infrastructure
4. **Pipeline 1 (positive selection)** -- builds on simulation infrastructure, adds HyPhy, Selectome empirical data
5. **Pipeline 3 (HMMER detection)** -- independent of simulation, tests data download pipeline
6. **GUIDANCE2 integration** -- add GUIDANCE2 as a method across all pipelines (may be done incrementally during steps 2-5)
7. **`figures.py`** -- publication figure generation from result JSONs
8. **`Containerfile.downstream`** -- once all pipelines work locally, containerise
9. **`__main__.py`** -- unified entry point with `--all`, `--figures`, `--quick`

---

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| INDELible build fails on ARM64 | Provide pre-compiled binary or use Docker x86 emulation |
| HyPhy FUBAR is too slow for 810 x 11 = 8,910 runs | Use `--max-sims` to limit; FUBAR is ~30s per small dataset so 12h is feasible |
| GUIDANCE2 installation is fragile (Perl deps) | Pin version in container; test in CI; fall back to GUIDANCE v1 if needed |
| GUIDANCE2 is extremely slow | Expected -- this is a key result; use `--max-sims` for dev; full run is overnight |
| Pfam FTP is unreliable | Cache aggressively; provide fallback mirror URLs |
| Swiss-Prot is large (~250MB) | Download once, cache; consider UniRef50 subset if too large |
| hmmsearch against full Swiss-Prot is slow | Pre-index with `hmmpress`; ~20s per family is manageable for 200 families |
| Confidence masking removes too many columns | Test multiple thresholds; report column retention alongside accuracy |
| Ensemble doesn't help for easy cases | Expected -- focus analysis on high-divergence / high-indel conditions where it matters |
| Pfam seed alignment is unbeatable | Likely for well-curated families -- the win may be on divergent/poorly-curated families, or on speed (automated vs manual curation) |
| Selectome data format changes | Pin dataset version; store downloaded data with provenance |
| Calibration is poor | This would be a negative result worth reporting honestly; may indicate need for calibration post-processing (Platt scaling) |
| Statistical tests show non-significant differences | Report effect sizes alongside p-values; some conditions may show strong effects even if global average is modest |
