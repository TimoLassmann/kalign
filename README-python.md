# Kalign Python Package

Python bindings for [Kalign](https://github.com/TimoLassmann/kalign), a fast multiple sequence alignment program for biological sequences (DNA, RNA, protein).

## Installation

```bash
pip install kalign-python
```

Optional dependencies for ecosystem integration:

```bash
pip install kalign-python[biopython]    # Biopython integration (fmt="biopython", I/O helpers)
pip install kalign-python[skbio]        # scikit-bio integration (fmt="skbio")
pip install kalign-python[io]           # I/O helpers (requires Biopython)
pip install kalign-python[analysis]     # pandas + matplotlib for downstream analysis
pip install kalign-python[all]          # all of the above
```

## Quick Start

```python
import kalign

sequences = [
    "ATCGATCGATCG",
    "ATCGTCGATCG",
    "ATCGATCATCG"
]

# Default mode (= fast): a single run
aligned = kalign.align(sequences)

# Accurate mode: a 5-member ensemble with per-column confidence
aligned = kalign.align(sequences, mode="accurate")
```

## Modes

Kalign 4 has two modes, each with a protein and a DNA/RNA setting, exposed
identically through the Python API and the `kalign` command-line tool. The
sequence type is detected from the input.

| Mode | Python | CLI | Description |
|------|--------|-----|-------------|
| **fast** | `mode="fast"`, or omit | `kalign` (or `--mode fast`) | Single run, fastest. The default. `"default"` is an alias. |
| **accurate** | `mode="accurate"` | `kalign --mode accurate` | 5-member ensemble with per-column confidence, slower. At most 5,000 sequences. `"recall"` is an alias. |

The four settings were fitted on reference alignments inside grouped nested
cross-validation: fast by multi-objective search (NSGA-III); accurate by an
NSGA-III ensemble search followed by a small grid over guide-tree noise and
consensus support. Each shipped setting is the final fit on all training data,
then scored on held-out data that the fit did not use. The fit code, the records and
the settings are in [presets/](presets/) (`presets/recipe`, `presets/docs`,
`presets/results`).

Explicit parameters always override mode defaults:

```python
# Accurate mode with a custom gap-open penalty
aligned = kalign.align(sequences, mode="accurate", gap_open=8.0)
```

Mode constants are also available: `kalign.MODE_FAST`, `kalign.MODE_DEFAULT`,
`kalign.MODE_RECALL`, `kalign.MODE_ACCURATE`.

## Core API

### `kalign.align()`

```python
aligned = kalign.align(
    sequences,              # list of str
    mode=None,              # "fast", "default", "recall", "accurate", or None (= default)
    seq_type="auto",        # "auto", "protein", "dna", "rna"
    gap_open=None,          # positive float, or None for defaults
    gap_extend=None,        # positive float, or None for defaults
    terminal_gap_extend=None,
    n_threads=None,         # int, or None for global default
    fmt="plain",            # "plain", "biopython", "skbio"
    ids=None,               # list of str (for biopython/skbio output)
    fill_in=None,           # ensemble modes only: None/"off", "realign", "project", "attach"
    fill_anchor_frac=None,  # fill-in anchor column threshold (0-1)
    fill_support=None,      # "attach" vote threshold (int >= 1)
)
```

Returns a list of aligned strings (default), a `Bio.Align.MultipleSeqAlignment` (`fmt="biopython"`), or a `skbio.TabularMSA` (`fmt="skbio"`).

For the ensemble modes (`mode="accurate"` or `mode="recall"`) with `fmt="biopython"`, per-residue confidence is attached as HMMER-style PP `letter_annotations["posterior_probability"]`.

### `kalign.align_from_file()`

Align sequences directly from a FASTA, MSF, or Clustal file:

```python
result = kalign.align_from_file("sequences.fasta", seq_type="protein")
for name, seq in zip(result.names, result.sequences):
    print(f"{name}: {seq}")
```

Returns an `AlignedSequences` object with `.names`, `.sequences`, and optional confidence fields (see below).

Additional parameters for advanced use:

```python
result = kalign.align_from_file(
    "input.fasta",
    mode="accurate",             # or "fast", "default", "recall"
    seq_type="protein",          # "auto", "protein", "dna", "rna"
    gap_open=None,               # override the preset's gap penalties
    n_threads=4,
)
```

### `kalign.compare()`

Score a test alignment against a reference using the Sum-of-Pairs (SP) score:

```python
score = kalign.compare("reference.msf", "test.fasta")
print(f"SP score: {score:.1f}")  # 0 (no match) to 100 (identical)
```

### `kalign.compare_detailed()`

Detailed alignment comparison returning POAR recall/precision/F1/TC:

```python
scores = kalign.compare_detailed("reference.msf", "test.fasta")
print(f"F1: {scores['f1']:.3f}, TC: {scores['tc']:.3f}")
```

### `kalign.write_alignment()`

Write aligned sequences to a file:

```python
kalign.write_alignment(aligned, "output.fasta", format="fasta", ids=ids)
```

Supported formats: `fasta`, `clustal`, `stockholm`, `phylip` (non-FASTA formats require Biopython).

## Ensemble Alignment & Confidence Scores

Ensemble mode runs multiple alignments with varied parameters and combines results via POAR (Pairs of Aligned Residues) consensus. Use `mode="accurate"` (5-member ensemble; `mode="recall"` is an alias). The number of runs and their parameters come from the mode preset; the Python API does not take an explicit run count.

```python
import kalign

# Accurate mode: 5-member ensemble
result = kalign.align_from_file("proteins.fasta", mode="accurate")

# Per-column confidence: average agreement across ensemble runs [0..1]
print(result.column_confidence[:10])   # e.g. [0.93, 0.87, 1.0, ...]

# Per-residue confidence: per-sequence, per-position agreement [0..1]
print(result.residue_confidence[0][:10])  # e.g. [0.95, 0.90, 1.0, ...]
```

### Stockholm Output with Confidence

Write confidence annotations in Stockholm format (`#=GR PP` per-residue, `#=GC PP_cons` per-column):

```python
result = kalign.align_from_file("input.fasta", mode="accurate")
kalign.write_alignment(
    result.sequences, "output.sto", format="stockholm",
    ids=result.names,
    column_confidence=result.column_confidence,
    residue_confidence=result.residue_confidence,
)
```

Confidence uses HMMER-style PP encoding: `*`=95%+, `9`=85-95%, ..., `0`=0-5%, `.`=gap.

### Biopython Per-Residue PP

When using `fmt="biopython"` with an ensemble mode, per-residue confidence is attached as `letter_annotations["posterior_probability"]`:

```python
aln = kalign.align(seqs, mode="accurate", fmt="biopython", ids=ids)
print(aln[0].letter_annotations["posterior_probability"])  # e.g. "998.76*9..."
```

## Threading

```python
import kalign

kalign.set_num_threads(4)        # set global default
n = kalign.get_num_threads()     # query current default

# or override per call
aligned = kalign.align(sequences, n_threads=8)
```

Thread settings are thread-local, so different threads can use different defaults.

## Utilities (`kalign.utils`)

Requires only NumPy (installed automatically):

```python
import kalign

aligned = kalign.align(sequences)

arr = kalign.utils.to_array(aligned)                          # numpy array
stats = kalign.utils.alignment_stats(aligned)                 # dict with gap_fraction, conservation, identity
consensus = kalign.utils.consensus_sequence(aligned, threshold=0.7)
matrix = kalign.utils.pairwise_identity_matrix(aligned)       # numpy array
trimmed = kalign.utils.remove_gap_columns(aligned)
region = kalign.utils.trim_alignment(aligned, start=2, end=10)
```

## Biopython Integration

Requires `pip install kalign-python[biopython]`.

```python
import kalign

# Return a Biopython MultipleSeqAlignment
aln = kalign.align(sequences, fmt="biopython", ids=["s1", "s2", "s3"])
print(aln.get_alignment_length())

# Write in various formats via Biopython
from Bio import AlignIO
AlignIO.write(aln, "output.clustal", "clustal")
```

### I/O helpers (`kalign.io`)

```python
sequences = kalign.io.read_fasta("input.fasta")
sequences, ids = kalign.io.read_sequences("input.fasta")

aligned = kalign.align(sequences)
kalign.io.write_fasta(aligned, "output.fasta", ids=ids)
kalign.io.write_clustal(aligned, "output.aln", ids=ids)
kalign.io.write_stockholm(aligned, "output.sto", ids=ids)
kalign.io.write_phylip(aligned, "output.phy", ids=ids)
```

## scikit-bio Integration

Requires `pip install kalign-python[skbio]`.

```python
import kalign

# Returns a TabularMSA of DNA, RNA, or Protein depending on seq_type
aln = kalign.align(sequences, seq_type="dna", fmt="skbio")
print(type(aln))  # <class 'skbio.alignment._tabular_msa.TabularMSA'>
```

## Constants

### Sequence Types

The same four types as the command line's `--type`. Anything but `"auto"`
overrides detection and selects the protein or the DNA/RNA presets (`"dna"`
and `"rna"` use the same presets); a `UserWarning` is issued when it differs
from the detected type.

| String | Constant | Description |
|--------|----------|-------------|
| `"auto"` | `kalign.AUTO` | Detect from the input (default) |
| `"protein"` | `kalign.PROTEIN` | Protein presets |
| `"dna"` | `kalign.DNA` | DNA/RNA presets |
| `"rna"` | `kalign.RNA` | DNA/RNA presets |

### Refinement Modes

| String | Constant | Description |
|--------|----------|-------------|
| `"none"` | `kalign.REFINE_NONE` | No refinement (default) |
| `"all"` | `kalign.REFINE_ALL` | Refine all columns |
| `"confident"` | `kalign.REFINE_CONFIDENT` | Refine only confident columns |
| `"inline"` | `kalign.REFINE_INLINE` | Inline refinement (disables parallelism) |

## Command-line Interface

```bash
# Modes
kalign-py -i sequences.fasta -o aligned.fasta                    # default mode
kalign-py --mode fast -i sequences.fasta -o aligned.fasta        # fast mode
kalign-py --mode accurate -i sequences.fasta -o aligned.fasta    # accurate mode

# I/O options
kalign-py -i sequences.fasta -o - --format clustal               # stdout
cat input.fa | kalign-py -i - -o aligned.fasta                   # stdin
kalign-py -i sequences.fasta -o aligned.fasta --type protein     # explicit type

kalign-py --version
```

## Development

```bash
git clone https://github.com/TimoLassmann/kalign.git
cd kalign
uv pip install -e .
uv run pytest tests/python/ -v
```

Requirements: Python 3.9+, CMake 3.18+, C++11 compiler, NumPy.

## Citation

If you use Kalign in your research, please cite:

> Lassmann, T. (2020). Kalign 3: multiple sequence alignment of large data sets.
> *Bioinformatics*, 36(6), 1928-1929.
> [doi:10.1093/bioinformatics/btz795](https://doi.org/10.1093/bioinformatics/btz795)

## License

Apache License, Version 2.0. See [COPYING](COPYING).
