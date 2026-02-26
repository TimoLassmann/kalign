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

# Default mode — consistency anchors + VSM (best general-purpose)
aligned = kalign.align(sequences)

# Fast mode — no consistency, fastest
aligned = kalign.align(sequences, mode="fast")

# Precise mode — ensemble + realign, highest precision
aligned = kalign.align(sequences, mode="precise")
```

## Modes

Kalign v3.5 provides three named modes that package the best configurations:

| Mode | Python | CLI | Description |
|------|--------|-----|-------------|
| **default** | `mode="default"` or omit | `kalign` | Consistency anchors + VSM. Best general-purpose. |
| **fast** | `mode="fast"` | `kalign --fast` | VSM only. Fastest, similar to kalign v3.4. |
| **precise** | `mode="precise"` | `kalign --precise` | Ensemble(3) + VSM + realign. Highest precision. |

Explicit parameters always override mode defaults:

```python
# Fast base + 5 ensemble runs
aligned = kalign.align(sequences, mode="fast", ensemble=5)

# Precise base + custom gap penalty
aligned = kalign.align(sequences, mode="precise", gap_open=8.0)
```

Mode constants are also available: `kalign.MODE_DEFAULT`, `kalign.MODE_FAST`, `kalign.MODE_PRECISE`.

## Core API

### `kalign.align()`

```python
aligned = kalign.align(
    sequences,              # list of str
    mode=None,              # "default", "fast", "precise", or None (= default)
    seq_type="auto",        # "auto", "dna", "rna", "protein", "divergent", "internal"
    gap_open=None,          # positive float, or None for defaults
    gap_extend=None,        # positive float, or None for defaults
    terminal_gap_extend=None,
    n_threads=None,         # int, or None for global default
    refine="none",          # refinement mode: "none", "all", "confident", "inline"
    ensemble=0,             # number of ensemble runs (0 = off, try 3-5)
    min_support=0,          # explicit consensus threshold (0 = auto)
    fmt="plain",            # "plain", "biopython", "skbio"
    ids=None,               # list of str (for biopython/skbio output)
)
```

Returns a list of aligned strings (default), a `Bio.Align.MultipleSeqAlignment` (`fmt="biopython"`), or a `skbio.TabularMSA` (`fmt="skbio"`).

When `ensemble > 0` and `fmt="biopython"`, per-residue confidence is attached as HMMER-style PP `letter_annotations["posterior_probability"]`.

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
    mode="precise",              # or "default", "fast"
    ensemble=5,                  # override: 5 runs instead of mode default (3)
    min_support=0,               # consensus threshold (0 = auto)
    save_poar="consensus.poar",  # save POAR table for re-thresholding
    # load_poar="consensus.poar",  # OR load pre-computed POAR
    refine="none",               # refinement mode
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

Ensemble mode runs multiple alignments with varied parameters and combines results via POAR (Pairs of Aligned Residues) consensus. The simplest way to use it is `mode="precise"` (ensemble=3 + realign). For more control, set `ensemble` directly.

```python
import kalign

# Precise mode: ensemble(3) + realign — highest precision
result = kalign.align_from_file("proteins.fasta", mode="precise")

# Or: explicit 5 ensemble runs
result = kalign.align_from_file("proteins.fasta", ensemble=5)

# Per-column confidence: average agreement across ensemble runs [0..1]
print(result.column_confidence[:10])   # e.g. [0.93, 0.87, 1.0, ...]

# Per-residue confidence: per-sequence, per-position agreement [0..1]
print(result.residue_confidence[0][:10])  # e.g. [0.95, 0.90, 1.0, ...]
```

### POAR Save/Load

Save the POAR consensus table to avoid re-running the ensemble when experimenting with thresholds:

```python
# First run: compute ensemble and save POAR
kalign.align_file_to_file("input.fa", "output.fa", ensemble=5,
                          save_poar="consensus.poar")

# Later: instant re-threshold from saved POAR (no re-alignment)
kalign.align_file_to_file("input.fa", "output2.fa",
                          load_poar="consensus.poar", min_support=3)
```

### Stockholm Output with Confidence

Write confidence annotations in Stockholm format (`#=GR PP` per-residue, `#=GC PP_cons` per-column):

```python
result = kalign.align_from_file("input.fasta", ensemble=5)
kalign.write_alignment(
    result.sequences, "output.sto", format="stockholm",
    ids=result.names,
    column_confidence=result.column_confidence,
    residue_confidence=result.residue_confidence,
)
```

Confidence uses HMMER-style PP encoding: `*`=95%+, `9`=85-95%, ..., `0`=0-5%, `.`=gap.

### Biopython Per-Residue PP

When using `fmt="biopython"` with ensemble, per-residue confidence is attached as `letter_annotations["posterior_probability"]`:

```python
aln = kalign.align(seqs, ensemble=3, fmt="biopython", ids=ids)
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

| String | Constant | Description |
|--------|----------|-------------|
| `"auto"` | `kalign.AUTO` | Auto-detect (default) |
| `"dna"` | `kalign.DNA` | DNA sequences |
| `"rna"` | `kalign.RNA` | RNA sequences |
| `"protein"` | `kalign.PROTEIN` | Protein sequences |
| `"divergent"` | `kalign.PROTEIN_DIVERGENT` | Divergent protein sequences |
| `"internal"` | `kalign.DNA_INTERNAL` | DNA with internal gap preference |

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
kalign-py --fast -i sequences.fasta -o aligned.fasta             # fast mode
kalign-py --precise -i sequences.fasta -o aligned.fasta          # precise mode

# I/O options
kalign-py -i sequences.fasta -o - --format clustal               # stdout
cat input.fa | kalign-py -i - -o aligned.fasta                   # stdin
kalign-py -i sequences.fasta -o aligned.fasta --type protein     # explicit type

# Ensemble with POAR save/load
kalign-py -i input.fa -o output.fa --ensemble 5 --save-poar consensus.poar
kalign-py -i input.fa -o output.fa --load-poar consensus.poar --min-support 3

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
