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

aligned = kalign.align(sequences, seq_type="dna")
for seq in aligned:
    print(seq)
```

## Core API

### `kalign.align()`

```python
aligned = kalign.align(
    sequences,              # list of str
    seq_type="auto",        # "auto", "dna", "rna", "protein", "divergent", "internal"
    gap_open=None,          # positive float, or None for defaults
    gap_extend=None,        # positive float, or None for defaults
    terminal_gap_extend=None,
    n_threads=None,         # int, or None for global default
    fmt="plain",            # "plain", "biopython", "skbio"
    ids=None,               # list of str (for biopython/skbio output)
)
```

Returns a list of aligned strings (default), a `Bio.Align.MultipleSeqAlignment` (`fmt="biopython"`), or a `skbio.TabularMSA` (`fmt="skbio"`).

### `kalign.align_from_file()`

Align sequences directly from a FASTA, MSF, or Clustal file:

```python
result = kalign.align_from_file("sequences.fasta", seq_type="protein")
for name, seq in zip(result.names, result.sequences):
    print(f"{name}: {seq}")
```

Returns an `AlignedSequences` named tuple with `.names` and `.sequences`.

### `kalign.compare()`

Score a test alignment against a reference using the Sum-of-Pairs (SP) score:

```python
score = kalign.compare("reference.msf", "test.fasta")
print(f"SP score: {score:.1f}")  # 0 (no match) to 100 (identical)
```

### `kalign.write_alignment()`

Write aligned sequences to a file:

```python
kalign.write_alignment(aligned, "output.fasta", format="fasta", ids=ids)
```

Supported formats: `fasta`, `clustal`, `stockholm`, `phylip` (non-FASTA formats require Biopython).

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

## Sequence Types

| String | Constant | Description |
|--------|----------|-------------|
| `"auto"` | `kalign.AUTO` | Auto-detect (default) |
| `"dna"` | `kalign.DNA` | DNA sequences |
| `"rna"` | `kalign.RNA` | RNA sequences |
| `"protein"` | `kalign.PROTEIN` | Protein sequences |
| `"divergent"` | `kalign.PROTEIN_DIVERGENT` | Divergent protein sequences |
| `"internal"` | `kalign.DNA_INTERNAL` | DNA with internal gap preference |

## Command-line Interface

```bash
kalign-py -i sequences.fasta -o aligned.fasta --format fasta --type protein
kalign-py -i sequences.fasta -o - --format clustal   # stdout
cat input.fa | kalign-py -i - -o aligned.fasta        # stdin
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

GNU General Public License v3.0 or later. See [COPYING](COPYING).
