# Kalign Python Package

Python bindings for [Kalign](https://github.com/TimoLassmann/kalign), a fast multiple sequence alignment program for biological sequences (DNA, RNA, protein).

## Installation

```bash
pip install kalign-python
```

Python 3.11 or later is required. On x86-64 the wheels are built with AVX2
and do not check the CPU at run time, so they need a CPU with AVX2 (Intel
Haswell, 2013, or later; AMD Excavator, 2015, or later). On older CPUs build
from source with `-DENABLE_AVX2=OFF`. The Apple silicon (arm64) wheels and Arm builds from
source do not use AVX2. Windows is not supported.

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
in the Python API and the `kalign` command-line tool with the same settings
and the same results. The sequence type is detected from the input. The one
difference is the default number of threads (see [Threading](#threading)).

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

Options that change a single alignment are only allowed with fast mode, and
options for the ensemble only with accurate mode. Anything else raises a
`ValueError`; nothing is changed silently.

- Fast mode only: `gap_open`, `gap_extend`, `terminal_gap_extend`. They
  override the fitted values. With `mode="accurate"` they raise
  `ValueError`, because the accurate ensemble uses its fitted settings.
- Accurate mode only: `fill_in`, `fill_anchor_frac`, `fill_support`. With
  fast mode they raise `ValueError`.
- Allowed in both modes: `seq_type`, `n_threads`, `fmt`, `ids`.

```python
# Fast mode with a custom gap-open penalty
aligned = kalign.align(sequences, mode="fast", gap_open=8.0)

# Raises ValueError: the accurate ensemble uses its fitted settings
kalign.align(sequences, mode="accurate", gap_open=8.0)
```

Mode constants are also available: `kalign.MODE_FAST`, `kalign.MODE_DEFAULT`,
`kalign.MODE_RECALL`, `kalign.MODE_ACCURATE`.

## Core API

### `kalign.align()`

```python
aligned = kalign.align(
    sequences,              # list of str
    mode=None,              # "fast", "default", "recall", "accurate", or None (= fast)
    seq_type="auto",        # "auto", "protein", "dna", "rna"
    gap_open=None,          # float >= 0 (fast mode only), or None for the fitted value
    gap_extend=None,        # float >= 0 (fast mode only), or None for the fitted value
    terminal_gap_extend=None,  # float >= 0 (fast mode only), or None for the fitted value
    n_threads=None,         # int >= 1, or None for the default (1; see Threading)
    fmt="plain",            # "plain", "biopython", "skbio"
    ids=None,               # list of str (for biopython/skbio output)
    fill_in=None,           # accurate mode only: None/"off", "realign", "project", "attach"
    fill_anchor_frac=None,  # fill-in anchor column threshold (0-1)
    fill_support=None,      # "attach" vote threshold (int >= 1)
)
```

Returns a list of aligned strings (default), a `Bio.Align.MultipleSeqAlignment` (`fmt="biopython"`), or a `skbio.TabularMSA` (`fmt="skbio"`).

With `mode="accurate"` (or `mode="recall"`) and `fmt="biopython"`, per-residue confidence is attached as HMMER-style PP `letter_annotations["posterior_probability"]`.

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
    # gap_open=8.0,              # only with mode="fast"
    seq_type="protein",          # "auto", "protein", "dna", "rna"
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

Supported formats: `fasta`, `clustal`, `stockholm`, `phylip` (non-FASTA formats require Biopython). `column_confidence` and `residue_confidence` are written only in Stockholm format.

### Other functions

- `kalign.align_file_to_file(input_file, output_file, format="fasta", ...)`: align a file and write the result in one call, without building Python lists. Output formats: `fasta`, `msf`, `clu`. Takes the same `mode`, `seq_type`, `n_threads` and fast-mode gap options as `align()`.
- `kalign.add_to_alignment(existing, new_sequences, output, format="fasta", n_threads=None, seq_type="auto")`: add sequences from a file to an existing alignment file. The existing rows keep their gaps. Output formats: `fasta`, `msf`, `clu`.
- `kalign.mask_alignment(result, threshold=0.5, style="lowercase")`: return a copy of an `AlignedSequences` result from `align_from_file(..., mode="accurate")` in which residues in columns below the confidence threshold are lowercased (`style="lowercase"`) or replaced by `-` (`style="remove"`).
- `kalign.filter_alignment(result, threshold=0.5)`: return a copy of such a result without the columns below the confidence threshold.
- `kalign.write_confidence(path, result)`: write the per-column confidence of such a result to a text file, one value per line.
- `kalign.generate_test_sequences(n_seq, n_obs, dna, length, seed=42)`: simulate sequences, for tests and benchmarks.

## Accurate Mode & Confidence Scores

Accurate mode runs multiple alignments with varied parameters and combines results via POAR (Pairs of Aligned Residues) consensus. The 5-member ensemble and its parameters come from the fitted setting (`mode="recall"` is an alias); the Python API does not take a run count.

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

Confidence uses HMMER-style PP encoding: a confidence of 0.95 or more is `*`; below that the character is the first digit of the confidence, so `9`=0.90-0.95, `8`=0.80-0.90, ..., `1`=0.10-0.20, `0`=below 0.10; `.`=gap.

### Biopython Per-Residue PP

When using `fmt="biopython"` with `mode="accurate"`, per-residue confidence is attached as `letter_annotations["posterior_probability"]`:

```python
aln = kalign.align(seqs, mode="accurate", fmt="biopython", ids=ids)
print(aln[0].letter_annotations["posterior_probability"])  # e.g. "998.76*9..."
```

## Threading

```python
import kalign

kalign.set_num_threads(4)        # set the default
n = kalign.get_num_threads()     # query the default

# or override per call
aligned = kalign.align(sequences, n_threads=8)
```

The Python library and `kalign-py` use 1 thread by default. The `kalign`
command line instead defaults to the number of cores minus one, at most 16.
To use more threads from Python, pass `n_threads=`, call
`kalign.set_num_threads()`, or use `kalign-py -n`.

Kalign uses at most as many threads as there are cores. A larger request is
reduced to the core count, without a warning from Python. The result does not
depend on the number of threads.

`kalign.set_num_threads(n)` sets the default for the calling thread and for
every thread that has not set its own. A thread that has called it keeps its own
value.

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

## Development interface

The compiled module `kalign._core` also has functions that take per-run
arrays of parameters, and constants such as the `REFINE_*` values and
`PROTEIN_DIVERGENT`. These were used during development to fit the settings.
They are not part of the supported API and may change. Use `kalign.align()`
and the other functions described above. The `kalign.REFINE_*` constants
re-exported at the top level belong to this interface too.

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

kalign-py -i sequences.fasta -o aligned.fasta -n 4               # more threads (default 1)

kalign-py --version
```

## Development

```bash
git clone https://github.com/TimoLassmann/kalign.git
cd kalign
uv pip install -e .
uv run pytest tests/python/ -v
```

Requirements: Python 3.11+, CMake 3.18+, C++11 compiler, NumPy.

## Citation

If you use Kalign in your research, please cite:

> Lassmann, T. (2020). Kalign 3: multiple sequence alignment of large data sets.
> *Bioinformatics*, 36(6), 1928-1929.
> [doi:10.1093/bioinformatics/btz795](https://doi.org/10.1093/bioinformatics/btz795)

## License

Apache License, Version 2.0. See [COPYING](COPYING).
