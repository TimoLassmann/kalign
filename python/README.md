# Kalign Python Package

Python bindings for [Kalign](https://github.com/TimoLassmann/kalign), a fast multiple sequence alignment program for biological sequences.

## Installation

Install directly from PyPI:

```bash
pip install kalign
```

## Quick Start

```python
import kalign

# Align DNA sequences
sequences = [
    "ATCGATCGATCG",
    "ATCGTCGATCG", 
    "ATCGATCATCG"
]

aligned = kalign.align(sequences, seq_type="dna")
for seq in aligned:
    print(seq)
```

Output:
```
ATCGATCGATCG
ATCG-TCGATCG
ATCGATC-ATCG
```

## Usage

### Basic Alignment

The main function is `kalign.align()`:

```python
import kalign

# Auto-detect sequence type
aligned = kalign.align(sequences)

# Specify sequence type explicitly
aligned = kalign.align(sequences, seq_type="protein")

# Use custom gap penalties
aligned = kalign.align(
    sequences,
    seq_type="dna",
    gap_open=-10.0,
    gap_extend=-1.0,
    terminal_gap_extend=0.0
)

# Use multiple threads
aligned = kalign.align(sequences, n_threads=4)
```

### Sequence Types

Supported sequence types:

- `"auto"` or `kalign.AUTO` - Auto-detect (default)
- `"dna"` or `kalign.DNA` - DNA sequences
- `"rna"` or `kalign.RNA` - RNA sequences  
- `"protein"` or `kalign.PROTEIN` - Protein sequences
- `"divergent"` or `kalign.PROTEIN_DIVERGENT` - Divergent protein sequences
- `"internal"` or `kalign.DNA_INTERNAL` - DNA with internal gap preference

### File-based Alignment

Align sequences directly from files:

```python
import kalign

# Read sequences from file and align
aligned = kalign.align_from_file("sequences.fasta", seq_type="protein")

# With custom parameters
aligned = kalign.align_from_file(
    "sequences.fasta",
    seq_type="protein",
    gap_open=-10.0,
    n_threads=4
)
```

Supported input formats:
- FASTA
- MSF
- Clustal
- Aligned FASTA

### Advanced Usage

```python
import kalign

# Protein alignment with custom parameters
protein_sequences = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL"
]

aligned = kalign.align(
    protein_sequences,
    seq_type="protein",
    gap_open=-10.0,
    gap_extend=-1.0,
    terminal_gap_extend=0.0,
    n_threads=2
)

print(f"Aligned {len(aligned)} sequences")
print(f"Alignment length: {len(aligned[0])}")
```

## API Reference

### `kalign.align(sequences, seq_type="auto", gap_open=None, gap_extend=None, terminal_gap_extend=None, n_threads=1)`

Align a list of sequences.

**Parameters:**
- `sequences` (list of str): Sequences to align
- `seq_type` (str or int): Sequence type specification
- `gap_open` (float, optional): Gap opening penalty
- `gap_extend` (float, optional): Gap extension penalty  
- `terminal_gap_extend` (float, optional): Terminal gap extension penalty
- `n_threads` (int): Number of threads to use

**Returns:**
- `list of str`: Aligned sequences

### `kalign.align_from_file(input_file, seq_type="auto", gap_open=None, gap_extend=None, terminal_gap_extend=None, n_threads=1)`

Align sequences from a file.

**Parameters:**
- `input_file` (str): Path to input file
- Other parameters same as `align()`

**Returns:**
- `list of str`: Aligned sequences

## Performance

Kalign is optimized for speed and accuracy:

- **Multi-threading**: Use `n_threads` parameter for parallel processing
- **SIMD optimizations**: Automatic vectorization on supported hardware
- **Memory efficient**: Optimized memory usage for large alignments
- **Fast algorithms**: Uses advanced alignment algorithms including bit-parallel methods

## Examples

### DNA Alignment
```python
import kalign

dna_sequences = [
    "ATCGATCGATCGATCG",
    "ATCGATCGTCGATCG",
    "ATCGATCGATCATCG",
    "ATCGATCGAGATCG"
]

aligned = kalign.align(dna_sequences, seq_type="dna")
```

### RNA Alignment
```python
import kalign

rna_sequences = [
    "AUCGAUCGAUCGAUCG",
    "AUCGAUCGUCGAUCG", 
    "AUCGAUCGAUCAUCG"
]

aligned = kalign.align(rna_sequences, seq_type="rna")
```

### Protein Alignment
```python
import kalign

protein_sequences = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL"
]

aligned = kalign.align(protein_sequences, seq_type="protein")
```

## Development

### Building from Source

To build the Python package from source:

```bash
# Clone the repository
git clone https://github.com/TimoLassmann/kalign.git
cd kalign/python

# Install build dependencies
pip install build scikit-build-core pybind11

# Build the package
python -m build

# Install in development mode
pip install -e .
```

### Requirements

- Python 3.9+
- CMake 3.18+
- C++11 compatible compiler
- NumPy

## License

This package is distributed under the GNU General Public License v3.0 or later.
See the [LICENSE](../COPYING) file for details.

## Citation

If you use Kalign in your research, please cite:

Lassmann, Timo. "Kalign 3: multiple sequence alignment of large data sets." *Bioinformatics* 36.6 (2020): 1928-1929.

## Links

- [Kalign Homepage](https://github.com/TimoLassmann/kalign)
- [PyPI Package](https://pypi.org/project/kalign/)
- [Documentation](https://github.com/TimoLassmann/kalign/blob/main/README.md)
- [Issue Tracker](https://github.com/TimoLassmann/kalign/issues)