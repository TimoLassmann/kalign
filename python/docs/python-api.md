# Kalign Python API Reference

This document provides comprehensive documentation for the Kalign Python package, designed to make multiple sequence alignment intuitive and seamless to integrate into your bioinformatics workflows.

## Table of Contents

- [Core Functions](#core-functions)
- [Return Format Options](#return-format-options)
- [Threading Control](#threading-control)
- [I/O Operations](#io-operations)
- [Utility Functions](#utility-functions)
- [Constants and Enums](#constants-and-enums)
- [Error Handling](#error-handling)

## Core Functions

### `kalign.align()`

The main alignment function with comprehensive format and ecosystem support.

```python
def align(
    sequences: List[str],
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: Optional[int] = None,
    fmt: Literal["plain", "biopython", "skbio"] = "plain",
    ids: Optional[List[str]] = None
) -> Union[List[str], Bio.Align.MultipleSeqAlignment, skbio.TabularMSA]
```

**Parameters:**

- **`sequences`** *(List[str])* - Sequences to align as strings
- **`seq_type`** *(str | int, optional)* - Sequence type:
  - `"auto"` or `kalign.AUTO` - Auto-detect (default)
  - `"dna"` or `kalign.DNA` - DNA sequences
  - `"rna"` or `kalign.RNA` - RNA sequences
  - `"protein"` or `kalign.PROTEIN` - Protein sequences
  - `"divergent"` or `kalign.PROTEIN_DIVERGENT` - Divergent proteins
  - `"internal"` or `kalign.DNA_INTERNAL` - DNA with internal gaps
- **`gap_open`** *(float, optional)* - Gap opening penalty (default: auto)
- **`gap_extend`** *(float, optional)* - Gap extension penalty (default: auto)
- **`terminal_gap_extend`** *(float, optional)* - Terminal gap penalty (default: auto)
- **`n_threads`** *(int, optional)* - Thread count (default: global setting)
- **`fmt`** *(str, optional)* - Return format:
  - `"plain"` - List of aligned sequences (fastest)
  - `"biopython"` - `Bio.Align.MultipleSeqAlignment` object
  - `"skbio"` - `skbio.TabularMSA` object
- **`ids`** *(List[str], optional)* - Sequence identifiers

**Returns:**
- List of aligned sequences, Biopython object, or scikit-bio object

**Examples:**

```python
import kalign

# Basic alignment
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]
aligned = kalign.align(sequences)
print(aligned)
# ['ATCGATCGATCG', 'ATCG-TCGATCG', 'ATCGATC-ATCG']

# With custom parameters
aligned = kalign.align(
    sequences,
    seq_type="dna",
    gap_open=-10.0,
    gap_extend=-1.0,
    n_threads=4
)

# Biopython integration
aln_bp = kalign.align(
    sequences, 
    fmt="biopython",
    ids=["seq1", "seq2", "seq3"]
)
print(type(aln_bp))
# <class 'Bio.Align.MultipleSeqAlignment'>

# scikit-bio integration
aln_sk = kalign.align(sequences, fmt="skbio")
print(type(aln_sk))
# <class 'skbio.alignment._tabular_msa.TabularMSA'>
```

### `kalign.align_from_file()`

Align sequences directly from files.

```python
def align_from_file(
    input_file: str,
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: int = 1
) -> List[str]
```

**Supported formats:** FASTA, MSF, Clustal, aligned FASTA

**Example:**

```python
import kalign

# Align from FASTA file
aligned = kalign.align_from_file(
    "sequences.fasta",
    seq_type="protein",
    n_threads=4
)
```

### `kalign.write_alignment()`

Write aligned sequences to files.

```python
def write_alignment(
    sequences: List[str],
    output_file: str,
    format: str = "fasta"
) -> None
```

**Supported formats:** `"fasta"`, `"msf"`, `"clustal"`

**Example:**

```python
kalign.write_alignment(aligned, "output.fasta", format="fasta")
```

## Return Format Options

### Plain Format (Default)

```python
aligned = kalign.align(sequences)
# Returns: ['ATCGATCGATCG', 'ATCG-TCGATCG', 'ATCGATC-ATCG']
```

**Use when:**
- You need simple string manipulation
- Performance is critical
- Working with custom analysis pipelines

### Biopython Format

```python
aln_bp = kalign.align(sequences, fmt="biopython", ids=["s1", "s2", "s3"])

# Rich Biopython functionality
from Bio import AlignIO
AlignIO.write(aln_bp, "output.clustal", "clustal")
print(f"Length: {aln_bp.get_alignment_length()}")

# Access individual records
for record in aln_bp:
    print(f"{record.id}: {record.seq}")
```

**Use when:**
- Working with existing Biopython pipelines
- Need rich alignment manipulation features
- Exporting to multiple formats
- Phylogenetic analysis workflows

### scikit-bio Format

```python
aln_sk = kalign.align(sequences, fmt="skbio")

# scikit-bio functionality
aln_sk.write("output.fasta", format="fasta")
consensus = aln_sk.consensus()
conservation = aln_sk.conservation()

# Statistical analysis
diversity = aln_sk.distribution()
```

**Use when:**
- Statistical analysis of alignments
- Machine learning workflows
- Advanced sequence analysis
- Ecological/diversity studies

## Threading Control

### Global Thread Settings

Set default threading behavior for all alignment operations:

```python
import kalign

# Set global default
kalign.set_num_threads(8)

# All subsequent aligns use 8 threads
aligned1 = kalign.align(sequences1)  # Uses 8 threads
aligned2 = kalign.align(sequences2)  # Uses 8 threads

# Override for specific calls
aligned3 = kalign.align(sequences3, n_threads=16)  # Uses 16 threads

# Check current setting
print(f"Current default: {kalign.get_num_threads()}")
```

### Thread-Local Storage

Threading settings are thread-local for safety in multi-threaded applications:

```python
import threading
import kalign

def worker_function(thread_id):
    # Each thread can have its own default
    kalign.set_num_threads(thread_id * 2)
    
    # This thread's alignments use its specific setting
    aligned = kalign.align(sequences)
    
    return aligned

# Create multiple worker threads
threads = []
for i in range(1, 5):
    t = threading.Thread(target=worker_function, args=(i,))
    threads.append(t)
    t.start()
```

### Performance Guidelines

```python
import os
import kalign

# Automatic detection
cpu_count = os.cpu_count()
kalign.set_num_threads(cpu_count)

# Conservative approach for shared systems
kalign.set_num_threads(max(1, cpu_count // 2))

# High-performance computing
kalign.set_num_threads(min(16, cpu_count))  # Diminishing returns beyond 16
```

## I/O Operations

The `kalign.io` module provides convenient file operations with ecosystem integration.

### Reading Sequences

```python
import kalign

# Read sequences only
sequences = kalign.io.read_fasta("input.fasta")

# Read sequences with IDs
sequences, ids = kalign.io.read_sequences("input.fasta")

# Direct alignment from file
aligned = kalign.align(sequences)
```

### Writing Alignments

```python
# Write in different formats
kalign.io.write_fasta(aligned, "output.fasta", ids=ids)
kalign.io.write_clustal(aligned, "output.aln", ids=ids)
kalign.io.write_stockholm(aligned, "output.sto", ids=ids)
kalign.io.write_phylip(aligned, "output.phy", ids=ids)
```

### Complete Workflow Example

```python
import kalign

# Read → Align → Write pipeline
sequences, ids = kalign.io.read_sequences("input.fasta")

aligned = kalign.align(
    sequences,
    seq_type="protein",
    n_threads=4
)

# Write in multiple formats
kalign.io.write_fasta(aligned, "aligned.fasta", ids=ids)
kalign.io.write_clustal(aligned, "aligned.aln", ids=ids)
```

## Utility Functions

The `kalign.utils` module provides analysis and manipulation tools.

### Array Conversion

```python
import kalign
import numpy as np

aligned = kalign.align(sequences)

# Convert to NumPy array for analysis
arr = kalign.utils.to_array(aligned)
print(f"Shape: {arr.shape}")  # (n_sequences, alignment_length)

# Character-level analysis
gap_positions = np.where(arr == '-')
conserved_positions = []
for col in range(arr.shape[1]):
    if len(np.unique(arr[:, col])) == 1:
        conserved_positions.append(col)
```

### Alignment Statistics

```python
stats = kalign.utils.alignment_stats(aligned)

print(f"Alignment length: {stats['length']}")
print(f"Number of sequences: {stats['n_sequences']}")
print(f"Gap fraction: {stats['gap_fraction']:.2f}")
print(f"Conservation: {stats['conservation']:.2f}")
print(f"Average identity: {stats['identity']:.2f}")
```

### Consensus Generation

```python
# Generate consensus with different thresholds
consensus_50 = kalign.utils.consensus_sequence(aligned, threshold=0.5)
consensus_75 = kalign.utils.consensus_sequence(aligned, threshold=0.75)

print(f"50% consensus: {consensus_50}")
print(f"75% consensus: {consensus_75}")
```

### Pairwise Identity Analysis

```python
# Calculate pairwise identity matrix
identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)

print("Pairwise identities:")
print(identity_matrix)

# Find most similar sequences
max_identity = np.max(identity_matrix[np.triu_indices_from(identity_matrix, k=1)])
print(f"Highest pairwise identity: {max_identity:.2f}")
```

### Alignment Manipulation

```python
# Remove gap-only columns
trimmed = kalign.utils.remove_gap_columns(aligned)

# Remove columns with >50% gaps
filtered = kalign.utils.remove_gap_columns(aligned, threshold=0.5)

# Extract specific region
region = kalign.utils.trim_alignment(aligned, start=10, end=50)
```

## Constants and Enums

### Sequence Type Constants

```python
import kalign

# Use string constants (recommended)
aligned = kalign.align(sequences, seq_type="dna")

# Or use integer constants
aligned = kalign.align(sequences, seq_type=kalign.DNA)

# Available constants:
kalign.AUTO              # Automatic detection
kalign.DNA               # DNA sequences
kalign.RNA               # RNA sequences
kalign.PROTEIN           # Protein sequences
kalign.PROTEIN_DIVERGENT # Divergent protein sequences
kalign.DNA_INTERNAL      # DNA with internal gap preference
```

## Error Handling

### Common Exceptions

```python
import kalign

try:
    aligned = kalign.align(sequences)
except ValueError as e:
    # Input validation errors
    print(f"Invalid input: {e}")
except RuntimeError as e:
    # Alignment computation errors
    print(f"Alignment failed: {e}")
except ImportError as e:
    # Missing optional dependencies
    print(f"Missing dependency: {e}")
```

### Specific Error Cases

```python
# Empty sequence list
try:
    kalign.align([])
except ValueError:
    print("Cannot align empty sequence list")

# Invalid sequence type
try:
    kalign.align(sequences, seq_type="invalid")
except ValueError:
    print("Invalid sequence type")

# Missing ecosystem dependencies
try:
    kalign.align(sequences, fmt="biopython")
except ImportError:
    print("Install with: pip install kalign[biopython]")

# ID count mismatch
try:
    kalign.align(sequences, ids=["seq1"])  # Too few IDs
except ValueError:
    print("Number of IDs must match number of sequences")
```

### Robust Error Handling

```python
def safe_align(sequences, **kwargs):
    """Safely align sequences with comprehensive error handling."""
    try:
        if not sequences:
            raise ValueError("No sequences provided")
        
        if not all(isinstance(seq, str) for seq in sequences):
            raise ValueError("All sequences must be strings")
        
        if not all(seq.strip() for seq in sequences):
            raise ValueError("All sequences must be non-empty")
        
        return kalign.align(sequences, **kwargs)
        
    except ValueError as e:
        print(f"Input validation error: {e}")
        return None
    except RuntimeError as e:
        print(f"Alignment computation error: {e}")
        return None
    except ImportError as e:
        print(f"Dependency error: {e}")
        print("Consider installing optional dependencies:")
        print("  pip install kalign[biopython]")
        print("  pip install kalign[skbio]")
        return None

# Usage
aligned = safe_align(sequences, seq_type="protein", n_threads=4)
if aligned:
    print("Alignment successful!")
```

## Version Information

```python
import kalign

print(f"Kalign version: {kalign.__version__}")
print(f"Author: {kalign.__author__}")
print(f"Contact: {kalign.__email__}")
```

---

For more examples and use cases, see the [Quick Start Guide](python-quickstart.md) and [Ecosystem Integration Guide](python-ecosystem.md).