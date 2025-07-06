# Kalign Python Quick Start Guide

Get up and running with Kalign for Python in minutes! This guide covers everything from installation to advanced usage patterns.

## Table of Contents

- [Installation](#installation)
- [Your First Alignment](#your-first-alignment)
- [Working with Files](#working-with-files)
- [Sequence Types](#sequence-types)
- [Performance Tuning](#performance-tuning)
- [Ecosystem Integration](#ecosystem-integration)
- [Common Workflows](#common-workflows)
- [Next Steps](#next-steps)

## Installation

### Basic Installation

```bash
# Install from PyPI
pip install kalign
```

### With Ecosystem Support

```bash
# For Biopython integration
pip install kalign[biopython]

# For scikit-bio integration
pip install kalign[skbio]

# For everything
pip install kalign[all]
```

### Verify Installation

```python
import kalign
print(f"Kalign version: {kalign.__version__}")

# Quick test
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]
aligned = kalign.align(sequences)
print("✅ Kalign is working!")
```

## Your First Alignment

### DNA Sequences

```python
import kalign

# Define your sequences
dna_sequences = [
    "ATCGATCGATCGATCG",
    "ATCGATCGTCGATCG",
    "ATCGATCGATCATCG",
    "ATCGATCGAGATCG"
]

# Align them (auto-detects DNA)
aligned = kalign.align(dna_sequences)

# Print results
for i, seq in enumerate(aligned):
    print(f"Seq {i+1}: {seq}")
```

**Output:**
```
Seq 1: ATCGATCGATCGATCG
Seq 2: ATCGATCG-TCGATCG
Seq 3: ATCGATCGATC-ATCG
Seq 4: ATCGATCGA--GATCG
```

### Protein Sequences

```python
import kalign

protein_sequences = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQF",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP"
]

# Explicitly specify protein type for better performance
aligned = kalign.align(protein_sequences, seq_type="protein")

for i, seq in enumerate(aligned):
    print(f"Protein {i+1}: {seq}")
```

### RNA Sequences

```python
import kalign

rna_sequences = [
    "AUCGAUCGAUCGAUCG",
    "AUCGAUCGUCGAUCG",
    "AUCGAUCGAUCAUCG"
]

aligned = kalign.align(rna_sequences, seq_type="rna")

for i, seq in enumerate(aligned):
    print(f"RNA {i+1}: {seq}")
```

## Working with Files

### Reading from Files

```python
import kalign

# Method 1: Use align_from_file (simple)
aligned = kalign.align_from_file("sequences.fasta", seq_type="protein")

# Method 2: Read then align (more control)
sequences, ids = kalign.io.read_sequences("sequences.fasta")
aligned = kalign.align(sequences, seq_type="protein")

# Print with IDs
for seq_id, seq in zip(ids, aligned):
    print(f"{seq_id}: {seq}")
```

### Writing to Files

```python
import kalign

# Align sequences
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]
aligned = kalign.align(sequences)

# Method 1: Simple write
kalign.write_alignment(aligned, "output.fasta")

# Method 2: With custom IDs and format
ids = ["sequence_1", "sequence_2", "sequence_3"]
kalign.io.write_fasta(aligned, "output.fasta", ids=ids)
kalign.io.write_clustal(aligned, "output.aln", ids=ids)
```

### Complete File Workflow

```python
import kalign

def align_file_workflow(input_file, output_prefix):
    """Complete file-based alignment workflow."""
    
    # Read sequences
    print(f"Reading sequences from {input_file}...")
    sequences, ids = kalign.io.read_sequences(input_file)
    print(f"Found {len(sequences)} sequences")
    
    # Align with optimal settings
    print("Aligning sequences...")
    aligned = kalign.align(
        sequences,
        seq_type="auto",  # Auto-detect
        n_threads=4       # Use 4 threads
    )
    
    # Write in multiple formats
    print(f"Writing alignments to {output_prefix}.*")
    kalign.io.write_fasta(aligned, f"{output_prefix}.fasta", ids=ids)
    kalign.io.write_clustal(aligned, f"{output_prefix}.aln", ids=ids)
    
    print("✅ Complete!")
    return aligned, ids

# Usage
aligned, ids = align_file_workflow("input.fasta", "aligned")
```

## Sequence Types

### Auto-Detection (Recommended)

```python
import kalign

# Kalign automatically detects sequence type
mixed_sequences = ["ATCGATCG", "AUCGAUCG", "ACDEFGHIK"]
aligned = kalign.align(mixed_sequences)  # seq_type="auto" is default
```

### Explicit Types for Better Performance

```python
import kalign

# DNA sequences
dna_seqs = ["ATCGATCG", "ATCGTCG", "ATCGCG"]
aligned = kalign.align(dna_seqs, seq_type="dna")

# Protein sequences
protein_seqs = ["ACDEFGHIK", "ACDEFGH", "ACDEFGHIKL"]
aligned = kalign.align(protein_seqs, seq_type="protein")

# Divergent proteins (for highly divergent sequences)
divergent_seqs = ["ACDEFGHIK", "MNPQRSTVW", "ACDEFGHIKL"]
aligned = kalign.align(divergent_seqs, seq_type="divergent")
```

### Using Constants

```python
import kalign

# Using string constants (recommended)
aligned = kalign.align(sequences, seq_type="protein")

# Using integer constants
aligned = kalign.align(sequences, seq_type=kalign.PROTEIN)

# Available options:
# kalign.AUTO, kalign.DNA, kalign.RNA, kalign.PROTEIN, 
# kalign.PROTEIN_DIVERGENT, kalign.DNA_INTERNAL
```

## Performance Tuning

### Threading

```python
import kalign
import os

# Set global thread count
cpu_count = os.cpu_count()
kalign.set_num_threads(cpu_count)

# All subsequent alignments use all CPUs
aligned1 = kalign.align(sequences1)
aligned2 = kalign.align(sequences2)

# Override for specific alignment
aligned3 = kalign.align(sequences3, n_threads=1)  # Single-threaded

# Check current setting
print(f"Using {kalign.get_num_threads()} threads by default")
```

### Custom Gap Penalties

```python
import kalign

# Conservative alignment (fewer gaps)
aligned = kalign.align(
    sequences,
    gap_open=-15.0,      # Higher penalty for opening gaps
    gap_extend=-2.0      # Higher penalty for extending gaps
)

# Aggressive alignment (more gaps)
aligned = kalign.align(
    sequences,
    gap_open=-5.0,       # Lower penalty for opening gaps
    gap_extend=-0.5      # Lower penalty for extending gaps
)

# Custom terminal gap handling
aligned = kalign.align(
    sequences,
    gap_open=-10.0,
    gap_extend=-1.0,
    terminal_gap_extend=0.0  # No penalty for terminal gaps
)
```

### Performance Tips

```python
import kalign

def optimize_alignment(sequences, seq_type_hint=None):
    """Optimize alignment for performance."""
    
    # 1. Specify sequence type if known
    if seq_type_hint is None:
        # Auto-detect for first few sequences
        sample = sequences[:min(10, len(sequences))]
        seq_type_hint = "auto"
    
    # 2. Use appropriate thread count
    import os
    n_threads = min(16, os.cpu_count())  # Diminishing returns after 16
    
    # 3. For large datasets, consider chunking
    if len(sequences) > 10000:
        print(f"Large dataset ({len(sequences)} sequences)")
        print("Consider splitting into smaller chunks")
    
    # 4. Perform alignment
    aligned = kalign.align(
        sequences,
        seq_type=seq_type_hint,
        n_threads=n_threads
    )
    
    return aligned

# Usage
aligned = optimize_alignment(sequences, seq_type_hint="protein")
```

## Ecosystem Integration

### Biopython Integration

```python
import kalign
from Bio import AlignIO

# Align and get Biopython object
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]
ids = ["seq1", "seq2", "seq3"]

aln_bp = kalign.align(sequences, fmt="biopython", ids=ids)

# Use Biopython's rich functionality
print(f"Alignment length: {aln_bp.get_alignment_length()}")
print(f"Number of sequences: {len(aln_bp)}")

# Export to different formats
AlignIO.write(aln_bp, "output.clustal", "clustal")
AlignIO.write(aln_bp, "output.phylip", "phylip")

# Access individual sequences
for record in aln_bp:
    print(f"{record.id}: {record.seq}")
```

### scikit-bio Integration

```python
import kalign
import skbio

# Align and get scikit-bio object
aln_sk = kalign.align(sequences, fmt="skbio")

# Use scikit-bio's analysis tools
print(f"Consensus: {aln_sk.consensus()}")

# Calculate conservation
conservation = aln_sk.conservation()
print(f"Conservation scores: {conservation}")

# Export
aln_sk.write("output.fasta", format="fasta")

# Statistical analysis
diversity = aln_sk.distribution()
```

### Interoperability

```python
import kalign

# Start with plain alignment
aligned = kalign.align(sequences)

# Convert to ecosystem objects when needed
def convert_to_biopython(aligned_seqs, ids=None):
    """Convert plain alignment to Biopython."""
    if ids is None:
        ids = [f"seq{i}" for i in range(len(aligned_seqs))]
    
    return kalign.align(
        [seq.replace('-', '') for seq in aligned_seqs],  # Remove gaps
        fmt="biopython",
        ids=ids
    )

# Usage
ids = ["sequence_1", "sequence_2", "sequence_3"]
aln_bp = convert_to_biopython(aligned, ids)
```

## Common Workflows

### Workflow 1: Basic Sequence Analysis

```python
import kalign

def analyze_sequences(sequences, name="analysis"):
    """Complete sequence analysis workflow."""
    
    print(f"🧬 Analyzing {len(sequences)} sequences ({name})")
    
    # 1. Align sequences
    print("  Aligning...")
    aligned = kalign.align(sequences, seq_type="auto", n_threads=4)
    
    # 2. Calculate statistics
    print("  Calculating statistics...")
    stats = kalign.utils.alignment_stats(aligned)
    
    print(f"  ✅ Alignment length: {stats['length']}")
    print(f"  ✅ Gap fraction: {stats['gap_fraction']:.2%}")
    print(f"  ✅ Conservation: {stats['conservation']:.2%}")
    print(f"  ✅ Average identity: {stats['identity']:.2%}")
    
    # 3. Generate consensus
    consensus = kalign.utils.consensus_sequence(aligned, threshold=0.5)
    print(f"  ✅ Consensus: {consensus}")
    
    return aligned, stats, consensus

# Usage
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]
aligned, stats, consensus = analyze_sequences(sequences, "DNA test")
```

### Workflow 2: Comparative Analysis

```python
import kalign
import numpy as np

def compare_sequences(sequences, ids=None):
    """Compare sequences with pairwise identity analysis."""
    
    if ids is None:
        ids = [f"Seq{i+1}" for i in range(len(sequences))]
    
    # Align sequences
    aligned = kalign.align(sequences)
    
    # Calculate pairwise identities
    identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
    
    # Find most and least similar pairs
    n = len(sequences)
    max_identity = 0
    min_identity = 1
    max_pair = None
    min_pair = None
    
    for i in range(n):
        for j in range(i+1, n):
            identity = identity_matrix[i, j]
            if identity > max_identity:
                max_identity = identity
                max_pair = (ids[i], ids[j])
            if identity < min_identity:
                min_identity = identity
                min_pair = (ids[i], ids[j])
    
    print(f"Most similar: {max_pair[0]} vs {max_pair[1]} ({max_identity:.2%})")
    print(f"Least similar: {min_pair[0]} vs {min_pair[1]} ({min_identity:.2%})")
    
    return aligned, identity_matrix

# Usage
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG", "GCTAGCTAGCTA"]
ids = ["Human", "Mouse", "Rat", "Chicken"]
aligned, matrix = compare_sequences(sequences, ids)
```

### Workflow 3: File Processing Pipeline

```python
import kalign
import os
from pathlib import Path

def process_sequence_files(input_dir, output_dir, file_pattern="*.fasta"):
    """Process multiple sequence files in a directory."""
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Find all matching files
    files = list(input_path.glob(file_pattern))
    print(f"Found {len(files)} files to process")
    
    results = []
    
    for file_path in files:
        print(f"Processing {file_path.name}...")
        
        try:
            # Read sequences
            sequences, ids = kalign.io.read_sequences(str(file_path))
            
            # Align
            aligned = kalign.align(sequences, n_threads=4)
            
            # Write aligned sequences
            output_file = output_path / f"aligned_{file_path.name}"
            kalign.io.write_fasta(aligned, str(output_file), ids=ids)
            
            # Calculate stats
            stats = kalign.utils.alignment_stats(aligned)
            
            results.append({
                'file': file_path.name,
                'n_sequences': len(sequences),
                'alignment_length': stats['length'],
                'conservation': stats['conservation']
            })
            
            print(f"  ✅ {len(sequences)} sequences, length {stats['length']}")
            
        except Exception as e:
            print(f"  ❌ Error: {e}")
    
    # Summary report
    print(f"\n📊 Summary of {len(results)} successful alignments:")
    for result in results:
        print(f"  {result['file']}: {result['n_sequences']} seqs, "
              f"length {result['alignment_length']}, "
              f"conservation {result['conservation']:.2%}")
    
    return results

# Usage
# results = process_sequence_files("input_sequences", "aligned_sequences")
```

### Workflow 4: Quality Control

```python
import kalign

def quality_control_alignment(sequences, ids=None, min_identity=0.3):
    """Perform quality control on sequence alignment."""
    
    if ids is None:
        ids = [f"seq{i}" for i in range(len(sequences))]
    
    print(f"🔍 Quality control for {len(sequences)} sequences")
    
    # 1. Check sequence lengths
    lengths = [len(seq.replace('-', '')) for seq in sequences]
    mean_length = np.mean(lengths)
    std_length = np.std(lengths)
    
    print(f"  Original lengths: {mean_length:.1f} ± {std_length:.1f}")
    
    # 2. Align sequences
    aligned = kalign.align(sequences)
    
    # 3. Calculate statistics
    stats = kalign.utils.alignment_stats(aligned)
    
    # 4. Check for outliers
    identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
    
    outliers = []
    for i, seq_id in enumerate(ids):
        # Calculate average identity to all other sequences
        avg_identity = np.mean([identity_matrix[i, j] 
                               for j in range(len(ids)) if i != j])
        
        if avg_identity < min_identity:
            outliers.append((seq_id, avg_identity))
    
    # 5. Report results
    print(f"  ✅ Alignment length: {stats['length']}")
    print(f"  ✅ Gap fraction: {stats['gap_fraction']:.2%}")
    print(f"  ✅ Conservation: {stats['conservation']:.2%}")
    
    if outliers:
        print(f"  ⚠️  Potential outliers ({len(outliers)}):")
        for seq_id, identity in outliers:
            print(f"    {seq_id}: {identity:.2%} average identity")
    else:
        print(f"  ✅ No outliers detected (min identity: {min_identity:.2%})")
    
    return aligned, stats, outliers

# Usage
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG", "GGGGGGGGGGGG"]
ids = ["seq1", "seq2", "seq3", "outlier"]
aligned, stats, outliers = quality_control_alignment(sequences, ids)
```

## Next Steps

Now that you're familiar with the basics, explore these advanced topics:

1. **[Ecosystem Integration Guide](python-ecosystem.md)** - Deep dive into Biopython and scikit-bio integration
2. **[Performance Tuning Guide](python-performance.md)** - Optimize for large-scale alignments
3. **[API Reference](python-api.md)** - Complete function documentation
4. **[Troubleshooting Guide](python-troubleshooting.md)** - Common issues and solutions

### Community Resources

- **GitHub Repository**: [TimoLassmann/kalign](https://github.com/TimoLassmann/kalign)
- **Issues & Support**: [GitHub Issues](https://github.com/TimoLassmann/kalign/issues)
- **PyPI Package**: [kalign](https://pypi.org/project/kalign/)

### Getting Help

```python
import kalign

# Check version and info
print(f"Kalign version: {kalign.__version__}")
print(f"Author: {kalign.__author__}")
print(f"Contact: {kalign.__email__}")

# View help for specific functions
help(kalign.align)
help(kalign.utils.alignment_stats)
```

Happy aligning! 🧬✨