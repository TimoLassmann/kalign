# Kalign Python Examples

This directory contains comprehensive examples demonstrating various aspects of the Kalign Python package.

## Examples Overview

### 1. `basic_usage.py` - Essential Functionality
Demonstrates core Kalign features that every user should know:

- **Basic alignment** of DNA, RNA, and protein sequences
- **File input/output** operations
- **Multi-threading** for performance
- **Alignment analysis** and statistics
- **Error handling** best practices

**Run with:**
```bash
python basic_usage.py
```

**Key learning points:**
- How to align sequences with different types
- Using threading for performance
- Reading from and writing to files
- Calculating alignment statistics and pairwise identities

### 2. `ecosystem_integration.py` - Bioinformatics Ecosystem
Shows integration with popular bioinformatics tools:

- **Biopython** integration for rich alignment objects
- **scikit-bio** integration for statistical analysis
- **Pandas** for data analysis and comparison
- **Matplotlib** for visualization
- **Comprehensive workflows** combining multiple tools

**Run with:**
```bash
python ecosystem_integration.py
```

**Prerequisites:**
```bash
# Install optional dependencies for full functionality
pip install kalign[all]  # or individual packages:
pip install biopython scikit-bio pandas matplotlib
```

**Key learning points:**
- Converting between different object formats
- Using ecosystem tools for advanced analysis
- Creating visualizations of alignments
- Building complete analysis pipelines

### 3. `performance_benchmarks.py` - Performance Optimization
Comprehensive performance testing and optimization:

- **Thread scaling** benchmarks
- **Memory usage** analysis
- **Parameter optimization** testing
- **Scalability** analysis with different dataset sizes
- **System profiling** and recommendations

**Run with:**
```bash
python performance_benchmarks.py
```

**Optional dependencies:**
```bash
pip install psutil matplotlib  # For memory monitoring and plotting
```

**Key learning points:**
- Finding optimal thread counts for your system
- Understanding memory requirements
- Optimizing gap penalties for your data
- Scaling considerations for large datasets

## Quick Start

1. **Install Kalign with examples dependencies:**
```bash
pip install kalign[all]
pip install psutil  # For performance monitoring
```

2. **Run the basic usage example:**
```bash
python python-examples/basic_usage.py
```

3. **Try ecosystem integration:**
```bash
python python-examples/ecosystem_integration.py
```

4. **Benchmark your system:**
```bash
python python-examples/performance_benchmarks.py
```

## Example Outputs

### Basic Usage Output
```
🧬 Kalign Python Examples
==================================================
Example 1: Basic Sequence Alignment
==================================================
Input sequences:
  Seq 1: ATCGATCGATCGATCG
  Seq 2: ATCGATCGTCGATCG
  Seq 3: ATCGATCGATCATCG
  Seq 4: ATCGATCGAGATCG

Aligned sequences:
  Seq 1: ATCGATCGATCGATCG
  Seq 2: ATCGATCG-TCGATCG
  Seq 3: ATCGATCGATC-ATCG
  Seq 4: ATCGATCGA--GATCG

Alignment statistics:
  Length: 16
  Gap fraction: 18.75%
  Conservation: 75.00%
  Average identity: 81.67%
```

### Performance Benchmark Output
```
🔬 Kalign Performance Benchmarking Suite
📊 System Information:
   Platform: darwin
   CPU cores: 8
   Memory: 16.0 GB
   Kalign version: 3.4.1

Thread Scaling Benchmark
Testing 1 threads...
  Average: 0.245s (±0.003s)
  Speedup: 1.00x
  Efficiency: 100.0%

Testing 4 threads...
  Average: 0.089s (±0.002s)
  Speedup: 2.75x
  Efficiency: 68.8%

🎯 Recommended thread count: 4 (efficiency: 68.8%)
```

## Understanding the Examples

### Code Structure
Each example follows a consistent pattern:

1. **Import statements** with optional dependency checking
2. **Helper functions** for specific tasks
3. **Main examples** demonstrating key concepts
4. **Error handling** and graceful degradation
5. **Cleanup** of temporary files

### Best Practices Demonstrated

- **Input validation** before alignment
- **Memory-efficient** processing for large datasets
- **Thread optimization** for your hardware
- **Format conversion** between different objects
- **Comprehensive error handling**
- **Resource cleanup** and temporary file management

### Customization

The examples are designed to be easily customizable:

```python
# Modify test sequences
sequences = [
    "YOUR_SEQUENCE_1",
    "YOUR_SEQUENCE_2",
    "YOUR_SEQUENCE_3"
]

# Adjust parameters
aligned = kalign.align(
    sequences,
    seq_type="protein",    # Change sequence type
    gap_open=-12.0,        # Adjust gap penalties
    n_threads=8            # Set thread count
)

# Add custom analysis
stats = kalign.utils.alignment_stats(aligned)
consensus = kalign.utils.consensus_sequence(aligned, threshold=0.8)
```

## Troubleshooting

### Common Issues

1. **Import errors for optional dependencies:**
```bash
# Install missing packages
pip install biopython scikit-bio pandas matplotlib psutil
```

2. **Performance issues:**
```bash
# Check system resources
python -c "import psutil; print(f'CPU: {psutil.cpu_count()}, Memory: {psutil.virtual_memory().total/1024**3:.1f}GB')"
```

3. **Memory errors with large datasets:**
```python
# Use chunked processing
def chunked_alignment(sequences, chunk_size=1000):
    results = []
    for i in range(0, len(sequences), chunk_size):
        chunk = sequences[i:i+chunk_size]
        aligned_chunk = kalign.align(chunk)
        results.extend(aligned_chunk)
    return results
```

### Getting Help

- Check the main documentation: `python-docs/`
- Run the diagnostic script: `python -c "import kalign; help(kalign)"`
- Report issues: [GitHub Issues](https://github.com/TimoLassmann/kalign/issues)

## Advanced Usage

### Creating Custom Examples

Use these examples as templates for your own applications:

```python
#!/usr/bin/env python3
"""
Custom Kalign Application

Adapt this template for your specific use case.
"""

import kalign

def your_custom_analysis(sequences):
    """Your custom analysis pipeline."""
    
    # Step 1: Align sequences
    aligned = kalign.align(sequences, seq_type="auto", n_threads=4)
    
    # Step 2: Your custom analysis
    stats = kalign.utils.alignment_stats(aligned)
    
    # Step 3: Custom output
    print(f"Analysis complete: {stats['conservation']:.2%} conservation")
    
    return aligned, stats

# Your application logic here
if __name__ == "__main__":
    sequences = ["ATCG", "ATCG", "TACG"]
    aligned, stats = your_custom_analysis(sequences)
```

### Integration with Other Tools

The examples show how to integrate Kalign with:

- **Phylogenetic tools** (Bio.Phylo)
- **Statistical analysis** (scipy, numpy)
- **Data processing** (pandas)
- **Visualization** (matplotlib, seaborn)
- **Machine learning** (scikit-learn)
- **Web frameworks** (Flask, Django)

## Contributing

To contribute additional examples:

1. Follow the existing code style and structure
2. Include comprehensive error handling
3. Add informative comments and docstrings
4. Test with various input types and sizes
5. Update this README with your example

## License

These examples are distributed under the same license as Kalign (GPL-3.0-or-later).