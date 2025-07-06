# Kalign Python Troubleshooting Guide

This guide helps resolve common issues when using Kalign Python package, from installation problems to performance issues and ecosystem integration challenges.

## Table of Contents

- [Installation Issues](#installation-issues)
- [Runtime Errors](#runtime-errors)
- [Performance Problems](#performance-problems)
- [Ecosystem Integration Issues](#ecosystem-integration-issues)
- [Memory Issues](#memory-issues)
- [Threading Problems](#threading-problems)
- [File I/O Issues](#file-io-issues)
- [Debugging Tools](#debugging-tools)
- [Getting Help](#getting-help)

## Installation Issues

### Problem: Package Installation Fails

**Symptoms:**
```bash
pip install kalign
ERROR: Failed building wheel for kalign
```

**Solutions:**

1. **Update build tools:**
```bash
pip install --upgrade pip setuptools wheel
pip install --upgrade cmake scikit-build-core pybind11
pip install kalign
```

2. **Install system dependencies (Linux):**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install cmake build-essential

# CentOS/RHEL
sudo yum install cmake gcc gcc-c++ make
```

3. **Install system dependencies (macOS):**
```bash
# Install Xcode command line tools
xcode-select --install

# Or install via Homebrew
brew install cmake
```

4. **Windows-specific:**
```bash
# Install Visual Studio Build Tools
# Or use conda-forge
conda install -c conda-forge kalign
```

### Problem: CMake Not Found

**Symptoms:**
```
CMake must be installed to build kalign
```

**Solutions:**

1. **Install CMake:**
```bash
pip install cmake
# or
conda install cmake
# or system package manager
```

2. **Check CMake version:**
```bash
cmake --version
# Requires CMake 3.18+
```

3. **Manual CMake installation:**
- Download from https://cmake.org/download/
- Add to PATH

### Problem: Compilation Errors

**Symptoms:**
```
error: Microsoft Visual C++ 14.0 is required
```

**Solutions:**

1. **Windows - Install Visual Studio Build Tools:**
- Download Visual Studio Build Tools
- Install C++ build tools
- Or install full Visual Studio Community

2. **Use pre-compiled wheels:**
```bash
pip install --only-binary=kalign kalign
```

3. **Try conda-forge:**
```bash
conda install -c conda-forge kalign
```

## Runtime Errors

### Problem: ImportError on First Use

**Symptoms:**
```python
import kalign
# ImportError: DLL load failed (Windows)
# ImportError: cannot find shared library (Linux)
```

**Solutions:**

1. **Check installation:**
```python
import kalign
print(kalign.__version__)
print(kalign.__file__)
```

2. **Reinstall package:**
```bash
pip uninstall kalign
pip install kalign --no-cache-dir
```

3. **Check dependencies:**
```python
import numpy  # Should work
```

4. **Try different installation method:**
```bash
conda install -c conda-forge kalign
```

### Problem: Alignment Fails with ValueError

**Symptoms:**
```python
kalign.align(sequences)
# ValueError: All sequences must be non-empty
```

**Diagnostic Script:**
```python
import kalign

def diagnose_sequences(sequences):
    """Diagnose common sequence problems."""
    
    print(f"🔍 Diagnosing {len(sequences)} sequences...")
    
    if not sequences:
        print("❌ Empty sequence list")
        return False
    
    for i, seq in enumerate(sequences):
        if not isinstance(seq, str):
            print(f"❌ Sequence {i} is not a string: {type(seq)}")
            return False
        
        if not seq.strip():
            print(f"❌ Sequence {i} is empty or whitespace only")
            return False
        
        if len(seq) < 1:
            print(f"❌ Sequence {i} has length 0")
            return False
        
        # Check for unusual characters
        unusual_chars = set(seq) - set('ATCGURYSWKMBDHVN-acdefghiklmnpqrstvwy')
        if unusual_chars:
            print(f"⚠️  Sequence {i} contains unusual characters: {unusual_chars}")
    
    print("✅ Sequences appear valid")
    return True

# Example usage
sequences = ["ATCG", "", "TACG"]  # Contains empty sequence
if diagnose_sequences(sequences):
    aligned = kalign.align(sequences)
```

**Solutions:**

1. **Remove empty sequences:**
```python
sequences = [seq for seq in sequences if seq.strip()]
```

2. **Validate input:**
```python
def clean_sequences(sequences):
    """Clean and validate sequences."""
    cleaned = []
    for seq in sequences:
        if isinstance(seq, str) and seq.strip():
            # Remove whitespace and convert to uppercase
            clean_seq = seq.strip().upper()
            cleaned.append(clean_seq)
    return cleaned

sequences = clean_sequences(raw_sequences)
```

### Problem: RuntimeError During Alignment

**Symptoms:**
```python
kalign.align(sequences)
# RuntimeError: Alignment failed: ...
```

**Diagnostic Script:**
```python
import kalign

def safe_align(sequences, **kwargs):
    """Safely align sequences with detailed error reporting."""
    
    try:
        # Basic validation
        if not sequences:
            raise ValueError("No sequences provided")
        
        if len(sequences) < 2:
            raise ValueError("Need at least 2 sequences for alignment")
        
        # Check sequence lengths
        lengths = [len(seq) for seq in sequences]
        if max(lengths) / min(lengths) > 10:
            print(f"⚠️  Large length variation: {min(lengths)} - {max(lengths)}")
        
        # Check for extremely long sequences
        if max(lengths) > 50000:
            print(f"⚠️  Very long sequences detected: max {max(lengths)}")
            print("Consider splitting or using chunked processing")
        
        # Try alignment with error handling
        print(f"Aligning {len(sequences)} sequences...")
        result = kalign.align(sequences, **kwargs)
        
        print(f"✅ Alignment successful: length {len(result[0])}")
        return result
        
    except ValueError as e:
        print(f"❌ Input validation error: {e}")
        return None
    
    except RuntimeError as e:
        print(f"❌ Alignment runtime error: {e}")
        
        # Try with different parameters
        print("Trying with modified parameters...")
        try:
            modified_kwargs = kwargs.copy()
            modified_kwargs.update({
                'gap_open': -10.0,
                'gap_extend': -1.0,
                'n_threads': 1
            })
            result = kalign.align(sequences, **modified_kwargs)
            print("✅ Alignment successful with modified parameters")
            return result
        except:
            print("❌ Alignment failed even with modified parameters")
            return None
    
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return None

# Usage
result = safe_align(sequences, seq_type="auto", n_threads=4)
```

## Performance Problems

### Problem: Alignment is Very Slow

**Symptoms:**
- Long wait times for alignment
- High CPU usage but slow progress

**Diagnostic Script:**
```python
import kalign
import time
import psutil

def performance_diagnostic(sequences, **kwargs):
    """Diagnose performance issues."""
    
    print("🚀 Performance Diagnostic")
    print(f"Sequences: {len(sequences)}")
    print(f"Lengths: {[len(seq) for seq in sequences[:5]]}...")
    
    # Check system resources
    print(f"CPU cores: {psutil.cpu_count()}")
    print(f"Available memory: {psutil.virtual_memory().available / 1024**3:.1f} GB")
    
    # Test different thread counts
    thread_counts = [1, 2, 4] if len(sequences) > 10 else [1]
    
    for threads in thread_counts:
        print(f"\nTesting with {threads} threads...")
        start_time = time.time()
        
        try:
            result = kalign.align(sequences, n_threads=threads, **kwargs)
            end_time = time.time()
            
            duration = end_time - start_time
            print(f"✅ Completed in {duration:.2f}s")
            print(f"   Throughput: {len(sequences)/duration:.1f} sequences/second")
            
            if threads == 1:
                baseline_time = duration
            else:
                speedup = baseline_time / duration
                print(f"   Speedup: {speedup:.2f}x")
                
        except Exception as e:
            print(f"❌ Failed: {e}")

# Usage
performance_diagnostic(sequences, seq_type="auto")
```

**Solutions:**

1. **Optimize thread count:**
```python
import os
optimal_threads = min(8, os.cpu_count())
kalign.set_num_threads(optimal_threads)
```

2. **Process in chunks for large datasets:**
```python
def chunked_alignment(sequences, chunk_size=1000):
    """Process large sequence sets in chunks."""
    results = []
    
    for i in range(0, len(sequences), chunk_size):
        chunk = sequences[i:i+chunk_size]
        print(f"Processing chunk {i//chunk_size + 1}...")
        aligned_chunk = kalign.align(chunk)
        results.extend(aligned_chunk)
    
    return results
```

3. **Check sequence type specification:**
```python
# Explicit type is faster than auto-detection
aligned = kalign.align(sequences, seq_type="protein")  # vs "auto"
```

### Problem: Memory Usage Too High

**Symptoms:**
- Out of memory errors
- System becomes unresponsive

**Memory Monitoring Script:**
```python
import kalign
import psutil
import gc

def memory_aware_alignment(sequences, max_memory_gb=8, **kwargs):
    """Align sequences with memory monitoring."""
    
    def get_memory_mb():
        process = psutil.Process()
        return process.memory_info().rss / 1024**2
    
    initial_memory = get_memory_mb()
    max_memory_mb = max_memory_gb * 1024
    
    print(f"Initial memory: {initial_memory:.1f} MB")
    print(f"Memory limit: {max_memory_mb:.1f} MB")
    
    # Estimate memory requirements
    total_sequence_length = sum(len(seq) for seq in sequences)
    estimated_memory = (total_sequence_length * len(sequences) * 8) / 1024**2  # Rough estimate
    
    print(f"Estimated memory needed: {estimated_memory:.1f} MB")
    
    if estimated_memory > max_memory_mb * 0.8:
        print("⚠️  High memory usage expected, using chunked processing")
        return chunked_alignment(sequences, chunk_size=500)
    
    # Monitor memory during alignment
    gc.collect()  # Clean up before starting
    
    try:
        result = kalign.align(sequences, **kwargs)
        
        final_memory = get_memory_mb()
        memory_increase = final_memory - initial_memory
        
        print(f"Final memory: {final_memory:.1f} MB")
        print(f"Memory increase: {memory_increase:.1f} MB")
        
        return result
        
    except MemoryError:
        print("❌ Out of memory - trying chunked processing")
        return chunked_alignment(sequences, chunk_size=200)

# Usage
result = memory_aware_alignment(sequences, max_memory_gb=4)
```

## Ecosystem Integration Issues

### Problem: Biopython Import Errors

**Symptoms:**
```python
kalign.align(sequences, fmt="biopython")
# ImportError: Biopython not installed
```

**Solutions:**

1. **Install Biopython:**
```bash
pip install kalign[biopython]
# or
pip install biopython
```

2. **Verify installation:**
```python
try:
    import Bio
    print(f"Biopython version: {Bio.__version__}")
    
    # Test specific modules
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    print("✅ Biopython modules accessible")
    
except ImportError as e:
    print(f"❌ Biopython issue: {e}")
```

3. **Test integration:**
```python
import kalign

def test_biopython_integration():
    """Test Biopython integration."""
    
    sequences = ["ATCG", "ATCG", "TACG"]
    
    try:
        # Test Biopython format
        aln_bp = kalign.align(sequences, fmt="biopython")
        print(f"✅ Biopython format works: {type(aln_bp)}")
        
        # Test functionality
        print(f"   Length: {aln_bp.get_alignment_length()}")
        print(f"   Sequences: {len(aln_bp)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Biopython integration failed: {e}")
        return False

test_biopython_integration()
```

### Problem: scikit-bio Import Errors

**Symptoms:**
```python
kalign.align(sequences, fmt="skbio")
# ImportError: scikit-bio not installed
```

**Solutions:**

1. **Install scikit-bio:**
```bash
pip install kalign[skbio]
# or
conda install -c conda-forge scikit-bio
```

2. **Handle version compatibility:**
```python
def test_skbio_compatibility():
    """Test scikit-bio compatibility."""
    
    try:
        import skbio
        print(f"scikit-bio version: {skbio.__version__}")
        
        # Test specific classes
        from skbio import TabularMSA
        from skbio.sequence import DNA
        print("✅ scikit-bio modules accessible")
        
        # Test with Kalign
        sequences = ["ATCG", "ATCG", "TACG"]
        aln_sk = kalign.align(sequences, fmt="skbio")
        print(f"✅ scikit-bio format works: {type(aln_sk)}")
        
        return True
        
    except ImportError as e:
        print(f"❌ scikit-bio not available: {e}")
        return False
    
    except Exception as e:
        print(f"❌ scikit-bio integration failed: {e}")
        return False

test_skbio_compatibility()
```

## Memory Issues

### Problem: Memory Leaks

**Symptoms:**
- Memory usage grows over time
- Performance degrades with repeated alignments

**Memory Leak Detection:**
```python
import kalign
import psutil
import gc

def detect_memory_leaks(n_iterations=10):
    """Detect potential memory leaks."""
    
    def get_memory():
        process = psutil.Process()
        return process.memory_info().rss / 1024**2
    
    sequences = ["ATCGATCG" * 100] * 10  # Test sequences
    
    memory_readings = []
    
    for i in range(n_iterations):
        initial_memory = get_memory()
        
        # Perform alignment
        aligned = kalign.align(sequences)
        
        # Force cleanup
        del aligned
        gc.collect()
        
        final_memory = get_memory()
        memory_readings.append(final_memory)
        
        print(f"Iteration {i+1}: {final_memory:.1f} MB")
    
    # Check for memory growth trend
    if len(memory_readings) > 5:
        early_avg = sum(memory_readings[:3]) / 3
        late_avg = sum(memory_readings[-3:]) / 3
        growth = late_avg - early_avg
        
        if growth > 10:  # More than 10 MB growth
            print(f"⚠️  Potential memory leak detected: {growth:.1f} MB growth")
        else:
            print(f"✅ No significant memory growth: {growth:.1f} MB")
    
    return memory_readings

# Usage
memory_readings = detect_memory_leaks()
```

### Problem: Large Alignment Memory Errors

**Solutions:**

1. **Streaming alignment for huge datasets:**
```python
import kalign
import tempfile

def streaming_large_alignment(sequences, temp_dir=None):
    """Handle very large alignments with streaming."""
    
    if len(sequences) < 10000:
        # Normal processing for smaller sets
        return kalign.align(sequences)
    
    print(f"Large dataset detected ({len(sequences)} sequences)")
    print("Using streaming approach...")
    
    # Split into manageable chunks
    chunk_size = 1000
    temp_alignments = []
    
    with tempfile.TemporaryDirectory(dir=temp_dir) as temp_dir:
        # Process chunks
        for i in range(0, len(sequences), chunk_size):
            chunk = sequences[i:i+chunk_size]
            print(f"Processing chunk {i//chunk_size + 1}...")
            
            # Align chunk
            aligned_chunk = kalign.align(chunk)
            
            # Save to temporary file
            temp_file = f"{temp_dir}/chunk_{i//chunk_size}.fasta"
            with open(temp_file, 'w') as f:
                for j, seq in enumerate(aligned_chunk):
                    f.write(f">seq_{i+j}\n{seq}\n")
            
            temp_alignments.append(temp_file)
            
            # Clear memory
            del aligned_chunk
        
        # Combine results (simplified - real implementation would be more complex)
        print("Combining chunks...")
        final_alignment = []
        
        for temp_file in temp_alignments:
            with open(temp_file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if not line.startswith('>'):
                        final_alignment.append(line.strip())
        
        return final_alignment
```

## Threading Problems

### Problem: Threading Not Working

**Symptoms:**
- No performance improvement with more threads
- Single CPU core usage

**Threading Diagnostic:**
```python
import kalign
import threading
import time
import psutil

def test_threading():
    """Test if threading is working correctly."""
    
    sequences = ["ATCGATCG" * 200] * 20  # Large enough to benefit from threading
    
    print("Testing threading functionality...")
    
    # Test single thread
    start_time = time.time()
    kalign.align(sequences, n_threads=1)
    single_thread_time = time.time() - start_time
    
    # Test multiple threads
    start_time = time.time()
    kalign.align(sequences, n_threads=4)
    multi_thread_time = time.time() - start_time
    
    speedup = single_thread_time / multi_thread_time
    
    print(f"Single thread: {single_thread_time:.2f}s")
    print(f"Multi thread (4): {multi_thread_time:.2f}s")
    print(f"Speedup: {speedup:.2f}x")
    
    if speedup < 1.2:
        print("⚠️  Threading may not be working effectively")
        
        # Check CPU usage during alignment
        print("Monitoring CPU usage...")
        start_time = time.time()
        cpu_before = psutil.cpu_percent(interval=None)
        
        kalign.align(sequences, n_threads=4)
        
        cpu_after = psutil.cpu_percent(interval=1)
        print(f"CPU usage: {cpu_after}%")
        
        if cpu_after < 50:
            print("❌ Low CPU usage suggests threading issues")
        else:
            print("✅ CPU usage looks normal")
    else:
        print("✅ Threading appears to be working")

# Usage
test_threading()
```

**Solutions:**

1. **Check OpenMP support:**
```python
# Check if Kalign was compiled with OpenMP
import kalign
print("Kalign version:", kalign.__version__)

# Try setting environment variables
import os
os.environ['OMP_NUM_THREADS'] = '4'
```

2. **Verify thread setting:**
```python
import kalign

# Set global thread count
kalign.set_num_threads(4)
print(f"Thread count: {kalign.get_num_threads()}")

# Ensure n_threads parameter is used
aligned = kalign.align(sequences, n_threads=4)
```

## File I/O Issues

### Problem: File Reading Errors

**Symptoms:**
```python
kalign.align_from_file("sequences.fasta")
# FileNotFoundError or parsing errors
```

**File Diagnostic Script:**
```python
import kalign
import os

def diagnose_file_issues(file_path):
    """Diagnose file reading issues."""
    
    print(f"🔍 Diagnosing file: {file_path}")
    
    # Check file existence
    if not os.path.exists(file_path):
        print(f"❌ File does not exist: {file_path}")
        return False
    
    print(f"✅ File exists: {os.path.getsize(file_path)} bytes")
    
    # Check file permissions
    if not os.access(file_path, os.R_OK):
        print(f"❌ File not readable")
        return False
    
    print("✅ File is readable")
    
    # Check file format
    try:
        with open(file_path, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(5)]
        
        print("First few lines:")
        for i, line in enumerate(first_lines):
            print(f"  {i+1}: {line[:50]}")
        
        # Basic format check
        if any(line.startswith('>') for line in first_lines):
            print("✅ Appears to be FASTA format")
        else:
            print("⚠️  May not be FASTA format")
        
    except UnicodeDecodeError:
        print("❌ File encoding issue - not UTF-8 text")
        return False
    
    # Try reading with Kalign
    try:
        sequences, ids = kalign.io.read_sequences(file_path)
        print(f"✅ Successfully read {len(sequences)} sequences")
        
        if sequences:
            lengths = [len(seq) for seq in sequences]
            print(f"   Sequence lengths: {min(lengths)} - {max(lengths)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Kalign reading failed: {e}")
        return False

# Usage
success = diagnose_file_issues("sequences.fasta")
```

**Solutions:**

1. **Fix file format issues:**
```python
def fix_fasta_format(input_file, output_file):
    """Fix common FASTA format issues."""
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_header = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # Write previous sequence
                if current_header and current_sequence:
                    outfile.write(f"{current_header}\n")
                    outfile.write(f"{''.join(current_sequence)}\n")
                
                # Start new sequence
                current_header = line
                current_sequence = []
            
            elif line:  # Non-empty line
                # Remove any non-sequence characters
                clean_line = ''.join(c for c in line.upper() if c.isalpha() or c == '-')
                if clean_line:
                    current_sequence.append(clean_line)
        
        # Write last sequence
        if current_header and current_sequence:
            outfile.write(f"{current_header}\n")
            outfile.write(f"{''.join(current_sequence)}\n")
    
    print(f"Fixed FASTA written to {output_file}")

# Usage
fix_fasta_format("problematic.fasta", "fixed.fasta")
```

## Debugging Tools

### Comprehensive Diagnostic Script

```python
import kalign
import sys
import platform
import numpy as np

def run_full_diagnostic():
    """Run comprehensive diagnostic of Kalign installation and functionality."""
    
    print("🔧 Kalign Python Diagnostic Tool")
    print("=" * 50)
    
    # System information
    print("\n📊 System Information:")
    print(f"   Python version: {sys.version}")
    print(f"   Platform: {platform.platform()}")
    print(f"   Architecture: {platform.architecture()}")
    
    # Package versions
    print("\n📦 Package Information:")
    try:
        print(f"   Kalign version: {kalign.__version__}")
        print(f"   Kalign location: {kalign.__file__}")
    except Exception as e:
        print(f"   ❌ Kalign import error: {e}")
        return
    
    try:
        print(f"   NumPy version: {np.__version__}")
    except:
        print("   ❌ NumPy not available")
    
    # Test basic functionality
    print("\n🧪 Basic Functionality Tests:")
    
    # Test 1: Simple alignment
    try:
        test_sequences = ["ATCG", "ATCG", "TACG"]
        result = kalign.align(test_sequences)
        print(f"   ✅ Basic alignment: {len(result)} sequences")
    except Exception as e:
        print(f"   ❌ Basic alignment failed: {e}")
    
    # Test 2: Threading
    try:
        kalign.set_num_threads(2)
        current_threads = kalign.get_num_threads()
        print(f"   ✅ Threading control: {current_threads} threads")
    except Exception as e:
        print(f"   ❌ Threading failed: {e}")
    
    # Test 3: Sequence types
    try:
        for seq_type in ["auto", "dna", "rna", "protein"]:
            result = kalign.align(test_sequences, seq_type=seq_type)
            print(f"   ✅ Sequence type '{seq_type}': OK")
    except Exception as e:
        print(f"   ❌ Sequence type test failed: {e}")
    
    # Test 4: Utility functions
    try:
        aligned = kalign.align(test_sequences)
        stats = kalign.utils.alignment_stats(aligned)
        array = kalign.utils.to_array(aligned)
        consensus = kalign.utils.consensus_sequence(aligned)
        print(f"   ✅ Utility functions: OK")
    except Exception as e:
        print(f"   ❌ Utility functions failed: {e}")
    
    # Test 5: Ecosystem integration (optional)
    print("\n🌐 Ecosystem Integration Tests:")
    
    # Biopython
    try:
        aln_bp = kalign.align(test_sequences, fmt="biopython")
        print(f"   ✅ Biopython integration: {type(aln_bp)}")
    except ImportError:
        print("   ⚠️  Biopython not installed (optional)")
    except Exception as e:
        print(f"   ❌ Biopython integration failed: {e}")
    
    # scikit-bio
    try:
        aln_sk = kalign.align(test_sequences, fmt="skbio")
        print(f"   ✅ scikit-bio integration: {type(aln_sk)}")
    except ImportError:
        print("   ⚠️  scikit-bio not installed (optional)")
    except Exception as e:
        print(f"   ❌ scikit-bio integration failed: {e}")
    
    print("\n✅ Diagnostic complete!")
    print("\nIf you see any ❌ errors above, check the troubleshooting guide")
    print("or report the issue at: https://github.com/TimoLassmann/kalign/issues")

# Run diagnostic
if __name__ == "__main__":
    run_full_diagnostic()
```

### Performance Profiler

```python
import kalign
import cProfile
import pstats
import io

def profile_alignment(sequences, **kwargs):
    """Profile alignment performance."""
    
    print("🔍 Profiling alignment performance...")
    
    # Create profiler
    profiler = cProfile.Profile()
    
    # Profile alignment
    profiler.enable()
    result = kalign.align(sequences, **kwargs)
    profiler.disable()
    
    # Analyze results
    stream = io.StringIO()
    ps = pstats.Stats(profiler, stream=stream)
    ps.sort_stats('cumulative').print_stats(20)
    
    profile_output = stream.getvalue()
    print(profile_output)
    
    return result

# Usage
# sequences = ["ATCGATCG" * 100] * 10
# result = profile_alignment(sequences, n_threads=2)
```

## Getting Help

### Before Reporting Issues

1. **Run the diagnostic script:**
```python
# Copy and run the full diagnostic script above
run_full_diagnostic()
```

2. **Check version compatibility:**
```python
import kalign
print(f"Kalign version: {kalign.__version__}")
print(f"Python version: {sys.version}")
```

3. **Create minimal reproducible example:**
```python
import kalign

# Minimal example that reproduces the issue
sequences = ["ATCG", "ATCG"]  # Replace with your problematic sequences
try:
    result = kalign.align(sequences)
    print("Success")
except Exception as e:
    print(f"Error: {e}")
```

### Reporting Issues

When reporting issues, please include:

1. **System information:**
   - Operating system and version
   - Python version
   - Kalign version
   
2. **Installation method:**
   - pip, conda, source compilation
   - Any special installation steps
   
3. **Complete error message:**
   - Full traceback
   - Error context
   
4. **Minimal reproducible example:**
   - Smallest code that reproduces the issue
   - Sample data if needed
   
5. **Expected vs actual behavior:**
   - What you expected to happen
   - What actually happened

### Support Channels

- **GitHub Issues**: https://github.com/TimoLassmann/kalign/issues
- **Documentation**: Check the other documentation files in this directory
- **Email**: Contact information in `kalign.__email__`

### Common Issue Templates

**Performance Issue:**
```
**System**: [OS, Python version, CPU cores]
**Kalign version**: [version]
**Dataset size**: [number of sequences, average length]
**Current performance**: [time taken, memory usage]
**Expected performance**: [what you expected]
**Code**: [minimal example]
```

**Installation Issue:**
```
**System**: [OS and version]
**Installation method**: [pip/conda/source]
**Error message**: [complete error output]
**Build tools**: [cmake version, compiler version]
```

**Integration Issue:**
```
**Ecosystem**: [Biopython/scikit-bio/other]
**Versions**: [package versions]
**Error**: [complete traceback]
**Code**: [minimal example]
```

This troubleshooting guide covers the most common issues users encounter and provides systematic approaches to diagnose and resolve them. Remember to always start with the diagnostic script to get a complete picture of your system's status.