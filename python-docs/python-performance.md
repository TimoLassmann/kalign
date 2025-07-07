# Kalign Python Performance Tuning Guide

This guide provides comprehensive strategies for optimizing Kalign performance for different use cases, from small-scale interactive analyses to high-throughput production pipelines.

## Table of Contents

- [Understanding Performance](#understanding-performance)
- [Threading Optimization](#threading-optimization)
- [Memory Management](#memory-management)
- [Algorithm Parameter Tuning](#algorithm-parameter-tuning)
- [Batch Processing](#batch-processing)
- [Large-Scale Workflows](#large-scale-workflows)
- [Benchmarking Tools](#benchmarking-tools)
- [Hardware Considerations](#hardware-considerations)
- [Profiling and Debugging](#profiling-and-debugging)

## Understanding Performance

### Performance Factors

Kalign's performance depends on several key factors:

1. **Number of sequences** - Algorithmic complexity scales
2. **Sequence length** - Memory and computation increase
3. **Sequence similarity** - Affects alignment algorithm efficiency
4. **Thread count** - Parallel processing utilization
5. **Hardware** - CPU, memory, and cache performance
6. **Algorithm parameters** - Gap penalties affect computation

### Performance Characteristics

```python
import kalign
import time
import matplotlib.pyplot as plt
import numpy as np

def benchmark_scaling(max_sequences=20, seq_length=1000):
    """Benchmark how alignment time scales with number of sequences."""
    
    # Generate test sequences
    base_seq = "ATCG" * (seq_length // 4)
    
    results = []
    sequence_counts = range(2, max_sequences + 1, 2)
    
    for n_seqs in sequence_counts:
        print(f"Benchmarking {n_seqs} sequences...")
        
        # Create test sequences with some variation
        sequences = []
        for i in range(n_seqs):
            # Add some mutations
            seq = list(base_seq)
            n_mutations = seq_length // 20  # 5% mutations
            mutation_positions = np.random.choice(len(seq), n_mutations, replace=False)
            
            for pos in mutation_positions:
                seq[pos] = np.random.choice(['A', 'T', 'C', 'G'])
            
            sequences.append(''.join(seq))
        
        # Benchmark alignment
        start_time = time.time()
        aligned = kalign.align(sequences, n_threads=4)
        end_time = time.time()
        
        duration = end_time - start_time
        results.append({
            'n_sequences': n_seqs,
            'time_seconds': duration,
            'sequences_per_second': n_seqs / duration
        })
        
        print(f"  {n_seqs} sequences: {duration:.2f}s ({n_seqs/duration:.1f} seq/s)")
    
    # Plot results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    n_seqs = [r['n_sequences'] for r in results]
    times = [r['time_seconds'] for r in results]
    throughput = [r['sequences_per_second'] for r in results]
    
    ax1.plot(n_seqs, times, 'bo-')
    ax1.set_xlabel('Number of Sequences')
    ax1.set_ylabel('Alignment Time (seconds)')
    ax1.set_title('Scaling: Time vs Number of Sequences')
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(n_seqs, throughput, 'ro-')
    ax2.set_xlabel('Number of Sequences')
    ax2.set_ylabel('Throughput (sequences/second)')
    ax2.set_title('Scaling: Throughput vs Number of Sequences')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return results

# Example usage
# results = benchmark_scaling(max_sequences=16, seq_length=500)
```

## Threading Optimization

### Optimal Thread Count

```python
import kalign
import os
import time
import matplotlib.pyplot as plt

def find_optimal_threads(sequences, max_threads=None):
    """Find optimal thread count for given sequences."""
    
    if max_threads is None:
        max_threads = min(16, os.cpu_count())  # Diminishing returns after 16
    
    print(f"🔍 Finding optimal thread count (testing 1-{max_threads} threads)")
    
    results = []
    thread_counts = range(1, max_threads + 1)
    
    for n_threads in thread_counts:
        print(f"  Testing {n_threads} threads...", end=' ')
        
        # Run multiple iterations for accuracy
        times = []
        for _ in range(3):
            start_time = time.time()
            aligned = kalign.align(sequences, n_threads=n_threads)
            end_time = time.time()
            times.append(end_time - start_time)
        
        avg_time = np.mean(times)
        std_time = np.std(times)
        speedup = times[0] / avg_time if len(results) == 0 else results[0]['avg_time'] / avg_time
        
        results.append({
            'threads': n_threads,
            'avg_time': avg_time,
            'std_time': std_time,
            'speedup': speedup
        })
        
        print(f"{avg_time:.3f}s (±{std_time:.3f}s), speedup: {speedup:.2f}x")
    
    # Find optimal thread count (best speedup vs overhead trade-off)
    # Look for point where speedup gain is less than 10%
    optimal_threads = 1
    for i in range(1, len(results)):
        prev_speedup = results[i-1]['speedup']
        curr_speedup = results[i]['speedup']
        improvement = (curr_speedup - prev_speedup) / prev_speedup
        
        if improvement < 0.1:  # Less than 10% improvement
            optimal_threads = results[i-1]['threads']
            break
        optimal_threads = results[i]['threads']
    
    # Plot results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    threads = [r['threads'] for r in results]
    times = [r['avg_time'] for r in results]
    speedups = [r['speedup'] for r in results]
    
    # Time vs threads
    ax1.errorbar(threads, times, yerr=[r['std_time'] for r in results], 
                 marker='o', capsize=5)
    ax1.axvline(optimal_threads, color='red', linestyle='--', alpha=0.7,
                label=f'Optimal: {optimal_threads} threads')
    ax1.set_xlabel('Number of Threads')
    ax1.set_ylabel('Alignment Time (seconds)')
    ax1.set_title('Performance vs Thread Count')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Speedup vs threads
    ax2.plot(threads, speedups, 'go-')
    ax2.axvline(optimal_threads, color='red', linestyle='--', alpha=0.7,
                label=f'Optimal: {optimal_threads} threads')
    ax2.axhline(optimal_threads, color='gray', linestyle=':', alpha=0.5, 
                label='Linear speedup')
    ax2.set_xlabel('Number of Threads')
    ax2.set_ylabel('Speedup (vs 1 thread)')
    ax2.set_title('Speedup vs Thread Count')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    print(f"\n🎯 Optimal thread count: {optimal_threads}")
    print(f"   Best time: {results[optimal_threads-1]['avg_time']:.3f}s")
    print(f"   Speedup: {results[optimal_threads-1]['speedup']:.2f}x")
    
    return optimal_threads, results

# Example usage
def create_test_sequences(n_seqs=10, length=2000):
    """Create test sequences for benchmarking."""
    base = "ATCGATCGATCG" * (length // 12)
    sequences = []
    
    for i in range(n_seqs):
        seq = list(base[:length])
        # Add 5% mutations
        n_mutations = length // 20
        positions = np.random.choice(length, n_mutations, replace=False)
        for pos in positions:
            seq[pos] = np.random.choice(['A', 'T', 'C', 'G'])
        sequences.append(''.join(seq))
    
    return sequences

# test_sequences = create_test_sequences(n_seqs=8, length=1000)
# optimal_threads, thread_results = find_optimal_threads(test_sequences)
```

### Thread Pool Management

```python
import kalign
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time

class KalignThreadPool:
    """Thread pool manager for Kalign operations."""
    
    def __init__(self, max_workers=None):
        self.max_workers = max_workers or os.cpu_count()
        self.executor = None
    
    def __enter__(self):
        self.executor = ThreadPoolExecutor(max_workers=self.max_workers)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.executor:
            self.executor.shutdown(wait=True)
    
    def align_multiple(self, sequence_groups, **kwargs):
        """Align multiple sequence groups in parallel."""
        
        def align_group(sequences):
            return kalign.align(sequences, **kwargs)
        
        if self.executor is None:
            raise RuntimeError("Thread pool not initialized. Use with statement.")
        
        futures = [self.executor.submit(align_group, seqs) for seqs in sequence_groups]
        results = [future.result() for future in futures]
        
        return results
    
    def batch_align_sequences(self, all_sequences, batch_size=100, **kwargs):
        """Align sequences in batches for memory efficiency."""
        
        batches = [all_sequences[i:i+batch_size] 
                  for i in range(0, len(all_sequences), batch_size)]
        
        print(f"Processing {len(all_sequences)} sequences in {len(batches)} batches")
        
        results = []
        for i, batch in enumerate(batches):
            print(f"Processing batch {i+1}/{len(batches)}...")
            aligned_batch = kalign.align(batch, **kwargs)
            results.extend(aligned_batch)
        
        return results

# Example usage
def parallel_alignment_example():
    """Example of parallel alignment processing."""
    
    # Create multiple sequence groups
    sequence_groups = []
    for group in range(5):
        group_sequences = [
            f"ATCGATCG{'A' * group}TCGATCG",
            f"ATCGATCG{'T' * group}TCGATCG", 
            f"ATCGATCG{'C' * group}TCGATCG"
        ]
        sequence_groups.append(group_sequences)
    
    print(f"Aligning {len(sequence_groups)} groups with thread pool...")
    
    start_time = time.time()
    
    with KalignThreadPool(max_workers=4) as pool:
        alignments = pool.align_multiple(sequence_groups, n_threads=1)
    
    end_time = time.time()
    
    print(f"Completed in {end_time - start_time:.2f} seconds")
    print(f"Processed {len(alignments)} alignments")
    
    return alignments

# alignments = parallel_alignment_example()
```

## Memory Management

### Memory-Efficient Processing

```python
import kalign
import psutil
import gc
import sys

def memory_efficient_alignment(sequences, chunk_size=1000, monitor_memory=True):
    """Perform memory-efficient alignment for large sequence sets."""
    
    def get_memory_usage():
        """Get current memory usage in MB."""
        process = psutil.Process()
        return process.memory_info().rss / 1024 / 1024
    
    if monitor_memory:
        initial_memory = get_memory_usage()
        print(f"Initial memory usage: {initial_memory:.1f} MB")
    
    total_sequences = len(sequences)
    aligned_results = []
    
    # Process in chunks
    for start_idx in range(0, total_sequences, chunk_size):
        end_idx = min(start_idx + chunk_size, total_sequences)
        chunk = sequences[start_idx:end_idx]
        
        if monitor_memory:
            chunk_start_memory = get_memory_usage()
        
        print(f"Processing chunk {start_idx//chunk_size + 1}: "
              f"sequences {start_idx+1}-{end_idx}")
        
        # Align chunk
        chunk_aligned = kalign.align(chunk, n_threads=4)
        aligned_results.extend(chunk_aligned)
        
        if monitor_memory:
            chunk_end_memory = get_memory_usage()
            memory_increase = chunk_end_memory - chunk_start_memory
            print(f"  Memory increase: {memory_increase:.1f} MB")
            
            # Force garbage collection to free memory
            gc.collect()
            post_gc_memory = get_memory_usage()
            memory_freed = chunk_end_memory - post_gc_memory
            print(f"  Memory freed by GC: {memory_freed:.1f} MB")
    
    if monitor_memory:
        final_memory = get_memory_usage()
        total_memory_increase = final_memory - initial_memory
        print(f"\nFinal memory usage: {final_memory:.1f} MB")
        print(f"Total memory increase: {total_memory_increase:.1f} MB")
    
    return aligned_results

def stream_alignment_from_file(input_file, output_file, chunk_size=1000):
    """Stream alignment processing for very large files."""
    
    print(f"Streaming alignment from {input_file} to {output_file}")
    
    # Read sequences in chunks
    sequences_buffer = []
    ids_buffer = []
    processed_count = 0
    
    with open(output_file, 'w') as output_handle:
        # Read sequences (simplified - would use proper FASTA parser)
        try:
            sequences, ids = kalign.io.read_sequences(input_file)
            total_sequences = len(sequences)
            
            print(f"Processing {total_sequences} sequences in chunks of {chunk_size}")
            
            for i in range(0, total_sequences, chunk_size):
                chunk_sequences = sequences[i:i+chunk_size]
                chunk_ids = ids[i:i+chunk_size]
                
                print(f"Aligning chunk {i//chunk_size + 1}...")
                
                # Align chunk
                aligned_chunk = kalign.align(chunk_sequences)
                
                # Write aligned sequences immediately
                for seq_id, aligned_seq in zip(chunk_ids, aligned_chunk):
                    output_handle.write(f">{seq_id}\n{aligned_seq}\n")
                
                processed_count += len(chunk_sequences)
                print(f"  Processed {processed_count}/{total_sequences} sequences")
                
                # Clear memory
                del aligned_chunk
                gc.collect()
        
        except Exception as e:
            print(f"Error processing file: {e}")
            return False
    
    print(f"✅ Stream processing complete: {processed_count} sequences")
    return True

# Example: Monitor memory usage
def memory_usage_example():
    """Example of memory usage monitoring."""
    
    # Create large sequence set
    print("Creating large sequence set...")
    base_seq = "ATCGATCGATCG" * 100  # 1200 bp sequences
    large_sequence_set = []
    
    for i in range(500):  # 500 sequences
        seq = list(base_seq)
        # Add mutations
        n_mutations = len(seq) // 50  # 2% mutations
        positions = np.random.choice(len(seq), n_mutations, replace=False)
        for pos in positions:
            seq[pos] = np.random.choice(['A', 'T', 'C', 'G'])
        large_sequence_set.append(''.join(seq))
    
    print(f"Created {len(large_sequence_set)} sequences of length {len(large_sequence_set[0])}")
    
    # Process with memory monitoring
    aligned = memory_efficient_alignment(
        large_sequence_set, 
        chunk_size=100, 
        monitor_memory=True
    )
    
    print(f"Final result: {len(aligned)} aligned sequences")
    return aligned

# aligned = memory_usage_example()
```

## Algorithm Parameter Tuning

### Gap Penalty Optimization

```python
import kalign
import numpy as np
import matplotlib.pyplot as plt

def optimize_gap_penalties(sequences, penalty_ranges=None):
    """Find optimal gap penalties for given sequences."""
    
    if penalty_ranges is None:
        penalty_ranges = {
            'gap_open': np.arange(-20, -2, 2),
            'gap_extend': np.arange(-5, -0.5, 0.5)
        }
    
    print("🔧 Optimizing gap penalties...")
    print(f"Testing {len(penalty_ranges['gap_open'])} × {len(penalty_ranges['gap_extend'])} combinations")
    
    results = []
    best_score = -np.inf
    best_params = None
    
    for gap_open in penalty_ranges['gap_open']:
        for gap_extend in penalty_ranges['gap_extend']:
            print(f"Testing gap_open={gap_open}, gap_extend={gap_extend}...", end=' ')
            
            try:
                # Align with current parameters
                aligned = kalign.align(
                    sequences,
                    gap_open=gap_open,
                    gap_extend=gap_extend
                )
                
                # Calculate quality metrics
                stats = kalign.utils.alignment_stats(aligned)
                identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
                
                # Scoring function (customize based on your needs)
                # Higher conservation and identity = better
                # Lower gap fraction = better (usually)
                score = (stats['conservation'] * 0.4 + 
                        stats['identity'] * 0.4 + 
                        (1 - stats['gap_fraction']) * 0.2)
                
                results.append({
                    'gap_open': gap_open,
                    'gap_extend': gap_extend,
                    'score': score,
                    'conservation': stats['conservation'],
                    'identity': stats['identity'],
                    'gap_fraction': stats['gap_fraction'],
                    'alignment_length': stats['length']
                })
                
                if score > best_score:
                    best_score = score
                    best_params = (gap_open, gap_extend)
                
                print(f"score={score:.3f}")
                
            except Exception as e:
                print(f"failed ({e})")
                continue
    
    if not results:
        print("❌ No successful alignments")
        return None, None
    
    # Visualize results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Create heatmaps
    gap_opens = sorted(set(r['gap_open'] for r in results))
    gap_extends = sorted(set(r['gap_extend'] for r in results))
    
    # Score heatmap
    score_matrix = np.zeros((len(gap_extends), len(gap_opens)))
    for r in results:
        i = gap_extends.index(r['gap_extend'])
        j = gap_opens.index(r['gap_open'])
        score_matrix[i, j] = r['score']
    
    im1 = axes[0, 0].imshow(score_matrix, cmap='RdYlGn', aspect='auto')
    axes[0, 0].set_title('Alignment Quality Score')
    axes[0, 0].set_xlabel('Gap Open Penalty')
    axes[0, 0].set_ylabel('Gap Extend Penalty')
    axes[0, 0].set_xticks(range(len(gap_opens)))
    axes[0, 0].set_xticklabels([f'{x:.0f}' for x in gap_opens])
    axes[0, 0].set_yticks(range(len(gap_extends)))
    axes[0, 0].set_yticklabels([f'{x:.1f}' for x in gap_extends])
    plt.colorbar(im1, ax=axes[0, 0])
    
    # Mark best parameters
    best_i = gap_extends.index(best_params[1])
    best_j = gap_opens.index(best_params[0])
    axes[0, 0].plot(best_j, best_i, 'r*', markersize=15, markeredgecolor='black')
    
    # Conservation heatmap
    conservation_matrix = np.zeros((len(gap_extends), len(gap_opens)))
    for r in results:
        i = gap_extends.index(r['gap_extend'])
        j = gap_opens.index(r['gap_open'])
        conservation_matrix[i, j] = r['conservation']
    
    im2 = axes[0, 1].imshow(conservation_matrix, cmap='Blues', aspect='auto')
    axes[0, 1].set_title('Conservation')
    axes[0, 1].set_xlabel('Gap Open Penalty')
    axes[0, 1].set_ylabel('Gap Extend Penalty')
    axes[0, 1].set_xticks(range(len(gap_opens)))
    axes[0, 1].set_xticklabels([f'{x:.0f}' for x in gap_opens])
    axes[0, 1].set_yticks(range(len(gap_extends)))
    axes[0, 1].set_yticklabels([f'{x:.1f}' for x in gap_extends])
    plt.colorbar(im2, ax=axes[0, 1])
    
    # Gap fraction heatmap
    gap_matrix = np.zeros((len(gap_extends), len(gap_opens)))
    for r in results:
        i = gap_extends.index(r['gap_extend'])
        j = gap_opens.index(r['gap_open'])
        gap_matrix[i, j] = r['gap_fraction']
    
    im3 = axes[1, 0].imshow(gap_matrix, cmap='Reds', aspect='auto')
    axes[1, 0].set_title('Gap Fraction')
    axes[1, 0].set_xlabel('Gap Open Penalty')
    axes[1, 0].set_ylabel('Gap Extend Penalty')
    axes[1, 0].set_xticks(range(len(gap_opens)))
    axes[1, 0].set_xticklabels([f'{x:.0f}' for x in gap_opens])
    axes[1, 0].set_yticks(range(len(gap_extends)))
    axes[1, 0].set_yticklabels([f'{x:.1f}' for x in gap_extends])
    plt.colorbar(im3, ax=axes[1, 0])
    
    # Parameter correlation plot
    scores = [r['score'] for r in results]
    gap_open_values = [r['gap_open'] for r in results]
    
    scatter = axes[1, 1].scatter(gap_open_values, scores, 
                                c=[r['gap_extend'] for r in results], 
                                cmap='viridis', alpha=0.7)
    axes[1, 1].set_xlabel('Gap Open Penalty')
    axes[1, 1].set_ylabel('Alignment Quality Score')
    axes[1, 1].set_title('Score vs Gap Open (colored by Gap Extend)')
    plt.colorbar(scatter, ax=axes[1, 1], label='Gap Extend Penalty')
    
    # Mark best point
    best_result = next(r for r in results if 
                      r['gap_open'] == best_params[0] and 
                      r['gap_extend'] == best_params[1])
    axes[1, 1].plot(best_params[0], best_result['score'], 
                   'r*', markersize=15, markeredgecolor='black')
    
    plt.tight_layout()
    plt.show()
    
    print(f"\n🎯 Optimal parameters:")
    print(f"   Gap open: {best_params[0]}")
    print(f"   Gap extend: {best_params[1]}")
    print(f"   Score: {best_score:.3f}")
    print(f"   Conservation: {best_result['conservation']:.3f}")
    print(f"   Identity: {best_result['identity']:.3f}")
    print(f"   Gap fraction: {best_result['gap_fraction']:.3f}")
    
    return best_params, results

# Example usage
test_sequences = [
    "ATCGATCGATCGATCGATCG",
    "ATCGATCGTCGATCGATCG",
    "ATCGATCGATCATCGATCG",
    "ATCGATCGAGCTCGATCG"
]

# best_params, param_results = optimize_gap_penalties(test_sequences)
```

## Batch Processing

### High-Throughput Pipeline

```python
import kalign
import os
import json
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

class KalignBatchProcessor:
    """High-performance batch processor for Kalign."""
    
    def __init__(self, n_workers=None, temp_dir="./kalign_temp"):
        self.n_workers = n_workers or max(1, mp.cpu_count() - 1)
        self.temp_dir = Path(temp_dir)
        self.temp_dir.mkdir(exist_ok=True)
        
        # Performance tracking
        self.stats = {
            'total_files': 0,
            'successful_alignments': 0,
            'failed_alignments': 0,
            'total_sequences': 0,
            'total_time': 0
        }
    
    def process_file(self, input_file, output_dir=None, **kalign_kwargs):
        """Process a single file."""
        
        input_path = Path(input_file)
        if output_dir is None:
            output_dir = input_path.parent / "aligned"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        output_file = output_dir / f"aligned_{input_path.name}"
        
        try:
            start_time = time.time()
            
            # Read sequences
            sequences, ids = kalign.io.read_sequences(str(input_path))
            
            if len(sequences) < 2:
                return {
                    'input_file': str(input_path),
                    'output_file': None,
                    'status': 'skipped',
                    'message': 'Less than 2 sequences',
                    'n_sequences': len(sequences),
                    'time': 0
                }
            
            # Align sequences
            aligned = kalign.align(sequences, **kalign_kwargs)
            
            # Write output
            kalign.io.write_fasta(aligned, str(output_file), ids=ids)
            
            end_time = time.time()
            
            return {
                'input_file': str(input_path),
                'output_file': str(output_file),
                'status': 'success',
                'message': 'Alignment completed',
                'n_sequences': len(sequences),
                'alignment_length': len(aligned[0]),
                'time': end_time - start_time
            }
            
        except Exception as e:
            return {
                'input_file': str(input_path),
                'output_file': None,
                'status': 'error',
                'message': str(e),
                'n_sequences': 0,
                'time': 0
            }
    
    def process_directory(self, input_dir, output_dir=None, 
                         file_pattern="*.fasta", **kalign_kwargs):
        """Process all files in a directory."""
        
        input_path = Path(input_dir)
        if output_dir is None:
            output_dir = input_path / "aligned"
        
        # Find input files
        input_files = list(input_path.glob(file_pattern))
        
        if not input_files:
            print(f"❌ No files found matching pattern '{file_pattern}' in {input_dir}")
            return []
        
        print(f"🚀 Processing {len(input_files)} files with {self.n_workers} workers")
        
        start_time = time.time()
        results = []
        
        # Process files in parallel
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit all jobs
            futures = []
            for input_file in input_files:
                future = executor.submit(
                    self.process_file, 
                    input_file, 
                    output_dir, 
                    **kalign_kwargs
                )
                futures.append(future)
            
            # Collect results
            for i, future in enumerate(futures):
                try:
                    result = future.result(timeout=300)  # 5-minute timeout
                    results.append(result)
                    
                    # Update progress
                    if (i + 1) % 10 == 0 or i == len(futures) - 1:
                        success_count = sum(1 for r in results if r['status'] == 'success')
                        print(f"  Progress: {i+1}/{len(futures)} files, "
                              f"{success_count} successful")
                    
                except Exception as e:
                    results.append({
                        'input_file': str(input_files[i]),
                        'output_file': None,
                        'status': 'timeout',
                        'message': str(e),
                        'n_sequences': 0,
                        'time': 0
                    })
        
        end_time = time.time()
        
        # Update statistics
        self.stats['total_files'] += len(input_files)
        self.stats['successful_alignments'] += sum(1 for r in results if r['status'] == 'success')
        self.stats['failed_alignments'] += sum(1 for r in results if r['status'] != 'success')
        self.stats['total_sequences'] += sum(r.get('n_sequences', 0) for r in results)
        self.stats['total_time'] += end_time - start_time
        
        # Print summary
        self.print_summary(results, end_time - start_time)
        
        return results
    
    def print_summary(self, results, total_time):
        """Print processing summary."""
        
        successful = [r for r in results if r['status'] == 'success']
        failed = [r for r in results if r['status'] != 'success']
        
        print(f"\n📊 Batch Processing Summary:")
        print(f"   Total files: {len(results)}")
        print(f"   Successful: {len(successful)}")
        print(f"   Failed: {len(failed)}")
        print(f"   Total time: {total_time:.2f} seconds")
        print(f"   Average time per file: {total_time/len(results):.2f} seconds")
        
        if successful:
            total_sequences = sum(r['n_sequences'] for r in successful)
            total_alignment_time = sum(r['time'] for r in successful)
            print(f"   Total sequences aligned: {total_sequences}")
            print(f"   Average sequences per second: {total_sequences/total_alignment_time:.1f}")
        
        if failed:
            print(f"\n❌ Failed files:")
            for result in failed[:5]:  # Show first 5 failures
                print(f"   {result['input_file']}: {result['message']}")
            if len(failed) > 5:
                print(f"   ... and {len(failed) - 5} more")
    
    def save_report(self, results, output_file="kalign_batch_report.json"):
        """Save detailed processing report."""
        
        report = {
            'processing_stats': self.stats,
            'timestamp': time.time(),
            'results': results
        }
        
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"📄 Report saved to {output_file}")

# Example usage
def batch_processing_example():
    """Example of batch processing workflow."""
    
    # Create test files (in real usage, you'd have existing files)
    test_dir = Path("./test_sequences")
    test_dir.mkdir(exist_ok=True)
    
    # Create sample files
    for i in range(5):
        sequences = [
            f"ATCGATCG{'A' * i}TCGATCG",
            f"ATCGATCG{'T' * i}TCGATCG",
            f"ATCGATCG{'C' * i}TCGATCG"
        ]
        
        with open(test_dir / f"sequences_{i}.fasta", 'w') as f:
            for j, seq in enumerate(sequences):
                f.write(f">seq_{i}_{j}\n{seq}\n")
    
    # Process with batch processor
    processor = KalignBatchProcessor(n_workers=2)
    
    results = processor.process_directory(
        test_dir,
        output_dir=test_dir / "aligned",
        seq_type="dna",
        n_threads=1  # Use 1 thread per process since we're using multiple processes
    )
    
    # Save report
    processor.save_report(results, "batch_report.json")
    
    return results

# results = batch_processing_example()
```

## Large-Scale Workflows

### Distributed Processing

```python
import kalign
import pickle
import redis
from celery import Celery
import uuid

# Configure Celery for distributed processing
app = Celery('kalign_tasks', broker='redis://localhost:6379')

@app.task
def align_sequences_task(sequences, task_id=None, **kwargs):
    """Celery task for sequence alignment."""
    
    try:
        result = kalign.align(sequences, **kwargs)
        
        return {
            'task_id': task_id,
            'status': 'success',
            'result': result,
            'n_sequences': len(sequences),
            'alignment_length': len(result[0]) if result else 0
        }
    
    except Exception as e:
        return {
            'task_id': task_id,
            'status': 'error',
            'error': str(e),
            'n_sequences': len(sequences)
        }

class DistributedKalign:
    """Distributed Kalign processing using Celery."""
    
    def __init__(self, redis_host='localhost', redis_port=6379):
        self.redis_client = redis.Redis(host=redis_host, port=redis_port, decode_responses=True)
        self.pending_tasks = {}
    
    def submit_alignment(self, sequences, **kwargs):
        """Submit alignment task to distributed queue."""
        
        task_id = str(uuid.uuid4())
        
        # Submit task
        async_result = align_sequences_task.delay(sequences, task_id=task_id, **kwargs)
        
        # Track task
        self.pending_tasks[task_id] = {
            'async_result': async_result,
            'submitted_time': time.time(),
            'n_sequences': len(sequences)
        }
        
        print(f"Submitted task {task_id} with {len(sequences)} sequences")
        return task_id
    
    def get_result(self, task_id, timeout=None):
        """Get result for a specific task."""
        
        if task_id not in self.pending_tasks:
            raise ValueError(f"Task {task_id} not found")
        
        async_result = self.pending_tasks[task_id]['async_result']
        
        try:
            result = async_result.get(timeout=timeout)
            del self.pending_tasks[task_id]  # Clean up
            return result
        
        except Exception as e:
            return {
                'task_id': task_id,
                'status': 'error',
                'error': str(e)
            }
    
    def get_all_results(self, timeout=300):
        """Get results for all pending tasks."""
        
        results = []
        completed_tasks = []
        
        for task_id in self.pending_tasks:
            try:
                result = self.get_result(task_id, timeout=1)  # Short timeout
                results.append(result)
                completed_tasks.append(task_id)
            
            except:
                continue  # Task still pending
        
        print(f"Retrieved {len(results)} completed tasks")
        return results
    
    def submit_batch(self, sequence_groups, **kwargs):
        """Submit multiple alignment tasks."""
        
        task_ids = []
        
        for i, sequences in enumerate(sequence_groups):
            task_id = self.submit_alignment(sequences, **kwargs)
            task_ids.append(task_id)
        
        print(f"Submitted {len(task_ids)} batch tasks")
        return task_ids
    
    def wait_for_completion(self, task_ids=None, timeout=300, poll_interval=5):
        """Wait for tasks to complete."""
        
        if task_ids is None:
            task_ids = list(self.pending_tasks.keys())
        
        start_time = time.time()
        results = []
        
        while task_ids and (time.time() - start_time) < timeout:
            completed_this_round = []
            
            for task_id in task_ids:
                try:
                    result = self.get_result(task_id, timeout=1)
                    results.append(result)
                    completed_this_round.append(task_id)
                    print(f"Task {task_id} completed: {result['status']}")
                
                except:
                    continue  # Still pending
            
            # Remove completed tasks
            for task_id in completed_this_round:
                task_ids.remove(task_id)
            
            if task_ids:
                print(f"Waiting for {len(task_ids)} tasks...")
                time.sleep(poll_interval)
        
        if task_ids:
            print(f"Timeout: {len(task_ids)} tasks still pending")
        
        return results

# Example distributed workflow
def distributed_alignment_example():
    """Example of distributed alignment processing."""
    
    # Initialize distributed processor
    distributed = DistributedKalign()
    
    # Create sequence groups for processing
    sequence_groups = []
    for i in range(10):
        group = [
            f"ATCGATCG{'A' * (i % 5)}TCGATCG",
            f"ATCGATCG{'T' * (i % 5)}TCGATCG",
            f"ATCGATCG{'C' * (i % 5)}TCGATCG",
            f"ATCGATCG{'G' * (i % 5)}TCGATCG"
        ]
        sequence_groups.append(group)
    
    # Submit batch
    task_ids = distributed.submit_batch(
        sequence_groups,
        seq_type="dna",
        n_threads=2
    )
    
    # Wait for completion
    results = distributed.wait_for_completion(task_ids, timeout=600)
    
    # Analyze results
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'error']
    
    print(f"\nDistributed processing complete:")
    print(f"  Successful: {len(successful)}")
    print(f"  Failed: {len(failed)}")
    
    return results

# To run distributed processing:
# 1. Start Redis: redis-server
# 2. Start Celery workers: celery -A your_module worker --loglevel=info
# 3. Run: results = distributed_alignment_example()
```

This performance tuning guide provides comprehensive strategies for optimizing Kalign across different scales and use cases. From finding optimal thread counts to implementing distributed processing systems, users can now extract maximum performance from Kalign for their specific needs.

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"content": "Create comprehensive Python API documentation (docs/python-api.md)", "status": "completed", "priority": "high", "id": "1"}, {"content": "Create quick start guide with examples (docs/python-quickstart.md)", "status": "completed", "priority": "high", "id": "2"}, {"content": "Create ecosystem integration guide (docs/python-ecosystem.md)", "status": "completed", "priority": "high", "id": "3"}, {"content": "Create performance tuning guide (docs/python-performance.md)", "status": "completed", "priority": "medium", "id": "4"}, {"content": "Create troubleshooting guide (docs/python-troubleshooting.md)", "status": "in_progress", "priority": "medium", "id": "5"}, {"content": "Update main README with Python section", "status": "pending", "priority": "high", "id": "6"}, {"content": "Create example scripts directory", "status": "pending", "priority": "medium", "id": "7"}]