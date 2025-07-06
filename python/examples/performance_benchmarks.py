#!/usr/bin/env python3
"""
Performance Benchmarking Examples

Demonstrates performance optimization techniques and benchmarking tools.
"""

import kalign
import time
import os
import sys
import numpy as np

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    print("⚠️  psutil not available - memory monitoring disabled")

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("⚠️  matplotlib not available - plotting disabled")


def generate_test_sequences(n_sequences, seq_length, mutation_rate=0.05):
    """Generate test sequences for benchmarking."""
    
    # Create base sequence
    base_seq = "ATCG" * (seq_length // 4)
    if len(base_seq) < seq_length:
        base_seq += "ATCG"[:seq_length - len(base_seq)]
    
    sequences = []
    np.random.seed(42)  # For reproducible results
    
    for i in range(n_sequences):
        seq = list(base_seq)
        
        # Add mutations
        n_mutations = int(seq_length * mutation_rate)
        mutation_positions = np.random.choice(seq_length, n_mutations, replace=False)
        
        for pos in mutation_positions:
            seq[pos] = np.random.choice(['A', 'T', 'C', 'G'])
        
        sequences.append(''.join(seq))
    
    return sequences


def benchmark_thread_scaling():
    """Benchmark how performance scales with thread count."""
    print("=" * 60)
    print("Thread Scaling Benchmark")
    print("=" * 60)
    
    # Create test data
    sequences = generate_test_sequences(n_sequences=15, seq_length=1000)
    print(f"Test data: {len(sequences)} sequences of {len(sequences[0])} bp")
    
    # Test different thread counts
    max_threads = min(16, os.cpu_count())
    thread_counts = [1, 2, 4, 8] if max_threads >= 8 else list(range(1, max_threads + 1))
    
    results = []
    
    for n_threads in thread_counts:
        print(f"\nTesting {n_threads} threads...")
        
        # Run multiple iterations for accuracy
        times = []
        for iteration in range(3):
            print(f"  Iteration {iteration + 1}/3...", end=" ")
            
            start_time = time.time()
            aligned = kalign.align(sequences, n_threads=n_threads)
            end_time = time.time()
            
            duration = end_time - start_time
            times.append(duration)
            print(f"{duration:.3f}s")
        
        avg_time = np.mean(times)
        std_time = np.std(times)
        
        # Calculate speedup vs single thread
        if n_threads == 1:
            baseline_time = avg_time
            speedup = 1.0
        else:
            speedup = baseline_time / avg_time
        
        efficiency = speedup / n_threads
        
        results.append({
            'threads': n_threads,
            'avg_time': avg_time,
            'std_time': std_time,
            'speedup': speedup,
            'efficiency': efficiency
        })
        
        print(f"  Average: {avg_time:.3f}s (±{std_time:.3f}s)")
        print(f"  Speedup: {speedup:.2f}x")
        print(f"  Efficiency: {efficiency:.2%}")
    
    # Display results table
    print(f"\n{'Threads':<8} {'Time (s)':<10} {'Speedup':<8} {'Efficiency':<10}")
    print("-" * 36)
    for r in results:
        print(f"{r['threads']:<8} {r['avg_time']:<10.3f} {r['speedup']:<8.2f} {r['efficiency']:<10.2%}")
    
    # Find optimal thread count
    optimal = max(results, key=lambda x: x['efficiency'] if x['efficiency'] > 0.7 else 0)
    print(f"\n🎯 Recommended thread count: {optimal['threads']} "
          f"(efficiency: {optimal['efficiency']:.1%})")
    
    # Plot results if matplotlib available
    if MATPLOTLIB_AVAILABLE:
        plot_thread_scaling(results)
    
    return results


def plot_thread_scaling(results):
    """Plot thread scaling results."""
    
    threads = [r['threads'] for r in results]
    times = [r['avg_time'] for r in results]
    speedups = [r['speedup'] for r in results]
    efficiencies = [r['efficiency'] for r in results]
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot 1: Time vs threads
    ax1.plot(threads, times, 'bo-', linewidth=2, markersize=8)
    ax1.set_xlabel('Number of Threads')
    ax1.set_ylabel('Alignment Time (seconds)')
    ax1.set_title('Performance vs Thread Count')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Speedup vs threads
    ax2.plot(threads, speedups, 'go-', linewidth=2, markersize=8, label='Actual')
    ax2.plot(threads, threads, 'r--', alpha=0.7, label='Linear speedup')
    ax2.set_xlabel('Number of Threads')
    ax2.set_ylabel('Speedup (vs 1 thread)')
    ax2.set_title('Speedup vs Thread Count')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Efficiency vs threads
    ax3.plot(threads, [e * 100 for e in efficiencies], 'mo-', linewidth=2, markersize=8)
    ax3.axhline(100, color='r', linestyle='--', alpha=0.7, label='100% efficiency')
    ax3.set_xlabel('Number of Threads')
    ax3.set_ylabel('Efficiency (%)')
    ax3.set_title('Efficiency vs Thread Count')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('thread_scaling_benchmark.png', dpi=150, bbox_inches='tight')
    print("📊 Saved thread scaling plot: thread_scaling_benchmark.png")
    plt.close()


def benchmark_sequence_scaling():
    """Benchmark how performance scales with number of sequences."""
    print("\n" + "=" * 60)
    print("Sequence Count Scaling Benchmark")
    print("=" * 60)
    
    # Test different sequence counts
    seq_counts = [5, 10, 15, 20, 25, 30]
    seq_length = 500
    
    print(f"Testing sequence counts: {seq_counts}")
    print(f"Sequence length: {seq_length} bp")
    
    results = []
    
    for n_seqs in seq_counts:
        print(f"\nTesting {n_seqs} sequences...")
        
        # Generate test sequences
        sequences = generate_test_sequences(n_seqs, seq_length)
        
        # Benchmark alignment
        start_time = time.time()
        aligned = kalign.align(sequences, n_threads=4)
        end_time = time.time()
        
        duration = end_time - start_time
        throughput = n_seqs / duration
        
        results.append({
            'n_sequences': n_seqs,
            'time': duration,
            'throughput': throughput,
            'alignment_length': len(aligned[0])
        })
        
        print(f"  Time: {duration:.3f}s")
        print(f"  Throughput: {throughput:.1f} sequences/second")
        print(f"  Alignment length: {len(aligned[0])}")
    
    # Display results
    print(f"\n{'Sequences':<10} {'Time (s)':<10} {'Throughput':<15} {'Length':<8}")
    print("-" * 43)
    for r in results:
        print(f"{r['n_sequences']:<10} {r['time']:<10.3f} "
              f"{r['throughput']:<15.1f} {r['alignment_length']:<8}")
    
    # Analyze scaling
    if len(results) >= 3:
        # Simple complexity analysis
        times = [r['time'] for r in results]
        seq_counts_array = np.array([r['n_sequences'] for r in results])
        
        # Fit to n^2 model (typical for MSA)
        coeffs = np.polyfit(seq_counts_array**2, times, 1)
        r_squared = np.corrcoef(seq_counts_array**2, times)[0, 1]**2
        
        print(f"\n📈 Complexity analysis:")
        print(f"   Correlation with O(n²): R² = {r_squared:.3f}")
        if r_squared > 0.8:
            print("   ✅ Performance scales approximately quadratically")
        else:
            print("   ⚠️  Performance scaling is non-standard")
    
    return results


def benchmark_memory_usage():
    """Benchmark memory usage during alignment."""
    print("\n" + "=" * 60)
    print("Memory Usage Benchmark")
    print("=" * 60)
    
    if not PSUTIL_AVAILABLE:
        print("❌ psutil not available - skipping memory benchmark")
        return
    
    def get_memory_mb():
        process = psutil.Process()
        return process.memory_info().rss / 1024**2
    
    # Test different dataset sizes
    test_cases = [
        (10, 500),   # Small
        (20, 1000),  # Medium
        (30, 1500),  # Large
    ]
    
    results = []
    
    for n_seqs, seq_length in test_cases:
        print(f"\nTesting {n_seqs} sequences of {seq_length} bp...")
        
        # Measure initial memory
        initial_memory = get_memory_mb()
        
        # Generate sequences
        sequences = generate_test_sequences(n_seqs, seq_length)
        after_generation = get_memory_mb()
        
        # Perform alignment
        start_time = time.time()
        aligned = kalign.align(sequences, n_threads=4)
        end_time = time.time()
        after_alignment = get_memory_mb()
        
        # Calculate memory usage
        generation_memory = after_generation - initial_memory
        alignment_memory = after_alignment - after_generation
        total_memory = after_alignment - initial_memory
        
        duration = end_time - start_time
        
        results.append({
            'n_sequences': n_seqs,
            'seq_length': seq_length,
            'time': duration,
            'generation_memory': generation_memory,
            'alignment_memory': alignment_memory,
            'total_memory': total_memory,
            'peak_memory': after_alignment
        })
        
        print(f"  Time: {duration:.3f}s")
        print(f"  Generation memory: {generation_memory:.1f} MB")
        print(f"  Alignment memory: {alignment_memory:.1f} MB")
        print(f"  Total memory: {total_memory:.1f} MB")
        print(f"  Peak memory: {after_alignment:.1f} MB")
    
    # Display memory efficiency
    print(f"\n{'Case':<15} {'Time (s)':<10} {'Total Mem (MB)':<15} {'Mem/Seq (MB)':<15}")
    print("-" * 55)
    for i, r in enumerate(results):
        case_name = f"{r['n_sequences']}x{r['seq_length']}"
        mem_per_seq = r['total_memory'] / r['n_sequences']
        print(f"{case_name:<15} {r['time']:<10.3f} {r['total_memory']:<15.1f} {mem_per_seq:<15.2f}")
    
    return results


def benchmark_parameter_optimization():
    """Benchmark different parameter settings."""
    print("\n" + "=" * 60)
    print("Parameter Optimization Benchmark")
    print("=" * 60)
    
    # Test sequences
    sequences = generate_test_sequences(n_sequences=10, seq_length=800)
    print(f"Test data: {len(sequences)} sequences of {len(sequences[0])} bp")
    
    # Parameter combinations to test
    parameter_sets = [
        {"name": "Default", "gap_open": -1.0, "gap_extend": -1.0},
        {"name": "Conservative", "gap_open": -15.0, "gap_extend": -2.0},
        {"name": "Aggressive", "gap_open": -5.0, "gap_extend": -0.5},
        {"name": "Custom", "gap_open": -10.0, "gap_extend": -1.0},
    ]
    
    results = []
    
    for params in parameter_sets:
        print(f"\nTesting {params['name']} parameters...")
        print(f"  Gap open: {params['gap_open']}")
        print(f"  Gap extend: {params['gap_extend']}")
        
        # Benchmark alignment
        start_time = time.time()
        aligned = kalign.align(
            sequences,
            gap_open=params['gap_open'],
            gap_extend=params['gap_extend'],
            n_threads=4
        )
        end_time = time.time()
        
        # Calculate quality metrics
        stats = kalign.utils.alignment_stats(aligned)
        identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
        mean_identity = np.mean(identity_matrix[np.triu_indices_from(identity_matrix, k=1)])
        
        duration = end_time - start_time
        
        # Quality score (customize based on your needs)
        quality_score = (stats['conservation'] * 0.4 + 
                        mean_identity * 0.4 + 
                        (1 - stats['gap_fraction']) * 0.2)
        
        results.append({
            'name': params['name'],
            'gap_open': params['gap_open'],
            'gap_extend': params['gap_extend'],
            'time': duration,
            'alignment_length': stats['length'],
            'gap_fraction': stats['gap_fraction'],
            'conservation': stats['conservation'],
            'mean_identity': mean_identity,
            'quality_score': quality_score
        })
        
        print(f"  Time: {duration:.3f}s")
        print(f"  Length: {stats['length']}")
        print(f"  Gap fraction: {stats['gap_fraction']:.2%}")
        print(f"  Conservation: {stats['conservation']:.3f}")
        print(f"  Mean identity: {mean_identity:.3f}")
        print(f"  Quality score: {quality_score:.3f}")
    
    # Find best parameters
    best_quality = max(results, key=lambda x: x['quality_score'])
    best_speed = min(results, key=lambda x: x['time'])
    
    print(f"\n🏆 Best quality: {best_quality['name']} "
          f"(score: {best_quality['quality_score']:.3f})")
    print(f"🚀 Fastest: {best_speed['name']} "
          f"(time: {best_speed['time']:.3f}s)")
    
    return results


def run_comprehensive_benchmark():
    """Run all benchmarks and generate summary report."""
    print("🔬 Kalign Performance Benchmarking Suite")
    print("This will test various performance aspects of Kalign")
    
    # System information
    print(f"\n📊 System Information:")
    print(f"   Platform: {sys.platform}")
    print(f"   CPU cores: {os.cpu_count()}")
    if PSUTIL_AVAILABLE:
        print(f"   Memory: {psutil.virtual_memory().total / 1024**3:.1f} GB")
    print(f"   Kalign version: {kalign.__version__}")
    
    # Run benchmarks
    start_time = time.time()
    
    thread_results = benchmark_thread_scaling()
    sequence_results = benchmark_sequence_scaling()
    
    if PSUTIL_AVAILABLE:
        memory_results = benchmark_memory_usage()
    else:
        memory_results = None
    
    parameter_results = benchmark_parameter_optimization()
    
    total_time = time.time() - start_time
    
    # Generate summary report
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY REPORT")
    print("=" * 60)
    
    print(f"\nTotal benchmark time: {total_time:.1f} seconds")
    
    # Threading summary
    optimal_threads = max(thread_results, key=lambda x: x['efficiency'])
    print(f"\nOptimal threading: {optimal_threads['threads']} threads "
          f"({optimal_threads['efficiency']:.1%} efficiency)")
    
    # Performance summary
    largest_test = max(sequence_results, key=lambda x: x['n_sequences'])
    print(f"Largest test: {largest_test['n_sequences']} sequences in "
          f"{largest_test['time']:.2f}s "
          f"({largest_test['throughput']:.1f} seq/s)")
    
    # Memory summary
    if memory_results:
        max_memory = max(memory_results, key=lambda x: x['peak_memory'])
        print(f"Peak memory usage: {max_memory['peak_memory']:.1f} MB "
              f"({max_memory['n_sequences']} sequences)")
    
    # Parameter summary
    best_params = max(parameter_results, key=lambda x: x['quality_score'])
    print(f"Best parameters: {best_params['name']} "
          f"(quality score: {best_params['quality_score']:.3f})")
    
    print(f"\n✅ Benchmark suite completed successfully!")
    
    # Clean up plot files
    import glob
    for plot_file in glob.glob("*_benchmark.png"):
        if os.path.exists(plot_file):
            os.remove(plot_file)


def main():
    """Run performance benchmarks."""
    try:
        run_comprehensive_benchmark()
    except KeyboardInterrupt:
        print("\n⚠️  Benchmark interrupted by user")
    except Exception as e:
        print(f"\n❌ Benchmark failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()