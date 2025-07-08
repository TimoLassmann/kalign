#!/usr/bin/env python3

import subprocess
import time
import os
import sys
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich.columns import Columns
from rich import box

console = Console()


def build_with_thresholds(aln_threshold, kmeans_threshold):
    """Build kalign with specific threshold values"""
    build_dir = f"build_test_{aln_threshold}_{kmeans_threshold}"

    # Create build directory
    os.makedirs(build_dir, exist_ok=True)

    # Configure with specific thresholds
    cmd = [
        "cmake",
        f"-DKALIGN_ALN_SERIAL_THRESHOLD={aln_threshold}",
        f"-DKALIGN_KMEANS_UPGMA_THRESHOLD={kmeans_threshold}",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_PYTHON_MODULE=ON",
        "-DUSE_OPENMP=ON",  # Explicitly enable OpenMP
        "/Users/timo/code/kalign",
    ]

    try:
        result = subprocess.run(cmd, cwd=build_dir, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  CMake failed: {result.stderr}")
            return None

        result = subprocess.run(
            ["make", "-j8"], cwd=build_dir, capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"  Make failed: {result.stderr}")
            return None

        # Install the Python module
        python_dir = "/Users/timo/code/kalign/python"
        if os.path.exists(python_dir):
            result = subprocess.run(
                ["uv", "pip", "install", "-e", "."],
                cwd=python_dir,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                print(f"  Install failed: {result.stderr}")
                return None

        return build_dir

    except Exception as e:
        print(f"  Exception: {e}")
        return None


def run_benchmark(threads, problem_size="medium"):
    """Run benchmark with specific thread count and problem size"""
    python_dir = "/Users/timo/code/kalign/python"

    # Define different problem sizes - reduced sequence counts for speed
    if problem_size == "small":
        seqs_def = '["ATCGATCGATCGATCG" * 25] * 500'  # 500 seqs, 400bp each
    elif problem_size == "medium":
        # Medium: 1000 protein sequences, ~400aa each
        seqs_def = """
import random
random.seed(42)  # For reproducible sequences
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
base_protein = "".join(random.choices(amino_acids, k=400))  # 400 amino acids
seqs = []
for i in range(1000):  # 1000 sequences
    protein = list(base_protein)
    # Add 5% mutations for variation
    n_mutations = len(protein) // 20
    for _ in range(n_mutations):
        pos = random.randint(0, len(protein) - 1)
        protein[pos] = random.choice(amino_acids)
    seqs.append(''.join(protein))
seqs"""
    elif problem_size == "large":
        # Large: 2000 protein sequences, ~400aa each
        seqs_def = """
import random
random.seed(43)  # Different seed for variation
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
base_protein = "".join(random.choices(amino_acids, k=400))  # 400 amino acids
seqs = []
for i in range(2000):  # 2000 sequences
    protein = list(base_protein)
    # Add 5% mutations for variation
    n_mutations = len(protein) // 20
    for _ in range(n_mutations):
        pos = random.randint(0, len(protein) - 1)
        protein[pos] = random.choice(amino_acids)
    seqs.append(''.join(protein))
seqs"""
    elif problem_size == "long_dna":
        # NEW: Long DNA sequences test - fewer sequences but much longer (10kb each)
        seqs_def = """
import random
random.seed(42)  # For reproducible sequences
base_seq = "ATCGATCGATCGATCG" * 625  # 10,000 bp base sequence
seqs = []
for i in range(100):  # Only 100 sequences but 10kb each
    seq = list(base_seq)
    # Add 1% mutations for variation
    n_mutations = len(seq) // 100
    for _ in range(n_mutations):
        pos = random.randint(0, len(seq) - 1)
        seq[pos] = random.choice(['A', 'T', 'C', 'G'])
    seqs.append(''.join(seq))
seqs"""
    else:
        seqs_def = '["ATCGATCGATCGATCG" * 150] * 2000'  # 2000 seqs, 2400bp each

    # Run performance test
    if problem_size == "long_dna":
        test_script = f"""
import kalign
import time
import os
import random

# Set thread count multiple ways to ensure it takes effect
os.environ['OMP_NUM_THREADS'] = '{threads}'
os.environ['KALIGN_NUM_THREADS'] = '{threads}'

# Problem: {problem_size} - Generate long sequences
{seqs_def}

# Verify threading setup and sequence info
print(f"OMP_NUM_THREADS = {{os.environ.get('OMP_NUM_THREADS', 'not set')}}")
print(f"Generated {{len(seqs)}} sequences of {{len(seqs[0])}} bp each")

start = time.time()
result = kalign.align(seqs, n_threads={threads})
end = time.time()

print(f"Threads: {threads}, Size: {problem_size}, Time: {{end - start:.3f}}s")
"""
    elif problem_size in ["medium", "large"]:
        # Protein sequence test
        test_script = f"""
import kalign
import time
import os
import random

# Set thread count multiple ways to ensure it takes effect
os.environ['OMP_NUM_THREADS'] = '{threads}'
os.environ['KALIGN_NUM_THREADS'] = '{threads}'

# Problem: {problem_size} - Generate protein sequences
{seqs_def}

# Verify threading setup and sequence info
print(f"OMP_NUM_THREADS = {{os.environ.get('OMP_NUM_THREADS', 'not set')}}")
print(f"Generated {{len(seqs)}} protein sequences of {{len(seqs[0])}} aa each")

start = time.time()
result = kalign.align(seqs, seq_type="protein", n_threads={threads})
end = time.time()

print(f"Threads: {threads}, Size: {problem_size}, Time: {{end - start:.3f}}s")
"""
    else:
        # Default case (small)
        test_script = f"""
import kalign
import time
import os

# Set thread count multiple ways to ensure it takes effect
os.environ['OMP_NUM_THREADS'] = '{threads}'
os.environ['KALIGN_NUM_THREADS'] = '{threads}'

# Problem: {problem_size}
seqs = {seqs_def}

# Verify threading setup
print(f"OMP_NUM_THREADS = {{os.environ.get('OMP_NUM_THREADS', 'not set')}}")

start = time.time()
result = kalign.align(seqs, n_threads={threads})
end = time.time()

print(f"Threads: {threads}, Size: {problem_size}, Time: {{end - start:.3f}}s")
"""

    try:
        result = subprocess.run(
            ["uv", "run", "python", "-c", test_script],
            capture_output=True,
            text=True,
            cwd=python_dir,
        )

        if result.returncode == 0:
            # Extract timing from output
            for line in result.stdout.split("\n"):
                if "Time:" in line:
                    parts = line.split()
                    time_idx = parts.index("Time:") + 1
                    time_str = parts[time_idx].rstrip("s")
                    return float(time_str)
        else:
            print(f"  Benchmark failed: {result.stderr}")

    except Exception as e:
        print(f"  Exception: {e}")

    return None


def main():
    # Show configuration explanation
    console.print(
        Panel.fit(
            """
[bold blue]Kalign Threading Threshold Test (Optimized)[/bold blue]

[bold]What are these thresholds?[/bold]
• [yellow]aln_threshold[/yellow]: Below this many alignment positions, use serial instead of parallel processing
• [yellow]kmeans_threshold[/yellow]: Below this many sequences, use UPGMA instead of parallel k-means clustering

[bold]Test configurations:[/bold]
• [green]Original[/green]: aln=500, kmeans=100 (default values)
• [blue]Moderate[/blue]: aln=250, kmeans=50 (balanced parallelization)
• [cyan]Low[/cyan]: aln=50, kmeans=10 (more aggressive parallelization)
• [magenta]Always parallel[/magenta]: aln=0, kmeans=0 (maximum parallelization)

[bold]Test datasets (optimized for speed & coverage):[/bold]
• [blue]Medium[/blue]: 1,000 protein sequences × 400aa each
• [blue]Large[/blue]: 2,000 protein sequences × 400aa each  
• [yellow]Long DNA[/yellow]: 100 DNA sequences × 10,000bp each (NEW!)

[bold]Goal:[/bold] Find which threshold settings give best performance scaling beyond 8 threads.
""",
            box=box.ROUNDED,
        )
    )

    # Test different threshold combinations
    threshold_combinations = [
        (500, 100, "Original", "green"),
        (250, 50, "Moderate", "blue"),
        (50, 10, "Low", "cyan"),
        (0, 0, "Always parallel", "magenta"),
    ]

    # Test different thread counts
    thread_counts = [1, 2, 4, 8, 16]

    # Test different problem sizes - reduced counts for speed + new long DNA test
    problem_sizes = ["medium", "large", "long_dna"]

    all_results = []

    for aln_thresh, kmeans_thresh, label, color in threshold_combinations:
        console.print(f"\n[{color}]{'='*60}[/{color}]")
        console.print(
            f"[{color}]Testing {label} (aln={aln_thresh}, kmeans={kmeans_thresh})[/{color}]"
        )
        console.print(f"[{color}]{'='*60}[/{color}]")

        # Build once for this threshold combination
        with console.status(f"[{color}]Building {label} configuration...[/{color}]"):
            build_dir = build_with_thresholds(aln_thresh, kmeans_thresh)

        if build_dir is None:
            console.print(
                f"[red]Failed to build with thresholds {aln_thresh}, {kmeans_thresh}[/red]"
            )
            continue

        for problem_size in problem_sizes:
            console.print(f"\n[bold]Problem size: {problem_size}[/bold]")
            if problem_size == "medium":
                console.print("  🧬 1,000 protein sequences × 400aa each")
            elif problem_size == "large":
                console.print("  🧬 2,000 protein sequences × 400aa each")
            elif problem_size == "long_dna":
                console.print("  🧬 100 DNA sequences × 10,000bp each (long DNA test)")
            else:
                console.print(f"  📊 Problem size: {problem_size}")

            size_results = []
            for threads in thread_counts:
                with console.status(f"Testing {threads} threads..."):
                    runtime = run_benchmark(threads, problem_size)

                if runtime is not None:
                    size_results.append((threads, runtime))
                    all_results.append(
                        (
                            label,
                            aln_thresh,
                            kmeans_thresh,
                            problem_size,
                            threads,
                            runtime,
                            color,
                        )
                    )
                    console.print(f"  ✅ {threads} threads: {runtime:.3f}s")
                else:
                    console.print(f"  ❌ {threads} threads: FAILED")

            # Show scaling table for this problem size
            if size_results:
                table = Table(
                    title=f"Scaling Results: {problem_size} ({label})", box=box.SIMPLE
                )
                table.add_column("Threads", style="cyan", no_wrap=True)
                table.add_column("Time (s)", style="yellow")
                table.add_column("Speedup", style="green")
                table.add_column("Efficiency", style="blue")

                baseline = size_results[0][1]  # 1-thread time
                for threads, runtime in size_results:
                    speedup = baseline / runtime
                    efficiency = speedup / threads * 100

                    # Color code efficiency
                    if efficiency >= 80:
                        eff_style = "green"
                    elif efficiency >= 60:
                        eff_style = "yellow"
                    else:
                        eff_style = "red"

                    table.add_row(
                        str(threads),
                        f"{runtime:.3f}",
                        f"{speedup:.2f}x",
                        f"[{eff_style}]{efficiency:.1f}%[/{eff_style}]",
                    )

                console.print(table)

    # Overall comparison
    console.print("\n" + "=" * 80)
    console.print("[bold blue]PERFORMANCE COMPARISON[/bold blue]")
    console.print("=" * 80)

    # Create comparison table
    comparison_table = Table(title="Best Performance Summary", box=box.HEAVY_EDGE)
    comparison_table.add_column("Configuration", style="bold")
    comparison_table.add_column("Problem Size")
    comparison_table.add_column("Best Threads", style="cyan")
    comparison_table.add_column("Time (s)", style="yellow")
    comparison_table.add_column("Max Speedup", style="green")

    # Group by configuration and problem size
    configs = {}
    for (
        label,
        aln_thresh,
        kmeans_thresh,
        problem_size,
        threads,
        runtime,
        color,
    ) in all_results:
        key = (label, problem_size, color)
        if key not in configs:
            configs[key] = []
        configs[key].append((threads, runtime))

    # Find best speedup for each config
    for (label, problem_size, color), results in configs.items():
        results.sort()  # Sort by thread count
        baseline = results[0][1]  # 1-thread time

        # Find best speedup
        best_speedup = 0
        best_threads = 1
        best_time = baseline
        for threads, runtime in results:
            speedup = baseline / runtime
            if speedup > best_speedup:
                best_speedup = speedup
                best_threads = threads
                best_time = runtime

        comparison_table.add_row(
            f"[{color}]{label}[/{color}]",
            problem_size.title(),
            str(best_threads),
            f"{best_time:.3f}",
            f"{best_speedup:.2f}x",
        )

    console.print(comparison_table)

    # Threshold impact analysis
    console.print(f"\n[bold blue]THRESHOLD IMPACT AT HIGH THREAD COUNTS[/bold blue]")

    high_thread_results = {}
    for (
        label,
        aln_thresh,
        kmeans_thresh,
        problem_size,
        threads,
        runtime,
        color,
    ) in all_results:
        if threads >= 8:  # Look at 8 and 16 thread results
            key = (problem_size, threads)
            if key not in high_thread_results:
                high_thread_results[key] = {}
            high_thread_results[key][label] = (runtime, color)

    for (problem_size, threads), config_times in high_thread_results.items():
        if len(config_times) > 1:
            impact_table = Table(
                title=f"{problem_size.title()} Problem - {threads} Threads",
                box=box.SIMPLE,
            )
            impact_table.add_column("Configuration", style="bold")
            impact_table.add_column("Time (s)", style="yellow")
            impact_table.add_column("vs Fastest", style="red")

            sorted_configs = sorted(config_times.items(), key=lambda x: x[1][0])
            fastest_time = sorted_configs[0][1][0]

            for config, (time_val, color) in sorted_configs:
                if time_val == fastest_time:
                    vs_fastest = "🏆 Fastest"
                    vs_style = "green"
                else:
                    slowdown = (time_val / fastest_time - 1) * 100
                    vs_fastest = f"+{slowdown:.1f}% slower"
                    vs_style = "red"

                impact_table.add_row(
                    f"[{color}]{config}[/{color}]",
                    f"{time_val:.3f}",
                    f"[{vs_style}]{vs_fastest}[/{vs_style}]",
                )

            console.print(impact_table)

    # Final recommendations
    console.print(
        Panel.fit(
            """
[bold green]💡 RECOMMENDATIONS[/bold green]

Based on the results above:
1. Look for configurations that maintain good efficiency (>60%) at 16 threads
2. Compare 8-thread vs 16-thread performance to see if lower thresholds help scaling
3. The configuration with best 16-thread performance should be your optimal setting
4. [yellow]Pay special attention to Long DNA results[/yellow] - this tests alignment-heavy workloads

[yellow]Key insights:[/yellow] 
• If "Low" or "Always parallel" show better 16-thread performance than "Original", 
  the default thresholds are limiting your parallelization at high thread counts
• [blue]Long DNA test[/blue] reveals alignment threshold effects (aln_threshold)
• [blue]Protein tests (Medium/Large)[/blue] reveal clustering threshold effects (kmeans_threshold)
• [green]Different sequence types[/green] (protein vs DNA) may show different threshold sensitivities
""",
            box=box.ROUNDED,
            style="blue",
        )
    )


if __name__ == "__main__":
    main()
