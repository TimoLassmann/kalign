"""Optuna optimization of threadpool parallelization parameters.

Runs separate optimizations at each thread count (1, 2, 4, 8, 16, 32, 64)
to find optimal thresholds AND where thread scaling stops helping.

Usage:
    uv run python benchmarks/optimize_parallel.py --n-trials 100 --fresh
    uv run python benchmarks/optimize_parallel.py --n-trials 5 --dry-run
"""

import argparse
import faulthandler
import json
import os
import signal
import sys
import time
from pathlib import Path

import optuna

import kalign

# Enable faulthandler so segfaults print a traceback
faulthandler.enable()
if hasattr(signal, "SIGUSR1"):
    faulthandler.register(signal.SIGUSR1)

# ---------------------------------------------------------------------------
# Dataset generation
# ---------------------------------------------------------------------------

DSSIM_SPECS = [
    # (name, n_seq, length, dna, n_obs)
    ("prot_200", 200, 300, False, 50),
    ("prot_500", 500, 300, False, 50),
    ("prot_1000", 1000, 200, False, 50),
    ("dna_200", 200, 500, True, 50),
    ("dna_500", 500, 1000, True, 50),
    ("dna_1000", 1000, 500, True, 50),
    ("dna_viral", 200, 5000, True, 30),
    ("dna_viral_lg", 500, 3000, True, 30),
]

THREAD_COUNTS = [1, 2, 4, 8, 16, 32, 64]


def generate_dssim_datasets(cache_dir: Path) -> list[tuple[str, list[str], str]]:
    """Generate DSSim datasets, caching to disk."""
    datasets = []
    cache_dir.mkdir(parents=True, exist_ok=True)

    for name, n_seq, length, dna, n_obs in DSSIM_SPECS:
        cache_file = cache_dir / f"{name}.fasta"
        seq_type = "dna" if dna else "protein"

        if cache_file.exists():
            seqs = []
            current: list[str] = []
            for line in cache_file.read_text().splitlines():
                if line.startswith(">"):
                    if current:
                        seqs.append("".join(current))
                        current = []
                else:
                    current.append(line.strip())
            if current:
                seqs.append("".join(current))
            print(f"  {name}: {len(seqs)} seqs (cached)")
        else:
            print(f"  {name}: generating {n_seq} seqs, len={length}, "
                  f"{'DNA' if dna else 'protein'}...")
            seqs = kalign.generate_test_sequences(
                n_seq=n_seq, n_obs=n_obs, dna=dna, length=length, seed=42
            )
            with open(cache_file, "w") as f:
                for i, s in enumerate(seqs):
                    f.write(f">seq{i}\n{s}\n")
            print(f"  {name}: {len(seqs)} seqs generated")

        datasets.append((name, seqs, seq_type))

    return datasets


def load_balifam100() -> list[tuple[str, list[str], str]]:
    """Load BaliFam100 unaligned sequences."""
    try:
        try:
            from benchmarks.datasets import BALIFAM_DIR, balifam_download
        except ImportError:
            from datasets import BALIFAM_DIR, balifam_download

        if not BALIFAM_DIR.exists() or not any(BALIFAM_DIR.iterdir()):
            print("  Downloading BaliFam100...")
            balifam_download()

        in_dir = BALIFAM_DIR / "balifam100" / "in"
        if not in_dir.exists():
            in_dir = BALIFAM_DIR / "in"
        if not in_dir.exists():
            print("  WARNING: BaliFam100 in/ not found, skipping")
            return []

        datasets = []
        for fasta in sorted(in_dir.glob("*")):
            if not fasta.is_file():
                continue
            seqs = []
            current: list[str] = []
            for line in fasta.read_text().splitlines():
                if line.startswith(">"):
                    if current:
                        seqs.append("".join(current))
                        current = []
                else:
                    current.append(line.strip())
            if current:
                seqs.append("".join(current))
            if len(seqs) >= 10:
                datasets.append((f"bf100_{fasta.stem}", seqs, "protein"))

        print(f"  BaliFam100: {len(datasets)} families loaded")
        return datasets

    except Exception as e:
        print(f"  WARNING: Could not load BaliFam100: {e}")
        return []


# ---------------------------------------------------------------------------
# Benchmark runner
# ---------------------------------------------------------------------------


def time_alignment(seqs: list[str], seq_type: str, mode: str) -> float:
    start = time.perf_counter()
    kalign.align(seqs, seq_type=seq_type, mode=mode)
    return time.perf_counter() - start


def run_benchmark_suite(
    datasets: list[tuple[str, list[str], str]],
    mode: str,
    trial: optuna.trial.Trial | None = None,
) -> float:
    total = 0.0
    for i, (_name, seqs, seq_type) in enumerate(datasets):
        t = time_alignment(seqs, seq_type, mode)
        total += t
        if trial is not None:
            trial.report(total, i)
            if trial.should_prune():
                raise optuna.TrialPruned()
    return total


# ---------------------------------------------------------------------------
# Per-thread-count optimization
# ---------------------------------------------------------------------------


def optimize_for_thread_count(
    nt: int,
    datasets: list[tuple[str, list[str], str]],
    mode: str,
    n_trials: int,
    db_path: Path,
) -> dict:
    """Run Optuna optimization at a fixed thread count. Returns best result dict."""

    n_repeats = 3

    def objective(trial: optuna.trial.Trial) -> float:
        aln_serial = trial.suggest_int("aln_serial_threshold", 50, 1000, step=10)
        kmeans_upgma = trial.suggest_int("kmeans_upgma_threshold", 10, 150, step=5)
        dist_min = trial.suggest_int("dist_min_seqs", 0, 100, step=5)
        pfor_chunk = trial.suggest_int("pfor_min_chunk", 1, 32)

        sys.stderr.write(
            f"[{nt}T trial {trial.number}] aln={aln_serial} km={kmeans_upgma} "
            f"dist={dist_min} pfor={pfor_chunk}\n"
        )
        sys.stderr.flush()

        kalign.set_parallel_config(
            aln_serial_threshold=aln_serial,
            kmeans_upgma_threshold=kmeans_upgma,
            dist_min_seqs=dist_min,
            pfor_min_chunk=pfor_chunk,
        )
        kalign.set_num_threads(nt)

        try:
            times = []
            for r in range(n_repeats):
                t = run_benchmark_suite(datasets, mode, trial if r == 0 else None)
                times.append(t)
            times.sort()
            return times[len(times) // 2]
        except Exception:
            return float("inf")

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    storage = f"sqlite:///{db_path}"
    study_name = f"parallel_opt_{mode}_{nt}t"

    study = optuna.create_study(
        study_name=study_name,
        storage=storage,
        load_if_exists=True,
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=42),
        pruner=optuna.pruners.MedianPruner(
            n_startup_trials=5,
            n_warmup_steps=len(datasets) // 3,
        ),
    )

    existing = len(study.trials)
    remaining = max(0, n_trials - existing)

    if existing > 0 and remaining > 0:
        print(f"    Resuming ({existing} done, {remaining} remaining)")

    if remaining <= 0:
        print(f"    Already complete ({existing} trials)")
    else:
        t_start = time.perf_counter()
        n_done = [0]

        def callback(study: optuna.study.Study, trial: optuna.trial.FrozenTrial):
            n_done[0] += 1
            if trial.state == optuna.trial.TrialState.COMPLETE:
                best = study.best_trial
                is_new = trial.number == best.number
                elapsed = time.perf_counter() - t_start
                p = trial.params
                bp = best.params
                print(
                    f"    [{n_done[0]:>3}/{remaining}] "
                    f"{trial.value:>7.1f}s  "
                    f"aln={p['aln_serial_threshold']:<4} "
                    f"km={p['kmeans_upgma_threshold']:<4} "
                    f"dist={p['dist_min_seqs']:<4} "
                    f"pfor={p['pfor_min_chunk']:<3} "
                    f"| best={best.value:.1f}s "
                    f"(aln={bp['aln_serial_threshold']} km={bp['kmeans_upgma_threshold']} "
                    f"dist={bp['dist_min_seqs']} pfor={bp['pfor_min_chunk']})  "
                    f"({elapsed:.0f}s)"
                    f"{'  ***' if is_new else ''}"
                )
            elif trial.state == optuna.trial.TrialState.PRUNED:
                print(f"    [{n_done[0]:>3}/{remaining}] pruned")

        study.optimize(objective, n_trials=remaining, callbacks=[callback])

    best = study.best_trial
    return {
        "n_threads": nt,
        "best_time": best.value,
        "params": best.params,
        "trial_number": best.number,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Optimize threadpool parameters per thread count"
    )
    parser.add_argument(
        "--n-trials", type=int, default=100,
        help="Optuna trials per thread count (default: 100)"
    )
    parser.add_argument(
        "--mode", default="fast",
        choices=["fast", "default", "accurate"],
        help="Kalign mode preset (default: fast)"
    )
    parser.add_argument(
        "--no-balifam", action="store_true",
        help="Skip BaliFam100 datasets"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Just run baselines at each thread count, don't optimize"
    )
    parser.add_argument(
        "--cache-dir", type=str,
        default=str(Path(__file__).parent / "data" / "parallel_opt"),
        help="Directory for cached datasets and study DB"
    )
    parser.add_argument(
        "--fresh", action="store_true",
        help="Delete previous studies and start fresh"
    )
    parser.add_argument(
        "--thread-counts", type=str, default=None,
        help="Comma-separated thread counts (default: 1,2,4,8,16,32,64)"
    )
    args = parser.parse_args()

    if args.thread_counts:
        thread_counts = [int(x) for x in args.thread_counts.split(",")]
    else:
        max_cpu = os.cpu_count() or 64
        thread_counts = [t for t in THREAD_COUNTS if t <= max_cpu]

    cache_dir = Path(args.cache_dir)
    db_path = cache_dir / "parallel_opt.db"

    print(f"{'='*72}")
    print(f"  Threadpool Parameter Optimization")
    print(f"  Mode: {args.mode}  |  Trials/thread-count: {args.n_trials}")
    print(f"  Thread counts: {thread_counts}")
    print(f"{'='*72}")
    print()

    # Generate/load datasets
    print("Preparing datasets...")
    datasets = generate_dssim_datasets(cache_dir)

    if not args.no_balifam:
        bf = load_balifam100()
        datasets.extend(bf)

    print(f"\nTotal: {len(datasets)} benchmark cases")
    print()

    # Warm up
    print("Warming up...")
    warmup_seqs = kalign.generate_test_sequences(
        n_seq=50, n_obs=20, dna=False, length=100, seed=99
    )
    kalign.align(warmup_seqs, seq_type="protein", mode=args.mode)
    print()

    # Baseline at each thread count
    print("Running baselines (default params at each thread count)...")
    baselines = {}
    kalign.set_parallel_config()
    for nt in thread_counts:
        kalign.set_num_threads(nt)
        t = run_benchmark_suite(datasets, args.mode)
        baselines[nt] = t
        print(f"  {nt:>2} threads: {t:.2f}s")
    print()

    if args.dry_run:
        print("Dry run complete.")
        return

    # Fresh start
    if args.fresh and db_path.exists():
        print(f"--fresh: deleting {db_path}")
        db_path.unlink()
        print()

    # Optimize at each thread count
    results = []
    for nt in thread_counts:
        print(f"--- Optimizing for {nt} threads ({args.n_trials} trials) ---")
        result = optimize_for_thread_count(
            nt, datasets, args.mode, args.n_trials, db_path
        )
        result["baseline_time"] = baselines[nt]
        results.append(result)
        print()

    # Summary table
    print()
    print(f"{'='*90}")
    print(f"  SCALING & OPTIMIZATION RESULTS")
    print(f"{'='*90}")
    print()
    print(f"  {'Threads':>7}  {'Baseline':>10}  {'Optimized':>10}  {'Speedup':>8}  "
          f"{'aln_ser':>8}  {'km_upgma':>8}  {'dist_min':>8}  {'pfor_ch':>8}")
    print(f"  {'-'*7}  {'-'*10}  {'-'*10}  {'-'*8}  "
          f"{'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    best_overall = None
    for r in results:
        nt = r["n_threads"]
        bl = r["baseline_time"]
        opt = r["best_time"]
        speedup = bl / opt if opt > 0 else 0
        p = r["params"]
        is_best = best_overall is None or opt < best_overall["best_time"]
        if is_best:
            best_overall = r
        marker = "  <-- fastest" if is_best else ""
        print(
            f"  {nt:>7}  {bl:>9.2f}s  {opt:>9.2f}s  {speedup:>7.2f}x  "
            f"{p['aln_serial_threshold']:>8}  "
            f"{p['kmeans_upgma_threshold']:>8}  "
            f"{p['dist_min_seqs']:>8}  "
            f"{p['pfor_min_chunk']:>8}"
            f"{marker}"
        )

    print()
    print(f"  Fastest overall: {best_overall['n_threads']} threads, "
          f"{best_overall['best_time']:.2f}s")
    print()

    # Scaling analysis
    single = None
    for r in results:
        if r["n_threads"] == 1:
            single = r["best_time"]
            break
    if single:
        print("  Thread scaling (optimized):")
        for r in results:
            par_speedup = single / r["best_time"] if r["best_time"] > 0 else 0
            efficiency = par_speedup / r["n_threads"] * 100 if r["n_threads"] > 0 else 0
            bar = "#" * int(par_speedup * 2)
            print(f"  {r['n_threads']:>3}T: {par_speedup:>5.1f}x  "
                  f"({efficiency:>4.0f}% eff)  {bar}")
        print()

    # Save results
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out = {
        "mode": args.mode,
        "n_trials_per_thread_count": args.n_trials,
        "n_datasets": len(datasets),
        "results": results,
    }
    out_file = results_dir / "parallel_opt.json"
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Results saved to {out_file}")


if __name__ == "__main__":
    main()
