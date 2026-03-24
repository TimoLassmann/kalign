"""Benchmark kalign mode presets (fast/default/recall/accurate) vs external tools.

Compares all 4 kalign NSGA-III optimized mode presets against ClustalO,
MAFFT, and MUSCLE. Reports quality metrics and wall time.
Each timing measurement is repeated N times and the median is reported.

Usage (inside container):
    python -m benchmarks.bench_modes --threads 16
    python -m benchmarks.bench_modes --threads 8 --runs 3 --output results.json
"""

import argparse
import json
import statistics
import time
from pathlib import Path

from .datasets import get_cases, download_dataset
from .scoring import run_case


KALIGN_MODES = ["fast", "default", "recall", "accurate"]


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark kalign modes vs external tools",
        prog="python -m benchmarks.bench_modes",
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Threads for all tools (default: 1)",
    )
    parser.add_argument(
        "--runs", type=int, default=3,
        help="Timing repeats per method/case; report median (default: 3)",
    )
    parser.add_argument(
        "--output", type=str, default="",
        help="Save results to JSON file",
    )
    parser.add_argument(
        "--dataset", default="balibase",
        help="Dataset to benchmark (default: balibase)",
    )
    args = parser.parse_args()

    download_dataset(args.dataset)
    cases = get_cases(args.dataset)
    if not cases:
        print(f"No cases found for '{args.dataset}'.")
        return

    print(f"Dataset: {args.dataset} ({len(cases)} cases)")
    print(f"Threads: {args.threads}, Timing repeats: {args.runs}")
    print()

    # Build method list: 4 kalign modes + 3 external tools
    methods = []
    for mode in KALIGN_MODES:
        methods.append((f"kalign {mode}", dict(method="python_api", mode=mode)))
    for tool in ["clustalo", "mafft", "muscle"]:
        methods.append((tool, dict(method=tool)))

    all_results = {}
    for label, kwargs in methods:
        print(f"Running {label}...", flush=True)

        # First run: collect quality metrics
        quality_results = []
        timing_runs = {c.family: [] for c in cases}
        t0 = time.perf_counter()

        for i, case in enumerate(cases):
            r = run_case(case, n_threads=args.threads, **kwargs)
            quality_results.append(r)
            timing_runs[case.family].append(r.wall_time)
            if r.error:
                print(f"  [{i+1}/{len(cases)}] {r.family}: ERROR {r.error}")

        # Additional timing repeats
        for rep in range(1, args.runs):
            for case in cases:
                r = run_case(case, n_threads=args.threads, **kwargs)
                if not r.error:
                    timing_runs[case.family].append(r.wall_time)

        elapsed = time.perf_counter() - t0
        ok = [r for r in quality_results if not r.error]

        if ok:
            # Median wall time per case, then sum
            median_times = []
            for r in ok:
                times = timing_runs[r.family]
                median_times.append(statistics.median(times) if times else r.wall_time)
            total_median_time = sum(median_times)

            rec = statistics.mean([r.recall for r in ok])
            pre = statistics.mean([r.precision for r in ok])
            f1 = statistics.mean([r.f1 for r in ok])
            tc = statistics.mean([r.tc for r in ok])

            all_results[label] = {
                "n_cases": len(ok),
                "recall": round(rec, 4),
                "precision": round(pre, 4),
                "f1": round(f1, 4),
                "tc": round(tc, 4),
                "median_total_time": round(total_median_time, 2),
                "elapsed": round(elapsed, 1),
            }
            print(f"  {len(ok)} cases, median total {total_median_time:.1f}s "
                  f"(elapsed {elapsed:.1f}s)")
        else:
            all_results[label] = {"n_cases": 0, "error": "all failed"}
            print(f"  All cases failed")
        print()

    # Summary table
    print(f"{'Method':<20} {'N':>4} {'Recall':>8} {'Prec':>8} {'F1':>8} "
          f"{'TC':>8} {'Time(s)':>10}")
    print("-" * 72)
    for label, _ in methods:
        s = all_results.get(label, {})
        if s.get("n_cases", 0) == 0:
            print(f"{label:<20} {'0':>4} {'n/a':>8} {'n/a':>8} {'n/a':>8} "
                  f"{'n/a':>8} {'n/a':>10}")
            continue
        print(f"{label:<20} {s['n_cases']:>4} {s['recall']:>8.3f} {s['precision']:>8.3f} "
              f"{s['f1']:>8.3f} {s['tc']:>8.3f} {s['median_total_time']:>10.1f}")

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as f:
            json.dump(all_results, f, indent=2)
        print(f"\nResults saved to {out}")


if __name__ == "__main__":
    main()
