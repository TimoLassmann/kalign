"""Main benchmark orchestrator and CLI."""

import argparse
import json
import statistics
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

from .datasets import download_dataset, get_cases, DATASETS
from .scoring import AlignmentResult, EXTERNAL_TOOLS, run_case


def _run_one(args):
    """Worker function for parallel execution."""
    case, method, binary, n_threads, mode = args
    return run_case(case, method=method, binary=binary, n_threads=n_threads,
                    mode=mode)


def _result_label(r) -> str:
    """Format a concise label for verbose output."""
    if r.method in EXTERNAL_TOOLS:
        return r.method
    return f"kalign {r.refine}"


def run_benchmark(
    dataset: str = "balibase",
    methods: Optional[List[str]] = None,
    modes: Optional[List[str]] = None,
    max_cases: int = 0,
    binary: str = "kalign",
    n_threads: int = 1,
    verbose: bool = False,
    parallel: int = 1,
) -> List[AlignmentResult]:
    """Run benchmark suite and return results."""
    if methods is None:
        methods = ["cli"]
    if modes is None:
        modes = ["default"]

    cases = get_cases(dataset, max_cases=max_cases if max_cases > 0 else None)

    if not cases:
        print(f"No benchmark cases found for dataset '{dataset}'.")
        print("Try running with --download-only first.")
        return []

    print(f"Running {len(cases)} cases from '{dataset}'")
    print(f"  Methods: {methods}")
    print(f"  Modes:   {modes}")
    if parallel > 1:
        print(f"  Workers: {parallel}")
    print()

    # Build work items
    work = []
    for case in cases:
        for method in methods:
            if method in EXTERNAL_TOOLS:
                work.append((case, method, binary, n_threads, "default"))
            else:
                for mode in modes:
                    work.append((case, method, binary, n_threads, mode))

    total = len(work)

    if parallel <= 1:
        results = []
        for i, item in enumerate(work):
            result = _run_one(item)
            results.append(result)
            if verbose:
                label = _result_label(result)
                if result.error:
                    print(f"[{i+1}/{total}] {result.family:<12} {label:<25} ERROR: {result.error}")
                else:
                    print(f"[{i+1}/{total}] {result.family:<12} {label:<25} SP={result.recall:.3f}  TC={result.tc:.3f}  F1={result.f1:.3f}  {result.wall_time:.1f}s")
    else:
        results = []
        done = 0
        with ProcessPoolExecutor(max_workers=parallel) as pool:
            futures = {pool.submit(_run_one, item): i for i, item in enumerate(work)}
            indexed_results = [None] * total
            for future in as_completed(futures):
                idx = futures[future]
                result = future.result()
                indexed_results[idx] = result
                done += 1
                if verbose:
                    label = _result_label(result)
                    if result.error:
                        print(f"[{done}/{total}] {result.family:<12} {label:<25} ERROR: {result.error}")
                    else:
                        print(f"[{done}/{total}] {result.family:<12} {label:<25} SP={result.recall:.3f}  TC={result.tc:.3f}  F1={result.f1:.3f}  {result.wall_time:.1f}s")
        results = [r for r in indexed_results if r is not None]

    return results


def print_summary(results: List[AlignmentResult]) -> None:
    """Print aggregate summary of benchmark results."""
    by_group = {}
    for r in results:
        if r.error:
            continue
        if r.method in EXTERNAL_TOOLS:
            key = r.method
        else:
            key = f"kalign {r.refine}"
        by_group.setdefault(key, []).append(r)

    print(f"\n{'Method':<24} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8} {'N':>5}")
    print("-" * 75)

    for group, group_results in sorted(by_group.items()):
        recalls = [r.recall for r in group_results]
        precisions = [r.precision for r in group_results]
        f1s = [r.f1 for r in group_results]
        tcs = [r.tc for r in group_results]
        total_time = sum(r.wall_time for r in group_results)

        print(f"{group:<24} {statistics.mean(recalls):>8.3f} {statistics.mean(precisions):>8.3f} "
              f"{statistics.mean(f1s):>8.3f} {statistics.mean(tcs):>8.3f} "
              f"{total_time:>7.0f}s {len(group_results):>5}")

    # Per-category breakdown
    categories = sorted({r.dataset for r in results if not r.error})
    if len(categories) > 1:
        for cat in categories:
            cat_results = [r for r in results if r.dataset == cat and not r.error]
            if not cat_results:
                continue
            cat_groups = {}
            for r in cat_results:
                key = r.method if r.method in EXTERNAL_TOOLS else f"kalign {r.refine}"
                cat_groups.setdefault(key, []).append(r)

            cat_name = cat.replace("balibase_", "")
            n = len(next(iter(cat_groups.values())))
            print(f"\n--- {cat_name} ({n} cases) ---")
            print(f"{'Method':<24} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
            print("-" * 60)
            for group, gr in sorted(cat_groups.items()):
                print(f"{group:<24} {statistics.mean(r.recall for r in gr):>8.3f} "
                      f"{statistics.mean(r.precision for r in gr):>8.3f} "
                      f"{statistics.mean(r.f1 for r in gr):>8.3f} "
                      f"{statistics.mean(r.tc for r in gr):>8.3f}")

    errors = [r for r in results if r.error]
    if errors:
        print(f"\n{len(errors)} error(s):")
        for r in errors:
            print(f"  {r.family} ({r.method}): {r.error}")


def save_results(results: List[AlignmentResult], path: str) -> None:
    """Save results as JSON."""
    data = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "results": [r.to_dict() for r in results],
        "summary": {},
    }

    by_group = {}
    for r in results:
        if r.error:
            continue
        key = r.method if r.method in EXTERNAL_TOOLS else f"kalign_{r.refine}"
        by_group.setdefault(key, []).append(r)

    for group, group_results in by_group.items():
        recalls = [r.recall for r in group_results]
        precisions = [r.precision for r in group_results]
        f1s = [r.f1 for r in group_results]
        tcs = [r.tc for r in group_results]
        data["summary"][group] = {
            "n_cases": len(group_results),
            "recall_mean": statistics.mean(recalls),
            "precision_mean": statistics.mean(precisions),
            "f1_mean": statistics.mean(f1s),
            "tc_mean": statistics.mean(tcs),
            "total_time": sum(r.wall_time for r in group_results),
        }

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\nResults saved to {path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Kalign alignment benchmark suite",
        prog="python -m benchmarks",
    )
    parser.add_argument(
        "--dataset",
        default="balibase",
        choices=list(DATASETS.keys()) + ["all"],
        help="Which dataset to benchmark (default: balibase)",
    )
    parser.add_argument(
        "--method",
        nargs="+",
        default=["cli"],
        choices=["python_api", "cli", "clustalo", "mafft", "muscle"],
        help="Alignment method(s) to test (default: cli)",
    )
    parser.add_argument(
        "--mode",
        nargs="+",
        default=["default"],
        choices=["fast", "default", "recall", "accurate"],
        help="Kalign mode preset(s) to test (default: default)",
    )
    parser.add_argument(
        "--max-cases",
        type=int,
        default=0,
        help="Limit number of test cases (0 = all)",
    )
    parser.add_argument(
        "--binary",
        default="kalign",
        help="Path to kalign binary for CLI method (default: kalign)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads per alignment (default: 1)",
    )
    parser.add_argument(
        "--output",
        default="",
        help="Output JSON file for results",
    )
    parser.add_argument(
        "-j", "--parallel",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Only download datasets, don't run benchmarks",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output",
    )

    args = parser.parse_args()

    if args.download_only:
        download_dataset(args.dataset)
        print("Download complete.")
        return

    results = run_benchmark(
        dataset=args.dataset,
        methods=args.method,
        modes=args.mode,
        max_cases=args.max_cases,
        binary=args.binary,
        n_threads=args.threads,
        verbose=args.verbose,
        parallel=args.parallel,
    )

    if results:
        print_summary(results)
        if args.output:
            save_results(results, args.output)
