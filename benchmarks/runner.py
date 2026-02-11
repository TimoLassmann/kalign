"""Main benchmark orchestrator and CLI."""

import argparse
import json
import statistics
import time
from pathlib import Path
from typing import List, Optional

from .datasets import download_dataset, get_cases, DATASETS
from .scoring import AlignmentResult, run_case


def run_benchmark(
    dataset: str = "balibase",
    methods: Optional[List[str]] = None,
    refine_modes: Optional[List[str]] = None,
    max_cases: int = 0,
    binary: str = "kalign",
    n_threads: int = 1,
    verbose: bool = False,
    adaptive_budget: bool = False,
) -> List[AlignmentResult]:
    """Run benchmark suite and return results."""
    if methods is None:
        methods = ["python_api"]
    if refine_modes is None:
        refine_modes = ["none"]

    cases = get_cases(dataset, max_cases=max_cases if max_cases > 0 else None)

    if not cases:
        print(f"No benchmark cases found for dataset '{dataset}'.")
        print("Try running with --download-only first.")
        return []

    print(f"Running {len(cases)} cases from '{dataset}' with methods: {methods}, refine: {refine_modes}")
    print()

    results = []
    total = len(cases) * len(methods) * len(refine_modes)
    done = 0

    for case in cases:
        for method in methods:
            for refine in refine_modes:
                done += 1
                if verbose:
                    print(f"[{done}/{total}] {case.family} (refine={refine}) ... ", end="", flush=True)

                result = run_case(case, method=method, binary=binary, n_threads=n_threads, refine=refine, adaptive_budget=adaptive_budget)
                results.append(result)

                if verbose:
                    if result.error:
                        print(f"ERROR: {result.error}")
                    else:
                        print(f"SP={result.sp_score:.1f}  time={result.wall_time:.2f}s")

    return results


def print_summary(results: List[AlignmentResult]) -> None:
    """Print aggregate summary of benchmark results."""
    by_group = {}
    for r in results:
        if r.error:
            continue
        key = f"{r.method} refine={r.refine}"
        by_group.setdefault(key, []).append(r)

    for group, group_results in sorted(by_group.items()):
        scores = [r.sp_score for r in group_results]
        times = [r.wall_time for r in group_results]

        print(f"\n--- {group} ({len(group_results)} cases) ---")
        print(f"  SP score:  mean={statistics.mean(scores):.1f}  "
              f"median={statistics.median(scores):.1f}  "
              f"min={min(scores):.1f}  max={max(scores):.1f}")
        print(f"  Time (s):  total={sum(times):.1f}  "
              f"mean={statistics.mean(times):.2f}  "
              f"max={max(times):.2f}")

    errors = [r for r in results if r.error]
    if errors:
        print(f"\n{len(errors)} error(s):")
        for r in errors:
            print(f"  {r.family} ({r.method} refine={r.refine}): {r.error}")


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
        key = f"{r.method}_refine={r.refine}"
        by_group.setdefault(key, []).append(r)

    for group, group_results in by_group.items():
        scores = [r.sp_score for r in group_results]
        data["summary"][group] = {
            "n_cases": len(group_results),
            "sp_mean": statistics.mean(scores),
            "sp_median": statistics.median(scores),
            "sp_min": min(scores),
            "sp_max": max(scores),
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
        default=["python_api"],
        choices=["python_api", "cli"],
        help="Alignment method(s) to test (default: python_api)",
    )
    parser.add_argument(
        "--max-cases",
        type=int,
        default=0,
        help="Limit number of test cases (0 = all)",
    )
    parser.add_argument(
        "--binary",
        default="build/src/kalign",
        help="Path to C-compiled kalign binary for CLI method (default: build/src/kalign)",
    )
    parser.add_argument(
        "--refine",
        nargs="+",
        default=["none"],
        choices=["none", "all", "confident"],
        help="Refinement mode(s) to test (default: none)",
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
        "--adaptive-budget",
        action="store_true",
        help="Scale trial count by uncertainty",
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
        refine_modes=args.refine,
        max_cases=args.max_cases,
        binary=args.binary,
        n_threads=args.threads,
        verbose=args.verbose,
        adaptive_budget=args.adaptive_budget,
    )

    if results:
        print_summary(results)
        if args.output:
            save_results(results, args.output)
