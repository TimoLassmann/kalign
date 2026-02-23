"""Sweep dist_scale values across BAliBASE to find optimal gap scaling.

Usage:
    uv run python -m benchmarks.sweep_dist_scale
    uv run python -m benchmarks.sweep_dist_scale --dist-scales 0.0 0.3 0.5 0.7 1.0
    uv run python -m benchmarks.sweep_dist_scale --categories RV11 RV12
    uv run python -m benchmarks.sweep_dist_scale --csv results.csv
"""

import argparse
import json
import statistics
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional

import kalign

from .datasets import balibase_cases, balibase_is_available, balibase_download
from .scoring import score_alignment_detailed


RESULTS_DIR = Path(__file__).parent / "results"


def align_and_score(case, dist_scale: float, n_threads: int = 1) -> dict:
    """Align one case with given dist_scale, return scores."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / f"{case.family}_aln.fasta"
        start = time.perf_counter()
        kalign.align_file_to_file(
            str(case.unaligned),
            str(output),
            format="fasta",
            seq_type=case.seq_type,
            n_threads=n_threads,
            dist_scale=dist_scale,
        )
        wall_time = time.perf_counter() - start
        detailed = score_alignment_detailed(case.reference, output)
        return {
            "family": case.family,
            "dataset": case.dataset,
            "dist_scale": dist_scale,
            "recall": detailed["recall"],
            "precision": detailed["precision"],
            "f1": detailed["f1"],
            "tc": detailed["tc"],
            "wall_time": wall_time,
        }


def run_sweep(
    dist_scales: List[float],
    categories: Optional[List[str]] = None,
    n_threads: int = 1,
    verbose: bool = False,
) -> List[dict]:
    """Run sweep across all dist_scale values and BAliBASE cases."""
    if not balibase_is_available():
        print("Downloading BAliBASE...")
        balibase_download()

    cases = balibase_cases()
    if categories:
        cats = set(c.upper() for c in categories)
        cases = [c for c in cases if any(cat in c.dataset.upper() for cat in cats)]

    if not cases:
        print("No cases found!")
        return []

    total = len(cases) * len(dist_scales)
    print(f"Sweeping {len(dist_scales)} dist_scale values across {len(cases)} cases ({total} alignments)")
    print(f"dist_scale values: {dist_scales}")
    print()

    results = []
    done = 0
    for ds in dist_scales:
        for case in cases:
            done += 1
            try:
                r = align_and_score(case, ds, n_threads)
                results.append(r)
                if verbose:
                    print(f"  [{done}/{total}] {case.family} ds={ds:.2f} → "
                          f"F1={r['f1']:.3f} prec={r['precision']:.3f} rec={r['recall']:.3f} TC={r['tc']:.3f}")
            except Exception as e:
                print(f"  [{done}/{total}] {case.family} ds={ds:.2f} → ERROR: {e}")
                results.append({
                    "family": case.family,
                    "dataset": case.dataset,
                    "dist_scale": ds,
                    "recall": 0.0,
                    "precision": 0.0,
                    "f1": 0.0,
                    "tc": 0.0,
                    "wall_time": 0.0,
                    "error": str(e),
                })

    return results


def summarize(results: List[dict]) -> None:
    """Print summary table grouped by (category, dist_scale)."""
    # Group by (dataset, dist_scale)
    groups: Dict[tuple, List[dict]] = {}
    for r in results:
        if "error" in r:
            continue
        key = (r["dataset"], r["dist_scale"])
        groups.setdefault(key, []).append(r)

    # Get all categories and dist_scales
    all_cats = sorted(set(r["dataset"] for r in results))
    all_ds = sorted(set(r["dist_scale"] for r in results))

    # Print per-category tables
    print("\n" + "=" * 90)
    print("RESULTS BY CATEGORY")
    print("=" * 90)

    for cat in all_cats:
        print(f"\n--- {cat} ---")
        print(f"{'dist_scale':>10}  {'F1':>7}  {'Prec':>7}  {'Rec':>7}  {'TC':>7}  {'n':>3}")
        print("-" * 50)
        for ds in all_ds:
            key = (cat, ds)
            if key not in groups:
                continue
            rows = groups[key]
            n = len(rows)
            f1 = statistics.mean(r["f1"] for r in rows)
            prec = statistics.mean(r["precision"] for r in rows)
            rec = statistics.mean(r["recall"] for r in rows)
            tc = statistics.mean(r["tc"] for r in rows)
            print(f"{ds:>10.2f}  {f1:>7.3f}  {prec:>7.3f}  {rec:>7.3f}  {tc:>7.3f}  {n:>3}")

    # Print overall summary
    print(f"\n--- OVERALL (all categories) ---")
    print(f"{'dist_scale':>10}  {'F1':>7}  {'Prec':>7}  {'Rec':>7}  {'TC':>7}  {'n':>3}")
    print("-" * 50)
    for ds in all_ds:
        all_rows = [r for r in results if r["dist_scale"] == ds and "error" not in r]
        if not all_rows:
            continue
        n = len(all_rows)
        f1 = statistics.mean(r["f1"] for r in all_rows)
        prec = statistics.mean(r["precision"] for r in all_rows)
        rec = statistics.mean(r["recall"] for r in all_rows)
        tc = statistics.mean(r["tc"] for r in all_rows)
        print(f"{ds:>10.2f}  {f1:>7.3f}  {prec:>7.3f}  {rec:>7.3f}  {tc:>7.3f}  {n:>3}")

    # Find best dist_scale per category
    print(f"\n--- BEST dist_scale PER CATEGORY (by F1) ---")
    for cat in all_cats:
        best_ds = None
        best_f1 = -1.0
        for ds in all_ds:
            key = (cat, ds)
            if key not in groups:
                continue
            f1 = statistics.mean(r["f1"] for r in groups[key])
            if f1 > best_f1:
                best_f1 = f1
                best_ds = ds
        baseline_key = (cat, 0.0)
        baseline_f1 = statistics.mean(r["f1"] for r in groups.get(baseline_key, [{"f1": 0.0}]))
        delta = best_f1 - baseline_f1 if baseline_key in groups else 0.0
        sign = "+" if delta >= 0 else ""
        print(f"  {cat}: best dist_scale={best_ds:.2f}  F1={best_f1:.3f}  ({sign}{delta:.3f} vs baseline)")


def write_csv(results: List[dict], path: str) -> None:
    """Write results to CSV."""
    import csv
    fields = ["family", "dataset", "dist_scale", "recall", "precision", "f1", "tc", "wall_time"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    print(f"\nWrote {len(results)} rows to {path}")


def main():
    parser = argparse.ArgumentParser(description="Sweep dist_scale across BAliBASE")
    parser.add_argument(
        "--dist-scales", nargs="+", type=float,
        default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0],
        help="dist_scale values to test (default: 0.0 0.1 0.2 0.3 0.4 0.5 0.7 1.0)",
    )
    parser.add_argument(
        "--categories", nargs="+",
        help="BAliBASE categories to test (default: all). e.g. RV11 RV12 RV20",
    )
    parser.add_argument("--threads", type=int, default=1, help="Threads per alignment")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show per-case results")
    parser.add_argument("--csv", type=str, help="Write CSV output to this path")
    parser.add_argument("--json", type=str, help="Write JSON output to this path")

    args = parser.parse_args()

    results = run_sweep(
        dist_scales=args.dist_scales,
        categories=args.categories,
        n_threads=args.threads,
        verbose=args.verbose,
    )

    summarize(results)

    if args.csv:
        write_csv(results, args.csv)

    if args.json:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        json_path = args.json
        with open(json_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Wrote JSON to {json_path}")


if __name__ == "__main__":
    main()
