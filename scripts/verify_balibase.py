#!/usr/bin/env python3
"""BAliBASE verification: generate checksums or compare against baseline.

Usage:
    # Generate baseline checksums (run with CURRENT code, 1 thread)
    uv run python scripts/verify_balibase.py baseline --output baseline_checksums.json

    # Compare against baseline (run AFTER code change + rebuild)
    uv run python scripts/verify_balibase.py compare --baseline baseline_checksums.json

    # Timing benchmark (20 largest cases, multiple thread counts)
    uv run python scripts/verify_balibase.py timing --threads 1,4,8,16

    # Quick test (10 cases only)
    uv run python scripts/verify_balibase.py baseline --output baseline.json --max-cases 10
"""

import argparse
import hashlib
import json
import os
import sys
import tempfile
import time
from pathlib import Path

import kalign


BB_ROOT = Path(__file__).parent.parent / "benchmarks" / "data" / "downloads" / "bb3_release"
MODES = ["fast", "default", "accurate"]


def get_cases(max_cases=0):
    """Discover all BAliBASE .tfa files, sorted by size (largest first for timing)."""
    cases = []
    for rv_dir in sorted(BB_ROOT.iterdir()):
        if not rv_dir.is_dir() or not rv_dir.name.startswith("RV"):
            continue
        for tfa in sorted(rv_dir.glob("*.tfa")):
            msf = tfa.with_suffix(".msf")
            if msf.exists():
                cases.append({
                    "name": tfa.stem,
                    "category": rv_dir.name,
                    "unaligned": str(tfa),
                    "reference": str(msf),
                    "size": tfa.stat().st_size,
                })
    # Sort by size descending (largest first — useful for timing)
    cases.sort(key=lambda c: -c["size"])
    if max_cases > 0:
        cases = cases[:max_cases]
    return cases


def sha256_file(path):
    """Compute SHA256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def run_alignment(case, mode, n_threads, output_path):
    """Run one alignment, return (sha256, elapsed_seconds)."""
    t0 = time.perf_counter()
    kalign.align_file_to_file(
        case["unaligned"], output_path,
        mode=mode, n_threads=n_threads, format="fasta",
    )
    elapsed = time.perf_counter() - t0
    checksum = sha256_file(output_path)
    return checksum, elapsed


def cmd_baseline(args):
    """Generate baseline checksums for all cases × all modes."""
    cases = get_cases(args.max_cases)
    modes = args.modes.split(",") if args.modes else MODES
    n_threads = 1  # baseline always single-threaded for determinism

    print(f"Generating baseline: {len(cases)} cases × {len(modes)} modes, {n_threads} thread")
    print(f"BAliBASE root: {BB_ROOT}")

    results = {}
    total = len(cases) * len(modes)
    done = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        for case in cases:
            for mode in modes:
                out_path = os.path.join(tmpdir, f"{case['name']}_{mode}.fa")
                checksum, elapsed = run_alignment(case, mode, n_threads, out_path)

                key = f"{case['name']}:{mode}"
                results[key] = {
                    "checksum": checksum,
                    "time": round(elapsed, 3),
                    "category": case["category"],
                }

                done += 1
                if args.verbose:
                    print(f"  [{done}/{total}] {key}: {checksum[:12]}... ({elapsed:.2f}s)")
                elif done % 50 == 0:
                    print(f"  [{done}/{total}]...")

    # Save
    output = {
        "n_cases": len(cases),
        "modes": modes,
        "n_threads": n_threads,
        "checksums": results,
    }
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nBaseline saved to {args.output}: {len(results)} checksums")


def cmd_compare(args):
    """Compare current code output against baseline checksums."""
    with open(args.baseline) as f:
        baseline = json.load(f)

    cases_by_name = {c["name"]: c for c in get_cases()}
    modes = baseline["modes"]
    n_threads = args.threads  # can test with >1 thread

    checksums = baseline["checksums"]
    print(f"Comparing against baseline: {len(checksums)} checksums, {n_threads} threads")

    mismatches = []
    matches = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        total = len(checksums)
        done = 0

        for key, bdata in checksums.items():
            name, mode = key.split(":")
            if name not in cases_by_name:
                print(f"  SKIP {key}: case not found")
                continue

            case = cases_by_name[name]
            out_path = os.path.join(tmpdir, f"{name}_{mode}.fa")
            checksum, elapsed = run_alignment(case, mode, n_threads, out_path)

            done += 1
            if checksum == bdata["checksum"]:
                matches += 1
                if args.verbose:
                    print(f"  [{done}/{total}] {key}: OK ({elapsed:.2f}s)")
            else:
                mismatches.append({
                    "key": key,
                    "baseline": bdata["checksum"],
                    "current": checksum,
                })
                print(f"  [{done}/{total}] {key}: MISMATCH!")
                print(f"    baseline: {bdata['checksum'][:16]}...")
                print(f"    current:  {checksum[:16]}...")

            if not args.verbose and done % 50 == 0:
                print(f"  [{done}/{total}] ({matches} match, {len(mismatches)} mismatch)...")

    print(f"\n{'='*60}")
    print(f"Results: {matches} identical, {len(mismatches)} mismatched (of {done})")

    if mismatches:
        print(f"\nMISMATCHES:")
        for m in mismatches:
            print(f"  {m['key']}")
        return 1
    else:
        print("ALL IDENTICAL — verification passed.")
        return 0


def cmd_timing(args):
    """Timing benchmark on largest cases."""
    n_cases = args.max_cases or 20
    cases = get_cases(n_cases)
    modes = args.modes.split(",") if args.modes else MODES
    thread_counts = [int(x) for x in args.threads.split(",")]

    print(f"Timing benchmark: {len(cases)} cases × {len(modes)} modes × {len(thread_counts)} thread configs")

    results = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        for mode in modes:
            print(f"\n=== Mode: {mode} ===")
            for n_threads in thread_counts:
                total_time = 0.0
                for case in cases:
                    out_path = os.path.join(tmpdir, f"{case['name']}_{mode}_{n_threads}.fa")
                    _, elapsed = run_alignment(case, mode, n_threads, out_path)
                    total_time += elapsed

                key = f"{mode}:{n_threads}"
                results[key] = round(total_time, 2)

                baseline_key = f"{mode}:{thread_counts[0]}"
                baseline_time = results.get(baseline_key, total_time)
                speedup = baseline_time / total_time if total_time > 0 else 0

                print(f"  {n_threads:2d} threads: {total_time:7.1f}s total  (speedup: {speedup:.2f}x)")

    # Print summary table
    print(f"\n{'='*60}")
    print(f"{'Mode':<12s}", end="")
    for t in thread_counts:
        print(f"  {t:>2d}T", end="")
    print("  speedup")
    print("-" * (12 + 5 * len(thread_counts) + 10))
    for mode in modes:
        print(f"{mode:<12s}", end="")
        baseline = results.get(f"{mode}:{thread_counts[0]}", 1)
        for t in thread_counts:
            val = results.get(f"{mode}:{t}", 0)
            print(f"  {val:>4.0f}", end="")
        last = results.get(f"{mode}:{thread_counts[-1]}", 1)
        print(f"  {baseline/last:.2f}x")


def main():
    parser = argparse.ArgumentParser(description="BAliBASE verification and timing")
    sub = parser.add_subparsers(dest="command", required=True)

    p_base = sub.add_parser("baseline", help="Generate baseline checksums")
    p_base.add_argument("--output", required=True, help="Output JSON file")
    p_base.add_argument("--modes", default=None, help="Comma-separated modes (default: fast,default,accurate)")
    p_base.add_argument("--max-cases", type=int, default=0, help="Limit number of cases (0=all)")
    p_base.add_argument("-v", "--verbose", action="store_true")

    p_cmp = sub.add_parser("compare", help="Compare against baseline")
    p_cmp.add_argument("--baseline", required=True, help="Baseline JSON file")
    p_cmp.add_argument("--threads", type=int, default=1, help="Thread count for comparison run")
    p_cmp.add_argument("-v", "--verbose", action="store_true")

    p_time = sub.add_parser("timing", help="Timing benchmark")
    p_time.add_argument("--threads", default="1,4,8,16", help="Comma-separated thread counts")
    p_time.add_argument("--modes", default=None, help="Comma-separated modes")
    p_time.add_argument("--max-cases", type=int, default=20, help="Number of largest cases")

    args = parser.parse_args()

    if not BB_ROOT.exists():
        print(f"ERROR: BAliBASE not found at {BB_ROOT}")
        print("Run: uv run python -m benchmarks --download-only --dataset balibase")
        return 1

    if args.command == "baseline":
        return cmd_baseline(args)
    elif args.command == "compare":
        return cmd_compare(args)
    elif args.command == "timing":
        return cmd_timing(args)


if __name__ == "__main__":
    sys.exit(main() or 0)
