#!/usr/bin/env python3
"""Sweep tgpe values on BAliBASE to find optimal terminal gap penalty.

Runs the full BAliBASE dataset (218 protein cases) with multiple tgpe values.
Default protein tgpe=1.0; this tests values from 0.0 to 2.0.
"""

import statistics
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kalign
from benchmarks.datasets import get_cases
from benchmarks.scoring import parse_balibase_xml


def _score_case(reference, test_output):
    """Score with BAliBASE XML core block mask."""
    xml_path = reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(reference), str(test_output), column_mask=mask)
    return kalign.compare_detailed(str(reference), str(test_output), max_gap_frac=-1.0)


def _run_one(args):
    """Worker: align one case with given tgpe, return scores."""
    case, tgpe_value, label = args
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / "aln.fasta"
        t0 = time.perf_counter()
        kalign.align_file_to_file(
            str(case.unaligned),
            str(output),
            format="fasta",
            seq_type=case.seq_type,
            terminal_gap_extend=tgpe_value,
            n_threads=1,
        )
        wall = time.perf_counter() - t0
        scores = _score_case(case.reference, output)
        return {
            "family": case.family,
            "label": label,
            "tgpe": tgpe_value,
            "recall": scores["recall"],
            "precision": scores["precision"],
            "f1": scores["f1"],
            "tc": scores["tc"],
            "wall_time": wall,
        }


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-j", "--parallel", type=int, default=8)
    parser.add_argument("--max-cases", type=int, default=0)
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases > 0 else None)
    print(f"Running {len(cases)} BAliBASE cases")

    # Sweep: default (1.0), plus values around it
    configs = [
        (-1.0, "default"),     # internal default = 1.0
        (0.25, "tgpe=0.25"),
        (0.50, "tgpe=0.50"),
        (0.75, "tgpe=0.75"),
        (1.00, "tgpe=1.00"),   # same as default, sanity check
        (1.25, "tgpe=1.25"),
        (1.50, "tgpe=1.50"),
        (2.00, "tgpe=2.00"),
    ]

    work = []
    for case in cases:
        for tgpe_val, label in configs:
            work.append((case, tgpe_val, label))

    print(f"Total work items: {len(work)} ({len(cases)} cases x {len(configs)} configs)")
    results = []
    done = 0
    total = len(work)

    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(_run_one, item): i for i, item in enumerate(work)}
        for future in as_completed(futures):
            r = future.result()
            results.append(r)
            done += 1
            if done % 100 == 0:
                print(f"  [{done}/{total}]")

    # Group results
    by_label = {}
    for r in results:
        by_label.setdefault(r["label"], []).append(r)

    # Print summary table
    labels = [c[1] for c in configs]
    print()
    print(f"{'Config':<12} {'SP(recall)':>10} {'Precision':>10} {'F1':>10} {'TC':>10} {'Time(s)':>10}")
    print("-" * 64)
    for label in labels:
        rs = by_label[label]
        print(
            f"{label:<12} "
            f"{statistics.mean(r['recall'] for r in rs):>10.4f} "
            f"{statistics.mean(r['precision'] for r in rs):>10.4f} "
            f"{statistics.mean(r['f1'] for r in rs):>10.4f} "
            f"{statistics.mean(r['tc'] for r in rs):>10.4f} "
            f"{sum(r['wall_time'] for r in rs):>10.1f}"
        )

    # Per-case comparison: default vs each config
    default_by_fam = {r["family"]: r for r in by_label["default"]}

    print()
    print("Per-case F1 deltas vs default:")
    print(f"{'Config':<12} {'Improved':>8} {'Worsened':>8} {'Unchanged':>9} {'Mean dF1':>10} {'Median dF1':>10}")
    print("-" * 60)
    for label in labels:
        if label == "default":
            continue
        config_by_fam = {r["family"]: r for r in by_label[label]}
        diffs = []
        for fam in default_by_fam:
            delta = config_by_fam[fam]["f1"] - default_by_fam[fam]["f1"]
            diffs.append(delta)
        improved = sum(1 for d in diffs if d > 0.001)
        worsened = sum(1 for d in diffs if d < -0.001)
        unchanged = len(diffs) - improved - worsened
        print(
            f"{label:<12} "
            f"{improved:>8d} "
            f"{worsened:>8d} "
            f"{unchanged:>9d} "
            f"{statistics.mean(diffs):>+10.4f} "
            f"{statistics.median(diffs):>+10.4f}"
        )

    # Per-RV breakdown for best non-default config
    print()
    print("Per-reference subset (RV) breakdown:")
    rv_map = {"BB1": "RV11", "BB2": "RV12/20", "BB3": "RV30", "BB4": "RV40", "BB5": "RV50"}
    for label in labels:
        if label == "default":
            continue
        config_by_fam = {r["family"]: r for r in by_label[label]}
        by_rv = {}
        for fam in default_by_fam:
            prefix = fam[:3]
            rv = rv_map.get(prefix, prefix)
            delta = config_by_fam[fam]["f1"] - default_by_fam[fam]["f1"]
            by_rv.setdefault(rv, []).append(delta)
        rv_summary = {rv: statistics.mean(ds) for rv, ds in sorted(by_rv.items())}
        parts = "  ".join(f"{rv}:{m:+.3f}" for rv, m in rv_summary.items())
        print(f"  {label:<12} {parts}")


if __name__ == "__main__":
    main()
