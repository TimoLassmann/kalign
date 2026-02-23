#!/usr/bin/env python3
"""Comparison: tgpo terminal_dist_scale sweep on BAliBASE."""

import argparse
import statistics
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import kalign
from benchmarks.datasets import get_cases
from benchmarks.scoring import parse_balibase_xml


def _score_case(reference, test_output):
    xml_path = reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(reference), str(test_output), column_mask=mask)
    return kalign.compare_detailed(str(reference), str(test_output), max_gap_frac=-1.0)


def _run_one(args):
    case_unaligned, case_ref, case_family, kwargs, label = args
    with tempfile.TemporaryDirectory() as tmpdir:
        out = Path(tmpdir) / "aln.fasta"
        kalign.align_file_to_file(str(case_unaligned), str(out), **kwargs)
        d = _score_case(case_ref, out)
        return label, case_family, d["recall"], d["precision"], d["f1"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--workers", type=int, default=12)
    args = parser.parse_args()

    cases = get_cases("balibase")

    # Sweep terminal_dist_scale values
    # tds=0 means no distance-dependent tgpo reduction (baseline-equivalent)
    # Higher tds = more aggressive reduction for divergent sequences
    configs = [
        (dict(terminal_dist_scale=0.0), "tds=0.0 (baseline)"),
        (dict(), "tds=2.0 (default)"),
        (dict(terminal_dist_scale=4.0), "tds=4.0"),
        (dict(terminal_dist_scale=8.0), "tds=8.0"),
    ]

    tasks = []
    for c in cases:
        for kwargs, label in configs:
            tasks.append((c.unaligned, c.reference, c.family, kwargs, label))

    results = {}
    per_case = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_run_one, t): t for t in tasks}
        done = 0
        for fut in as_completed(futs):
            done += 1
            label, family, rec, prec, f1 = fut.result()
            results.setdefault(label, []).append(f1)
            per_case.setdefault(family, {})[label] = f1
            if done % max(1, len(tasks) // 10) == 0:
                print(f"  {done}/{len(tasks)} done...")

    # Summary
    print()
    labels = [l for _, l in configs]
    for label in labels:
        vals = results[label]
        print(f"{label:>20s}: F1={statistics.mean(vals):.4f}  (n={len(vals)})")

    # Per-case delta: default vs baseline
    baseline_label = labels[0]
    default_label = labels[1]
    print(f"\nPer-case delta ({default_label} - {baseline_label}):")
    deltas = []
    for fam in sorted(per_case):
        b = per_case[fam].get(baseline_label, 0)
        d = per_case[fam].get(default_label, 0)
        deltas.append((d - b, fam, b, d))
    deltas.sort()
    print("\nBiggest losses:")
    for delta, fam, b, d in deltas[:5]:
        print(f"  {fam}: {b:.4f} -> {d:.4f}  ({delta:+.4f})")
    print("\nBiggest gains:")
    for delta, fam, b, d in deltas[-5:]:
        print(f"  {fam}: {b:.4f} -> {d:.4f}  ({delta:+.4f})")
