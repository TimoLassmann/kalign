"""Sweep dist_scale on BAliBASE protein to see if distance-dependent
gap scaling improves recall on divergent sequences.

Usage:
    uv run python -m benchmarks.protein_dist_scale_sweep -j 4
"""

import argparse
import json
import statistics
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kalign
from .datasets import get_cases
from .scoring import parse_balibase_xml


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(case.reference), str(output_path), column_mask=mask)
    return kalign.compare_detailed(str(case.reference), str(output_path))


def _run_one(args):
    case, dist_scale = args
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type, dist_scale=dist_scale,
        )
        sp = kalign.compare(str(case.reference), tmp_path)
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "dist_scale": dist_scale, "sp_score": sp,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "dist_scale": dist_scale, "sp_score": 0,
            "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description="dist_scale sweep on BAliBASE")
    parser.add_argument("-j", "--parallel", type=int, default=4)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--categories", nargs="*", default=None)
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)
    if args.categories:
        cats = [c.upper() for c in args.categories]
        cases = [c for c in cases if any(cat in c.dataset.upper() for cat in cats)]

    ds_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]

    tasks = [(case, ds) for case in cases for ds in ds_values]
    print(f"{len(cases)} cases x {len(ds_values)} dist_scale values = {len(tasks)} tasks")

    t0 = time.perf_counter()
    results = []
    done = 0
    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(_run_one, t): t for t in tasks}
        for f in as_completed(futures):
            done += 1
            r = f.result()
            results.append(r)
            if done % 200 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"All done in {elapsed:.0f}s\n")

    # Overall summary
    groups = defaultdict(list)
    for r in results:
        if "error" not in r:
            groups[r["dist_scale"]].append(r)

    print(f"{'ds':>6} {'SP':>8} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 50)
    for ds in sorted(groups.keys()):
        entries = groups[ds]
        sp = statistics.mean(r["sp_score"] for r in entries)
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        print(f"{ds:>6.1f} {sp:>8.1f} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")

    # Per category
    all_cats = sorted({r["dataset"] for r in results if "error" not in r})
    for cat in all_cats:
        cat_results = [r for r in results if r["dataset"] == cat and "error" not in r]
        cat_groups = defaultdict(list)
        for r in cat_results:
            cat_groups[r["dist_scale"]].append(r)

        cat_label = cat.replace("balibase_", "")
        n = len(cat_results) // len(cat_groups)
        print(f"\n=== {cat_label} ({n} cases) ===")
        print(f"{'ds':>6} {'SP':>8} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 50)
        for ds in sorted(cat_groups.keys()):
            entries = cat_groups[ds]
            sp = statistics.mean(r["sp_score"] for r in entries)
            rec = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            print(f"{ds:>6.1f} {sp:>8.1f} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")

    # Save
    out = Path("benchmarks/results/protein_dist_scale_sweep.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
