"""Sweep pseudocount values for seq_weights profile rebalancing.

Runs baseline (0.0) and several pseudocount values on BAliBASE
to find the optimal rebalancing strength.
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
    case, pseudo = args
    label = f"p={pseudo:.0f}" if pseudo > 0 else "baseline"
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        t0 = time.perf_counter()
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type,
            vsm_amax=-1.0,
            refine="confident",
            seq_weights=pseudo,
        )
        wall = time.perf_counter() - t0
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "pseudo": pseudo, "method": label,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
            "wall_time": wall,
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "pseudo": pseudo, "method": label,
            "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "wall_time": 0, "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _print_table(groups, methods, label="Overall"):
    n = len(next(iter(groups.values())))
    print(f"\n=== {label} ({n} cases) ===")
    print(f"{'Method':>12} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
    print("-" * 58)
    for m in methods:
        entries = groups[m]
        valid = [r for r in entries if "error" not in r]
        if not valid:
            print(f"{m:>12}  (all errors)")
            continue
        rec = statistics.mean(r["recall"] for r in valid)
        prec = statistics.mean(r["precision"] for r in valid)
        f1 = statistics.mean(r["f1"] for r in valid)
        tc = statistics.mean(r["tc"] for r in valid)
        total_t = sum(r["wall_time"] for r in valid)
        print(f"{m:>12} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {total_t:>7.1f}s")


def main():
    parser = argparse.ArgumentParser(description="Sweep seq_weights pseudocount on BAliBASE")
    parser.add_argument("-j", "--parallel", type=int, default=8)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--pseudos", nargs="*", type=float,
                        default=[0.0, 1.0, 3.0, 5.0, 10.0, 20.0, 50.0],
                        help="Pseudocount values to test (0=baseline)")
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)
    pseudos = args.pseudos

    print(f"{len(cases)} BAliBASE cases")
    print(f"Pseudocount values: {pseudos}")

    method_labels = [f"p={p:.0f}" if p > 0 else "baseline" for p in pseudos]

    tasks = [(case, p) for case in cases for p in pseudos]
    print(f"{len(pseudos)} methods x {len(cases)} cases = {len(tasks)} tasks")

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
    print(f"All done in {elapsed:.0f}s")

    # Overall summary
    groups = defaultdict(list)
    for r in results:
        groups[r["method"]].append(r)
    _print_table(groups, method_labels, "Overall")

    # Per-category
    all_cats = sorted({r["dataset"].replace("balibase_", "") for r in results})
    for cat in all_cats:
        cat_groups = defaultdict(list)
        for r in results:
            if cat in r["dataset"]:
                cat_groups[r["method"]].append(r)
        _print_table(cat_groups, method_labels, cat)

    # Save
    out = Path("benchmarks/results/seq_weights_sweep.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
