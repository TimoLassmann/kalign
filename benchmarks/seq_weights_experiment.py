"""Compare BAliBASE results with and without sequence weighting.

Runs the fastest good mode (vsm_amax=2.0, refine=confident, no ensemble)
with seq_weights=False vs seq_weights=True.
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
    case, use_weights = args
    label = "weights" if use_weights else "baseline"
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        t0 = time.perf_counter()
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type,
            vsm_amax=-1.0,          # use C default (2.0 for protein)
            refine="confident",
            seq_weights=use_weights,
        )
        wall = time.perf_counter() - t0
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "method": label,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
            "wall_time": wall,
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "method": label,
            "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "wall_time": 0, "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _print_table(groups, methods, label="Overall"):
    print(f"\n=== {label} ({len(next(iter(groups.values())))} cases) ===")
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
    parser = argparse.ArgumentParser(description="Seq weights A/B on BAliBASE")
    parser.add_argument("-j", "--parallel", type=int, default=8)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--categories", nargs="*", default=None,
                        help="BAliBASE categories (e.g. RV11 RV12)")
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)

    if args.categories:
        cats = [c.upper() for c in args.categories]
        cases = [c for c in cases if any(cat in c.dataset.upper() for cat in cats)]

    print(f"{len(cases)} BAliBASE cases")

    methods_bool = [False, True]
    methods_label = ["baseline", "weights"]

    tasks = [(case, w) for case in cases for w in methods_bool]
    print(f"2 methods x {len(cases)} cases = {len(tasks)} tasks")

    t0 = time.perf_counter()
    results = []
    done = 0
    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(_run_one, t): t for t in tasks}
        for f in as_completed(futures):
            done += 1
            r = f.result()
            results.append(r)
            if done % 100 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"All done in {elapsed:.0f}s")

    # Overall summary
    groups = defaultdict(list)
    for r in results:
        groups[r["method"]].append(r)
    _print_table(groups, methods_label, "Overall")

    # Per-category
    all_cats = sorted({r["dataset"].replace("balibase_", "") for r in results})
    for cat in all_cats:
        cat_groups = defaultdict(list)
        for r in results:
            if cat in r["dataset"]:
                cat_groups[r["method"]].append(r)
        _print_table(cat_groups, methods_label, cat)

    # Per-family deltas (biggest wins/losses)
    families = sorted({r["family"] for r in results})
    deltas = []
    for fam in families:
        fam_results = {r["method"]: r for r in results if r["family"] == fam}
        if "baseline" in fam_results and "weights" in fam_results:
            b = fam_results["baseline"]
            w = fam_results["weights"]
            if "error" not in b and "error" not in w:
                deltas.append({
                    "family": fam,
                    "dataset": b["dataset"],
                    "delta_f1": w["f1"] - b["f1"],
                    "delta_tc": w["tc"] - b["tc"],
                    "f1_base": b["f1"],
                    "f1_wt": w["f1"],
                    "tc_base": b["tc"],
                    "tc_wt": w["tc"],
                })

    deltas.sort(key=lambda d: d["delta_f1"])
    print(f"\n=== Biggest F1 regressions (weights vs baseline) ===")
    print(f"{'Family':>12} {'Cat':>6} {'F1_base':>8} {'F1_wt':>8} {'delta':>8}")
    for d in deltas[:10]:
        cat = d["dataset"].replace("balibase_", "")
        print(f"{d['family']:>12} {cat:>6} {d['f1_base']:>8.3f} {d['f1_wt']:>8.3f} {d['delta_f1']:>+8.3f}")

    print(f"\n=== Biggest F1 improvements (weights vs baseline) ===")
    print(f"{'Family':>12} {'Cat':>6} {'F1_base':>8} {'F1_wt':>8} {'delta':>8}")
    for d in deltas[-10:]:
        cat = d["dataset"].replace("balibase_", "")
        print(f"{d['family']:>12} {cat:>6} {d['f1_base']:>8.3f} {d['f1_wt']:>8.3f} {d['delta_f1']:>+8.3f}")

    # Save
    out = Path("benchmarks/results/seq_weights_experiment.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
