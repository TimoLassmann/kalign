"""Test VSM combined with distance-dependent gap scaling.

Sweeps vsm_amax and dist_scale together to find optimal combination.
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
    case, amax, ds = args
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type, vsm_amax=amax, dist_scale=ds,
        )
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "vsm_amax": amax, "dist_scale": ds,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "vsm_amax": amax, "dist_scale": ds,
            "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description="VSM + dist_scale combined sweep")
    parser.add_argument("-j", "--parallel", type=int, default=4)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--categories", nargs="*", default=None)
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)

    if args.categories:
        cats = [c.upper() for c in args.categories]
        cases = [c for c in cases if any(cat in c.dataset.upper() for cat in cats)]

    print(f"{len(cases)} BAliBASE cases")

    # Configs: (vsm_amax, dist_scale)
    configs = [
        (0.0, 0.0),   # baseline
        (0.5, 0.0),   # VSM only
        (0.8, 0.0),
        (1.0, 0.0),
        (0.0, 0.3),   # dist_scale only
        (0.0, 0.5),
        (0.5, 0.3),   # combined
        (0.5, 0.5),
        (0.8, 0.3),
        (0.8, 0.5),
        (1.0, 0.3),
        (1.0, 0.5),
    ]

    tasks = [(case, amax, ds) for case in cases for amax, ds in configs]
    print(f"{len(configs)} configs x {len(cases)} cases = {len(tasks)} tasks")

    t0 = time.perf_counter()
    results = []
    done = 0
    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(_run_one, t): t for t in tasks}
        for f in as_completed(futures):
            done += 1
            r = f.result()
            results.append(r)
            if done % 500 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"All done in {elapsed:.0f}s\n")

    # Summarize
    groups = defaultdict(list)
    for r in results:
        groups[(r["vsm_amax"], r["dist_scale"])].append(r)

    print(f"{'amax':>6} {'ds':>6} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 50)
    for (amax, ds) in sorted(groups.keys()):
        entries = groups[(amax, ds)]
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        print(f"{amax:>6.1f} {ds:>6.1f} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")

    # Per category
    all_cats = sorted({r["dataset"].replace("balibase_", "") for r in results})
    for cat in all_cats:
        cat_results = [r for r in results if cat in r["dataset"]]
        cat_groups = defaultdict(list)
        for r in cat_results:
            cat_groups[(r["vsm_amax"], r["dist_scale"])].append(r)

        print(f"\n=== {cat} ===")
        print(f"{'amax':>6} {'ds':>6} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 50)
        for (amax, ds) in sorted(cat_groups.keys()):
            entries = cat_groups[(amax, ds)]
            rec = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            print(f"{amax:>6.1f} {ds:>6.1f} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")

    # Save
    out = Path("benchmarks/results/vsm_combined_sweep.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
