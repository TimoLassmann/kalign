"""Measure cumulative/additive effect of kalign improvements on BAliBASE.

Tests individual features and their combinations to determine additivity:
  Individual features (from baseline):
    1. baseline:        vsm_amax=0, no extras
    2. +vsm:            vsm_amax=2.0 only
    3. +ref:            refine=confident only (vsm=0)
    4. +re1:            realign=1 only (vsm=0)
  Stacking (additive):
    5. +vsm+ref:        vsm + refine
    6. +vsm+re1:        vsm + realign=1 (no refine)
    7. +vsm+ref+re1:    vsm + refine + realign=1
  Ensemble path:
    8. +ens3:           ensemble=3 (uses vsm, internal post-refine)
    9. +ens5:           ensemble=5
   10. +ens5+adapt:     ensemble=5 + adaptive_budget

Usage:
    uv run python -m benchmarks.combined_improvements -j 4
    uv run python -m benchmarks.combined_improvements -j 4 --categories RV11 RV12
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


CONFIGS = [
    # --- Individual features (from baseline) ---
    {
        "name": "baseline",
        "desc": "old defaults (vsm=0, no extras)",
        "vsm_amax": 0.0,
        "dist_scale": 0.0,
        "refine": "none",
        "ensemble": 0,
        "adaptive_budget": False,
        "realign": 0,
    },
    {
        "name": "+vsm",
        "desc": "vsm_amax=2.0 only",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "none",
        "ensemble": 0,
        "adaptive_budget": False,
        "realign": 0,
    },
    {
        "name": "+ref",
        "desc": "refine=confident only (vsm=0)",
        "vsm_amax": 0.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 0,
        "adaptive_budget": False,
        "realign": 0,
    },
    {
        "name": "+re1",
        "desc": "realign=1 only (vsm=0)",
        "vsm_amax": 0.0,
        "dist_scale": 0.0,
        "refine": "none",
        "ensemble": 0,
        "adaptive_budget": False,
        "realign": 1,
    },
    # --- Stacking (additive combinations) ---
    {
        "name": "+vsm+ref",
        "desc": "vsm + refine",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 0,
        "adaptive_budget": False,
        "realign": 0,
    },
    {
        "name": "+vsm+re1",
        "desc": "vsm + realign=1 (no refine)",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "none",
        "ensemble": 0,
        "adaptive_budget": False,
        "realign": 1,
    },
    {
        "name": "+vsm+ref+re1",
        "desc": "vsm + refine + realign=1 (full stack)",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 0,
        "adaptive_budget": False,
        "realign": 1,
    },
    # --- Ensemble path ---
    {
        "name": "+ens3",
        "desc": "ensemble=3 (uses vsm, internal post-refine)",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 3,
        "adaptive_budget": False,
        "realign": 0,
    },
    {
        "name": "+ens5",
        "desc": "ensemble=5 (uses vsm, internal post-refine)",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 5,
        "adaptive_budget": False,
        "realign": 0,
    },
    {
        "name": "+ens5+adapt",
        "desc": "ensemble=5 + adaptive_budget",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 5,
        "adaptive_budget": True,
        "realign": 0,
    },
    # --- Ensemble + post-realign ---
    {
        "name": "+ens3+re1",
        "desc": "ensemble=3 + post-ensemble realign=1",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 3,
        "adaptive_budget": False,
        "realign": 1,
    },
    {
        "name": "+ens5+re1",
        "desc": "ensemble=5 + post-ensemble realign=1",
        "vsm_amax": 2.0,
        "dist_scale": 0.0,
        "refine": "confident",
        "ensemble": 5,
        "adaptive_budget": False,
        "realign": 1,
    },
]


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(case.reference), str(output_path), column_mask=mask)
    return kalign.compare_detailed(str(case.reference), str(output_path))


def _run_one(args):
    case, config = args
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        t0 = time.perf_counter()
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type,
            vsm_amax=config["vsm_amax"],
            dist_scale=config["dist_scale"],
            refine=config["refine"],
            ensemble=config["ensemble"],
            adaptive_budget=config["adaptive_budget"],
            realign=config.get("realign", 0),
        )
        wall_time = time.perf_counter() - t0
        sp_score = kalign.compare(str(case.reference), tmp_path)
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "sp_score": sp_score,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
            "wall_time": wall_time,
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "sp_score": 0, "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "wall_time": 0, "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description="Combined improvements sweep on BAliBASE")
    parser.add_argument("-j", "--parallel", type=int, default=4)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--categories", nargs="*", default=None,
                        help="BAliBASE categories (e.g. RV11 RV12)")
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)
    if args.categories:
        cats = [c.upper() for c in args.categories]
        cases = [c for c in cases if any(cat in c.dataset.upper() for cat in cats)]

    print(f"{len(cases)} BAliBASE cases x {len(CONFIGS)} configs = {len(cases) * len(CONFIGS)} tasks")

    tasks = [(case, config) for case in cases for config in CONFIGS]

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
    print(f"All done in {elapsed:.0f}s\n")

    # Overall summary
    print("=" * 90)
    print("OVERALL")
    print("=" * 90)
    _print_config_table(results)

    # Per-category
    all_cats = sorted({r["dataset"] for r in results})
    for cat in all_cats:
        cat_results = [r for r in results if r["dataset"] == cat]
        n_cases = len(cat_results) // len(CONFIGS)
        cat_label = cat.replace("balibase_", "")
        print(f"\n{'=' * 90}")
        print(f"{cat_label} ({n_cases} cases)")
        print(f"{'=' * 90}")
        _print_config_table(cat_results)

    # Save
    out = Path("benchmarks/results/combined_improvements.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump({"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "results": results}, f, indent=2)
    print(f"\nSaved to {out}")


def _print_config_table(results):
    config_names = [c["name"] for c in CONFIGS]
    groups = defaultdict(list)
    for r in results:
        if "error" not in r:
            groups[r["config"]].append(r)

    print(f"{'Config':<16} {'SP':>8} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time(s)':>9}  {'dSP':>6} {'dF1':>6}")
    print("-" * 92)
    baseline_sp = None
    baseline_f1 = None
    for name in config_names:
        entries = groups.get(name, [])
        if not entries:
            continue
        sp = statistics.mean(r["sp_score"] for r in entries)
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        total_time = sum(r["wall_time"] for r in entries)
        if baseline_f1 is None:
            baseline_sp = sp
            baseline_f1 = f1
            dsp = ""
            df1 = ""
        else:
            dsp = f"{sp - baseline_sp:+.1f}"
            df1 = f"{f1 - baseline_f1:+.3f}"
        print(f"{name:<16} {sp:>8.1f} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {total_time:>9.1f}  {dsp:>6} {df1:>6}")


if __name__ == "__main__":
    main()
