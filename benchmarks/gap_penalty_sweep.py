"""2D gap penalty grid search with PFASUM43 on BAliBASE.

Sweep gpo and gpe to find the best combination. No adaptive tricks,
no per-category tuning — just find what works best overall.

Usage:
    uv run python -m benchmarks.gap_penalty_sweep -j 8
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


# Grid: gpo x gpe
GPO_VALUES = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0]
GPE_VALUES = [1.0, 1.5, 2.0, 2.5]


def _build_configs():
    configs = []
    for gpo in GPO_VALUES:
        for gpe in GPE_VALUES:
            configs.append({
                "name": f"gpo={gpo:.1f}_gpe={gpe:.1f}",
                "gpo": gpo, "gpe": gpe, "tgpe": 1.0,
            })
    return configs


CONFIGS = _build_configs()


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
            seq_type="pfasum43",
            gap_open=config["gpo"],
            gap_extend=config["gpe"],
            terminal_gap_extend=config["tgpe"],
            vsm_amax=-1.0,
            refine="confident",
        )
        wall_time = time.perf_counter() - t0
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "gpo": config["gpo"], "gpe": config["gpe"],
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
            "wall_time": wall_time,
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "gpo": config["gpo"], "gpe": config["gpe"],
            "f1": 0, "recall": 0, "precision": 0, "tc": 0,
            "wall_time": 0, "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description="Gap penalty sweep with PFASUM43")
    parser.add_argument("-j", "--parallel", type=int, default=4)
    args = parser.parse_args()

    cases = get_cases("balibase")
    print(f"{len(cases)} cases x {len(CONFIGS)} gap combos = {len(cases) * len(CONFIGS)} tasks")

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
            if done % 500 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"All done in {elapsed:.0f}s\n")

    # Build 2D table
    print("=" * 100)
    print("PFASUM43 + confident refine — F1 by gap open (rows) x gap extend (cols)")
    print("=" * 100)

    grid = defaultdict(dict)
    for r in results:
        if "error" not in r:
            key = (r["gpo"], r["gpe"])
            if key not in grid:
                grid[key] = {"f1s": [], "recs": [], "precs": [], "tcs": []}
            grid[key]["f1s"].append(r["f1"])
            grid[key]["recs"].append(r["recall"])
            grid[key]["precs"].append(r["precision"])
            grid[key]["tcs"].append(r["tc"])

    # Print F1 grid
    print(f"\n{'gpo\\gpe':>10}", end="")
    for gpe in GPE_VALUES:
        print(f"  {gpe:>7.1f}", end="")
    print()
    print("-" * (10 + 9 * len(GPE_VALUES)))

    best_f1 = 0
    best_params = None
    for gpo in GPO_VALUES:
        print(f"{gpo:>10.1f}", end="")
        for gpe in GPE_VALUES:
            key = (gpo, gpe)
            if key in grid:
                f1 = statistics.mean(grid[key]["f1s"])
                print(f"  {f1:>7.4f}", end="")
                if f1 > best_f1:
                    best_f1 = f1
                    best_params = key
            else:
                print(f"  {'—':>7}", end="")
        print()

    print(f"\nBest: gpo={best_params[0]}, gpe={best_params[1]} → F1={best_f1:.4f}")

    # Print recall grid
    print(f"\n{'Recall':>10}", end="")
    for gpe in GPE_VALUES:
        print(f"  {gpe:>7.1f}", end="")
    print()
    print("-" * (10 + 9 * len(GPE_VALUES)))
    for gpo in GPO_VALUES:
        print(f"{gpo:>10.1f}", end="")
        for gpe in GPE_VALUES:
            key = (gpo, gpe)
            if key in grid:
                rec = statistics.mean(grid[key]["recs"])
                print(f"  {rec:>7.4f}", end="")
            else:
                print(f"  {'—':>7}", end="")
        print()

    # Print precision grid
    print(f"\n{'Precision':>10}", end="")
    for gpe in GPE_VALUES:
        print(f"  {gpe:>7.1f}", end="")
    print()
    print("-" * (10 + 9 * len(GPE_VALUES)))
    for gpo in GPO_VALUES:
        print(f"{gpo:>10.1f}", end="")
        for gpe in GPE_VALUES:
            key = (gpo, gpe)
            if key in grid:
                prec = statistics.mean(grid[key]["precs"])
                print(f"  {prec:>7.4f}", end="")
            else:
                print(f"  {'—':>7}", end="")
        print()

    # Print TC grid
    print(f"\n{'TC':>10}", end="")
    for gpe in GPE_VALUES:
        print(f"  {gpe:>7.1f}", end="")
    print()
    print("-" * (10 + 9 * len(GPE_VALUES)))
    for gpo in GPO_VALUES:
        print(f"{gpo:>10.1f}", end="")
        for gpe in GPE_VALUES:
            key = (gpo, gpe)
            if key in grid:
                tc = statistics.mean(grid[key]["tcs"])
                print(f"  {tc:>7.4f}", end="")
            else:
                print(f"  {'—':>7}", end="")
        print()

    # Per-category for the best params
    best_gpo, best_gpe = best_params
    best_name = f"gpo={best_gpo:.1f}_gpe={best_gpe:.1f}"
    print(f"\n{'=' * 100}")
    print(f"Per-category breakdown for best: {best_name}")
    print(f"{'=' * 100}")
    print(f"{'Category':<12} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 45)
    all_cats = sorted({r["dataset"] for r in results})
    for cat in all_cats:
        rows = [r for r in results if r["dataset"] == cat
                and r["config"] == best_name and "error" not in r]
        if rows:
            cat_label = cat.replace("balibase_", "")
            rec = statistics.mean(r["recall"] for r in rows)
            prec = statistics.mean(r["precision"] for r in rows)
            f1 = statistics.mean(r["f1"] for r in rows)
            tc = statistics.mean(r["tc"] for r in rows)
            print(f"{cat_label:<12} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    # Save
    out = Path("benchmarks/results/gap_penalty_sweep.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump({"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "results": results}, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
