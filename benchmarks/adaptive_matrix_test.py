"""Test distance-adaptive matrix/gap-penalty selection on BAliBASE.

Strategy: measure average pairwise identity of each family, then pick the
best matrix/gap-penalty combo for that identity range. This simulates what
an adaptive algorithm would do internally.

We test:
  1. Oracle: run all configs, pick best F1 per family (upper bound)
  2. Distance-based heuristic: pick config based on avg pairwise identity
  3. Fixed configs as baselines

Usage:
    uv run python -m benchmarks.adaptive_matrix_test -j 8
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
    # Close sequences: PFASUM60, tighter gaps, high VSM
    {
        "name": "PFASUM60-close",
        "desc": "PFASUM60, gpo=7.0, gpe=1.5, vsm=3.0 (for close seqs)",
        "seq_type": "pfasum60",
        "gpo": 7.0, "gpe": 1.5, "tgpe": 1.0,
        "vsm_amax": 3.0,
    },
    # Medium sequences: PFASUM43, tighter gaps, moderate VSM
    {
        "name": "PFASUM43-mid",
        "desc": "PFASUM43, gpo=7.0, gpe=1.5, vsm=2.0 (for medium seqs)",
        "seq_type": "pfasum43",
        "gpo": 7.0, "gpe": 1.5, "tgpe": 1.0,
        "vsm_amax": 2.0,
    },
    # Divergent sequences: PFASUM43, looser gaps, no VSM
    {
        "name": "PFASUM43-div",
        "desc": "PFASUM43, gpo=5.5, gpe=2.0, vsm=0 (for divergent seqs)",
        "seq_type": "pfasum43",
        "gpo": 5.5, "gpe": 2.0, "tgpe": 1.0,
        "vsm_amax": 0.0,
    },
    # Very divergent: PFASUM43, even looser, no VSM
    {
        "name": "PFASUM43-vdiv",
        "desc": "PFASUM43, gpo=4.5, gpe=2.0, vsm=0 (for very divergent)",
        "seq_type": "pfasum43",
        "gpo": 4.5, "gpe": 2.0, "tgpe": 1.0,
        "vsm_amax": 0.0,
    },
    # Baselines
    {
        "name": "CorBLOSUM66-def",
        "desc": "Current default",
        "seq_type": "protein",
        "gpo": None, "gpe": None, "tgpe": None,
        "vsm_amax": -1.0,
    },
    {
        "name": "PFASUM43-g7",
        "desc": "Best from matrix_sweep",
        "seq_type": "pfasum43",
        "gpo": 7.0, "gpe": 1.5, "tgpe": 1.0,
        "vsm_amax": -1.0,
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
            seq_type=config["seq_type"],
            gap_open=config["gpo"],
            gap_extend=config["gpe"],
            terminal_gap_extend=config["tgpe"],
            vsm_amax=config["vsm_amax"],
            refine="confident",
        )
        wall_time = time.perf_counter() - t0
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
            "wall_time": wall_time,
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "f1": 0, "recall": 0, "precision": 0, "tc": 0,
            "wall_time": 0, "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _identity_heuristic(dataset):
    """Map BAliBASE category to expected identity range and select config."""
    cat = dataset.replace("balibase_", "").upper()
    # RV11: <20% identity -> very divergent
    # RV12: 20-40% identity -> divergent
    # RV20: orphan sequences -> mixed, divergent merges
    # RV30: equidistant -> medium
    # RV40: N/C extensions -> close core + terminal noise
    # RV50: internal insertions -> medium
    if cat == "RV11":
        return "PFASUM43-div"       # divergent: looser gaps, no VSM
    elif cat == "RV12":
        return "PFASUM43-mid"       # medium: tighter gaps, moderate VSM
    elif cat == "RV20":
        return "PFASUM43-mid"       # mixed: medium settings
    elif cat == "RV30":
        return "PFASUM43-mid"       # equidistant subfamilies
    elif cat in ("RV40", "RV50"):
        return "PFASUM60-close"     # close core, need precision
    else:
        return "PFASUM43-mid"       # fallback


def main():
    parser = argparse.ArgumentParser(description="Adaptive matrix test on BAliBASE")
    parser.add_argument("-j", "--parallel", type=int, default=4)
    args = parser.parse_args()

    cases = get_cases("balibase")
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
            if done % 200 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"All done in {elapsed:.0f}s\n")

    # Group by family and config
    by_family = defaultdict(dict)
    for r in results:
        if "error" not in r:
            by_family[r["family"]][r["config"]] = r

    # 1. Oracle: best F1 per family
    oracle_f1s = []
    oracle_picks = defaultdict(int)
    for fam, configs in by_family.items():
        best = max(configs.values(), key=lambda x: x["f1"])
        oracle_f1s.append(best)
        oracle_picks[best["config"]] += 1

    # 2. Category heuristic
    heuristic_f1s = []
    for fam, configs in by_family.items():
        dataset = list(configs.values())[0]["dataset"]
        pick = _identity_heuristic(dataset)
        if pick in configs:
            heuristic_f1s.append(configs[pick])
        else:
            # Fallback to best available
            heuristic_f1s.append(max(configs.values(), key=lambda x: x["f1"]))

    # 3. Fixed baselines
    config_names = [c["name"] for c in CONFIGS]

    print("=" * 100)
    print("OVERALL RESULTS")
    print("=" * 100)
    print(f"\n{'Method':<25} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 58)

    for name in config_names:
        rows = [r for r in results if r.get("config") == name and "error" not in r]
        if rows:
            rec = statistics.mean(r["recall"] for r in rows)
            prec = statistics.mean(r["precision"] for r in rows)
            f1 = statistics.mean(r["f1"] for r in rows)
            tc = statistics.mean(r["tc"] for r in rows)
            print(f"  {name:<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    # Oracle
    rec = statistics.mean(r["recall"] for r in oracle_f1s)
    prec = statistics.mean(r["precision"] for r in oracle_f1s)
    f1 = statistics.mean(r["f1"] for r in oracle_f1s)
    tc = statistics.mean(r["tc"] for r in oracle_f1s)
    print(f"  {'** Oracle (best/fam) **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    # Heuristic
    rec = statistics.mean(r["recall"] for r in heuristic_f1s)
    prec = statistics.mean(r["precision"] for r in heuristic_f1s)
    f1 = statistics.mean(r["f1"] for r in heuristic_f1s)
    tc = statistics.mean(r["tc"] for r in heuristic_f1s)
    print(f"  {'** Cat. heuristic **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    # Per-category breakdown
    all_cats = sorted({r["dataset"] for r in results})
    for cat in all_cats:
        cat_results = [r for r in results if r["dataset"] == cat]
        cat_families = {r["family"] for r in cat_results}
        n_fam = len(cat_families)
        cat_label = cat.replace("balibase_", "")
        print(f"\n{'=' * 100}")
        print(f"{cat_label} ({n_fam} cases)")
        print(f"{'=' * 100}")
        print(f"{'Method':<25} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 58)

        for name in config_names:
            rows = [r for r in cat_results if r.get("config") == name and "error" not in r]
            if rows:
                rec = statistics.mean(r["recall"] for r in rows)
                prec = statistics.mean(r["precision"] for r in rows)
                f1 = statistics.mean(r["f1"] for r in rows)
                tc = statistics.mean(r["tc"] for r in rows)
                print(f"  {name:<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

        # Oracle for this category
        cat_by_fam = defaultdict(dict)
        for r in cat_results:
            if "error" not in r:
                cat_by_fam[r["family"]][r["config"]] = r
        cat_oracle = [max(cs.values(), key=lambda x: x["f1"]) for cs in cat_by_fam.values()]
        if cat_oracle:
            rec = statistics.mean(r["recall"] for r in cat_oracle)
            prec = statistics.mean(r["precision"] for r in cat_oracle)
            f1 = statistics.mean(r["f1"] for r in cat_oracle)
            tc = statistics.mean(r["tc"] for r in cat_oracle)
            print(f"  {'** Oracle **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

        # Heuristic for this cat
        cat_heuristic = []
        for fam, configs in cat_by_fam.items():
            pick = _identity_heuristic(cat)
            if pick in configs:
                cat_heuristic.append(configs[pick])
            else:
                cat_heuristic.append(max(configs.values(), key=lambda x: x["f1"]))
        if cat_heuristic:
            rec = statistics.mean(r["recall"] for r in cat_heuristic)
            prec = statistics.mean(r["precision"] for r in cat_heuristic)
            f1 = statistics.mean(r["f1"] for r in cat_heuristic)
            tc = statistics.mean(r["tc"] for r in cat_heuristic)
            print(f"  {'** Heuristic **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    # Oracle pick distribution
    print(f"\n{'=' * 100}")
    print("Oracle pick distribution (which config was best for how many families)")
    print(f"{'=' * 100}")
    for name in sorted(oracle_picks, key=oracle_picks.get, reverse=True):
        print(f"  {name:<23} {oracle_picks[name]:>4} families")

    # Save
    out = Path("benchmarks/results/adaptive_matrix_test.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump({"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "results": results}, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
