#!/usr/bin/env python3
"""Evaluate selected configs from an NSGA-III checkpoint on full BAliBASE.

Reports true SP (recall), precision, F1, TC, and wall time — not CV estimates.
"""

import pickle
import sys
import time
from pathlib import Path


# Ensure the benchmarks package is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmarks.optimize_unified import (
    decode_unified_params,
    evaluate_unified,
    set_active_profile,
)
from benchmarks.datasets import balibase_cases, balibase_download, balibase_is_available


def main():
    checkpoint_path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/gen_checkpoint.pkl"

    with open(checkpoint_path, "rb") as f:
        ckpt = pickle.load(f)

    X = ckpt["pop_X"]
    F = ckpt["pop_F"]
    max_runs = ckpt["max_runs"]
    f1_cv = -F[:, 0]
    tc_cv = -F[:, 1]

    # Selected configs to evaluate
    selected = {
        122: "Hyper-fast bare (n=1,no refine)",
        67:  "Fast+refine (n=1,ref=CONF)",
    }

    # Set up BAliBASE
    profile = ckpt.get("profile", "protein")
    set_active_profile(profile)

    if not balibase_is_available():
        print("Downloading BAliBASE...", flush=True)
        balibase_download()

    cases = balibase_cases()
    print(f"Loaded {len(cases)} BAliBASE cases (profile: {profile})\n")

    # Print header
    print(f"{'Idx':>4}  {'Label':<25} {'CV_F1':>6} {'CV_TC':>6} "
          f"{'Recall':>7} {'Prec':>7} {'F1':>7} {'TC':>7} {'Time':>8}")
    print("-" * 95)

    results = {}
    for idx, label in selected.items():
        x = X[idx]
        params = decode_unified_params(x, max_runs)

        print(f"[{idx:3d}]  {label:<25} {f1_cv[idx]:.4f} {tc_cv[idx]:.4f} ",
              end="", flush=True)

        start = time.perf_counter()
        result = evaluate_unified(params, cases, n_threads=1, quiet=False)
        elapsed = time.perf_counter() - start

        print(f"{result['recall']:.4f}  {result['precision']:.4f}  "
              f"{result['f1']:.4f}  {result['tc']:.4f}  {elapsed:7.1f}s")

        results[idx] = {
            "label": label,
            "cv_f1": f1_cv[idx],
            "cv_tc": tc_cv[idx],
            "params": params,
            **result,
        }

    # Summary table
    print("\n" + "=" * 95)
    print("FULL BENCHMARK RESULTS vs COMPETITORS")
    print("=" * 95)
    print(f"\n{'Method':<30} {'Recall(SP)':>10} {'Prec':>7} {'F1':>7} {'TC':>7} {'Time':>8}")
    print("-" * 75)

    for idx, label in selected.items():
        r = results[idx]
        print(f"kalign [{idx}] {label:<20} {r['recall']:>10.4f} {r['precision']:>7.4f} "
              f"{r['f1']:>7.4f} {r['tc']:>7.4f} {r['wall_time']:>7.1f}s")

    print("-" * 75)
    print(f"{'kalign fast (shipped)':<30} {'0.809':>10} {'0.663':>7} {'0.723':>7} {'0.482':>7} {'10':>7}s")
    print(f"{'kalign default (shipped)':<30} {'0.816':>10} {'0.758':>7} {'0.780':>7} {'0.490':>7} {'101':>7}s")
    print(f"{'kalign accurate (shipped)':<30} {'0.837':>10} {'0.719':>7} {'0.769':>7} {'0.518':>7} {'262':>7}s")
    print(f"{'ClustalO':<30} {'0.840':>10} {'0.710':>7} {'0.764':>7} {'0.559':>7}")
    print(f"{'MAFFT':<30} {'0.867':>10} {'0.715':>7} {'0.778':>7} {'0.590':>7}")
    print(f"{'MUSCLE':<30} {'0.870':>10} {'0.721':>7} {'0.783':>7} {'0.581':>7}")

    # Per-category breakdown for each config
    for idx, label in selected.items():
        r = results[idx]
        if "per_category" in r and r["per_category"]:
            print(f"\n--- [{idx}] {label} per-category ---")
            for cat in sorted(r["per_category"]):
                c = r["per_category"][cat]
                print(f"  {cat:<20} Recall={c['recall']:.4f} Prec={c['precision']:.4f} "
                      f"F1={c['f1']:.4f} TC={c['tc']:.4f} (n={c['n']})")


if __name__ == "__main__":
    main()
