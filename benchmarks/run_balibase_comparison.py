"""Run BAliBASE quality benchmark across all 4 kalign mode presets."""

import json
from pathlib import Path

import kalign
from benchmarks.datasets import get_cases, download_dataset
from benchmarks.scoring import run_case

download_dataset("balibase")
cases = get_cases("balibase")
print(f"Loaded {len(cases)} BAliBASE cases")

modes = ["fast", "default", "recall", "accurate"]

all_results = {}

for mode in modes:
    print(f"\n{'=' * 60}")
    print(f"  mode={mode}")
    print(f"{'=' * 60}")

    recalls, precs, f1s, tcs, times = [], [], [], [], []
    errors = 0
    for i, case in enumerate(cases):
        r = run_case(case, method="python_api", n_threads=8, mode=mode)
        if r.error:
            errors += 1
            continue
        recalls.append(r.recall)
        precs.append(r.precision)
        f1s.append(r.f1)
        tcs.append(r.tc)
        times.append(r.wall_time)
        if (i + 1) % 50 == 0:
            n = len(f1s)
            print(f"  [{i+1}/{len(cases)}] F1={sum(f1s)/n:.4f} so far...")

    n = len(f1s)
    avg_r = sum(recalls) / n if n else 0
    avg_p = sum(precs) / n if n else 0
    avg_f = sum(f1s) / n if n else 0
    avg_t = sum(tcs) / n if n else 0
    total_t = sum(times)

    all_results[mode] = {
        "mode": mode,
        "n_cases": n,
        "errors": errors,
        "recall": round(avg_r, 4),
        "precision": round(avg_p, 4),
        "f1": round(avg_f, 4),
        "tc": round(avg_t, 4),
        "total_time": round(total_t, 1),
    }
    print(f"  Recall={avg_r:.4f}  Prec={avg_p:.4f}  "
          f"F1={avg_f:.4f}  TC={avg_t:.4f}  "
          f"Time={total_t:.1f}s  ({errors} errors)")

# Summary
print()
print("=" * 72)
print("  RESULTS")
print("=" * 72)
print()
print(f"  {'Mode':<12} {'N':>4} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
print(f"  {'-'*12} {'-'*4} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
for mode in modes:
    s = all_results[mode]
    print(f"  {s['mode']:<12} {s['n_cases']:>4} {s['recall']:>8.4f} {s['precision']:>8.4f} "
          f"{s['f1']:>8.4f} {s['tc']:>8.4f} {s['total_time']:>7.1f}s")

out_file = Path("benchmarks/results/balibase_modes.json")
out_file.parent.mkdir(exist_ok=True)
with open(out_file, "w") as f:
    json.dump(all_results, f, indent=2)
print(f"\nSaved to {out_file}")
