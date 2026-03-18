"""Benchmark kalign mode presets (fast/default/accurate) vs external tools on BAliBASE."""

import statistics
import time
from pathlib import Path

from .datasets import get_cases
from .scoring import run_case, EXTERNAL_TOOLS


def main():
    cases = get_cases("balibase")
    if not cases:
        print("No BAliBASE cases found. Run: uv run python -m benchmarks --download-only")
        return

    print(f"BAliBASE: {len(cases)} cases\n")

    methods = [
        ("kalign fast",     dict(method="python_api", mode="fast")),
        ("kalign default",  dict(method="python_api", mode="default")),
        ("kalign accurate", dict(method="python_api", mode="accurate")),
        ("clustalo",        dict(method="clustalo")),
        ("mafft",           dict(method="mafft")),
        ("muscle",          dict(method="muscle")),
    ]

    all_results = {}
    for label, kwargs in methods:
        print(f"Running {label}...", flush=True)
        results = []
        t0 = time.perf_counter()
        for i, case in enumerate(cases):
            r = run_case(case, n_threads=1, **kwargs)
            results.append(r)
            if r.error:
                print(f"  [{i+1}/{len(cases)}] {r.family}: ERROR {r.error}")
        elapsed = time.perf_counter() - t0

        ok = [r for r in results if not r.error]
        all_results[label] = ok
        if ok:
            print(f"  {len(ok)} cases, total {elapsed:.1f}s")
        print()

    # Summary table
    print(f"{'Method':<20} {'N':>4} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time(s)':>10}")
    print("-" * 72)
    for label, _ in methods:
        ok = all_results.get(label, [])
        if not ok:
            print(f"{label:<20} {'0':>4} {'n/a':>8} {'n/a':>8} {'n/a':>8} {'n/a':>8} {'n/a':>10}")
            continue
        rec = statistics.mean([r.recall for r in ok])
        pre = statistics.mean([r.precision for r in ok])
        f1 = statistics.mean([r.f1 for r in ok])
        tc = statistics.mean([r.tc for r in ok])
        t = sum(r.wall_time for r in ok)
        print(f"{label:<20} {len(ok):>4} {rec:>8.3f} {pre:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>10.1f}")


if __name__ == "__main__":
    main()
