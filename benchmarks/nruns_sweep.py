"""Sweep N ensemble runs x proportional support thresholds.

Tests whether more runs with proportional thresholds (e.g. 75% of N)
can achieve better precision-recall tradeoff than fewer runs with unanimity.
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
    case, n_runs, ms = args
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type, ensemble=n_runs, min_support=ms,
        )
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "n_runs": n_runs, "min_support": ms,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "n_runs": n_runs, "min_support": ms,
            "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description="Sweep N runs x support thresholds")
    parser.add_argument("-j", "--parallel", type=int, default=16)
    parser.add_argument("--max-cases", type=int, default=0)
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)
    print(f"{len(cases)} BAliBASE cases")

    # Configs: (n_runs, min_support)
    configs = []
    # N=4: all thresholds
    for ms in range(1, 5):
        configs.append((4, ms))
    # N=8: all thresholds
    for ms in range(1, 9):
        configs.append((8, ms))

    tasks = [(case, nr, ms) for case in cases for nr, ms in configs]
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
        groups[(r["n_runs"], r["min_support"])].append(r)

    print(f"{'N':>4} {'s':>4} {'s/N':>6} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 52)
    for (nr, ms) in sorted(groups.keys()):
        entries = groups[(nr, ms)]
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        print(f"{nr:>4} {ms:>4} {ms/nr:>6.0%} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")

    # Save
    out = Path("benchmarks/results/mumsa_nruns_sweep.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
