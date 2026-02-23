"""Compare baseline vs anchor consistency on BAliBASE."""

import statistics
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kalign
from benchmarks.datasets import get_cases
from benchmarks.scoring import parse_balibase_xml


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(
            str(case.reference), str(output_path), column_mask=mask
        )
    return kalign.compare_detailed(str(case.reference), str(output_path))


def _run_one(args):
    case, config_name, kwargs = args
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        t0 = time.perf_counter()
        kalign.align_file_to_file(
            str(case.unaligned),
            tmp_path,
            format="fasta",
            seq_type=case.seq_type,
            **kwargs,
        )
        wall_time = time.perf_counter() - t0
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family,
            "dataset": case.dataset,
            "config": config_name,
            **scores,
            "wall_time": wall_time,
        }
    except Exception as e:
        return {
            "family": case.family,
            "dataset": case.dataset,
            "config": config_name,
            "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    cases = get_cases("balibase")
    print(f"BAliBASE: {len(cases)} cases")

    configs = {
        "baseline":       {"refine": "none", "vsm_amax": 2.0},
        "consistency=5":  {"refine": "none", "vsm_amax": 2.0, "consistency": 5},
        "consistency=10": {"refine": "none", "vsm_amax": 2.0, "consistency": 10},
        "consistency=20": {"refine": "none", "vsm_amax": 2.0, "consistency": 20},
    }

    tasks = [(c, name, kw) for c in cases for name, kw in configs.items()]
    print(f"{len(configs)} configs x {len(cases)} cases = {len(tasks)} tasks\n")

    results = []
    done = 0
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = {pool.submit(_run_one, t): t for t in tasks}
        for f in as_completed(futures):
            done += 1
            r = f.result()
            results.append(r)
            if done % 50 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\nDone in {elapsed:.0f}s\n")

    # Summary by config
    print(f"{'Config':<20} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
    print("-" * 72)

    groups = defaultdict(list)
    for r in results:
        if "error" not in r:
            groups[r["config"]].append(r)

    baseline_f1 = None
    baseline_tc = None
    for name in configs:
        entries = groups.get(name, [])
        if not entries:
            print(f"{name:<20} (no results)")
            continue
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        t = sum(r["wall_time"] for r in entries)

        if baseline_f1 is None:
            baseline_f1 = f1
            baseline_tc = tc
            delta = ""
        else:
            delta = f"  F1{f1-baseline_f1:+.3f} TC{tc-baseline_tc:+.3f}"

        print(f"{name:<20} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s{delta}")

    # Per-category breakdown for best config
    print(f"\n{'Category':<12} {'baseline F1':>12} {'c=10 F1':>12} {'delta':>8}")
    print("-" * 50)
    by_cat = defaultdict(lambda: defaultdict(list))
    for r in results:
        if "error" not in r:
            cat = r["dataset"].replace("balibase_", "")
            by_cat[cat][r["config"]].append(r["f1"])
    for cat in sorted(by_cat):
        b = statistics.mean(by_cat[cat].get("baseline", [0]))
        c = statistics.mean(by_cat[cat].get("consistency=10", [0]))
        print(f"{cat:<12} {b:>12.3f} {c:>12.3f} {c-b:>+8.3f}")


if __name__ == "__main__":
    main()
