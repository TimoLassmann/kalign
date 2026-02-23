"""Grid search: consistency anchors x weight x vsm_amax on BAliBASE."""

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

    # Grid search parameters
    consistency_vals = [0, 5, 10, 15, 20]
    weight_vals = [0.5, 1.0, 2.0, 4.0, 8.0]
    vsm_vals = [0.0, 1.0, 2.0, 3.0]

    # Base config: refine=confident, seq_weights=1.0
    configs = {}
    for c in consistency_vals:
        for w in weight_vals:
            for v in vsm_vals:
                if c == 0 and w != weight_vals[0]:
                    continue  # skip redundant weight combos when no consistency
                name = f"c{c}_w{w:.1f}_v{v:.1f}"
                kw = {
                    "refine": "confident",
                    "seq_weights": 1.0,
                    "vsm_amax": v,
                }
                if c > 0:
                    kw["consistency"] = c
                    kw["consistency_weight"] = w
                configs[name] = kw

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
            if done % 200 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\nDone in {elapsed:.0f}s\n")

    # Aggregate results
    groups = defaultdict(list)
    for r in results:
        if "error" not in r:
            groups[r["config"]].append(r)

    # Sort by F1 descending
    summary = []
    for name in configs:
        entries = groups.get(name, [])
        if not entries:
            continue
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        t = sum(r["wall_time"] for r in entries)
        summary.append((name, rec, prec, f1, tc, t))

    summary.sort(key=lambda x: -x[3])  # sort by F1

    print(f"{'Config':<22} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
    print("-" * 76)
    for name, rec, prec, f1, tc, t in summary[:30]:
        print(f"{name:<22} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s")

    # Per-category breakdown for top 5 configs
    print(f"\n--- Per-category breakdown (top 5 configs) ---")
    top5 = [s[0] for s in summary[:5]]
    by_cat = defaultdict(lambda: defaultdict(list))
    for r in results:
        if "error" not in r and r["config"] in top5:
            cat = r["dataset"].replace("balibase_", "")
            by_cat[cat][r["config"]].append(r["f1"])

    header = f"{'Category':<12}" + "".join(f" {n:>20}" for n in top5)
    print(header)
    print("-" * len(header))
    for cat in sorted(by_cat):
        row = f"{cat:<12}"
        for name in top5:
            vals = by_cat[cat].get(name, [])
            if vals:
                row += f" {statistics.mean(vals):>20.3f}"
            else:
                row += f" {'—':>20}"
        print(row)

    # Baseline comparison
    baseline_name = "c0_w0.5_v2.0"
    baseline = dict((s[0], s) for s in summary)
    if baseline_name in baseline:
        _, _, _, bf1, btc, _ = baseline[baseline_name]
        print(f"\nBaseline ({baseline_name}): F1={bf1:.3f}, TC={btc:.3f}")
        print(f"\n{'Config':<22} {'F1':>8} {'dF1':>8} {'TC':>8} {'dTC':>8}")
        print("-" * 56)
        for name, rec, prec, f1, tc, t in summary[:15]:
            print(f"{name:<22} {f1:>8.3f} {f1-bf1:>+8.3f} {tc:>8.3f} {tc-btc:>+8.3f}")


if __name__ == "__main__":
    main()
