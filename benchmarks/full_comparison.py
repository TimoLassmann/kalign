"""Full BAliBASE comparison: all kalign modes + external tools (cached).

Shows ALL metrics: SP(=Recall), Precision, F1, TC, Time — overall and per-category.
"""

import json
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


def _print_table(title, display_order, summary_data, configs):
    """Print one table with ALL metrics: SP, Prec, F1, TC, Time."""
    print(f"\n{'=' * 100}")
    print(f"  {title}")
    print(f"{'=' * 100}")
    print(
        f"{'Method':<24} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8}"
        f" {'Time':>8} {'N':>4}"
    )
    print("-" * 100)
    for name in display_order:
        if name not in summary_data:
            continue
        sp, prec, f1, tc, t, n = summary_data[name]
        marker = " " if name in configs else "*"
        t_str = f"{t:.0f}s" if t is not None else "—"
        print(
            f"{marker}{name:<23} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f}"
            f" {tc:>8.3f} {t_str:>8} {n:>4}"
        )
    print("-" * 100)
    print("  SP = Sum-of-Pairs recall (BAliBASE convention)")
    print("  * = external tool (cached results)")


def main():
    cases = get_cases("balibase")
    print(f"BAliBASE: {len(cases)} cases\n")

    # Kalign configurations to test
    configs = {
        "baseline(vsm=0)": {"vsm_amax": 0.0},
        "+vsm": {"vsm_amax": 2.0},
        "+vsm+ref": {"vsm_amax": 2.0, "refine": "confident"},
        "+vsm+ref+sw": {
            "vsm_amax": 2.0,
            "refine": "confident",
            "seq_weights": 1.0,
        },
        "+vsm+ref+sw+c5": {
            "vsm_amax": 2.0,
            "refine": "confident",
            "seq_weights": 1.0,
            "consistency": 5,
            "consistency_weight": 2.0,
        },
        "ens3+vsm": {
            "vsm_amax": 2.0,
            "ensemble": 3,
        },
        "ens3+vsm+ref": {
            "vsm_amax": 2.0,
            "ensemble": 3,
            "refine": "confident",
        },
        "ens3+vsm+ref+ra1": {
            "vsm_amax": 2.0,
            "ensemble": 3,
            "refine": "confident",
            "realign": 1,
        },
        "ens3+vsm+ref+ra1+c5": {
            "vsm_amax": 2.0,
            "ensemble": 3,
            "refine": "confident",
            "realign": 1,
            "consistency": 5,
            "consistency_weight": 2.0,
        },
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
            if done % 200 == 0:
                elapsed = time.perf_counter() - t0
                print(f"  {done}/{len(tasks)} ({elapsed:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\nKalign runs done in {elapsed:.0f}s\n")

    # Load cached external results
    ext_path = Path("benchmarks/data/external_balibase.json")
    if ext_path.exists():
        ext_data = json.load(open(ext_path))
        for r in ext_data:
            r["config"] = r.pop("method")
        results.extend(ext_data)
        ext_methods = set(r["config"] for r in ext_data)
        print(f"Loaded {len(ext_data)} cached external results: {ext_methods}\n")

    # Group by config
    groups = defaultdict(list)
    for r in results:
        if "error" not in r:
            groups[r["config"]].append(r)

    # Display order: kalign configs first, then external sorted
    display_order = list(configs.keys()) + sorted(
        k for k in groups if k not in configs
    )

    # Compute overall summary: (sp, prec, f1, tc, total_time, count)
    overall = {}
    for name in display_order:
        entries = groups.get(name, [])
        if not entries:
            continue
        sp = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        t = sum(r["wall_time"] for r in entries)
        overall[name] = (sp, prec, f1, tc, t, len(entries))

    # Print overall table
    _print_table("OVERALL (218 cases)", display_order, overall, configs)

    # Group results by category
    categories = ["RV11", "RV12", "RV20", "RV30", "RV40", "RV50"]
    by_cat = defaultdict(lambda: defaultdict(list))
    for r in results:
        if "error" not in r:
            cat = r["dataset"].replace("balibase_", "")
            by_cat[cat][r["config"]].append(r)

    # Per-category tables — each with ALL metrics
    for cat in categories:
        cat_summary = {}
        for name in display_order:
            entries = by_cat[cat].get(name, [])
            if not entries:
                continue
            sp = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            t = sum(r["wall_time"] for r in entries)
            cat_summary[name] = (sp, prec, f1, tc, t, len(entries))
        _print_table(f"{cat} ({len(by_cat[cat].get(display_order[0], []))} cases)", display_order, cat_summary, configs)


if __name__ == "__main__":
    main()
