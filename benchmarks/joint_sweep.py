"""Joint parameter sweep for single-run kalign on BAliBASE.

Sweeps all key parameters together to find the true optimum,
avoiding sequential tuning bias.
"""

import statistics
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product
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
    print(f"BAliBASE: {len(cases)} cases\n")

    # === Joint parameter grid ===
    # Key interactions to test:
    # 1. VSM x consistency (do they complement or overlap?)
    # 2. seq_weights x consistency (both affect profile building)
    # 3. refine x consistency (does consistency make refine useful again?)
    # 4. gap penalties (never tuned!)

    configs = {}

    # Baseline (original kalign3)
    configs["original"] = {"vsm_amax": 0.0}

    # --- Grid 1: core features without gap tuning ---
    vsm_vals = [0.0, 1.0, 2.0, 3.0]
    consistency_vals = [0, 5, 10]
    cweight_vals = [1.0, 2.0, 4.0]
    sw_vals = [0.0, 1.0]
    refine_vals = ["none", "confident"]

    for vsm, cons, sw, ref in product(vsm_vals, consistency_vals, sw_vals, refine_vals):
        # Skip redundant: consistency=0 doesn't need weight variations
        kw = {"vsm_amax": vsm}
        name_parts = [f"v{vsm:.0f}"]

        if cons > 0:
            # Only sweep weights for non-zero consistency
            for cw in cweight_vals:
                kw2 = dict(kw, consistency=cons, consistency_weight=cw)
                name2 = f"v{vsm:.0f}_c{cons}_w{cw:.0f}"
                if sw > 0:
                    kw2["seq_weights"] = sw
                    name2 += "_sw"
                if ref != "none":
                    kw2["refine"] = ref
                    name2 += "_ref"
                configs[name2] = kw2
        else:
            kw2 = dict(kw)
            name2 = f"v{vsm:.0f}"
            if sw > 0:
                kw2["seq_weights"] = sw
                name2 += "_sw"
            if ref != "none":
                kw2["refine"] = ref
                name2 += "_ref"
            configs[name2] = kw2

    # --- Grid 2: gap penalty sweep with best-guess features ---
    # Test gap penalties with vsm=2, consistency=5/w=2
    gpo_vals = [-3.0, -5.0, -8.0]   # default is around -6 for protein
    gpe_vals = [-0.5, -1.0, -2.0]   # default is around -0.8
    for gpo, gpe in product(gpo_vals, gpe_vals):
        name = f"v2_c5_w2_gpo{gpo:.0f}_gpe{gpe:.1f}"
        configs[name] = {
            "vsm_amax": 2.0,
            "consistency": 5,
            "consistency_weight": 2.0,
            "gap_open": gpo,
            "gap_extend": gpe,
        }

    # Deduplicate (some combos are identical)
    seen = {}
    deduped = {}
    for name, kw in configs.items():
        key = tuple(sorted(kw.items()))
        if key not in seen:
            seen[key] = name
            deduped[name] = kw
    configs = deduped

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
            if done % 500 == 0:
                elapsed = time.perf_counter() - t0
                print(f"  {done}/{len(tasks)} ({elapsed:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\nDone in {elapsed:.0f}s\n")

    # Aggregate
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
        sp = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        t = sum(r["wall_time"] for r in entries)
        summary.append((name, sp, prec, f1, tc, t, len(entries)))

    summary.sort(key=lambda x: -x[3])  # sort by F1

    # Overall top 30
    print("=" * 110)
    print("  TOP 30 BY F1")
    print("=" * 110)
    print(f"{'#':>3} {'Config':<32} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8} {'N':>4}")
    print("-" * 110)
    for i, (name, sp, prec, f1, tc, t, n) in enumerate(summary[:30]):
        print(f"{i+1:>3} {name:<32} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s {n:>4}")

    # Top 30 by TC
    summary_tc = sorted(summary, key=lambda x: -x[4])
    print(f"\n{'=' * 110}")
    print("  TOP 30 BY TC")
    print(f"{'=' * 110}")
    print(f"{'#':>3} {'Config':<32} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8} {'N':>4}")
    print("-" * 110)
    for i, (name, sp, prec, f1, tc, t, n) in enumerate(summary_tc[:30]):
        print(f"{i+1:>3} {name:<32} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s {n:>4}")

    # Per-category for top 10 F1 configs
    categories = ["RV11", "RV12", "RV20", "RV30", "RV40", "RV50"]
    by_cat = defaultdict(lambda: defaultdict(list))
    for r in results:
        if "error" not in r:
            cat = r["dataset"].replace("balibase_", "")
            by_cat[cat][r["config"]].append(r)

    top10 = [s[0] for s in summary[:10]]
    for metric_name, metric_key in [("F1", "f1"), ("TC", "tc")]:
        print(f"\n{'=' * 140}")
        print(f"  Per-category {metric_name} (top 10 configs by overall F1)")
        print(f"{'=' * 140}")
        header = f"{'Config':<32}" + "".join(f" {cat:>10}" for cat in categories) + f" {'Overall':>10}"
        print(header)
        print("-" * len(header))
        for name in top10:
            row = f"{name:<32}"
            for cat in categories:
                vals = by_cat[cat].get(name, [])
                if vals:
                    v = statistics.mean(r[metric_key] for r in vals)
                    row += f" {v:>10.3f}"
                else:
                    row += f" {'—':>10}"
            overall_f1 = next(s[3] for s in summary if s[0] == name)
            overall_tc = next(s[4] for s in summary if s[0] == name)
            val = overall_f1 if metric_key == "f1" else overall_tc
            row += f" {val:>10.3f}"
            print(row)

    # Gap penalty analysis
    gap_configs = [(n, sp, pr, f1, tc, t, cnt) for n, sp, pr, f1, tc, t, cnt in summary if "gpo" in n]
    if gap_configs:
        gap_configs.sort(key=lambda x: -x[3])
        print(f"\n{'=' * 110}")
        print("  GAP PENALTY SWEEP (with vsm=2, consistency=5, weight=2)")
        print(f"{'=' * 110}")
        print(f"{'Config':<32} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
        print("-" * 80)
        for name, sp, prec, f1, tc, t, n in gap_configs:
            print(f"{name:<32} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s")


if __name__ == "__main__":
    main()
