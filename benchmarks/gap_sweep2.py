"""Proper gap penalty sweep around PFASUM43 defaults.

PFASUM43 defaults: gpo=5.5, gpe=2.0, tgpe=1.0
Substitution scores range: roughly [-5, 13]
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

    configs = {}

    # Reference: defaults (gpo=5.5, gpe=2.0, tgpe=1.0)
    base = {"vsm_amax": 2.0, "consistency": 5, "consistency_weight": 2.0}
    configs["default(5.5/2.0/1.0)"] = dict(base)

    # GPO sweep (fix gpe=2.0, tgpe=1.0)
    for gpo in [2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0]:
        name = f"gpo{gpo:.1f}"
        configs[name] = dict(base, gap_open=gpo, gap_extend=2.0, terminal_gap_extend=1.0)

    # GPE sweep (fix gpo=5.5, tgpe=1.0)
    for gpe in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
        name = f"gpe{gpe:.1f}"
        configs[name] = dict(base, gap_open=5.5, gap_extend=gpe, terminal_gap_extend=1.0)

    # TGPE sweep (fix gpo=5.5, gpe=2.0)
    for tgpe in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]:
        name = f"tgpe{tgpe:.1f}"
        configs[name] = dict(base, gap_open=5.5, gap_extend=2.0, terminal_gap_extend=tgpe)

    # Joint GPO x GPE grid around default
    for gpo, gpe in product([3.0, 4.0, 5.0, 5.5, 6.5, 8.0],
                             [1.0, 1.5, 2.0, 2.5, 3.0]):
        name = f"j_gpo{gpo:.1f}_gpe{gpe:.1f}"
        configs[name] = dict(base, gap_open=gpo, gap_extend=gpe, terminal_gap_extend=1.0)

    # Same without consistency (to check interaction)
    base_nocons = {"vsm_amax": 2.0}
    configs["nocons_default"] = dict(base_nocons)
    for gpo, gpe in product([3.0, 4.0, 5.5, 7.0],
                             [1.0, 2.0, 3.0]):
        name = f"nocons_gpo{gpo:.1f}_gpe{gpe:.1f}"
        configs[name] = dict(base_nocons, gap_open=gpo, gap_extend=gpe, terminal_gap_extend=1.0)

    # Deduplicate
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

    # With consistency - by F1
    cons = [(n, sp, pr, f1, tc, t, cnt) for n, sp, pr, f1, tc, t, cnt in summary if "nocons" not in n]
    cons.sort(key=lambda x: -x[3])
    print("=" * 110)
    print("  WITH CONSISTENCY — sorted by F1")
    print("=" * 110)
    print(f"{'#':>3} {'Config':<28} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
    print("-" * 85)
    for i, (name, sp, prec, f1, tc, t, n) in enumerate(cons):
        print(f"{i+1:>3} {name:<28} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s")

    # With consistency - by TC
    cons_tc = sorted(cons, key=lambda x: -x[4])
    print(f"\n{'=' * 110}")
    print("  WITH CONSISTENCY — sorted by TC")
    print(f"{'=' * 110}")
    print(f"{'#':>3} {'Config':<28} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
    print("-" * 85)
    for i, (name, sp, prec, f1, tc, t, n) in enumerate(cons_tc[:20]):
        print(f"{i+1:>3} {name:<28} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s")

    # Without consistency
    nocons = [(n, sp, pr, f1, tc, t, cnt) for n, sp, pr, f1, tc, t, cnt in summary if "nocons" in n]
    nocons.sort(key=lambda x: -x[3])
    print(f"\n{'=' * 110}")
    print("  WITHOUT CONSISTENCY — sorted by F1")
    print(f"{'=' * 110}")
    print(f"{'#':>3} {'Config':<28} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8}")
    print("-" * 85)
    for i, (name, sp, prec, f1, tc, t, n) in enumerate(nocons):
        print(f"{i+1:>3} {name:<28} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s")


if __name__ == "__main__":
    main()
