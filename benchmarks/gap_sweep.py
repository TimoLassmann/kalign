"""Gap penalty sweep for single-run kalign on BAliBASE.

Gap penalties in kalign are POSITIVE numbers (internally negated).
The C default for PFASUM43 is gpo=55, gpe=5.5, tgpe=5.5
(in scaled integer units embedded in the profile vectors).

But the user-facing API uses smaller float values.
Let's first discover what the actual default scale is.
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

    # Reference: best config with default gaps
    configs["default_gaps"] = {
        "vsm_amax": 2.0,
        "consistency": 5,
        "consistency_weight": 2.0,
    }

    # Gap penalty sweep with POSITIVE values
    # PFASUM43 defaults from aln_param.c: gpo=55, gpe=5.5, tgpe=5.5
    # But those are in profile-embedded units. The user API gpo/gpe are
    # applied via `if(gpo >= 0.0) ap->gpo = gpo;` so same scale.
    # Let's sweep around the defaults.
    gpo_vals = [20, 35, 55, 75, 100]
    gpe_vals = [2.0, 4.0, 5.5, 8.0, 12.0]
    tgpe_vals = [5.5]  # keep tgpe fixed at default for now

    for gpo, gpe in product(gpo_vals, gpe_vals):
        name = f"gpo{gpo}_gpe{gpe:.1f}"
        configs[name] = {
            "vsm_amax": 2.0,
            "consistency": 5,
            "consistency_weight": 2.0,
            "gap_open": float(gpo),
            "gap_extend": float(gpe),
        }

    # Also try without consistency to see if gap penalties interact
    for gpo, gpe in product([35, 55, 75], [4.0, 5.5, 8.0]):
        name = f"nocons_gpo{gpo}_gpe{gpe:.1f}"
        configs[name] = {
            "vsm_amax": 2.0,
            "gap_open": float(gpo),
            "gap_extend": float(gpe),
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
            if done % 500 == 0:
                elapsed = time.perf_counter() - t0
                print(f"  {done}/{len(tasks)} ({elapsed:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\nDone in {elapsed:.0f}s\n")

    # Check for errors
    errors = [r for r in results if "error" in r]
    if errors:
        print(f"WARNING: {len(errors)} errors")
        for e in errors[:5]:
            print(f"  {e['config']}: {e['error']}")
        print()

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

    # With consistency
    cons_configs = [(n, sp, pr, f1, tc, t, cnt) for n, sp, pr, f1, tc, t, cnt in summary if "nocons" not in n]
    cons_configs.sort(key=lambda x: -x[3])
    print("=" * 110)
    print("  WITH CONSISTENCY (vsm=2, c=5, w=2) — sorted by F1")
    print("=" * 110)
    print(f"{'#':>3} {'Config':<24} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8} {'N':>4}")
    print("-" * 110)
    for i, (name, sp, prec, f1, tc, t, n) in enumerate(cons_configs):
        print(f"{i+1:>3} {name:<24} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s {n:>4}")

    # Without consistency
    nocons_configs = [(n, sp, pr, f1, tc, t, cnt) for n, sp, pr, f1, tc, t, cnt in summary if "nocons" in n]
    nocons_configs.sort(key=lambda x: -x[3])
    print(f"\n{'=' * 110}")
    print("  WITHOUT CONSISTENCY (vsm=2 only) — sorted by F1")
    print(f"{'=' * 110}")
    print(f"{'#':>3} {'Config':<24} {'SP':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>8} {'N':>4}")
    print("-" * 110)
    for i, (name, sp, prec, f1, tc, t, n) in enumerate(nocons_configs):
        print(f"{i+1:>3} {name:<24} {sp:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {t:>7.0f}s {n:>4}")


if __name__ == "__main__":
    main()
