"""Sweep vsm_amax and dist_scale on BRAliBASE (RNA).

RNA substitution matrix scores are in the hundreds (base ~283-383), so
vsm_amax values need to be much larger than for protein (where scores
are -4 to 13) to have proportional effect.

dist_scale controls distance-dependent gap penalty scaling:
  scale = max(0.3, min(1.0, 1.0 - dist_scale * avg_divergence))
  So dist_scale=0.5 with avg_div=1.0 gives scale=0.5 (halved gap penalties).

Usage:
    uv run python -m benchmarks.rna_param_sweep
    uv run python -m benchmarks.rna_param_sweep -j 4
    uv run python -m benchmarks.rna_param_sweep --sweep dist_scale
    uv run python -m benchmarks.rna_param_sweep --sweep both
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
from .scoring import _fasta_ref_has_gaps


def _score_case(case, output_path):
    """Score alignment, skipping gapless references."""
    if case.reference.suffix in ('.fa', '.fasta') and not _fasta_ref_has_gaps(case.reference):
        raise RuntimeError("Reference alignment has no gaps — skipping")
    return kalign.compare_detailed(str(case.reference), str(output_path), max_gap_frac=-1.0)


def _run_one(args):
    case, vsm_amax, dist_scale = args
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type, vsm_amax=vsm_amax, dist_scale=dist_scale,
        )
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "vsm_amax": vsm_amax, "dist_scale": dist_scale,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "vsm_amax": vsm_amax, "dist_scale": dist_scale,
            "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def print_summary(results, param_key):
    """Print summary table grouped by parameter value."""
    groups = defaultdict(list)
    for r in results:
        if "error" not in r:
            groups[r[param_key]].append(r)

    print(f"\n{'':>8} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}  n")
    print("-" * 52)
    for val in sorted(groups.keys()):
        entries = groups[val]
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        print(f"{val:>8.1f} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}  {len(entries)}")

    # Per RNA family
    all_families = sorted({r["dataset"] for r in results if "error" not in r})
    for fam in all_families:
        fam_results = [r for r in results if r["dataset"] == fam and "error" not in r]
        fam_groups = defaultdict(list)
        for r in fam_results:
            fam_groups[r[param_key]].append(r)

        print(f"\n=== {fam} ({len(fam_results) // len(fam_groups)} cases) ===")
        print(f"{'':>8} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 44)
        for val in sorted(fam_groups.keys()):
            entries = fam_groups[val]
            rec = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            print(f"{val:>8.1f} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")


def main():
    parser = argparse.ArgumentParser(description="RNA parameter sweep on BRAliBASE")
    parser.add_argument("-j", "--parallel", type=int, default=4)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--sweep", choices=["vsm_amax", "dist_scale", "both"],
                        default="both", help="Which parameter(s) to sweep")
    args = parser.parse_args()

    cases = get_cases("bralibase", max_cases=args.max_cases if args.max_cases else None)
    print(f"{len(cases)} BRAliBASE cases")

    # Build parameter grid
    # RNA matrix scores are ~283-383. For proportional effect similar to protein:
    #   protein amax=2.0 on score range [-4,13] is ~12% of max
    #   RNA equivalent: 0.12 * 383 ≈ 46
    # But gap penalties are also much larger (gpo=217, gpe=39.4)
    # Test a wide range to see what happens
    vsm_values = [0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]
    ds_values = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]

    tasks = []
    if args.sweep in ("vsm_amax", "both"):
        for case in cases:
            for amax in vsm_values:
                tasks.append((case, amax, 0.0))  # sweep vsm_amax, dist_scale=0
    if args.sweep in ("dist_scale", "both"):
        for case in cases:
            for ds in ds_values:
                if ds == 0.0 and args.sweep == "both":
                    continue  # already covered by vsm sweep baseline
                tasks.append((case, 0.0, ds))  # sweep dist_scale, vsm_amax=0 (RNA default)

    print(f"{len(tasks)} total tasks")

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
    print(f"All done in {elapsed:.0f}s")

    # Separate results by sweep type
    vsm_results = [r for r in results if r["dist_scale"] == 0.0]
    ds_results = [r for r in results if r["vsm_amax"] == 0.0]

    if vsm_results:
        print("\n" + "=" * 60)
        print("VSM_AMAX SWEEP (dist_scale=0.0)")
        print("=" * 60)
        print_summary(vsm_results, "vsm_amax")

    if ds_results:
        print("\n" + "=" * 60)
        print("DIST_SCALE SWEEP (vsm_amax=0.0)")
        print("=" * 60)
        print_summary(ds_results, "dist_scale")

    # Save
    out = Path("benchmarks/results/rna_param_sweep.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
