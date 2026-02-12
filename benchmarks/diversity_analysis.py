"""Measure diversity and quality of ensemble parameter presets.

For each preset, runs kalign independently (no ensemble, no tree noise)
with the scaled gap penalties. Measures SP score against reference and
pairwise agreement between presets.

Usage:
    uv run python -m benchmarks.diversity_analysis [--max-cases 50] [-j 8]
"""

import argparse
import tempfile
import statistics
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kalign

from .datasets import get_cases

# Same presets as ensemble.c (gap scaling factors)
PRESETS = [
    {"name": "default",          "gpo": 1.0, "gpe": 1.0, "tgpe": 1.0},
    {"name": "low-open/hi-ext",  "gpo": 0.5, "gpe": 1.5, "tgpe": 0.8},
    {"name": "hi-open/low-ext",  "gpo": 1.5, "gpe": 0.5, "tgpe": 1.2},
    {"name": "relaxed-all",      "gpo": 0.7, "gpe": 0.7, "tgpe": 0.5},
    {"name": "strict-all",       "gpo": 1.4, "gpe": 1.4, "tgpe": 1.5},
    {"name": "open-relax",       "gpo": 0.8, "gpe": 1.2, "tgpe": 1.0},
    {"name": "strict-open",      "gpo": 1.3, "gpe": 0.8, "tgpe": 0.7},
    {"name": "relax-open-only",  "gpo": 0.6, "gpe": 1.0, "tgpe": 1.3},
    {"name": "relax-ext-only",   "gpo": 1.0, "gpe": 0.6, "tgpe": 0.6},
    {"name": "v-strict-open",    "gpo": 1.8, "gpe": 1.0, "tgpe": 1.0},
    {"name": "v-strict-ext",     "gpo": 1.0, "gpe": 1.8, "tgpe": 1.8},
    {"name": "v-relaxed-all",    "gpo": 0.4, "gpe": 0.4, "tgpe": 0.3},
]

# Default protein gap penalties (from aln_param.c CorBLOSUM66)
BASE_GPO = 5.5
BASE_GPE = 2.0
BASE_TGPE = 1.0


def _run_one_preset(args):
    """Run one preset on one case. Returns (case_family, preset_idx, sp_score)."""
    case, preset_idx, preset = args
    gpo = BASE_GPO * preset["gpo"]
    gpe = BASE_GPE * preset["gpe"]
    tgpe = BASE_TGPE * preset["tgpe"]

    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / "aln.fasta"
        kalign.align_file_to_file(
            str(case.unaligned),
            str(output),
            format="fasta",
            seq_type=case.seq_type,
            n_threads=1,
            gap_open=gpo,
            gap_extend=gpe,
            terminal_gap_extend=tgpe,
        )
        sp = kalign.compare(str(case.reference), str(output))
    return case.family, preset_idx, sp


def main():
    parser = argparse.ArgumentParser(description="Analyze ensemble preset diversity")
    parser.add_argument("--max-cases", type=int, default=50)
    parser.add_argument("-j", "--parallel", type=int, default=1)
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases)
    n_presets = len(PRESETS)
    print(f"Running {n_presets} presets on {len(cases)} cases...")

    # Build work items: every (case, preset) combination
    work = []
    for case in cases:
        for pi, preset in enumerate(PRESETS):
            work.append((case, pi, preset))

    # Run all
    results = {}  # (family, preset_idx) -> sp_score
    done = 0
    total = len(work)

    if args.parallel <= 1:
        for item in work:
            fam, pi, sp = _run_one_preset(item)
            results[(fam, pi)] = sp
            done += 1
            if done % n_presets == 0:
                print(f"  [{done}/{total}] {fam} done", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.parallel) as pool:
            futures = {pool.submit(_run_one_preset, item): item for item in work}
            for future in as_completed(futures):
                fam, pi, sp = future.result()
                results[(fam, pi)] = sp
                done += 1
                if done % (n_presets * 2) == 0:
                    print(f"  [{done}/{total}] ...", flush=True)

    families = [c.family for c in cases]

    # === Analysis 1: Per-preset quality ===
    print(f"\n{'Preset':<20s} {'Mean SP':>8s} {'Med SP':>8s} {'vs Def':>8s} {'Wins':>5s} {'Loss':>5s}")
    print("-" * 60)
    default_scores = {fam: results[(fam, 0)] for fam in families}

    for pi, preset in enumerate(PRESETS):
        scores = [results[(fam, pi)] for fam in families]
        deltas = [results[(fam, pi)] - default_scores[fam] for fam in families]
        wins = sum(1 for d in deltas if d > 0.5)
        losses = sum(1 for d in deltas if d < -0.5)
        mean_delta = statistics.mean(deltas)
        print(f"{preset['name']:<20s} {statistics.mean(scores):8.2f} {statistics.median(scores):8.2f} "
              f"{mean_delta:+8.2f} {wins:5d} {losses:5d}")

    # === Analysis 2: Pairwise agreement between presets ===
    # For each pair of presets, count how often they produce the same SP score
    # (within 0.5 points = effectively same alignment)
    print(f"\n=== Pairwise agreement (% of cases with <0.5 SP difference) ===")
    print(f"{'':>20s}", end="")
    for pi in range(n_presets):
        print(f" {pi:>4d}", end="")
    print()

    for pi in range(n_presets):
        print(f"{PRESETS[pi]['name']:>20s}", end="")
        for pj in range(n_presets):
            agree = sum(1 for fam in families
                        if abs(results[(fam, pi)] - results[(fam, pj)]) < 0.5)
            pct = 100 * agree / len(families)
            print(f" {pct:4.0f}", end="")
        print()

    # === Analysis 3: "Oracle" best — what if we always picked the best preset? ===
    oracle_scores = []
    for fam in families:
        best = max(results[(fam, pi)] for pi in range(n_presets))
        oracle_scores.append(best)

    print(f"\nOracle (best of {n_presets} presets per case):")
    print(f"  Mean SP: {statistics.mean(oracle_scores):.2f}  "
          f"(default: {statistics.mean(default_scores.values()):.2f}, "
          f"delta: {statistics.mean(oracle_scores) - statistics.mean(default_scores.values()):+.2f})")

    # === Analysis 4: Per-preset "unique wins" ===
    # Cases where this preset is the ONLY one that beats the default
    print(f"\n{'Preset':<20s} {'Unique wins':>12s} {'Best count':>12s}")
    print("-" * 48)
    for pi, preset in enumerate(PRESETS):
        unique_wins = 0
        best_count = 0
        for fam in families:
            sp_pi = results[(fam, pi)]
            is_best = all(sp_pi >= results[(fam, pj)] - 0.01 for pj in range(n_presets))
            if is_best:
                best_count += 1
            # Unique win: beats default by >0.5 and no other preset beats default by >0.5
            if sp_pi > default_scores[fam] + 0.5:
                others_win = any(
                    results[(fam, pj)] > default_scores[fam] + 0.5
                    for pj in range(n_presets) if pj != pi
                )
                if not others_win:
                    unique_wins += 1
        print(f"{preset['name']:<20s} {unique_wins:12d} {best_count:12d}")


if __name__ == "__main__":
    main()
