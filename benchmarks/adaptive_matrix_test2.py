"""Refined adaptive matrix test — moderate configs, actual distances.

Uses kalign to measure average pairwise identity per family, then tests
distance-based parameter switching.

Usage:
    uv run python -m benchmarks.adaptive_matrix_test2 -j 8
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


CONFIGS = [
    # A: Close (>50% identity): PFASUM43, tight gaps, full VSM
    {
        "name": "A-close",
        "desc": "PFASUM43, gpo=8.0, gpe=1.5, vsm=2.5",
        "seq_type": "pfasum43",
        "gpo": 8.0, "gpe": 1.5, "tgpe": 1.0,
        "vsm_amax": 2.5,
    },
    # B: Medium (25-50%): PFASUM43, moderate gaps, moderate VSM
    {
        "name": "B-medium",
        "desc": "PFASUM43, gpo=7.0, gpe=1.5, vsm=2.0",
        "seq_type": "pfasum43",
        "gpo": 7.0, "gpe": 1.5, "tgpe": 1.0,
        "vsm_amax": 2.0,
    },
    # C: Divergent (<25%): PFASUM43, default gaps, reduced VSM
    {
        "name": "C-divergent",
        "desc": "PFASUM43, gpo=5.5, gpe=2.0, vsm=1.0",
        "seq_type": "pfasum43",
        "gpo": 5.5, "gpe": 2.0, "tgpe": 1.0,
        "vsm_amax": 1.0,
    },
    # D: Very divergent (<15%): PFASUM43, loose, no VSM
    {
        "name": "D-vdivergent",
        "desc": "PFASUM43, gpo=5.0, gpe=2.0, vsm=0",
        "seq_type": "pfasum43",
        "gpo": 5.0, "gpe": 2.0, "tgpe": 1.0,
        "vsm_amax": 0.0,
    },
    # E: Baselines
    {
        "name": "CorBLOSUM66-def",
        "desc": "Current default",
        "seq_type": "protein",
        "gpo": None, "gpe": None, "tgpe": None,
        "vsm_amax": -1.0,
    },
    {
        "name": "PFASUM43-g7",
        "desc": "Best fixed from sweep",
        "seq_type": "pfasum43",
        "gpo": 7.0, "gpe": 1.5, "tgpe": 1.0,
        "vsm_amax": -1.0,
    },
    # Extra: testing gap penalties in the divergent regime
    {
        "name": "E-div-g6",
        "desc": "PFASUM43, gpo=6.0, gpe=1.5, vsm=0.5",
        "seq_type": "pfasum43",
        "gpo": 6.0, "gpe": 1.5, "tgpe": 1.0,
        "vsm_amax": 0.5,
    },
]


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(case.reference), str(output_path), column_mask=mask)
    return kalign.compare_detailed(str(case.reference), str(output_path))


def _measure_avg_identity(fasta_path):
    """Quick measurement: align with defaults, compute avg pairwise identity."""
    import re
    seqs = []
    current = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                current = []
            else:
                current.append(line)
    if current:
        seqs.append("".join(current))

    if len(seqs) < 2:
        return 1.0

    # Sample up to 100 pairs
    import random
    pairs = []
    n = len(seqs)
    if n * (n - 1) // 2 <= 100:
        for i in range(n):
            for j in range(i + 1, n):
                pairs.append((i, j))
    else:
        rng = random.Random(42)
        seen = set()
        while len(pairs) < 100:
            i = rng.randint(0, n - 1)
            j = rng.randint(0, n - 1)
            if i != j and (i, j) not in seen:
                seen.add((i, j))
                pairs.append((i, j))

    # Compute k-mer identity (fast approximation using 3-mer Jaccard)
    def kmer_set(seq, k=3):
        return set(seq[i:i+k] for i in range(len(seq) - k + 1))

    identities = []
    for i, j in pairs:
        s1 = kmer_set(seqs[i])
        s2 = kmer_set(seqs[j])
        if not s1 or not s2:
            continue
        jaccard = len(s1 & s2) / len(s1 | s2)
        identities.append(jaccard)

    return statistics.mean(identities) if identities else 0.5


def _run_one(args):
    case, config = args
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        t0 = time.perf_counter()
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=config["seq_type"],
            gap_open=config["gpo"],
            gap_extend=config["gpe"],
            terminal_gap_extend=config["tgpe"],
            vsm_amax=config["vsm_amax"],
            refine="confident",
        )
        wall_time = time.perf_counter() - t0
        scores = _score_case(case, tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
            "wall_time": wall_time,
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "config": config["name"],
            "f1": 0, "recall": 0, "precision": 0, "tc": 0,
            "wall_time": 0, "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description="Adaptive matrix test v2")
    parser.add_argument("-j", "--parallel", type=int, default=4)
    args = parser.parse_args()

    cases = get_cases("balibase")

    # Measure identity for each family
    print("Measuring pairwise identity for each family...")
    family_identity = {}
    for case in cases:
        if case.family not in family_identity:
            ident = _measure_avg_identity(str(case.unaligned))
            family_identity[case.family] = ident

    # Show distribution
    print(f"\nIdentity distribution ({len(family_identity)} families):")
    bins = [(0, 0.15, "very div"), (0.15, 0.25, "divergent"),
            (0.25, 0.40, "medium"), (0.40, 0.60, "close"), (0.60, 1.01, "v.close")]
    for lo, hi, label in bins:
        count = sum(1 for v in family_identity.values() if lo <= v < hi)
        print(f"  {label:>10} ({lo:.0%}-{hi:.0%}): {count} families")

    print(f"\n{len(cases)} cases x {len(CONFIGS)} configs = {len(cases) * len(CONFIGS)} tasks")

    tasks = [(case, config) for case in cases for config in CONFIGS]

    t0 = time.perf_counter()
    results = []
    done = 0
    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(_run_one, t): t for t in tasks}
        for f in as_completed(futures):
            done += 1
            r = f.result()
            results.append(r)
            if done % 200 == 0:
                print(f"  {done}/{len(tasks)} ({time.perf_counter()-t0:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"All done in {elapsed:.0f}s\n")

    by_family = defaultdict(dict)
    for r in results:
        if "error" not in r:
            by_family[r["family"]][r["config"]] = r

    # Distance-based heuristic using measured identity
    def pick_config(identity):
        if identity >= 0.40:
            return "A-close"
        elif identity >= 0.25:
            return "B-medium"
        elif identity >= 0.15:
            return "E-div-g6"
        else:
            return "D-vdivergent"

    heuristic_results = []
    for fam, configs in by_family.items():
        ident = family_identity.get(fam, 0.5)
        pick = pick_config(ident)
        if pick in configs:
            heuristic_results.append(configs[pick])
        else:
            heuristic_results.append(configs["PFASUM43-g7"])

    # Oracle
    oracle_results = []
    oracle_picks = defaultdict(int)
    for fam, configs in by_family.items():
        best = max(configs.values(), key=lambda x: x["f1"])
        oracle_results.append(best)
        oracle_picks[best["config"]] += 1

    config_names = [c["name"] for c in CONFIGS]

    # Print overall
    print("=" * 100)
    print("OVERALL")
    print("=" * 100)
    print(f"{'Method':<25} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 58)
    for name in config_names:
        rows = [r for r in results if r.get("config") == name and "error" not in r]
        if rows:
            f1 = statistics.mean(r["f1"] for r in rows)
            rec = statistics.mean(r["recall"] for r in rows)
            prec = statistics.mean(r["precision"] for r in rows)
            tc = statistics.mean(r["tc"] for r in rows)
            print(f"  {name:<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    f1 = statistics.mean(r["f1"] for r in oracle_results)
    rec = statistics.mean(r["recall"] for r in oracle_results)
    prec = statistics.mean(r["precision"] for r in oracle_results)
    tc = statistics.mean(r["tc"] for r in oracle_results)
    print(f"  {'** Oracle **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    f1 = statistics.mean(r["f1"] for r in heuristic_results)
    rec = statistics.mean(r["recall"] for r in heuristic_results)
    prec = statistics.mean(r["precision"] for r in heuristic_results)
    tc = statistics.mean(r["tc"] for r in heuristic_results)
    print(f"  {'** Heuristic **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    # Per-category
    all_cats = sorted({r["dataset"] for r in results})
    for cat in all_cats:
        cat_results = [r for r in results if r["dataset"] == cat]
        n = len({r["family"] for r in cat_results})
        cat_label = cat.replace("balibase_", "")
        print(f"\n{'=' * 100}")
        print(f"{cat_label} ({n} cases)")
        print(f"{'=' * 100}")
        print(f"{'Method':<25} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 58)
        for name in config_names:
            rows = [r for r in cat_results if r.get("config") == name and "error" not in r]
            if rows:
                f1 = statistics.mean(r["f1"] for r in rows)
                rec = statistics.mean(r["recall"] for r in rows)
                prec = statistics.mean(r["precision"] for r in rows)
                tc = statistics.mean(r["tc"] for r in rows)
                print(f"  {name:<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

        # Oracle + heuristic for this cat
        cat_by_fam = defaultdict(dict)
        for r in cat_results:
            if "error" not in r:
                cat_by_fam[r["family"]][r["config"]] = r

        cat_oracle = [max(cs.values(), key=lambda x: x["f1"]) for cs in cat_by_fam.values()]
        f1 = statistics.mean(r["f1"] for r in cat_oracle)
        rec = statistics.mean(r["recall"] for r in cat_oracle)
        prec = statistics.mean(r["precision"] for r in cat_oracle)
        tc = statistics.mean(r["tc"] for r in cat_oracle)
        print(f"  {'** Oracle **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

        cat_heur = []
        for fam, configs in cat_by_fam.items():
            ident = family_identity.get(fam, 0.5)
            pick = pick_config(ident)
            if pick in configs:
                cat_heur.append(configs[pick])
            else:
                cat_heur.append(configs["PFASUM43-g7"])
        f1 = statistics.mean(r["f1"] for r in cat_heur)
        rec = statistics.mean(r["recall"] for r in cat_heur)
        prec = statistics.mean(r["precision"] for r in cat_heur)
        tc = statistics.mean(r["tc"] for r in cat_heur)
        print(f"  {'** Heuristic **':<23} {rec:8.3f} {prec:8.3f} {f1:8.3f} {tc:8.3f}")

    # Oracle picks
    print(f"\nOracle pick distribution:")
    for name in sorted(oracle_picks, key=oracle_picks.get, reverse=True):
        print(f"  {name:<23} {oracle_picks[name]:>4}")

    # Save
    out = Path("benchmarks/results/adaptive_matrix_test2.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump({
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "family_identity": family_identity,
            "results": results,
        }, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
