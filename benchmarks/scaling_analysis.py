#!/usr/bin/env python3
"""Analyze how alignment quality degrades with number of sequences.

Loads existing benchmark results and correlates F1/TC/precision/recall
with sequence count per BAliBASE case. Compares kalign modes vs
external tools (MAFFT, MUSCLE, ClustalO).
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

RESULTS_DIR = Path(__file__).parent / "results"
BB_DIR = Path(__file__).parent / "data" / "downloads" / "bb3_release"


def count_sequences(tfa_path: Path) -> tuple[int, float]:
    """Return (num_sequences, mean_ungapped_length) from a .tfa file."""
    nseq = 0
    lengths = []
    buf = []
    with open(tfa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if buf:
                    seq = "".join(buf).replace("-", "").replace(".", "")
                    lengths.append(len(seq))
                nseq += 1
                buf = []
            elif line:
                buf.append(line)
    if buf:
        seq = "".join(buf).replace("-", "").replace(".", "")
        lengths.append(len(seq))
    mean_len = sum(lengths) / len(lengths) if lengths else 0
    return nseq, mean_len


def get_case_metadata() -> dict:
    """Get nseq and mean_len for every BAliBASE case."""
    meta = {}
    for tfa in sorted(BB_DIR.rglob("BB*.tfa")):
        family = tfa.stem
        nseq, mean_len = count_sequences(tfa)
        # Get RV category from path
        for part in tfa.parts:
            if part.startswith("RV"):
                rv = part
                break
        else:
            rv = "unknown"
        meta[family] = {"nseq": nseq, "mean_len": mean_len, "rv": rv}
    return meta


def load_results() -> dict[str, list[dict]]:
    """Load results from multiple JSON files, keyed by method label."""
    results = defaultdict(list)

    # 1. full_comparison.json — external tools + ens5 kalign
    with open(RESULTS_DIR / "full_comparison.json") as f:
        data = json.load(f)
    for r in data["results"]:
        method = r["method"]
        if method in ("clustalo", "mafft", "muscle"):
            label = method
        else:
            ens = r.get("ensemble", 0)
            ref = r.get("refine", "none")
            label = f"kalign_ens{ens}_ref={ref}"
        results[label].append(r)

    # 2. combined_improvements.json — kalign ablation configs
    with open(RESULTS_DIR / "combined_improvements.json") as f:
        data = json.load(f)
    for r in data["results"]:
        label = f"kalign_{r['config']}"
        results[label].append(r)

    # 3. best_mode_seqweights.json — best modes with seq weights
    with open(RESULTS_DIR / "best_mode_seqweights.json") as f:
        data = json.load(f)
    for r in data:
        label = f"kalign_{r['method']}"
        results[label].append(r)

    return dict(results)


def bin_label(nseq: int) -> str:
    """Assign a sequence count bin."""
    if nseq <= 5:
        return "2-5"
    elif nseq <= 10:
        return "6-10"
    elif nseq <= 20:
        return "11-20"
    elif nseq <= 50:
        return "21-50"
    elif nseq <= 100:
        return "51-100"
    else:
        return "101+"


def analyze(methods_to_show=None):
    meta = get_case_metadata()
    all_results = load_results()

    if methods_to_show is None:
        # Key comparisons: baseline, best single, best ensemble, competitors
        methods_to_show = [
            "kalign_baseline",
            "kalign_+vsm+ref+re1",
            "kalign_ens3+vsm+ref+ra1 (old)",
            "mafft",
            "muscle",
            "clustalo",
        ]

    # Verify methods exist
    available = sorted(all_results.keys())
    for m in methods_to_show:
        if m not in all_results:
            print(f"WARNING: '{m}' not found. Available: {available}")

    bin_order = ["2-5", "6-10", "11-20", "21-50", "51-100", "101+"]

    # Per-case detailed output
    print("=" * 120)
    print("PER-CASE DETAILS (kalign baseline)")
    print("=" * 120)
    print(f"{'Family':<10} {'RV':<5} {'Nseq':>5} {'MeanLen':>8} {'F1':>6} {'TC':>6} {'Recall':>7} {'Prec':>6}")
    print("-" * 120)

    baseline_key = "kalign_baseline"
    if baseline_key in all_results:
        case_data = {r["family"]: r for r in all_results[baseline_key]}
        for family in sorted(meta.keys()):
            m = meta[family]
            if family in case_data:
                r = case_data[family]
                print(
                    f"{family:<10} {m['rv']:<5} {m['nseq']:>5} {m['mean_len']:>8.0f} "
                    f"{r.get('f1', 0):>6.3f} {r.get('tc', 0):>6.3f} "
                    f"{r.get('recall', 0):>7.3f} {r.get('precision', 0):>6.3f}"
                )

    # Binned summary for each method
    print()
    print("=" * 120)
    print("QUALITY BY SEQUENCE COUNT BINS")
    print("=" * 120)

    for metric in ["f1", "tc", "recall", "precision"]:
        print(f"\n--- {metric.upper()} ---")
        header = f"{'Method':<35}"
        for b in bin_order:
            header += f" {b:>8}"
        header += f" {'ALL':>8}"
        print(header)
        print("-" * len(header))

        for method_label in methods_to_show:
            if method_label not in all_results:
                continue
            # Bin the results
            bins = defaultdict(list)
            all_vals = []
            for r in all_results[method_label]:
                family = r["family"]
                if family not in meta:
                    continue
                val = r.get(metric)
                if val is None:
                    continue
                nseq = meta[family]["nseq"]
                b = bin_label(nseq)
                bins[b].append(val)
                all_vals.append(val)

            row = f"{method_label:<35}"
            for b in bin_order:
                if bins[b]:
                    mean = sum(bins[b]) / len(bins[b])
                    row += f" {mean:>7.3f}({len(bins[b]):>1})"  # just truncate for space
                else:
                    row += f" {'—':>8}"
            if all_vals:
                overall = sum(all_vals) / len(all_vals)
                row += f" {overall:>8.3f}"
            print(row)

    # Distribution of case sizes
    print()
    print("=" * 80)
    print("CASE SIZE DISTRIBUTION")
    print("=" * 80)
    size_bins = defaultdict(list)
    for family, m in meta.items():
        b = bin_label(m["nseq"])
        size_bins[b].append((family, m["nseq"], m["mean_len"], m["rv"]))

    for b in bin_order:
        cases = size_bins[b]
        if not cases:
            continue
        nseqs = [c[1] for c in cases]
        print(
            f"\n{b}: {len(cases)} cases, "
            f"nseq range [{min(nseqs)}-{max(nseqs)}], "
            f"mean nseq={sum(nseqs)/len(nseqs):.1f}"
        )
        # RV breakdown
        from collections import Counter
        rv_counts = Counter(c[3] for c in cases)
        print(f"  RV breakdown: {dict(sorted(rv_counts.items()))}")

    # Correlation: nseq vs F1 for each method
    print()
    print("=" * 80)
    print("CORRELATION: nseq vs F1 (Pearson r)")
    print("=" * 80)
    for method_label in methods_to_show:
        if method_label not in all_results:
            continue
        xs, ys = [], []
        for r in all_results[method_label]:
            family = r["family"]
            if family not in meta or r.get("f1") is None:
                continue
            xs.append(meta[family]["nseq"])
            ys.append(r["f1"])
        if len(xs) > 2:
            # Pearson correlation
            n = len(xs)
            mx = sum(xs) / n
            my = sum(ys) / n
            cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
            sx = (sum((x - mx) ** 2 for x in xs) / n) ** 0.5
            sy = (sum((y - my) ** 2 for y in ys) / n) ** 0.5
            r_val = cov / (sx * sy) if sx > 0 and sy > 0 else 0
            print(f"  {method_label:<35} r={r_val:+.3f}  (n={n})")

    # Delta analysis: where does kalign lose most vs MAFFT?
    print()
    print("=" * 80)
    print("F1 DELTA: kalign_baseline vs mafft (sorted by delta, worst first)")
    print("=" * 80)
    if "kalign_baseline" in all_results and "mafft" in all_results:
        kalign_f1 = {r["family"]: r.get("f1", 0) for r in all_results["kalign_baseline"]}
        mafft_f1 = {r["family"]: r.get("f1", 0) for r in all_results["mafft"]}

        deltas = []
        for family in kalign_f1:
            if family in mafft_f1 and family in meta:
                delta = kalign_f1[family] - mafft_f1[family]
                deltas.append((family, meta[family]["nseq"], meta[family]["rv"],
                               kalign_f1[family], mafft_f1[family], delta))

        deltas.sort(key=lambda x: x[5])
        print(f"{'Family':<10} {'Nseq':>5} {'RV':<5} {'kalign':>7} {'mafft':>7} {'delta':>7}")
        print("-" * 50)
        for family, nseq, rv, kf1, mf1, delta in deltas[:30]:
            print(f"{family:<10} {nseq:>5} {rv:<5} {kf1:>7.3f} {mf1:>7.3f} {delta:>+7.3f}")

        # Binned delta
        print(f"\n--- Mean F1 delta (kalign - mafft) by nseq bin ---")
        bin_deltas = defaultdict(list)
        for family, nseq, rv, kf1, mf1, delta in deltas:
            bin_deltas[bin_label(nseq)].append(delta)
        for b in bin_order:
            if bin_deltas[b]:
                md = sum(bin_deltas[b]) / len(bin_deltas[b])
                print(f"  {b:>8}: {md:+.3f}  (n={len(bin_deltas[b])})")


if __name__ == "__main__":
    analyze()
