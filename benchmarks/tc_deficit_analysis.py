#!/usr/bin/env python3
"""
TC (Total Column) deficit analysis for kalign vs MAFFT/MUSCLE/ClustalO on BAliBASE.

Investigates why kalign's TC scores are disproportionately lower than its SP scores
compared to other tools.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path
import math

RESULTS_DIR = Path(__file__).parent / "results"


def load_data():
    """Load and index all three JSON result files."""
    # 1. Kalign baseline (refine=confident is the default mode)
    with open(RESULTS_DIR / "balibase_baseline.json") as f:
        baseline_raw = json.load(f)["results"]

    kalign = {}
    for r in baseline_raw:
        if r["refine"] == "confident" and r["error"] is None:
            kalign[r["family"]] = r

    # 2. External tools from mumsa_precision.json
    with open(RESULTS_DIR / "mumsa_precision.json") as f:
        mumsa_raw = json.load(f)["results"]

    external = defaultdict(dict)  # tool -> family -> record
    for r in mumsa_raw:
        if r["method"] in ("mafft", "muscle", "clustalo"):
            external[r["method"]][r["family"]] = r

    # 3. Overgap analysis (has nseq, mean_seq_len, category, default_ratio)
    with open(RESULTS_DIR / "overgap_analysis.json") as f:
        overgap_raw = json.load(f)["results"]

    overgap = {}
    for r in overgap_raw:
        if r["error"] is None:
            overgap[r["family"]] = r

    return kalign, external, overgap


def category_from_dataset(dataset):
    """Extract category like RV11 from dataset string like balibase_RV11."""
    return dataset.replace("balibase_", "")


def print_section(title):
    print()
    print("=" * 75)
    print(f"  {title}")
    print("=" * 75)
    print()


def main():
    kalign, external, overgap = load_data()
    mafft = external["mafft"]
    muscle = external["muscle"]
    clustalo = external["clustalo"]

    # Common families across all datasets
    all_families = sorted(
        set(kalign.keys()) & set(mafft.keys()) & set(muscle.keys())
        & set(clustalo.keys()) & set(overgap.keys())
    )
    print(f"Total families in common across all datasets: {len(all_families)}")

    # =========================================================================
    # 1. Overall stats
    # =========================================================================
    print_section("1. OVERALL AVERAGES")

    tools = {
        "kalign": {f: kalign[f] for f in all_families},
        "mafft": {f: mafft[f] for f in all_families},
        "muscle": {f: muscle[f] for f in all_families},
        "clustalo": {f: clustalo[f] for f in all_families},
    }

    print(f"{'Tool':<12} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'TC/SP':>8} {'TC/F1':>8}")
    print("-" * 68)
    for name, fam_data in tools.items():
        vals = list(fam_data.values())
        n = len(vals)
        avg_recall = sum(v["recall"] for v in vals) / n
        avg_prec = sum(v["precision"] for v in vals) / n
        avg_f1 = sum(v["f1"] for v in vals) / n
        avg_tc = sum(v["tc"] for v in vals) / n
        # SP = recall for BAliBASE scoring
        avg_sp = avg_recall
        tc_sp = avg_tc / avg_sp if avg_sp > 0 else 0
        tc_f1 = avg_tc / avg_f1 if avg_f1 > 0 else 0
        print(f"{name:<12} {avg_recall:>8.3f} {avg_prec:>8.3f} {avg_f1:>8.3f} "
              f"{avg_tc:>8.3f} {tc_sp:>8.3f} {tc_f1:>8.3f}")

    # =========================================================================
    # 2a. TC/SP ratio per family, grouped by category
    # =========================================================================
    print_section("2a. TC/SP RATIO BY CATEGORY (kalign vs MAFFT)")

    categories = defaultdict(list)
    for f in all_families:
        cat = category_from_dataset(kalign[f]["dataset"])
        categories[cat].append(f)

    print(f"{'Category':<10} {'N':>4}  "
          f"{'kalign TC':>10} {'kalign SP':>10} {'kalign TC/SP':>12}  "
          f"{'MAFFT TC':>10} {'MAFFT SP':>10} {'MAFFT TC/SP':>12}")
    print("-" * 100)
    for cat in sorted(categories.keys()):
        fams = categories[cat]
        n = len(fams)
        k_tc = sum(kalign[f]["tc"] for f in fams) / n
        k_sp = sum(kalign[f]["recall"] for f in fams) / n
        m_tc = sum(mafft[f]["tc"] for f in fams) / n
        m_sp = sum(mafft[f]["recall"] for f in fams) / n
        k_ratio = k_tc / k_sp if k_sp > 0 else 0
        m_ratio = m_tc / m_sp if m_sp > 0 else 0
        print(f"{cat:<10} {n:>4}  {k_tc:>10.3f} {k_sp:>10.3f} {k_ratio:>12.3f}  "
              f"{m_tc:>10.3f} {m_sp:>10.3f} {m_ratio:>12.3f}")

    # =========================================================================
    # 2b. TC=0 counts
    # =========================================================================
    print_section("2b. FAMILIES WITH TC=0")

    for name, fam_data in tools.items():
        tc0 = sum(1 for f in all_families if fam_data[f]["tc"] == 0.0)
        print(f"{name:<12}: {tc0:>3} families with TC=0 ({100*tc0/len(all_families):.1f}%)")

    # =========================================================================
    # 2c,d. TC distribution
    # =========================================================================
    print_section("2c/d. TC VALUE DISTRIBUTION")

    thresholds = [0.8, 0.5, 0.2, 0.0]
    print(f"{'Tool':<12} {'TC=0':>8} {'TC>0':>8} {'TC>0.2':>8} {'TC>0.5':>8} {'TC>0.8':>8}")
    print("-" * 55)
    for name, fam_data in tools.items():
        tc_vals = [fam_data[f]["tc"] for f in all_families]
        n = len(tc_vals)
        tc0 = sum(1 for v in tc_vals if v == 0.0)
        gt0 = sum(1 for v in tc_vals if v > 0.0)
        gt02 = sum(1 for v in tc_vals if v > 0.2)
        gt05 = sum(1 for v in tc_vals if v > 0.5)
        gt08 = sum(1 for v in tc_vals if v > 0.8)
        print(f"{name:<12} {tc0:>8} {gt0:>8} {gt02:>8} {gt05:>8} {gt08:>8}")

    print()
    print("As fractions:")
    print(f"{'Tool':<12} {'TC=0':>8} {'TC>0':>8} {'TC>0.2':>8} {'TC>0.5':>8} {'TC>0.8':>8}")
    print("-" * 55)
    for name, fam_data in tools.items():
        tc_vals = [fam_data[f]["tc"] for f in all_families]
        n = len(tc_vals)
        tc0 = sum(1 for v in tc_vals if v == 0.0) / n
        gt0 = sum(1 for v in tc_vals if v > 0.0) / n
        gt02 = sum(1 for v in tc_vals if v > 0.2) / n
        gt05 = sum(1 for v in tc_vals if v > 0.5) / n
        gt08 = sum(1 for v in tc_vals if v > 0.8) / n
        print(f"{name:<12} {tc0:>8.3f} {gt0:>8.3f} {gt02:>8.3f} {gt05:>8.3f} {gt08:>8.3f}")

    # =========================================================================
    # 3. TC-killed families
    # =========================================================================
    print_section("3. TC-KILLED FAMILIES (kalign TC=0 but MAFFT TC>0)")

    tc_killed = []
    for f in all_families:
        if kalign[f]["tc"] == 0.0 and mafft[f]["tc"] > 0.0:
            og = overgap.get(f, {})
            tc_killed.append({
                "family": f,
                "category": category_from_dataset(kalign[f]["dataset"]),
                "nseq": og.get("nseq", "?"),
                "mean_len": og.get("mean_seq_len", "?"),
                "overgap_ratio": og.get("default_ratio", "?"),
                "kalign_sp": kalign[f]["recall"],
                "mafft_tc": mafft[f]["tc"],
                "mafft_sp": mafft[f]["recall"],
            })

    print(f"Total TC-killed families: {len(tc_killed)}")
    print()

    # Summarize by category
    killed_by_cat = defaultdict(list)
    for t in tc_killed:
        killed_by_cat[t["category"]].append(t)

    print(f"{'Category':<10} {'Count':>6} {'Avg nseq':>10} {'Avg len':>10} {'Avg overgap':>12} {'Avg k_SP':>10} {'Avg m_TC':>10}")
    print("-" * 75)
    for cat in sorted(killed_by_cat.keys()):
        items = killed_by_cat[cat]
        n = len(items)
        avg_nseq = sum(t["nseq"] for t in items if t["nseq"] != "?") / max(1, sum(1 for t in items if t["nseq"] != "?"))
        avg_len = sum(t["mean_len"] for t in items if t["mean_len"] != "?") / max(1, sum(1 for t in items if t["mean_len"] != "?"))
        avg_og = sum(t["overgap_ratio"] for t in items if t["overgap_ratio"] != "?") / max(1, sum(1 for t in items if t["overgap_ratio"] != "?"))
        avg_ksp = sum(t["kalign_sp"] for t in items) / n
        avg_mtc = sum(t["mafft_tc"] for t in items) / n
        print(f"{cat:<10} {n:>6} {avg_nseq:>10.1f} {avg_len:>10.1f} {avg_og:>12.3f} {avg_ksp:>10.3f} {avg_mtc:>10.3f}")

    print()
    print("Worst TC-killed families (highest MAFFT TC):")
    tc_killed_sorted = sorted(tc_killed, key=lambda x: -x["mafft_tc"])
    print(f"{'Family':<12} {'Cat':<8} {'nseq':>5} {'len':>8} {'overgap':>8} {'k_SP':>8} {'m_TC':>8} {'m_SP':>8}")
    print("-" * 75)
    for t in tc_killed_sorted[:20]:
        print(f"{t['family']:<12} {t['category']:<8} {t['nseq']:>5} {t['mean_len']:>8.0f} "
              f"{t['overgap_ratio']:>8.3f} {t['kalign_sp']:>8.3f} {t['mafft_tc']:>8.3f} {t['mafft_sp']:>8.3f}")

    # =========================================================================
    # 4. TC gap excluding TC=0 families
    # =========================================================================
    print_section("4. TC GAP WHERE BOTH kalign AND MAFFT HAVE TC>0")

    both_gt0 = [f for f in all_families if kalign[f]["tc"] > 0 and mafft[f]["tc"] > 0]
    print(f"Families where both kalign TC>0 and MAFFT TC>0: {len(both_gt0)} / {len(all_families)}")

    if both_gt0:
        k_tc_avg = sum(kalign[f]["tc"] for f in both_gt0) / len(both_gt0)
        m_tc_avg = sum(mafft[f]["tc"] for f in both_gt0) / len(both_gt0)
        mu_tc_avg = sum(muscle[f]["tc"] for f in both_gt0) / len(both_gt0)
        cl_tc_avg = sum(clustalo[f]["tc"] for f in both_gt0) / len(both_gt0)

        k_sp_avg = sum(kalign[f]["recall"] for f in both_gt0) / len(both_gt0)
        m_sp_avg = sum(mafft[f]["recall"] for f in both_gt0) / len(both_gt0)

        print(f"  kalign   mean TC = {k_tc_avg:.3f}   mean SP = {k_sp_avg:.3f}   TC/SP = {k_tc_avg/k_sp_avg:.3f}")
        print(f"  MAFFT    mean TC = {m_tc_avg:.3f}   mean SP = {m_sp_avg:.3f}   TC/SP = {m_tc_avg/m_sp_avg:.3f}")
        print(f"  MUSCLE   mean TC = {mu_tc_avg:.3f}")
        print(f"  ClustalO mean TC = {cl_tc_avg:.3f}")
        print(f"  TC gap (MAFFT - kalign) = {m_tc_avg - k_tc_avg:.3f}")

        # Per-family TC gap distribution
        tc_gaps = [mafft[f]["tc"] - kalign[f]["tc"] for f in both_gt0]
        tc_gaps_sorted = sorted(tc_gaps)
        print()
        print(f"  Per-family TC gap (MAFFT - kalign) distribution:")
        print(f"    Mean:   {sum(tc_gaps)/len(tc_gaps):>+.3f}")
        print(f"    Median: {tc_gaps_sorted[len(tc_gaps_sorted)//2]:>+.3f}")
        print(f"    Min:    {tc_gaps_sorted[0]:>+.3f}")
        print(f"    Max:    {tc_gaps_sorted[-1]:>+.3f}")
        print(f"    kalign better: {sum(1 for g in tc_gaps if g < 0)} families")
        print(f"    MAFFT better:  {sum(1 for g in tc_gaps if g > 0)} families")
        print(f"    Equal:         {sum(1 for g in tc_gaps if g == 0)} families")

    # =========================================================================
    # 5. Correlations for kalign families with TC>0
    # =========================================================================
    print_section("5. CORRELATIONS (kalign families with TC>0)")

    k_tc_gt0 = [f for f in all_families if kalign[f]["tc"] > 0]
    print(f"Families with kalign TC>0: {len(k_tc_gt0)}")

    def pearson_r(xs, ys):
        n = len(xs)
        if n < 3:
            return float('nan')
        mx = sum(xs) / n
        my = sum(ys) / n
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
        sy = math.sqrt(sum((y - my) ** 2 for y in ys))
        if sx == 0 or sy == 0:
            return float('nan')
        return cov / (sx * sy)

    # Overgap ratio vs TC
    og_vals = []
    tc_vals = []
    for f in k_tc_gt0:
        if f in overgap and overgap[f].get("default_ratio") is not None:
            og_vals.append(overgap[f]["default_ratio"])
            tc_vals.append(kalign[f]["tc"])
    r_og_tc = pearson_r(og_vals, tc_vals)
    print(f"  Overgap ratio vs TC:  r = {r_og_tc:+.3f}  (n={len(og_vals)})")

    # nseq vs TC
    nseq_vals = []
    tc_vals2 = []
    for f in k_tc_gt0:
        if f in overgap and overgap[f].get("nseq") is not None:
            nseq_vals.append(overgap[f]["nseq"])
            tc_vals2.append(kalign[f]["tc"])
    r_nseq_tc = pearson_r(nseq_vals, tc_vals2)
    print(f"  nseq vs TC:           r = {r_nseq_tc:+.3f}  (n={len(nseq_vals)})")

    # SP vs TC for kalign
    sp_k = [kalign[f]["recall"] for f in k_tc_gt0]
    tc_k = [kalign[f]["tc"] for f in k_tc_gt0]
    r_sp_tc_kalign = pearson_r(sp_k, tc_k)
    print(f"  SP vs TC (kalign):    r = {r_sp_tc_kalign:+.3f}  (n={len(sp_k)})")

    # SP vs TC for mafft (families where MAFFT TC>0)
    m_tc_gt0 = [f for f in all_families if mafft[f]["tc"] > 0]
    sp_m = [mafft[f]["recall"] for f in m_tc_gt0]
    tc_m = [mafft[f]["tc"] for f in m_tc_gt0]
    r_sp_tc_mafft = pearson_r(sp_m, tc_m)
    print(f"  SP vs TC (MAFFT):     r = {r_sp_tc_mafft:+.3f}  (n={len(sp_m)})")

    # =========================================================================
    # 6. Is the TC problem driven by TC=0 families or a broad deficit?
    # =========================================================================
    print_section("6. TC PROBLEM: FEW FAILURES vs BROAD DEFICIT")

    # Exclude families where ANY tool has TC=0
    no_zero = [f for f in all_families
                if kalign[f]["tc"] > 0 and mafft[f]["tc"] > 0
                and muscle[f]["tc"] > 0 and clustalo[f]["tc"] > 0]
    print(f"Families where ALL tools have TC>0: {len(no_zero)} / {len(all_families)}")

    if no_zero:
        print()
        print(f"{'Tool':<12} {'Mean TC':>10} {'Mean SP':>10} {'TC/SP':>10}")
        print("-" * 48)
        for name, fam_data in tools.items():
            avg_tc = sum(fam_data[f]["tc"] for f in no_zero) / len(no_zero)
            avg_sp = sum(fam_data[f]["recall"] for f in no_zero) / len(no_zero)
            ratio = avg_tc / avg_sp if avg_sp > 0 else 0
            print(f"{name:<12} {avg_tc:>10.3f} {avg_sp:>10.3f} {ratio:>10.3f}")

    # Decompose the overall TC gap into contributions
    print()
    print("Decomposition of kalign vs MAFFT TC gap:")
    overall_k_tc = sum(kalign[f]["tc"] for f in all_families) / len(all_families)
    overall_m_tc = sum(mafft[f]["tc"] for f in all_families) / len(all_families)
    overall_gap = overall_m_tc - overall_k_tc
    print(f"  Overall TC gap (MAFFT - kalign): {overall_gap:.3f}")

    # Contribution from families where kalign TC=0 but MAFFT TC>0
    k_tc0_m_gt0 = [f for f in all_families if kalign[f]["tc"] == 0 and mafft[f]["tc"] > 0]
    # Contribution from families where MAFFT TC=0 but kalign TC>0
    m_tc0_k_gt0 = [f for f in all_families if mafft[f]["tc"] == 0 and kalign[f]["tc"] > 0]
    # Both TC=0
    both_tc0 = [f for f in all_families if kalign[f]["tc"] == 0 and mafft[f]["tc"] == 0]
    # Both TC>0
    both_gt0_all = [f for f in all_families if kalign[f]["tc"] > 0 and mafft[f]["tc"] > 0]

    gap_from_k0_mgt0 = sum(mafft[f]["tc"] for f in k_tc0_m_gt0) / len(all_families)
    gap_from_m0_kgt0 = -sum(kalign[f]["tc"] for f in m_tc0_k_gt0) / len(all_families)
    gap_from_both_gt0 = sum(mafft[f]["tc"] - kalign[f]["tc"] for f in both_gt0_all) / len(all_families)
    # both_tc0 contributes 0

    print(f"  Families: kalign TC=0 & MAFFT TC>0: {len(k_tc0_m_gt0):>4} -> TC gap contribution: {gap_from_k0_mgt0:>+.3f}")
    print(f"  Families: MAFFT TC=0 & kalign TC>0: {len(m_tc0_k_gt0):>4} -> TC gap contribution: {gap_from_m0_kgt0:>+.3f}")
    print(f"  Families: both TC=0:                {len(both_tc0):>4} -> TC gap contribution: +0.000")
    print(f"  Families: both TC>0:                {len(both_gt0_all):>4} -> TC gap contribution: {gap_from_both_gt0:>+.3f}")
    print(f"  Sum of contributions:                      {gap_from_k0_mgt0 + gap_from_m0_kgt0 + gap_from_both_gt0:>+.3f}")
    print(f"  Actual gap:                                {overall_gap:>+.3f}")

    pct_from_killed = gap_from_k0_mgt0 / overall_gap * 100 if overall_gap > 0 else 0
    pct_from_broad = gap_from_both_gt0 / overall_gap * 100 if overall_gap > 0 else 0
    pct_from_m0 = gap_from_m0_kgt0 / overall_gap * 100 if overall_gap > 0 else 0
    print()
    print(f"  -> {pct_from_killed:.1f}% of TC gap from TC-killed families (kalign=0, MAFFT>0)")
    print(f"  -> {pct_from_broad:.1f}% of TC gap from broad deficit (both>0)")
    print(f"  -> {pct_from_m0:.1f}% offset from families where kalign>0, MAFFT=0")

    # =========================================================================
    # 7. TC/SP ratio comparison per family
    # =========================================================================
    print_section("7. TC/SP RATIO COMPARISON (kalign vs MAFFT)")

    # For families where both have TC>0 and SP>0
    tc_sp_ratios_k = []
    tc_sp_ratios_m = []
    ratio_diffs = []
    for f in both_gt0_all:
        k_sp = kalign[f]["recall"]
        m_sp = mafft[f]["recall"]
        if k_sp > 0 and m_sp > 0:
            kr = kalign[f]["tc"] / k_sp
            mr = mafft[f]["tc"] / m_sp
            tc_sp_ratios_k.append(kr)
            tc_sp_ratios_m.append(mr)
            ratio_diffs.append(mr - kr)

    if tc_sp_ratios_k:
        n = len(tc_sp_ratios_k)
        print(f"Families with both TC>0 and SP>0: {n}")
        print(f"  kalign mean TC/SP: {sum(tc_sp_ratios_k)/n:.3f}")
        print(f"  MAFFT  mean TC/SP: {sum(tc_sp_ratios_m)/n:.3f}")
        print(f"  Mean (MAFFT TC/SP - kalign TC/SP): {sum(ratio_diffs)/n:+.3f}")

        # By category
        print()
        print(f"{'Category':<10} {'N':>4}  {'kalign TC/SP':>14} {'MAFFT TC/SP':>14} {'Diff':>10}")
        print("-" * 60)
        cat_data = defaultdict(lambda: {"k": [], "m": []})
        for f in both_gt0_all:
            k_sp = kalign[f]["recall"]
            m_sp = mafft[f]["recall"]
            if k_sp > 0 and m_sp > 0:
                cat = category_from_dataset(kalign[f]["dataset"])
                cat_data[cat]["k"].append(kalign[f]["tc"] / k_sp)
                cat_data[cat]["m"].append(mafft[f]["tc"] / m_sp)
        for cat in sorted(cat_data.keys()):
            kvals = cat_data[cat]["k"]
            mvals = cat_data[cat]["m"]
            n = len(kvals)
            kavg = sum(kvals) / n
            mavg = sum(mvals) / n
            print(f"{cat:<10} {n:>4}  {kavg:>14.3f} {mavg:>14.3f} {mavg-kavg:>+10.3f}")

    # =========================================================================
    # Additional: Analyze the log-odds view
    # =========================================================================
    print_section("8. LOG-SPACE ANALYSIS: TC vs SP^(nseq_choose_2)")

    # If pair errors were independent across columns, TC ~ SP^(npairs) is too
    # extreme. Instead, let's look at log(TC)/log(SP) as an effective "pair count"
    # -- how many independent pair errors does it take to explain the observed TC?
    print("Effective pair-error multiplier: log(TC) / log(SP)")
    print("  If = 1: TC = SP (errors perfectly correlated across columns)")
    print("  If = npairs: TC = SP^npairs (errors fully independent)")
    print()

    eff_mult_k = []
    eff_mult_m = []
    eff_mult_data = []
    for f in both_gt0_all:
        k_sp = kalign[f]["recall"]
        k_tc = kalign[f]["tc"]
        m_sp = mafft[f]["recall"]
        m_tc = mafft[f]["tc"]
        nseq = overgap.get(f, {}).get("nseq", None)
        if nseq is None:
            continue
        npairs = nseq * (nseq - 1) / 2

        if 0 < k_sp < 1 and 0 < k_tc < 1:
            em_k = math.log(k_tc) / math.log(k_sp)
            eff_mult_k.append(em_k)
        else:
            em_k = None

        if 0 < m_sp < 1 and 0 < m_tc < 1:
            em_m = math.log(m_tc) / math.log(m_sp)
            eff_mult_m.append(em_m)
        else:
            em_m = None

        if em_k is not None and em_m is not None:
            eff_mult_data.append({
                "family": f,
                "category": category_from_dataset(kalign[f]["dataset"]),
                "nseq": nseq,
                "npairs": npairs,
                "em_kalign": em_k,
                "em_mafft": em_m,
            })

    if eff_mult_k:
        print(f"kalign effective multiplier: mean={sum(eff_mult_k)/len(eff_mult_k):.2f}, "
              f"median={sorted(eff_mult_k)[len(eff_mult_k)//2]:.2f} (n={len(eff_mult_k)})")
    if eff_mult_m:
        print(f"MAFFT  effective multiplier: mean={sum(eff_mult_m)/len(eff_mult_m):.2f}, "
              f"median={sorted(eff_mult_m)[len(eff_mult_m)//2]:.2f} (n={len(eff_mult_m)})")

    print()
    print("Interpretation: higher multiplier = errors more spread across columns")
    print("                (more columns have at least one error)")

    # By category
    if eff_mult_data:
        print()
        cat_em = defaultdict(lambda: {"k": [], "m": [], "npairs": []})
        for d in eff_mult_data:
            cat_em[d["category"]]["k"].append(d["em_kalign"])
            cat_em[d["category"]]["m"].append(d["em_mafft"])
            cat_em[d["category"]]["npairs"].append(d["npairs"])

        print(f"{'Category':<10} {'N':>4}  {'kalign EM':>12} {'MAFFT EM':>12} {'Avg npairs':>12}")
        print("-" * 58)
        for cat in sorted(cat_em.keys()):
            n = len(cat_em[cat]["k"])
            kavg = sum(cat_em[cat]["k"]) / n
            mavg = sum(cat_em[cat]["m"]) / n
            pavg = sum(cat_em[cat]["npairs"]) / n
            print(f"{cat:<10} {n:>4}  {kavg:>12.2f} {mavg:>12.2f} {pavg:>12.1f}")

    # =========================================================================
    # 9. Summary of key findings
    # =========================================================================
    print_section("9. KEY FINDINGS SUMMARY")

    n_total = len(all_families)
    k_tc0 = sum(1 for f in all_families if kalign[f]["tc"] == 0)
    m_tc0 = sum(1 for f in all_families if mafft[f]["tc"] == 0)
    k_only_tc0 = len(k_tc0_m_gt0)

    print(f"1. SCOPE: kalign has {k_tc0} families with TC=0 vs MAFFT's {m_tc0}")
    print(f"   - {k_only_tc0} families have kalign TC=0 but MAFFT TC>0 (TC-killed)")
    print()
    print(f"2. DECOMPOSITION of total TC gap ({overall_gap:.3f}):")
    print(f"   - TC-killed families contribute {gap_from_k0_mgt0:.3f} ({pct_from_killed:.1f}%)")
    print(f"   - Broad deficit across remaining families contributes {gap_from_both_gt0:.3f} ({pct_from_broad:.1f}%)")
    print()

    if no_zero:
        nz_k_tc = sum(kalign[f]["tc"] for f in no_zero) / len(no_zero)
        nz_m_tc = sum(mafft[f]["tc"] for f in no_zero) / len(no_zero)
        nz_mu_tc = sum(muscle[f]["tc"] for f in no_zero) / len(no_zero)
        nz_cl_tc = sum(clustalo[f]["tc"] for f in no_zero) / len(no_zero)
        print(f"3. EXCLUDING ALL TC=0 FAMILIES ({len(no_zero)} remain):")
        print(f"   kalign TC={nz_k_tc:.3f}  MAFFT TC={nz_m_tc:.3f}  MUSCLE TC={nz_mu_tc:.3f}  ClustalO TC={nz_cl_tc:.3f}")
        print(f"   Gap shrinks from {overall_gap:.3f} to {nz_m_tc - nz_k_tc:.3f}")
        remaining_gap = nz_m_tc - nz_k_tc
    else:
        remaining_gap = 0

    print()
    print(f"4. ANSWER: Is kalign's TC problem mainly:")
    print(f"   (a) A few families with TC=0 dragging down the average?  -> {pct_from_killed:.0f}%")
    print(f"   (b) A broad deficit across many families?                -> {pct_from_broad:.0f}%")
    if pct_from_killed > pct_from_broad:
        print(f"   VERDICT: Primarily driven by TC-killed families ({pct_from_killed:.0f}%)")
    elif pct_from_broad > pct_from_killed:
        print(f"   VERDICT: Primarily a broad deficit ({pct_from_broad:.0f}%)")
    else:
        print(f"   VERDICT: Roughly equal contribution from both")

    print()
    if eff_mult_k and eff_mult_m:
        em_k_mean = sum(eff_mult_k) / len(eff_mult_k)
        em_m_mean = sum(eff_mult_m) / len(eff_mult_m)
        print(f"5. ERROR DISTRIBUTION: kalign's errors are more spread across columns")
        print(f"   (effective multiplier: kalign={em_k_mean:.2f} vs MAFFT={em_m_mean:.2f})")
        print(f"   This means kalign's alignment errors contaminate more columns,")
        print(f"   killing TC even when SP (pairwise accuracy) is reasonable.")

    # =========================================================================
    # 10. Per-category detailed breakdown
    # =========================================================================
    print_section("10. PER-CATEGORY DETAILED BREAKDOWN")

    print(f"{'Cat':<8} {'N':>4} | {'k_TC=0':>7} {'m_TC=0':>7} {'k_only0':>8} | "
          f"{'k_TC':>7} {'m_TC':>7} {'gap':>7} | {'k_SP':>7} {'m_SP':>7}")
    print("-" * 90)
    for cat in sorted(categories.keys()):
        fams = categories[cat]
        n = len(fams)
        k0 = sum(1 for f in fams if kalign[f]["tc"] == 0)
        m0 = sum(1 for f in fams if mafft[f]["tc"] == 0)
        ko = sum(1 for f in fams if kalign[f]["tc"] == 0 and mafft[f]["tc"] > 0)
        ktc = sum(kalign[f]["tc"] for f in fams) / n
        mtc = sum(mafft[f]["tc"] for f in fams) / n
        ksp = sum(kalign[f]["recall"] for f in fams) / n
        msp = sum(mafft[f]["recall"] for f in fams) / n
        gap = mtc - ktc
        print(f"{cat:<8} {n:>4} | {k0:>7} {m0:>7} {ko:>8} | "
              f"{ktc:>7.3f} {mtc:>7.3f} {gap:>+7.3f} | {ksp:>7.3f} {msp:>7.3f}")

    # =========================================================================
    # 11. Nseq sensitivity analysis
    # =========================================================================
    print_section("11. NSEQ SENSITIVITY: TC DROPS FASTER WITH MORE SEQUENCES?")

    nseq_bins = [(2, 5), (6, 10), (11, 20), (21, 50), (51, 1000)]
    print(f"{'nseq range':<12} {'N':>4} | {'k_TC':>7} {'m_TC':>7} {'gap':>7} | "
          f"{'k_TC=0':>7} {'m_TC=0':>7} | {'k_TC/SP':>8} {'m_TC/SP':>8}")
    print("-" * 90)
    for lo, hi in nseq_bins:
        fams_in = [f for f in all_families
                   if f in overgap and lo <= overgap[f].get("nseq", 0) <= hi]
        if not fams_in:
            continue
        n = len(fams_in)
        ktc = sum(kalign[f]["tc"] for f in fams_in) / n
        mtc = sum(mafft[f]["tc"] for f in fams_in) / n
        ksp = sum(kalign[f]["recall"] for f in fams_in) / n
        msp = sum(mafft[f]["recall"] for f in fams_in) / n
        k0 = sum(1 for f in fams_in if kalign[f]["tc"] == 0)
        m0 = sum(1 for f in fams_in if mafft[f]["tc"] == 0)
        ktcsp = ktc / ksp if ksp > 0 else 0
        mtcsp = mtc / msp if msp > 0 else 0
        gap = mtc - ktc
        print(f"{lo:>3}-{hi:<7} {n:>4} | {ktc:>7.3f} {mtc:>7.3f} {gap:>+7.3f} | "
              f"{k0:>7} {m0:>7} | {ktcsp:>8.3f} {mtcsp:>8.3f}")


if __name__ == "__main__":
    main()
