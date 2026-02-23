"""Compare ProbMSA vs mafft vs true alignments at depth 4.

Runs inside the downstream container. Produces a diagnostic report
showing WHERE ProbMSA fails compared to mafft.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def parse_fasta(path):
    names, seqs = [], []
    name, seq = None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                names.append(name)
                seqs.append("".join(seq))
            name = line[1:].split()[0]
            seq = []
        elif line:
            seq.append(line)
    if name is not None:
        names.append(name)
        seqs.append("".join(seq))
    return names, seqs


def alignment_to_pairs(names, seqs):
    """Extract all aligned residue pairs (i.e. same column, both non-gap)."""
    pairs = set()
    n = len(seqs)
    alen = len(seqs[0])
    # Build residue index arrays
    res_idx = []
    for s in range(n):
        ri = []
        c = 0
        for col in range(alen):
            if seqs[s][col] != "-":
                ri.append(c)
                c += 1
            else:
                ri.append(-1)
        res_idx.append(ri)

    for col in range(alen):
        for i in range(n):
            if res_idx[i][col] < 0:
                continue
            for j in range(i + 1, n):
                if res_idx[j][col] < 0:
                    continue
                pairs.add((i, j, res_idx[i][col], res_idx[j][col]))
    return pairs


def gap_stats(names, seqs):
    """Compute gap statistics for an alignment."""
    n = len(seqs)
    alen = len(seqs[0]) if seqs else 0

    total_gaps = 0
    gap_lengths = []
    gaps_per_seq = []

    for s in range(n):
        seq = seqs[s]
        ngaps = seq.count("-")
        total_gaps += ngaps
        gaps_per_seq.append(ngaps)

        # Gap length distribution
        in_gap = False
        glen = 0
        for c in seq:
            if c == "-":
                glen += 1
                in_gap = True
            else:
                if in_gap:
                    gap_lengths.append(glen)
                    glen = 0
                    in_gap = False
        if in_gap:
            gap_lengths.append(glen)

    ungapped_len = sum(len(s.replace("-", "")) for s in seqs)

    return {
        "alen": alen,
        "total_gaps": total_gaps,
        "gap_frac": total_gaps / (n * alen) if alen > 0 else 0,
        "n_gap_blocks": len(gap_lengths),
        "mean_gap_len": np.mean(gap_lengths) if gap_lengths else 0,
        "median_gap_len": np.median(gap_lengths) if gap_lengths else 0,
        "max_gap_len": max(gap_lengths) if gap_lengths else 0,
        "mean_gaps_per_seq": np.mean(gaps_per_seq),
    }


def error_analysis(true_names, true_seqs, test_names, test_seqs):
    """Classify errors: false positive pairs (FP) and false negative pairs (FN)."""
    # Reorder test to match true name order
    test_map = {n: s for n, s in zip(test_names, test_seqs)}
    test_seqs_ordered = [test_map[n] for n in true_names]

    true_pairs = alignment_to_pairs(true_names, true_seqs)
    test_pairs = alignment_to_pairs(true_names, test_seqs_ordered)

    tp = true_pairs & test_pairs
    fp = test_pairs - true_pairs  # pairs in test but not in true
    fn = true_pairs - test_pairs  # pairs in true but not in test

    return {
        "tp": len(tp),
        "fp": len(fp),
        "fn": len(fn),
        "precision": len(tp) / len(test_pairs) if test_pairs else 1.0,
        "recall": len(tp) / len(true_pairs) if true_pairs else 1.0,
    }


def run_diagnostics():
    from benchmarks.downstream.simulation import (
        PROTEIN_GRID,
        generate_indelible_dataset,
        iter_simulation_params,
        random_birth_death_tree,
    )
    from benchmarks.downstream.utils import run_method

    data_dir = Path("benchmarks/data")
    sim_dir = data_dir / "diag_sims"
    sim_dir.mkdir(parents=True, exist_ok=True)

    # Generate depth-4 simulations
    all_params = list(iter_simulation_params(PROTEIN_GRID, "WAG"))
    depth4_params = [p for p in all_params if p["target_depth"] == 4.0]
    print(f"Generating {len(depth4_params)} depth-4 simulations...")

    sims = []
    for sp in depth4_params:
        sim_out = sim_dir / sp["sim_id"]
        tree = random_birth_death_tree(
            n_taxa=sp["n_taxa"],
            target_depth=sp["target_depth"],
            seed=sp.get("seed", 42),
        )
        ds = generate_indelible_dataset(
            tree=tree,
            model=sp["model"],
            seq_length=sp.get("seq_length", 300),
            indel_rate=sp["indel_rate"],
            indel_length_mean=sp.get("indel_length_mean", 2.0),
            output_dir=sim_out,
            seed=sp.get("seed", 42),
        )
        ds.params.update(sp)
        sims.append(ds)
    print(f"Generated {len(sims)} simulations.")

    # Collect per-case diagnostics
    results = []

    methods = ["kalign_probmsa", "mafft"]

    for si, sim in enumerate(sims):
        params = sim.params
        case_id = params.get("sim_id", str(si))
        work_dir = sim_dir / case_id / "diag"
        work_dir.mkdir(parents=True, exist_ok=True)

        true_names, true_seqs = parse_fasta(sim.true_alignment)
        true_gs = gap_stats(true_names, true_seqs)

        case = {
            "case_id": case_id,
            "n_taxa": params.get("n_taxa"),
            "indel_rate": params.get("indel_rate"),
            "indel_length_mean": params.get("indel_length_mean"),
            "replicate": params.get("replicate"),
            "true_alen": true_gs["alen"],
            "true_gap_frac": true_gs["gap_frac"],
            "true_n_gap_blocks": true_gs["n_gap_blocks"],
            "true_mean_gap_len": true_gs["mean_gap_len"],
        }

        for method in methods:
            try:
                aln = run_method(method, sim.unaligned, work_dir, seq_type="protein")
                gs = gap_stats(aln.names, aln.sequences)
                ea = error_analysis(true_names, true_seqs, aln.names, aln.sequences)

                case[f"{method}_alen"] = gs["alen"]
                case[f"{method}_gap_frac"] = gs["gap_frac"]
                case[f"{method}_n_gap_blocks"] = gs["n_gap_blocks"]
                case[f"{method}_mean_gap_len"] = gs["mean_gap_len"]
                case[f"{method}_max_gap_len"] = gs["max_gap_len"]
                case[f"{method}_precision"] = ea["precision"]
                case[f"{method}_recall"] = ea["recall"]
                case[f"{method}_tp"] = ea["tp"]
                case[f"{method}_fp"] = ea["fp"]
                case[f"{method}_fn"] = ea["fn"]
            except Exception as e:
                print(f"  FAILED {method} on {case_id}: {e}", file=sys.stderr)

        results.append(case)
        if (si + 1) % 10 == 0:
            print(f"  {si + 1}/{len(sims)} cases done")

    # ── Summary report ───────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("ALIGNMENT DIAGNOSTICS: ProbMSA vs MAFFT at depth 4.0")
    print("=" * 80)

    # Overall averages
    print("\n--- Overall Averages ---")
    print(f"{'Metric':<35} {'True':>10} {'ProbMSA':>10} {'MAFFT':>10}")
    print("-" * 70)

    avg = lambda key: np.mean([r[key] for r in results if key in r])

    print(f"{'Alignment length':<35} {avg('true_alen'):>10.1f} {avg('kalign_probmsa_alen'):>10.1f} {avg('mafft_alen'):>10.1f}")
    print(f"{'Gap fraction':<35} {avg('true_gap_frac'):>10.3f} {avg('kalign_probmsa_gap_frac'):>10.3f} {avg('mafft_gap_frac'):>10.3f}")
    print(f"{'# Gap blocks':<35} {avg('true_n_gap_blocks'):>10.1f} {avg('kalign_probmsa_n_gap_blocks'):>10.1f} {avg('mafft_n_gap_blocks'):>10.1f}")
    print(f"{'Mean gap length':<35} {avg('true_mean_gap_len'):>10.2f} {avg('kalign_probmsa_mean_gap_len'):>10.2f} {avg('mafft_mean_gap_len'):>10.2f}")
    print(f"{'Max gap length':<35} {'':>10} {avg('kalign_probmsa_max_gap_len'):>10.1f} {avg('mafft_max_gap_len'):>10.1f}")
    print(f"{'Precision':<35} {'':>10} {avg('kalign_probmsa_precision'):>10.3f} {avg('mafft_precision'):>10.3f}")
    print(f"{'Recall':<35} {'':>10} {avg('kalign_probmsa_recall'):>10.3f} {avg('mafft_recall'):>10.3f}")

    # Stratified by indel_rate
    print("\n--- Stratified by indel_rate ---")
    for ir in sorted(set(r["indel_rate"] for r in results)):
        subset = [r for r in results if r["indel_rate"] == ir]
        n = len(subset)
        pm_prec = np.mean([r["kalign_probmsa_precision"] for r in subset])
        pm_rec = np.mean([r["kalign_probmsa_recall"] for r in subset])
        mf_prec = np.mean([r["mafft_precision"] for r in subset])
        mf_rec = np.mean([r["mafft_recall"] for r in subset])
        pm_gf = np.mean([r["kalign_probmsa_gap_frac"] for r in subset])
        mf_gf = np.mean([r["mafft_gap_frac"] for r in subset])
        true_gf = np.mean([r["true_gap_frac"] for r in subset])
        pm_mgl = np.mean([r["kalign_probmsa_mean_gap_len"] for r in subset])
        mf_mgl = np.mean([r["mafft_mean_gap_len"] for r in subset])
        true_mgl = np.mean([r["true_mean_gap_len"] for r in subset])

        print(f"\n  indel_rate={ir} (n={n})")
        print(f"    {'':30} {'True':>8} {'ProbMSA':>8} {'MAFFT':>8}")
        print(f"    {'Gap fraction':<30} {true_gf:>8.3f} {pm_gf:>8.3f} {mf_gf:>8.3f}")
        print(f"    {'Mean gap length':<30} {true_mgl:>8.2f} {pm_mgl:>8.2f} {mf_mgl:>8.2f}")
        print(f"    {'Precision':<30} {'':>8} {pm_prec:>8.3f} {mf_prec:>8.3f}")
        print(f"    {'Recall':<30} {'':>8} {pm_rec:>8.3f} {mf_rec:>8.3f}")

    # Stratified by indel_length_mean
    print("\n--- Stratified by indel_length_mean ---")
    for ilm in sorted(set(r["indel_length_mean"] for r in results)):
        subset = [r for r in results if r["indel_length_mean"] == ilm]
        n = len(subset)
        pm_prec = np.mean([r["kalign_probmsa_precision"] for r in subset])
        pm_rec = np.mean([r["kalign_probmsa_recall"] for r in subset])
        mf_prec = np.mean([r["mafft_precision"] for r in subset])
        mf_rec = np.mean([r["mafft_recall"] for r in subset])
        pm_gf = np.mean([r["kalign_probmsa_gap_frac"] for r in subset])
        mf_gf = np.mean([r["mafft_gap_frac"] for r in subset])
        true_gf = np.mean([r["true_gap_frac"] for r in subset])

        print(f"\n  indel_length_mean={ilm} (n={n})")
        print(f"    {'':30} {'True':>8} {'ProbMSA':>8} {'MAFFT':>8}")
        print(f"    {'Gap fraction':<30} {true_gf:>8.3f} {pm_gf:>8.3f} {mf_gf:>8.3f}")
        print(f"    {'Precision':<30} {'':>8} {pm_prec:>8.3f} {mf_prec:>8.3f}")
        print(f"    {'Recall':<30} {'':>8} {pm_rec:>8.3f} {mf_rec:>8.3f}")

    # Stratified by n_taxa
    print("\n--- Stratified by n_taxa ---")
    for nt in sorted(set(r["n_taxa"] for r in results)):
        subset = [r for r in results if r["n_taxa"] == nt]
        n = len(subset)
        pm_prec = np.mean([r["kalign_probmsa_precision"] for r in subset])
        pm_rec = np.mean([r["kalign_probmsa_recall"] for r in subset])
        mf_prec = np.mean([r["mafft_precision"] for r in subset])
        mf_rec = np.mean([r["mafft_recall"] for r in subset])

        print(f"\n  n_taxa={nt} (n={n})")
        print(f"    {'Precision':<30} {'':>8} {pm_prec:>8.3f} {mf_prec:>8.3f}")
        print(f"    {'Recall':<30} {'':>8} {pm_rec:>8.3f} {mf_rec:>8.3f}")

    # Gap fraction error (over/under-gapping)
    print("\n--- Gap Fraction Error (test - true) ---")
    pm_gap_err = [r["kalign_probmsa_gap_frac"] - r["true_gap_frac"] for r in results]
    mf_gap_err = [r["mafft_gap_frac"] - r["true_gap_frac"] for r in results]
    print(f"  ProbMSA: mean={np.mean(pm_gap_err):+.4f}  std={np.std(pm_gap_err):.4f}  (>0 = over-gapping)")
    print(f"  MAFFT:   mean={np.mean(mf_gap_err):+.4f}  std={np.std(mf_gap_err):.4f}")

    # Alignment length error
    print("\n--- Alignment Length Error (test - true) ---")
    pm_alen_err = [r["kalign_probmsa_alen"] - r["true_alen"] for r in results if "kalign_probmsa_alen" in r]
    mf_alen_err = [r["mafft_alen"] - r["true_alen"] for r in results if "mafft_alen" in r]
    print(f"  ProbMSA: mean={np.mean(pm_alen_err):+.1f}  std={np.std(pm_alen_err):.1f}")
    print(f"  MAFFT:   mean={np.mean(mf_alen_err):+.1f}  std={np.std(mf_alen_err):.1f}")

    # FP vs FN analysis
    print("\n--- Error Type Breakdown (avg per case) ---")
    pm_fp = np.mean([r["kalign_probmsa_fp"] for r in results])
    pm_fn = np.mean([r["kalign_probmsa_fn"] for r in results])
    mf_fp = np.mean([r["mafft_fp"] for r in results])
    mf_fn = np.mean([r["mafft_fn"] for r in results])
    pm_tp = np.mean([r["kalign_probmsa_tp"] for r in results])
    mf_tp = np.mean([r["mafft_tp"] for r in results])

    print(f"  {'':20} {'ProbMSA':>10} {'MAFFT':>10}")
    print(f"  {'True positives':<20} {pm_tp:>10.0f} {mf_tp:>10.0f}")
    print(f"  {'False positives':<20} {pm_fp:>10.0f} {mf_fp:>10.0f}")
    print(f"  {'False negatives':<20} {pm_fn:>10.0f} {mf_fn:>10.0f}")
    print(f"  {'FP/(FP+FN) ratio':<20} {pm_fp/(pm_fp+pm_fn):>10.3f} {mf_fp/(mf_fp+mf_fn):>10.3f}")

    # Correlation: where ProbMSA improves over mafft
    print("\n--- Cases where ProbMSA beats MAFFT (by precision) ---")
    pm_wins = [r for r in results if r.get("kalign_probmsa_precision", 0) > r.get("mafft_precision", 0)]
    print(f"  ProbMSA wins: {len(pm_wins)}/{len(results)} cases")
    if pm_wins:
        print(f"  Avg indel_rate in ProbMSA wins: {np.mean([r['indel_rate'] for r in pm_wins]):.3f}")
        print(f"  Avg indel_rate in MAFFT wins:   {np.mean([r['indel_rate'] for r in results if r not in pm_wins]):.3f}")

    # Save raw results
    out_path = Path("benchmarks/results") / "alignment_diagnostics.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nRaw results saved to {out_path}")


if __name__ == "__main__":
    import os
    os.chdir("/kalign")
    sys.path.insert(0, "/kalign")
    run_diagnostics()
