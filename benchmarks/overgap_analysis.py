"""Over-gapping analysis for kalign on BAliBASE.

Measures alignment length ratios (kalign output / reference) to investigate
whether kalign produces longer-than-necessary alignments, and correlates
this with TC score deficits.
"""

import json
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kalign

from benchmarks.datasets import get_cases, BenchmarkCase

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_fasta_alignment_length(fasta_path: Path) -> int:
    """Read a FASTA alignment and return the alignment length (columns)."""
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # Accumulate the first sequence
                seq = []
                for line2 in f:
                    if line2.startswith(">"):
                        break
                    seq.append(line2.strip())
                return len("".join(seq))
    raise ValueError(f"No sequences found in {fasta_path}")


def read_msf_alignment_length(msf_path: Path) -> int:
    """Read an MSF alignment and return the alignment length from the header."""
    with open(msf_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("MSF:"):
                # "   MSF:   96  Type: P ..."
                parts = line.split()
                idx = parts.index("MSF:")
                return int(parts[idx + 1])
            if "MSF:" in line:
                # Handle "   MSF:   96  Type: P    Check: ..."
                for i, tok in enumerate(line.split()):
                    if tok == "MSF:":
                        return int(line.split()[i + 1])
    raise ValueError(f"Could not parse MSF length from {msf_path}")


def count_sequences_and_mean_length(fasta_path: Path) -> tuple:
    """Read unaligned FASTA, return (n_sequences, mean_sequence_length)."""
    seqs = []
    current = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if current:
                    seqs.append(len("".join(current)))
                current = []
            else:
                current.append(line.strip())
    if current:
        seqs.append(len("".join(current)))
    if not seqs:
        return 0, 0.0
    return len(seqs), sum(seqs) / len(seqs)


def get_ref_alignment_length(case: BenchmarkCase) -> int:
    """Get the reference alignment length."""
    ref = case.reference
    if ref.suffix == ".msf":
        return read_msf_alignment_length(ref)
    else:
        # FASTA format reference
        return read_fasta_alignment_length(ref)


# ---------------------------------------------------------------------------
# Worker function for parallel execution
# ---------------------------------------------------------------------------

def process_case(case_dict: dict) -> dict:
    """Process a single BAliBASE case. Takes a dict (for pickling) and returns results."""
    family = case_dict["family"]
    dataset = case_dict["dataset"]
    unaligned = Path(case_dict["unaligned"])
    reference = Path(case_dict["reference"])
    seq_type = case_dict["seq_type"]

    case = BenchmarkCase(
        family=family,
        dataset=dataset,
        unaligned=unaligned,
        reference=reference,
        seq_type=seq_type,
    )

    result = {
        "family": family,
        "dataset": dataset,
        "category": dataset.replace("balibase_", ""),
    }

    try:
        # Reference alignment length
        ref_len = get_ref_alignment_length(case)
        result["ref_len"] = ref_len

        # Sequence stats
        nseq, mean_seq_len = count_sequences_and_mean_length(unaligned)
        result["nseq"] = nseq
        result["mean_seq_len"] = mean_seq_len

        with tempfile.TemporaryDirectory() as tmpdir:
            # Default mode: vsm_amax=2.0, refine=confident
            out_default = Path(tmpdir) / "default.fasta"
            kalign.align_file_to_file(
                str(unaligned),
                str(out_default),
                format="fasta",
                seq_type=seq_type,
                vsm_amax=2.0,
                refine="confident",
                n_threads=1,
            )
            default_len = read_fasta_alignment_length(out_default)
            result["default_len"] = default_len
            result["default_ratio"] = default_len / ref_len if ref_len > 0 else float("nan")

            # Best mode: ensemble=3, vsm_amax=2.0, refine=confident, realign=1
            out_best = Path(tmpdir) / "best.fasta"
            kalign.align_file_to_file(
                str(unaligned),
                str(out_best),
                format="fasta",
                seq_type=seq_type,
                vsm_amax=2.0,
                refine="confident",
                ensemble=3,
                realign=1,
                n_threads=1,
            )
            best_len = read_fasta_alignment_length(out_best)
            result["best_len"] = best_len
            result["best_ratio"] = best_len / ref_len if ref_len > 0 else float("nan")

        result["error"] = None

    except Exception as e:
        result["error"] = str(e)

    return result


# ---------------------------------------------------------------------------
# Load existing TC scores for correlation
# ---------------------------------------------------------------------------

def load_tc_scores(results_dir: Path) -> dict:
    """Load TC scores from existing benchmark result files.

    Returns dict of {(family, method): tc_score}
    """
    tc_scores = {}

    # balibase_baseline.json has default+refine=confident scores
    baseline_path = results_dir / "balibase_baseline.json"
    if baseline_path.exists():
        with open(baseline_path) as f:
            data = json.load(f)
        for r in data.get("results", []):
            if r.get("refine") == "confident" and r.get("ensemble", 0) == 0:
                tc_scores[(r["family"], "default")] = r.get("tc", 0.0)
            elif r.get("refine") == "none" and r.get("ensemble", 0) == 0:
                tc_scores[(r["family"], "baseline_norefine")] = r.get("tc", 0.0)

    # balibase_best.json has best mode scores
    best_path = results_dir / "balibase_best.json"
    if best_path.exists():
        with open(best_path) as f:
            data = json.load(f)
        if isinstance(data, list):
            results = data
        else:
            results = data.get("results", [])
        for r in results:
            tc_scores[(r["family"], "best")] = r.get("tc", 0.0)

    # full_comparison.json has external tool scores
    comp_path = results_dir / "full_comparison.json"
    if comp_path.exists():
        with open(comp_path) as f:
            data = json.load(f)
        for r in data.get("results", []):
            method = r.get("method", "")
            if method in ("clustalo", "mafft", "muscle"):
                tc_scores[(r["family"], method)] = r.get("tc", 0.0)

    return tc_scores


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("BAliBASE Over-Gapping Analysis")
    print("=" * 72)
    print()

    # Get cases
    cases = get_cases("balibase")
    print(f"Found {len(cases)} BAliBASE cases")
    print()

    # Load TC scores
    results_dir = Path(__file__).parent / "results"
    tc_scores = load_tc_scores(results_dir)
    print(f"Loaded TC scores for {len(tc_scores)} (family, method) pairs")
    print()

    # Convert cases to dicts for pickling
    case_dicts = [
        {
            "family": c.family,
            "dataset": c.dataset,
            "unaligned": str(c.unaligned),
            "reference": str(c.reference),
            "seq_type": c.seq_type,
        }
        for c in cases
    ]

    # Process in parallel
    results = []
    n_workers = 4
    print(f"Processing {len(case_dicts)} cases with {n_workers} workers...")
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_case, cd): cd["family"] for cd in case_dicts}
        done_count = 0
        for future in as_completed(futures):
            done_count += 1
            family = futures[future]
            try:
                result = future.result()
                results.append(result)
                if done_count % 20 == 0 or done_count == len(case_dicts):
                    print(f"  [{done_count}/{len(case_dicts)}] done")
            except Exception as e:
                print(f"  ERROR on {family}: {e}")

    # Filter out errors
    errors = [r for r in results if r.get("error")]
    ok = [r for r in results if not r.get("error")]
    if errors:
        print(f"\n{len(errors)} cases had errors:")
        for e in errors[:5]:
            print(f"  {e['family']}: {e['error']}")
        if len(errors) > 5:
            print(f"  ... and {len(errors) - 5} more")

    print(f"\nSuccessfully processed {len(ok)} cases")
    print()

    # ---------------------------------------------------------------------------
    # Summary by category
    # ---------------------------------------------------------------------------
    print("=" * 72)
    print("ALIGNMENT LENGTH RATIOS BY CATEGORY (kalign_len / ref_len)")
    print("=" * 72)
    print()

    categories = sorted(set(r["category"] for r in ok))

    # Build per-category stats
    cat_stats = defaultdict(lambda: {
        "default_ratios": [], "best_ratios": [],
        "nseqs": [], "mean_seq_lens": [],
        "ref_lens": [], "default_lens": [], "best_lens": [],
    })

    for r in ok:
        cat = r["category"]
        cat_stats[cat]["default_ratios"].append(r["default_ratio"])
        cat_stats[cat]["best_ratios"].append(r["best_ratio"])
        cat_stats[cat]["nseqs"].append(r["nseq"])
        cat_stats[cat]["mean_seq_lens"].append(r["mean_seq_len"])
        cat_stats[cat]["ref_lens"].append(r["ref_len"])
        cat_stats[cat]["default_lens"].append(r["default_len"])
        cat_stats[cat]["best_lens"].append(r["best_len"])

    def mean(xs):
        return sum(xs) / len(xs) if xs else 0.0

    def median(xs):
        s = sorted(xs)
        n = len(s)
        if n == 0:
            return 0.0
        if n % 2 == 1:
            return s[n // 2]
        return (s[n // 2 - 1] + s[n // 2]) / 2.0

    # Table header
    print(f"{'Category':<10} {'N':>4}  {'Mean nseq':>9}  {'Mean seqlen':>11}  "
          f"{'Default ratio':>14}  {'Best ratio':>11}  "
          f"{'Mean ref_len':>12}  {'Mean def_len':>12}  {'Mean best_len':>13}")
    print("-" * 120)

    all_default_ratios = []
    all_best_ratios = []

    for cat in categories:
        s = cat_stats[cat]
        n = len(s["default_ratios"])
        dr = mean(s["default_ratios"])
        br = mean(s["best_ratios"])
        ns = mean(s["nseqs"])
        sl = mean(s["mean_seq_lens"])
        rl = mean(s["ref_lens"])
        dl = mean(s["default_lens"])
        bl = mean(s["best_lens"])

        all_default_ratios.extend(s["default_ratios"])
        all_best_ratios.extend(s["best_ratios"])

        print(f"{cat:<10} {n:>4}  {ns:>9.1f}  {sl:>11.1f}  "
              f"{dr:>14.4f}  {br:>11.4f}  "
              f"{rl:>12.1f}  {dl:>12.1f}  {bl:>13.1f}")

    print("-" * 120)
    print(f"{'OVERALL':<10} {len(ok):>4}  {mean([r['nseq'] for r in ok]):>9.1f}  "
          f"{mean([r['mean_seq_len'] for r in ok]):>11.1f}  "
          f"{mean(all_default_ratios):>14.4f}  {mean(all_best_ratios):>11.4f}  "
          f"{mean([r['ref_len'] for r in ok]):>12.1f}  "
          f"{mean([r['default_len'] for r in ok]):>12.1f}  "
          f"{mean([r['best_len'] for r in ok]):>13.1f}")
    print()

    # ---------------------------------------------------------------------------
    # Distribution of over-gap ratios
    # ---------------------------------------------------------------------------
    print("=" * 72)
    print("DISTRIBUTION OF OVER-GAP RATIOS (default mode)")
    print("=" * 72)
    print()

    for cat in categories:
        ratios = sorted(cat_stats[cat]["default_ratios"])
        n = len(ratios)
        under = sum(1 for r in ratios if r < 1.0)
        at = sum(1 for r in ratios if abs(r - 1.0) < 0.001)
        over_small = sum(1 for r in ratios if 1.001 <= r < 1.1)
        over_medium = sum(1 for r in ratios if 1.1 <= r < 1.3)
        over_large = sum(1 for r in ratios if r >= 1.3)

        print(f"  {cat}: n={n}, min={ratios[0]:.3f}, median={median(ratios):.3f}, "
              f"max={ratios[-1]:.3f}")
        print(f"    ratio<1.0: {under}  ratio~=1.0: {at}  "
              f"1.0-1.1: {over_small}  1.1-1.3: {over_medium}  >1.3: {over_large}")
    print()

    # ---------------------------------------------------------------------------
    # Correlation with TC scores
    # ---------------------------------------------------------------------------
    print("=" * 72)
    print("CORRELATION: OVER-GAP RATIO vs TC SCORE")
    print("=" * 72)
    print()

    # Gather paired (ratio, tc) for default mode
    paired_default = []
    paired_best = []
    for r in ok:
        fam = r["family"]
        tc_def = tc_scores.get((fam, "default"))
        tc_best = tc_scores.get((fam, "best"))
        if tc_def is not None:
            paired_default.append((r["default_ratio"], tc_def, fam, r["category"]))
        if tc_best is not None:
            paired_best.append((r["best_ratio"], tc_best, fam, r["category"]))

    def pearson_corr(pairs):
        """Compute Pearson correlation coefficient."""
        n = len(pairs)
        if n < 3:
            return float("nan")
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        mx = sum(xs) / n
        my = sum(ys) / n
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
        sx = (sum((x - mx) ** 2 for x in xs) / n) ** 0.5
        sy = (sum((y - my) ** 2 for y in ys) / n) ** 0.5
        if sx == 0 or sy == 0:
            return float("nan")
        return cov / (sx * sy)

    if paired_default:
        r_def = pearson_corr(paired_default)
        print(f"Default mode: n={len(paired_default)}, Pearson r = {r_def:.4f}")
        print(f"  (negative r means higher ratio => lower TC, i.e. over-gapping hurts TC)")
    else:
        print("No TC scores available for default mode")

    if paired_best:
        r_best = pearson_corr(paired_best)
        print(f"Best mode:    n={len(paired_best)}, Pearson r = {r_best:.4f}")
    else:
        print("No TC scores available for best mode")
    print()

    # Per-category correlation
    print("Per-category correlation (default mode):")
    for cat in categories:
        cat_pairs = [(ratio, tc, fam, c) for ratio, tc, fam, c in paired_default if c == cat]
        if len(cat_pairs) >= 3:
            r_cat = pearson_corr(cat_pairs)
            print(f"  {cat}: n={len(cat_pairs)}, r = {r_cat:.4f}")
        else:
            print(f"  {cat}: n={len(cat_pairs)}, too few for correlation")
    print()

    # ---------------------------------------------------------------------------
    # TC deficit vs external tools
    # ---------------------------------------------------------------------------
    print("=" * 72)
    print("TC DEFICIT: KALIGN vs EXTERNAL TOOLS")
    print("=" * 72)
    print()

    for ext_tool in ["clustalo", "mafft", "muscle"]:
        deficits = []
        for r in ok:
            fam = r["family"]
            tc_kal = tc_scores.get((fam, "default"))
            tc_ext = tc_scores.get((fam, ext_tool))
            if tc_kal is not None and tc_ext is not None:
                deficits.append({
                    "family": fam,
                    "category": r["category"],
                    "ratio": r["default_ratio"],
                    "tc_kalign": tc_kal,
                    "tc_external": tc_ext,
                    "tc_deficit": tc_kal - tc_ext,
                })

        if not deficits:
            print(f"  {ext_tool}: no paired data")
            continue

        cat_deficit = defaultdict(list)
        for d in deficits:
            cat_deficit[d["category"]].append(d)

        print(f"  vs {ext_tool} (negative deficit = kalign is worse):")
        print(f"    {'Category':<10} {'N':>4}  {'Mean TC deficit':>15}  "
              f"{'Mean ratio':>10}  {'Corr(ratio,deficit)':>20}")
        print(f"    {'-'*70}")

        all_deficits_pairs = []
        for cat in categories:
            items = cat_deficit.get(cat, [])
            if not items:
                continue
            mean_def = mean([d["tc_deficit"] for d in items])
            mean_rat = mean([d["ratio"] for d in items])
            pairs = [(d["ratio"], d["tc_deficit"], d["family"], d["category"]) for d in items]
            corr = pearson_corr(pairs) if len(pairs) >= 3 else float("nan")
            all_deficits_pairs.extend(pairs)
            print(f"    {cat:<10} {len(items):>4}  {mean_def:>15.4f}  "
                  f"{mean_rat:>10.4f}  {corr:>20.4f}")

        overall_corr = pearson_corr(all_deficits_pairs) if len(all_deficits_pairs) >= 3 else float("nan")
        print(f"    {'OVERALL':<10} {len(all_deficits_pairs):>4}  "
              f"{mean([d['tc_deficit'] for d in deficits]):>15.4f}  "
              f"{mean([d['ratio'] for d in deficits]):>10.4f}  "
              f"{overall_corr:>20.4f}")
        print()

    # ---------------------------------------------------------------------------
    # Worst over-gapped cases
    # ---------------------------------------------------------------------------
    print("=" * 72)
    print("TOP 20 MOST OVER-GAPPED CASES (default mode)")
    print("=" * 72)
    print()

    sorted_by_ratio = sorted(ok, key=lambda r: r["default_ratio"], reverse=True)
    print(f"{'Family':<12} {'Cat':<6} {'Ratio':>7} {'RefLen':>7} {'KalLen':>7} "
          f"{'Nseq':>5} {'MeanLen':>8} {'TC_kal':>7} {'TC_mafft':>8}")
    print("-" * 80)
    for r in sorted_by_ratio[:20]:
        fam = r["family"]
        tc_kal = tc_scores.get((fam, "default"), float("nan"))
        tc_mafft = tc_scores.get((fam, "mafft"), float("nan"))
        print(f"{fam:<12} {r['category']:<6} {r['default_ratio']:>7.3f} "
              f"{r['ref_len']:>7} {r['default_len']:>7} "
              f"{r['nseq']:>5} {r['mean_seq_len']:>8.1f} "
              f"{tc_kal:>7.3f} {tc_mafft:>8.3f}")
    print()

    # ---------------------------------------------------------------------------
    # Save raw data
    # ---------------------------------------------------------------------------
    out_path = Path(__file__).parent / "results" / "overgap_analysis.json"
    with open(out_path, "w") as f:
        json.dump({"results": ok, "errors": errors}, f, indent=2)
    print(f"Raw data saved to {out_path}")
    print()

    # ---------------------------------------------------------------------------
    # Key findings summary
    # ---------------------------------------------------------------------------
    print("=" * 72)
    print("KEY FINDINGS")
    print("=" * 72)
    print()

    overall_default = mean(all_default_ratios)
    overall_best = mean(all_best_ratios)
    print(f"1. Overall mean alignment length ratio:")
    print(f"   Default mode: {overall_default:.4f} (kalign alignments are "
          f"{(overall_default - 1) * 100:.1f}% longer than reference)")
    print(f"   Best mode:    {overall_best:.4f} (kalign alignments are "
          f"{(overall_best - 1) * 100:.1f}% longer than reference)")
    print()

    # Category with worst over-gapping
    worst_cat = max(categories, key=lambda c: mean(cat_stats[c]["default_ratios"]))
    worst_ratio = mean(cat_stats[worst_cat]["default_ratios"])
    print(f"2. Worst category for over-gapping: {worst_cat} "
          f"(mean ratio = {worst_ratio:.4f}, "
          f"{(worst_ratio - 1) * 100:.1f}% over reference)")
    print()

    if paired_default:
        r_def = pearson_corr(paired_default)
        print(f"3. Correlation between over-gap ratio and TC score: r = {r_def:.4f}")
        if r_def < -0.1:
            print("   => Over-gapping IS negatively correlated with TC score")
        elif r_def > 0.1:
            print("   => Surprisingly, longer alignments correlate with HIGHER TC")
        else:
            print("   => Weak or no correlation between alignment length and TC")
    print()

    # Does best mode reduce over-gapping?
    if all_default_ratios and all_best_ratios:
        improvement = mean(all_default_ratios) - mean(all_best_ratios)
        print(f"4. Best mode vs default mode:")
        if improvement > 0:
            print(f"   Best mode reduces over-gapping by {improvement:.4f} "
                  f"({improvement / mean(all_default_ratios) * 100:.1f}%)")
        else:
            print(f"   Best mode INCREASES alignment length by {-improvement:.4f}")
    print()


if __name__ == "__main__":
    main()
