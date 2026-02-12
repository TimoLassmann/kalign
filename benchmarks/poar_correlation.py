"""Measure correlation between POAR agreement scores and SP scores.

For each benchmark case, runs all presets independently, then:
1. Computes SP score against reference (ground truth quality)
2. Computes MUMSA-style POAR agreement score (how much each alignment
   agrees with the others)
3. Reports correlation between the two

This answers: "Can POAR scores reliably identify the best alignment?"

Usage:
    uv run python -m benchmarks.poar_correlation [--max-cases 50] [-j 8]
"""

import argparse
import tempfile
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kalign

from .datasets import get_cases

# Same presets as ensemble.c
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

BASE_GPO = 5.5
BASE_GPE = 2.0
BASE_TGPE = 1.0


def parse_fasta_alignment(path):
    """Parse a FASTA alignment file. Returns list of (name, seq) tuples."""
    seqs = []
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(parts)))
                name = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        seqs.append((name, "".join(parts)))
    return seqs


def extract_poars(seqs):
    """Extract all pairs of aligned residues from an alignment.

    Returns dict: (seq_i, seq_j) -> set of (res_pos_i, res_pos_j)
    """
    n = len(seqs)
    alnlen = len(seqs[0][1])

    # Precompute residue positions for each sequence
    res_pos = []  # res_pos[seq][col] = residue position or -1
    for s in range(n):
        positions = []
        pos = 0
        for col in range(alnlen):
            if seqs[s][1][col] != '-':
                positions.append(pos)
                pos += 1
            else:
                positions.append(-1)
        res_pos.append(positions)

    poars = {}
    for i in range(n - 1):
        for j in range(i + 1, n):
            pairs = set()
            for col in range(alnlen):
                ri = res_pos[i][col]
                rj = res_pos[j][col]
                if ri >= 0 and rj >= 0:
                    pairs.add((ri, rj))
            poars[(i, j)] = pairs
    return poars


def compute_poar_score(alignment_poars, all_poars_list, aln_idx):
    """Compute MUMSA-style POAR agreement score for one alignment.

    For each aligned pair in this alignment, count how many OTHER
    alignments also have that pair. Score = sum((support-1)/(m-1))
    where support includes self.

    This matches score_alignment_poar() in consensus_msa.c.
    """
    n_aln = len(all_poars_list)
    if n_aln <= 1:
        return 0.0

    denom = n_aln - 1
    total = 0.0

    for pair_key, pairs in alignment_poars.items():
        for poar in pairs:
            # Count how many alignments have this POAR
            support = sum(1 for k in range(n_aln)
                         if poar in all_poars_list[k].get(pair_key, set()))
            # support includes self; subtract 1 for other-agreement
            total += (support - 1) / denom

    return total


def _run_one_preset(args):
    """Run one preset on one case. Returns (family, preset_idx, sp_score, aln_path)."""
    case, preset_idx, preset, keep_dir = args
    gpo = BASE_GPO * preset["gpo"]
    gpe = BASE_GPE * preset["gpe"]
    tgpe = BASE_TGPE * preset["tgpe"]

    output = Path(keep_dir) / f"preset_{preset_idx}.fasta"
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
    return case.family, preset_idx, sp, str(output)


def main():
    parser = argparse.ArgumentParser(description="POAR vs SP score correlation")
    parser.add_argument("--max-cases", type=int, default=50)
    parser.add_argument("-j", "--parallel", type=int, default=1)
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases)
    n_presets = len(PRESETS)
    print(f"Running {n_presets} presets on {len(cases)} cases...")

    # For each case: run all presets, keep alignment files, compute POAR scores
    all_results = []  # list of (family, [(preset_idx, sp, poar_score), ...])

    for ci, case in enumerate(cases):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run all presets for this case
            if args.parallel <= 1:
                preset_results = []
                for pi, preset in enumerate(PRESETS):
                    item = (case, pi, preset, tmpdir)
                    fam, pi_out, sp, aln_path = _run_one_preset(item)
                    preset_results.append((pi_out, sp, aln_path))
            else:
                preset_results = [None] * n_presets
                work = [(case, pi, preset, tmpdir) for pi, preset in enumerate(PRESETS)]
                with ProcessPoolExecutor(max_workers=args.parallel) as pool:
                    futures = {pool.submit(_run_one_preset, item): item[1] for item in work}
                    for future in as_completed(futures):
                        fam, pi_out, sp, aln_path = future.result()
                        preset_results[pi_out] = (pi_out, sp, aln_path)

            # Parse all alignments and extract POARs
            all_poars = []  # one per preset
            for pi, sp, aln_path in preset_results:
                seqs = parse_fasta_alignment(aln_path)
                poars = extract_poars(seqs)
                all_poars.append(poars)

            # Compute POAR agreement score for each preset's alignment
            case_results = []
            for pi, sp, aln_path in preset_results:
                poar_score = compute_poar_score(all_poars[pi], all_poars, pi)
                case_results.append((pi, sp, poar_score))

            all_results.append((case.family, case_results))

        print(f"  [{ci+1}/{len(cases)}] {case.family} done", flush=True)

    # === Analysis 1: Per-case correlation ===
    print(f"\n{'Family':<12s} {'Pearson r':>10s}  {'POAR best':>10s} {'SP best':>10s} {'Match?':>7s}")
    print("-" * 55)

    per_case_correlations = []
    poar_selects_best = 0
    poar_selects_worst = 0
    total_cases = len(all_results)

    for family, results in all_results:
        sps = [r[1] for r in results]
        poars = [r[2] for r in results]

        # Which preset has best POAR score?
        best_poar_idx = max(range(n_presets), key=lambda k: poars[k])
        # Which preset has best SP score?
        best_sp_idx = max(range(n_presets), key=lambda k: sps[k])
        # Which has worst SP?
        worst_sp_idx = min(range(n_presets), key=lambda k: sps[k])

        match = best_poar_idx == best_sp_idx

        # Pearson correlation
        r = _pearson(sps, poars)
        per_case_correlations.append(r)

        if match:
            poar_selects_best += 1
        if best_poar_idx == worst_sp_idx:
            poar_selects_worst += 1

        print(f"{family:<12s} {r:10.3f}  {PRESETS[best_poar_idx]['name']:>10s} "
              f"{PRESETS[best_sp_idx]['name']:>10s} {'YES' if match else 'no':>7s}")

    print(f"\nSummary:")
    print(f"  Mean Pearson r: {statistics.mean(per_case_correlations):.3f}")
    print(f"  Median Pearson r: {statistics.median(per_case_correlations):.3f}")
    print(f"  POAR selects best SP: {poar_selects_best}/{total_cases} "
          f"({100*poar_selects_best/total_cases:.0f}%)")
    print(f"  POAR selects worst SP: {poar_selects_worst}/{total_cases} "
          f"({100*poar_selects_worst/total_cases:.0f}%)")

    # === Analysis 2: SP of POAR-selected vs default vs oracle ===
    print(f"\n=== Selection quality ===")
    default_sps = []
    poar_selected_sps = []
    oracle_sps = []

    for family, results in all_results:
        sps = [r[1] for r in results]
        poars = [r[2] for r in results]

        default_sps.append(sps[0])
        oracle_sps.append(max(sps))

        best_poar_idx = max(range(n_presets), key=lambda k: poars[k])
        poar_selected_sps.append(sps[best_poar_idx])

    print(f"  Default mean SP:       {statistics.mean(default_sps):.2f}")
    print(f"  POAR-selected mean SP: {statistics.mean(poar_selected_sps):.2f} "
          f"(delta: {statistics.mean(poar_selected_sps) - statistics.mean(default_sps):+.2f})")
    print(f"  Oracle mean SP:        {statistics.mean(oracle_sps):.2f} "
          f"(delta: {statistics.mean(oracle_sps) - statistics.mean(default_sps):+.2f})")
    print(f"  POAR captures {100*(statistics.mean(poar_selected_sps) - statistics.mean(default_sps)) / max(0.01, statistics.mean(oracle_sps) - statistics.mean(default_sps)):.1f}% of oracle gap")

    # === Analysis 3: Per-preset aggregated POAR vs SP ===
    print(f"\n{'Preset':<20s} {'Mean SP':>8s} {'Mean POAR':>10s} {'SP rank':>8s} {'POAR rank':>10s}")
    print("-" * 60)

    preset_mean_sp = []
    preset_mean_poar = []
    for pi in range(n_presets):
        sps = [results[pi][1] for _, results in all_results]
        poars = [results[pi][2] for _, results in all_results]
        preset_mean_sp.append(statistics.mean(sps))
        preset_mean_poar.append(statistics.mean(poars))

    sp_ranks = _rank_descending(preset_mean_sp)
    poar_ranks = _rank_descending(preset_mean_poar)

    for pi in range(n_presets):
        print(f"{PRESETS[pi]['name']:<20s} {preset_mean_sp[pi]:8.2f} {preset_mean_poar[pi]:10.1f} "
              f"{sp_ranks[pi]:8d} {poar_ranks[pi]:10d}")

    print(f"\nPreset-level Pearson r (mean SP vs mean POAR): "
          f"{_pearson(preset_mean_sp, preset_mean_poar):.3f}")


def _pearson(x, y):
    """Compute Pearson correlation coefficient."""
    n = len(x)
    if n < 2:
        return 0.0
    mx = sum(x) / n
    my = sum(y) / n
    sx = (sum((xi - mx) ** 2 for xi in x)) ** 0.5
    sy = (sum((yi - my) ** 2 for yi in y)) ** 0.5
    if sx == 0 or sy == 0:
        return 0.0
    return sum((xi - mx) * (yi - my) for xi, yi in zip(x, y)) / (sx * sy)


def _rank_descending(values):
    """Return 1-based ranks (1 = highest value)."""
    indexed = sorted(enumerate(values), key=lambda x: -x[1])
    ranks = [0] * len(values)
    for rank, (idx, _) in enumerate(indexed, 1):
        ranks[idx] = rank
    return ranks


if __name__ == "__main__":
    main()
