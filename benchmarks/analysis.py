"""RV11 alignment structure analysis.

Compares gap structure and alignment geometry between kalign, reference,
and external tools (mafft, muscle, clustalo) to understand HOW alignments
differ structurally — not just score differences.

Usage:
    # Locally (kalign only, external tools skipped if not installed):
    uv run python -m benchmarks.analysis

    # Inside container (has mafft/muscle/clustalo):
    python -m benchmarks.analysis

    # Specific dataset:
    python -m benchmarks.analysis --dataset balibase_RV11

    # Write CSV:
    python -m benchmarks.analysis --csv benchmarks/results/gap_analysis.csv
"""

import argparse
import csv
import json
import re
import statistics
import sys
import tempfile
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .datasets import balibase_cases, balibase_is_available

RESULTS_DIR = Path(__file__).parent / "results"


# ---------------------------------------------------------------------------
# MSF parser (reference alignments are in GCG MSF format)
# ---------------------------------------------------------------------------

def parse_msf(path: Path) -> Dict[str, str]:
    """Parse a GCG MSF file into {name: aligned_sequence}."""
    text = path.read_text()

    # Split at "//" separator
    parts = text.split("//")
    if len(parts) < 2:
        raise ValueError(f"No // separator found in {path}")

    body = parts[1]
    seqs: Dict[str, List[str]] = {}
    for line in body.splitlines():
        line = line.strip()
        if not line:
            continue
        tokens = line.split()
        if len(tokens) < 2:
            continue
        name = tokens[0]
        # Sequence characters (may contain dots for gaps)
        seq_parts = "".join(tokens[1:])
        seqs.setdefault(name, []).append(seq_parts)

    # Join blocks and normalise: dots → dashes, remove whitespace
    result = {}
    for name, blocks in seqs.items():
        seq = "".join(blocks).replace(".", "-").upper()
        result[name] = seq
    return result


# ---------------------------------------------------------------------------
# FASTA parser (kalign/tool outputs)
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> Dict[str, str]:
    """Parse a FASTA file into {name: sequence}."""
    seqs: Dict[str, str] = {}
    current = None
    parts: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(parts).upper()
            current = line[1:].split()[0]
            parts = []
        elif current is not None:
            parts.append(line)
    if current is not None:
        seqs[current] = "".join(parts).upper()
    return seqs


# ---------------------------------------------------------------------------
# Gap structure metrics
# ---------------------------------------------------------------------------

@dataclass
class GapStats:
    """Gap structure metrics for one alignment."""
    n_seqs: int
    alignment_length: int
    mean_seq_length: float  # unaligned (non-gap chars)
    expansion_factor: float  # alignment_length / mean_seq_length
    total_gaps: int
    gap_fraction: float  # total_gaps / (n_seqs * alignment_length)
    n_gap_blocks: int
    mean_gap_block_len: float
    mean_terminal_gap: float  # average leading+trailing gap per sequence
    mean_internal_gap: float  # average total internal gap chars per sequence
    n_gappy_columns: int  # columns where >50% of sequences have a gap
    gappy_column_fraction: float


def compute_gap_stats(seqs: Dict[str, str]) -> GapStats:
    """Compute gap structure metrics from aligned sequences."""
    sequences = list(seqs.values())
    n_seqs = len(sequences)
    if n_seqs == 0:
        return GapStats(0, 0, 0.0, 0.0, 0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0.0)

    aln_len = len(sequences[0])
    # Unaligned lengths (non-gap characters)
    ungapped_lens = [len(s.replace("-", "")) for s in sequences]
    mean_seq_len = statistics.mean(ungapped_lens)
    expansion = aln_len / mean_seq_len if mean_seq_len > 0 else 0.0

    total_gaps = sum(s.count("-") for s in sequences)
    total_chars = n_seqs * aln_len
    gap_frac = total_gaps / total_chars if total_chars > 0 else 0.0

    # Gap blocks and lengths
    all_block_lens: List[int] = []
    terminal_gaps: List[int] = []
    internal_gaps: List[int] = []

    for seq in sequences:
        # Find all gap blocks
        blocks = [(m.start(), m.end()) for m in re.finditer(r"-+", seq)]
        all_block_lens.extend(m.end() - m.start() for m in re.finditer(r"-+", seq))

        # Terminal: leading and trailing
        leading = len(seq) - len(seq.lstrip("-"))
        trailing = len(seq) - len(seq.rstrip("-"))
        terminal_gaps.append(leading + trailing)

        # Internal: everything that's not leading/trailing
        internal = sum(e - s for s, e in blocks)
        internal -= leading + trailing
        internal_gaps.append(max(0, internal))

    n_gap_blocks = len(all_block_lens)
    mean_block_len = statistics.mean(all_block_lens) if all_block_lens else 0.0
    mean_terminal = statistics.mean(terminal_gaps)
    mean_internal = statistics.mean(internal_gaps)

    # Gappy columns (>50% gaps)
    n_gappy = 0
    for col in range(aln_len):
        gaps_in_col = sum(1 for s in sequences if s[col] == "-")
        if gaps_in_col > n_seqs / 2:
            n_gappy += 1

    return GapStats(
        n_seqs=n_seqs,
        alignment_length=aln_len,
        mean_seq_length=mean_seq_len,
        expansion_factor=expansion,
        total_gaps=total_gaps,
        gap_fraction=gap_frac,
        n_gap_blocks=n_gap_blocks,
        mean_gap_block_len=mean_block_len,
        mean_terminal_gap=mean_terminal,
        mean_internal_gap=mean_internal,
        n_gappy_columns=n_gappy,
        gappy_column_fraction=n_gappy / aln_len if aln_len > 0 else 0.0,
    )


# ---------------------------------------------------------------------------
# Alignment generation helpers
# ---------------------------------------------------------------------------

def _align_kalign(unaligned: Path, output: Path, seq_type: str) -> None:
    """Run kalign via Python API."""
    import kalign
    kalign.align_file_to_file(
        str(unaligned), str(output), format="fasta", seq_type=seq_type,
    )


def _align_external(unaligned: Path, output: Path, tool: str) -> bool:
    """Run an external tool. Returns True if successful."""
    import shutil
    import subprocess

    if shutil.which(tool) is None:
        return False

    try:
        if tool == "mafft":
            with open(output, "w") as f:
                subprocess.run(
                    ["mafft", "--auto", str(unaligned)],
                    stdout=f, stderr=subprocess.PIPE, check=True,
                )
        elif tool == "clustalo":
            subprocess.run(
                ["clustalo", "-i", str(unaligned), "-o", str(output),
                 "--outfmt=fasta", "--force"],
                capture_output=True, check=True,
            )
        elif tool == "muscle":
            subprocess.run(
                ["muscle", "-align", str(unaligned), "-output", str(output)],
                capture_output=True, check=True,
            )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


# ---------------------------------------------------------------------------
# Per-case analysis row
# ---------------------------------------------------------------------------

@dataclass
class CaseRow:
    family: str
    method: str
    # Scores from results JSON (NaN if not available)
    recall: float
    precision: float
    f1: float
    tc: float
    # Gap stats
    alignment_length: int
    expansion_factor: float
    gap_fraction: float
    n_gap_blocks: int
    mean_gap_block_len: float
    mean_terminal_gap: float
    mean_internal_gap: float
    n_gappy_columns: int
    gappy_column_fraction: float


# ---------------------------------------------------------------------------
# Load frozen scores from full_comparison.json
# ---------------------------------------------------------------------------

def load_scores(json_path: Path, dataset_filter: str = "balibase_RV11") -> Dict[Tuple[str, str], dict]:
    """Load {(family, method_key): scores} from results JSON.

    method_key is 'kalign' for python_api/refine=none, or the external tool name.
    """
    with open(json_path) as f:
        data = json.load(f)

    scores: Dict[Tuple[str, str], dict] = {}
    for r in data["results"]:
        if dataset_filter and r["dataset"] != dataset_filter:
            continue

        if r["method"] == "python_api" and r["refine"] == "none":
            key = (r["family"], "kalign")
        elif r["method"] in ("clustalo", "mafft", "muscle"):
            key = (r["family"], r["method"])
        else:
            continue

        scores[key] = {
            "recall": r.get("recall", float("nan")),
            "precision": r.get("precision", float("nan")),
            "f1": r.get("f1", float("nan")),
            "tc": r.get("tc", float("nan")),
        }
    return scores


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyse_rv11(
    dataset_filter: str = "balibase_RV11",
    external_tools: Optional[List[str]] = None,
) -> List[CaseRow]:
    """Run gap analysis on all RV11 cases. Returns list of CaseRow."""
    if external_tools is None:
        external_tools = ["mafft", "muscle", "clustalo"]

    if not balibase_is_available():
        print("BAliBASE not downloaded. Run: uv run python -m benchmarks --download-only")
        return []

    cases = [c for c in balibase_cases() if c.dataset == dataset_filter]
    if not cases:
        print(f"No cases found for {dataset_filter}")
        return []

    # Load frozen scores
    json_path = RESULTS_DIR / "full_comparison.json"
    scores = load_scores(json_path, dataset_filter) if json_path.exists() else {}

    rows: List[CaseRow] = []

    for case in cases:
        print(f"  {case.family} ...", end="", flush=True)

        # --- Reference alignment ---
        ref_seqs = parse_msf(case.reference)
        ref_stats = compute_gap_stats(ref_seqs)
        sc = scores.get((case.family, "reference"), {})
        rows.append(CaseRow(
            family=case.family, method="reference",
            recall=1.0, precision=1.0, f1=1.0, tc=1.0,
            alignment_length=ref_stats.alignment_length,
            expansion_factor=ref_stats.expansion_factor,
            gap_fraction=ref_stats.gap_fraction,
            n_gap_blocks=ref_stats.n_gap_blocks,
            mean_gap_block_len=ref_stats.mean_gap_block_len,
            mean_terminal_gap=ref_stats.mean_terminal_gap,
            mean_internal_gap=ref_stats.mean_internal_gap,
            n_gappy_columns=ref_stats.n_gappy_columns,
            gappy_column_fraction=ref_stats.gappy_column_fraction,
        ))

        # --- Kalign ---
        with tempfile.TemporaryDirectory() as tmpdir:
            kalign_out = Path(tmpdir) / f"{case.family}_kalign.fa"
            _align_kalign(case.unaligned, kalign_out, case.seq_type)
            kalign_seqs = parse_fasta(kalign_out)
            kalign_stats = compute_gap_stats(kalign_seqs)

        sc = scores.get((case.family, "kalign"), {})
        rows.append(CaseRow(
            family=case.family, method="kalign",
            recall=sc.get("recall", float("nan")),
            precision=sc.get("precision", float("nan")),
            f1=sc.get("f1", float("nan")),
            tc=sc.get("tc", float("nan")),
            alignment_length=kalign_stats.alignment_length,
            expansion_factor=kalign_stats.expansion_factor,
            gap_fraction=kalign_stats.gap_fraction,
            n_gap_blocks=kalign_stats.n_gap_blocks,
            mean_gap_block_len=kalign_stats.mean_gap_block_len,
            mean_terminal_gap=kalign_stats.mean_terminal_gap,
            mean_internal_gap=kalign_stats.mean_internal_gap,
            n_gappy_columns=kalign_stats.n_gappy_columns,
            gappy_column_fraction=kalign_stats.gappy_column_fraction,
        ))

        # --- External tools ---
        for tool in external_tools:
            with tempfile.TemporaryDirectory() as tmpdir:
                tool_out = Path(tmpdir) / f"{case.family}_{tool}.fa"
                ok = _align_external(case.unaligned, tool_out, tool)
                if not ok:
                    continue
                tool_seqs = parse_fasta(tool_out)
                tool_stats = compute_gap_stats(tool_seqs)

            sc = scores.get((case.family, tool), {})
            rows.append(CaseRow(
                family=case.family, method=tool,
                recall=sc.get("recall", float("nan")),
                precision=sc.get("precision", float("nan")),
                f1=sc.get("f1", float("nan")),
                tc=sc.get("tc", float("nan")),
                alignment_length=tool_stats.alignment_length,
                expansion_factor=tool_stats.expansion_factor,
                gap_fraction=tool_stats.gap_fraction,
                n_gap_blocks=tool_stats.n_gap_blocks,
                mean_gap_block_len=tool_stats.mean_gap_block_len,
                mean_terminal_gap=tool_stats.mean_terminal_gap,
                mean_internal_gap=tool_stats.mean_internal_gap,
                n_gappy_columns=tool_stats.n_gappy_columns,
                gappy_column_fraction=tool_stats.gappy_column_fraction,
            ))

        print(" done")

    return rows


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def print_table(rows: List[CaseRow]) -> None:
    """Print a summary table to stdout, grouped by family."""
    if not rows:
        return

    families = sorted(set(r.family for r in rows))
    methods = sorted(set(r.method for r in rows))

    # Per-case comparison table
    hdr = f"{'Family':<10} {'Method':<10} {'Recall':>7} {'Prec':>7} {'F1':>7} {'AlnLen':>7} {'Expand':>7} {'GapFrac':>7} {'GapBlk':>7} {'MeanBL':>7} {'TermGap':>7} {'IntGap':>7} {'Gappy%':>7}"
    print("\n" + "=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))

    for fam in families:
        fam_rows = sorted(
            [r for r in rows if r.family == fam],
            key=lambda r: (r.method != "reference", r.method != "kalign", r.method),
        )
        for r in fam_rows:
            rec = f"{r.recall:.3f}" if r.recall == r.recall else "   n/a"
            pre = f"{r.precision:.3f}" if r.precision == r.precision else "   n/a"
            f1 = f"{r.f1:.3f}" if r.f1 == r.f1 else "   n/a"
            print(
                f"{r.family:<10} {r.method:<10} {rec:>7} {pre:>7} {f1:>7} "
                f"{r.alignment_length:>7} {r.expansion_factor:>7.2f} "
                f"{r.gap_fraction:>7.3f} {r.n_gap_blocks:>7} "
                f"{r.mean_gap_block_len:>7.1f} {r.mean_terminal_gap:>7.1f} "
                f"{r.mean_internal_gap:>7.1f} {r.gappy_column_fraction:>7.3f}"
            )
        print()

    # Aggregate summary by method
    print("=" * 80)
    print("AGGREGATE SUMMARY (means across all families)")
    print("-" * 80)
    fmt = "{:<10} {:>7} {:>7} {:>7} {:>8} {:>7} {:>7} {:>8} {:>8}"
    print(fmt.format("Method", "Recall", "Prec", "F1", "Expand", "GapFrac", "MeanBL", "TermGap", "IntGap"))
    print("-" * 80)

    for method in ["reference", "kalign"] + [m for m in methods if m not in ("reference", "kalign")]:
        method_rows = [r for r in rows if r.method == method]
        if not method_rows:
            continue
        def safe_mean(vals):
            clean = [v for v in vals if v == v]  # filter NaN
            return statistics.mean(clean) if clean else float("nan")

        rec = safe_mean([r.recall for r in method_rows])
        pre = safe_mean([r.precision for r in method_rows])
        f1 = safe_mean([r.f1 for r in method_rows])
        exp = statistics.mean([r.expansion_factor for r in method_rows])
        gf = statistics.mean([r.gap_fraction for r in method_rows])
        mbl = statistics.mean([r.mean_gap_block_len for r in method_rows])
        tg = statistics.mean([r.mean_terminal_gap for r in method_rows])
        ig = statistics.mean([r.mean_internal_gap for r in method_rows])

        rec_s = f"{rec:.3f}" if rec == rec else "   n/a"
        pre_s = f"{pre:.3f}" if pre == pre else "   n/a"
        f1_s = f"{f1:.3f}" if f1 == f1 else "   n/a"
        print(fmt.format(method, rec_s, pre_s, f1_s, f"{exp:.2f}", f"{gf:.3f}", f"{mbl:.1f}", f"{tg:.1f}", f"{ig:.1f}"))

    # Correlation analysis: precision vs expansion factor for kalign
    kalign_rows = [r for r in rows if r.method == "kalign" and r.precision == r.precision]
    if len(kalign_rows) >= 5:
        print("\n" + "=" * 80)
        print("CORRELATION: kalign precision vs gap metrics")
        print("-" * 80)

        # Simple Pearson correlation
        def pearson(xs, ys):
            n = len(xs)
            if n < 3:
                return float("nan")
            mx, my = statistics.mean(xs), statistics.mean(ys)
            num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            dx = sum((x - mx) ** 2 for x in xs) ** 0.5
            dy = sum((y - my) ** 2 for y in ys) ** 0.5
            return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")

        prec = [r.precision for r in kalign_rows]
        metrics = [
            ("expansion_factor", [r.expansion_factor for r in kalign_rows]),
            ("gap_fraction", [r.gap_fraction for r in kalign_rows]),
            ("mean_gap_block_len", [r.mean_gap_block_len for r in kalign_rows]),
            ("mean_terminal_gap", [r.mean_terminal_gap for r in kalign_rows]),
            ("mean_internal_gap", [r.mean_internal_gap for r in kalign_rows]),
            ("gappy_column_fraction", [r.gappy_column_fraction for r in kalign_rows]),
        ]
        for name, vals in metrics:
            r = pearson(prec, vals)
            print(f"  precision vs {name:<25}: r = {r:+.3f}")

        # Also: kalign expansion vs reference expansion
        ref_rows = {r.family: r for r in rows if r.method == "reference"}
        kalign_dict = {r.family: r for r in kalign_rows}
        common = sorted(set(ref_rows) & set(kalign_dict))
        if common:
            print(f"\n  Expansion factor comparison (kalign vs reference, n={len(common)}):")
            over = [f for f in common if kalign_dict[f].expansion_factor > ref_rows[f].expansion_factor * 1.05]
            under = [f for f in common if kalign_dict[f].expansion_factor < ref_rows[f].expansion_factor * 0.95]
            same = [f for f in common if f not in over and f not in under]
            print(f"    Over-expanded (>5%):  {len(over)}")
            print(f"    Under-expanded:       {len(under)}")
            print(f"    Similar (+/-5%):      {len(same)}")

            # Relative expansion ratio correlated with precision
            ratios = [kalign_dict[f].expansion_factor / max(ref_rows[f].expansion_factor, 0.01) for f in common]
            precs = [kalign_dict[f].precision for f in common]
            r_ratio = pearson(ratios, precs)
            print(f"    Correlation (expand_ratio vs precision): r = {r_ratio:+.3f}")

            if over:
                print(f"\n  Worst over-expanded cases (kalign expand / ref expand):")
                ranked = sorted(over, key=lambda f: kalign_dict[f].expansion_factor / max(ref_rows[f].expansion_factor, 0.01), reverse=True)
                for f in ranked[:10]:
                    ke = kalign_dict[f].expansion_factor
                    re = ref_rows[f].expansion_factor
                    kp = kalign_dict[f].precision
                    print(f"    {f}: kalign={ke:.2f} ref={re:.2f} ratio={ke/re:.2f} precision={kp:.3f}")

            if under:
                print(f"\n  Worst under-expanded cases (kalign expand / ref expand):")
                ranked = sorted(under, key=lambda f: kalign_dict[f].expansion_factor / max(ref_rows[f].expansion_factor, 0.01))
                for f in ranked[:10]:
                    ke = kalign_dict[f].expansion_factor
                    re = ref_rows[f].expansion_factor
                    kp = kalign_dict[f].precision
                    print(f"    {f}: kalign={ke:.2f} ref={re:.2f} ratio={ke/re:.2f} precision={kp:.3f}")


def write_csv(rows: List[CaseRow], path: str) -> None:
    """Write analysis results to CSV."""
    fieldnames = [f.name for f in fields(CaseRow)]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            d = {fn: getattr(row, fn) for fn in fieldnames}
            writer.writerow(d)
    print(f"\nCSV written to {path}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="RV11 alignment structure analysis",
        prog="python -m benchmarks.analysis",
    )
    parser.add_argument(
        "--dataset", default="balibase_RV11",
        help="Dataset filter (default: balibase_RV11)",
    )
    parser.add_argument(
        "--csv", default="",
        help="Write results to CSV file",
    )
    parser.add_argument(
        "--no-external", action="store_true",
        help="Skip external tools (mafft, muscle, clustalo)",
    )
    args = parser.parse_args()

    external = [] if args.no_external else ["mafft", "muscle", "clustalo"]

    print(f"Analysing {args.dataset} alignment structure...")
    rows = analyse_rv11(dataset_filter=args.dataset, external_tools=external)

    if not rows:
        sys.exit(1)

    print_table(rows)

    if args.csv:
        write_csv(rows, args.csv)


if __name__ == "__main__":
    main()
