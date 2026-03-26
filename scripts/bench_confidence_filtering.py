#!/usr/bin/env python3
"""Benchmark: SP/TC on confident columns only.

Scores BAliBASE alignments using only columns above a confidence threshold.
Shows that kalign's confidence scores identify reliable columns, and that
SP/TC on those columns is dramatically higher than on all columns.

Usage:
    uv run python scripts/bench_confidence_filtering.py
"""
import json
import os
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import kalign

BB_ROOT = Path(__file__).parent.parent / "benchmarks" / "data" / "downloads" / "bb3_release"


def parse_balibase_xml(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    colsco = root.find(".//column-score/colsco-data")
    if colsco is None or colsco.text is None:
        return None
    values = [int(v) for v in colsco.text.split()]
    return [1 if v == 1 else 0 for v in values]


def discover_cases():
    cases = []
    for rv_dir in sorted(BB_ROOT.iterdir()):
        if not rv_dir.is_dir() or not rv_dir.name.startswith("RV"):
            continue
        for tfa in sorted(rv_dir.glob("*.tfa")):
            if tfa.stem.startswith("BBS"):
                continue
            msf = tfa.with_suffix(".msf")
            xml = tfa.with_suffix(".xml")
            if msf.exists():
                cases.append({
                    "family": tfa.stem,
                    "category": rv_dir.name,
                    "unaligned": str(tfa),
                    "reference": str(msf),
                    "xml": str(xml) if xml.exists() else None,
                })
    return cases


def score_with_xml(ref_path, test_path, xml_path):
    """Score alignment using XML core block mask (standard BAliBASE scoring)."""
    if xml_path and Path(xml_path).exists():
        mask = parse_balibase_xml(xml_path)
        if mask:
            return kalign.compare_detailed(ref_path, test_path, column_mask=mask)
    return kalign.compare_detailed(ref_path, test_path, max_gap_frac=0.2)


def score_filtered(ref_path, test_sequences, test_names, conf, threshold, xml_path, tmpdir, family):
    """Score only confident columns: mask out low-confidence residues, then score."""
    # Write a filtered test alignment: low-confidence residues become gaps
    filtered_out = os.path.join(tmpdir, f"{family}_filt.fa")
    with open(filtered_out, 'w') as f:
        for name, seq in zip(test_names, test_sequences):
            filtered = []
            for col in range(len(seq)):
                if col < len(conf) and conf[col] >= threshold:
                    filtered.append(seq[col])
                else:
                    filtered.append('-')
            f.write(f">{name}\n{''.join(filtered)}\n")

    # Score with XML mask (normal BAliBASE scoring on the filtered alignment)
    return score_with_xml(ref_path, filtered_out, xml_path)


def main():
    if not BB_ROOT.exists():
        print(f"ERROR: BAliBASE not found at {BB_ROOT}")
        return 1

    cases = discover_cases()
    thresholds = [0.0, 0.3, 0.5, 0.7, 0.9]

    print(f"{len(cases)} BAliBASE cases, accurate mode")
    print()

    results_by_thresh = {t: {"recall": [], "precision": [], "f1": [], "tc": [], "n_cols": [], "n_total": []}
                         for t in thresholds}

    with tempfile.TemporaryDirectory() as tmpdir:
        for i, case in enumerate(cases):
            # Use align_from_file to get BOTH alignment and confidence in one call
            aln = kalign.align_from_file(case["unaligned"], mode="accurate")
            conf = aln.column_confidence

            # Write the alignment to file for scoring
            out = os.path.join(tmpdir, f"{case['family']}.fa")
            with open(out, 'w') as f:
                for name, seq in zip(aln.names, aln.sequences):
                    f.write(f">{name}\n{seq}\n")

            for thresh in thresholds:
                if thresh == 0.0:
                    # No confidence filtering — standard XML scoring
                    score = score_with_xml(case["reference"], out, case["xml"])
                else:
                    if conf is None:
                        continue
                    # Mask low-confidence residues to gaps, then score normally
                    score = score_filtered(
                        case["reference"], aln.sequences, aln.names,
                        conf, thresh, case["xml"], tmpdir, case["family"]
                    )

                if score is not None:
                    for k in ["recall", "precision", "f1", "tc"]:
                        results_by_thresh[thresh][k].append(score[k])
                    if conf is not None:
                        n_confident = sum(1 for c in conf if c >= thresh)
                        results_by_thresh[thresh]["n_cols"].append(n_confident)
                        results_by_thresh[thresh]["n_total"].append(len(conf))

            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(cases)}...")

    # Print results
    print()
    print(f"{'Threshold':>10s} {'Recall':>8s} {'Prec':>8s} {'F1':>8s} {'TC':>8s} {'Cols%':>7s}  (n)")
    print("=" * 60)
    for thresh in thresholds:
        r = results_by_thresh[thresh]
        n = len(r["recall"])
        if n == 0:
            continue
        avg_r = sum(r["recall"]) / n
        avg_p = sum(r["precision"]) / n
        avg_f = sum(r["f1"]) / n
        avg_t = sum(r["tc"]) / n

        if r["n_cols"] and r["n_total"]:
            avg_pct = sum(c / t * 100 for c, t in zip(r["n_cols"], r["n_total"])) / len(r["n_cols"])
        else:
            avg_pct = 100.0

        label = "all" if thresh == 0.0 else f">={thresh}"
        print(f"{label:>10s} {avg_r:8.3f} {avg_p:8.3f} {avg_f:8.3f} {avg_t:8.3f} {avg_pct:6.1f}%  ({n})")

    # Also print comparison to other tools on all columns
    print()
    print("Reference (all columns, from manuscript):")
    print(f"{'mafft':>10s} {'0.867':>8s} {'0.715':>8s} {'0.778':>8s} {'0.590':>8s}")
    print(f"{'muscle':>10s} {'0.870':>8s} {'0.721':>8s} {'0.783':>8s} {'0.581':>8s}")
    print(f"{'clustalo':>10s} {'0.840':>8s} {'0.710':>8s} {'0.764':>8s} {'0.559':>8s}")


if __name__ == "__main__":
    sys.exit(main() or 0)
