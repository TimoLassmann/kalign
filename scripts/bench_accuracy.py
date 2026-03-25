#!/usr/bin/env python3
"""Alignment accuracy benchmark — reuses scoring from the manuscript repo.

Discovers BAliBASE cases, runs kalign fast/default/recall/accurate,
scores with XML core block masks (same as manuscript pipeline).

Usage:
    uv run python scripts/bench_accuracy.py
    uv run python scripts/bench_accuracy.py --modes fast,accurate
"""
import argparse
import json
import logging
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import kalign

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

BB_ROOT = Path(__file__).parent.parent / "benchmarks" / "data" / "downloads" / "bb3_release"


def parse_balibase_xml(xml_path):
    """Parse BAliBASE XML to get core block column mask (from manuscript scoring.py)."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    colsco = root.find(".//column-score/colsco-data")
    if colsco is None or colsco.text is None:
        return None
    values = [int(v) for v in colsco.text.split()]
    return [1 if v == 1 else 0 for v in values]


def discover_cases():
    """Find all BAliBASE cases (excluding BBS supplement)."""
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


def score_case(ref_path, test_path, xml_path):
    """Score one alignment using XML core blocks if available."""
    if xml_path and Path(xml_path).exists():
        mask = parse_balibase_xml(xml_path)
        if mask:
            return kalign.compare_detailed(ref_path, test_path, column_mask=mask)
    return kalign.compare_detailed(ref_path, test_path, max_gap_frac=0.2)


def main():
    parser = argparse.ArgumentParser(description="BAliBASE accuracy benchmark")
    parser.add_argument("--modes", default="fast,default,recall,accurate")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", help="Save results JSON")
    args = parser.parse_args()

    if not BB_ROOT.exists():
        logger.error("BAliBASE not found at %s", BB_ROOT)
        return 1

    modes = args.modes.split(",")
    cases = discover_cases()
    logger.info("%d BAliBASE cases, modes: %s", len(cases), modes)

    results = []
    for mode in modes:
        with tempfile.TemporaryDirectory() as tmpdir:
            for i, case in enumerate(cases):
                out = str(Path(tmpdir) / f"{case['family']}.fa")
                kalign.align_file_to_file(
                    case["unaligned"], out,
                    mode=mode, n_threads=args.threads,
                )
                score = score_case(case["reference"], out, case["xml"])
                results.append({
                    "family": case["family"],
                    "category": case["category"],
                    "mode": mode,
                    **score,
                })
                if (i + 1) % 50 == 0:
                    logger.info("  %s: %d/%d", mode, i + 1, len(cases))
        logger.info("  %s: done", mode)

    # Aggregate and print
    print()
    print(f"{'Mode':12s} {'Recall':>8s} {'Prec':>8s} {'F1':>8s} {'TC':>8s}  (n)")
    print("=" * 52)
    for mode in modes:
        mc = [r for r in results if r["mode"] == mode]
        n = len(mc)
        r = sum(c["recall"] for c in mc) / n
        p = sum(c["precision"] for c in mc) / n
        f = sum(c["f1"] for c in mc) / n
        t = sum(c["tc"] for c in mc) / n
        print(f"{mode:12s} {r:8.3f} {p:8.3f} {f:8.3f} {t:8.3f}  ({n})")

    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        logger.info("Saved to %s", args.output)


if __name__ == "__main__":
    sys.exit(main() or 0)
