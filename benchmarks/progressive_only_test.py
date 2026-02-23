#!/usr/bin/env python3
"""Test whether competitors' TC advantage comes from iterative refinement.

Runs:
- mafft --auto (default, likely FFT-NS-2 = 2 tree iterations)
- mafft --retree 1 --maxiterate 0 (FFT-NS-1 = progressive only)
- muscle -align (default)
- clustalo (default)
- kalign (baseline)

If MAFFT's TC drops to kalign-level when refinement is disabled,
then iterative refinement is the key differentiator.
"""

import json
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import kalign

BB_DIR = Path("/kalign/benchmarks/data/downloads/bb3_release")
RESULTS_FILE = Path("/kalign/benchmarks/results/progressive_only_test.json")


def parse_balibase_xml(xml_path: Path) -> list[int]:
    """Parse BAliBASE XML annotation to extract core block column mask."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    colsco = root.find(".//column-score/colsco-data")
    if colsco is None or colsco.text is None:
        raise ValueError(f"No <colsco-data> element found in {xml_path}")
    values = [int(v) for v in colsco.text.split()]
    return [1 if v == 1 else 0 for v in values]


def count_sequences(tfa_path: Path) -> int:
    return sum(1 for line in open(tfa_path) if line.startswith(">"))


def find_cases():
    """Find all BAliBASE cases."""
    cases = []
    for tfa in sorted(BB_DIR.rglob("BB*.tfa")):
        family = tfa.stem
        if family.startswith("BBS"):
            continue
        msf = tfa.with_suffix(".msf")
        xml = tfa.with_suffix(".xml")
        if not msf.exists():
            continue
        rv = "unknown"
        for part in tfa.parts:
            if part.startswith("RV"):
                rv = part
                break
        cases.append({
            "family": family,
            "rv": rv,
            "unaligned": str(tfa),
            "reference": str(msf),
            "xml": str(xml) if xml.exists() else None,
            "nseq": count_sequences(tfa),
        })
    return cases


def run_alignment(method: str, input_fa: str, output_fa: str) -> float:
    """Run alignment, return wall time in seconds."""
    start = time.perf_counter()

    if method == "mafft_auto":
        with open(output_fa, "w") as f:
            r = subprocess.run(
                ["mafft", "--auto", "--quiet", input_fa],
                stdout=f, stderr=subprocess.PIPE, text=True,
            )
        if r.returncode != 0:
            raise RuntimeError(f"mafft --auto failed: {r.stderr[:200]}")

    elif method == "mafft_progressive":
        with open(output_fa, "w") as f:
            r = subprocess.run(
                ["mafft", "--retree", "1", "--maxiterate", "0", "--quiet", input_fa],
                stdout=f, stderr=subprocess.PIPE, text=True,
            )
        if r.returncode != 0:
            raise RuntimeError(f"mafft progressive failed: {r.stderr[:200]}")

    elif method == "muscle_default":
        r = subprocess.run(
            ["muscle", "-align", input_fa, "-output", output_fa],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            raise RuntimeError(f"muscle failed: {r.stderr[:200]}")

    elif method == "clustalo_default":
        r = subprocess.run(
            ["clustalo", "-i", input_fa, "-o", output_fa, "--force"],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            raise RuntimeError(f"clustalo failed: {r.stderr[:200]}")

    elif method == "kalign_baseline":
        r = subprocess.run(
            ["kalign", "-i", input_fa, "-o", output_fa, "--refine", "none"],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            raise RuntimeError(f"kalign failed: {r.stderr[:200]}")

    elif method == "kalign_refine":
        r = subprocess.run(
            ["kalign", "-i", input_fa, "-o", output_fa, "--refine", "confident"],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            raise RuntimeError(f"kalign refine failed: {r.stderr[:200]}")

    else:
        raise ValueError(f"Unknown method: {method}")

    return time.perf_counter() - start


def score_case(test_fa: str, ref_msf: str, xml_path: str | None) -> dict:
    """Score using kalign Python API with BAliBASE XML mask if available."""
    if xml_path and Path(xml_path).exists():
        mask = parse_balibase_xml(Path(xml_path))
        result = kalign.compare_detailed(ref_msf, test_fa, column_mask=mask)
    else:
        result = kalign.compare_detailed(ref_msf, test_fa)
    return {
        "recall": result["recall"],
        "precision": result["precision"],
        "f1": result["f1"],
        "tc": result["tc"],
    }


def main():
    import sys

    cases = find_cases()
    print(f"Found {len(cases)} BAliBASE cases", file=sys.stderr)

    methods = [
        "kalign_baseline",
        "kalign_refine",
        "mafft_auto",
        "mafft_progressive",
        "muscle_default",
        "clustalo_default",
    ]

    all_results = []

    for i, case in enumerate(cases):
        print(
            f"\r[{i+1}/{len(cases)}] {case['family']} (nseq={case['nseq']})...",
            end="", flush=True, file=sys.stderr,
        )

        for method in methods:
            output_fa = None
            try:
                with tempfile.NamedTemporaryFile(
                    suffix=".fa", delete=False, dir="/tmp"
                ) as tmp:
                    output_fa = tmp.name

                wall_time = run_alignment(method, case["unaligned"], output_fa)
                scores = score_case(output_fa, case["reference"], case["xml"])

                all_results.append({
                    "family": case["family"],
                    "rv": case["rv"],
                    "nseq": case["nseq"],
                    "method": method,
                    "wall_time": wall_time,
                    **scores,
                })
            except Exception as e:
                print(
                    f"\n  ERROR {method} on {case['family']}: {e}",
                    file=sys.stderr,
                )
                all_results.append({
                    "family": case["family"],
                    "rv": case["rv"],
                    "nseq": case["nseq"],
                    "method": method,
                    "error": str(e),
                })
            finally:
                if output_fa:
                    Path(output_fa).unlink(missing_ok=True)

    print(file=sys.stderr)

    RESULTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_FILE, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"Saved {len(all_results)} results to {RESULTS_FILE}", file=sys.stderr)

    print_summary(all_results)


def print_summary(all_results):
    methods = sorted(set(r["method"] for r in all_results if "error" not in r))

    by_method = defaultdict(list)
    for r in all_results:
        if "error" not in r:
            by_method[r["method"]].append(r)

    print("\n=== OVERALL SUMMARY ===")
    print(f"{'Method':<25} {'N':>4} {'F1':>7} {'TC':>7} {'Recall':>7} {'Prec':>7} {'Time':>7}")
    print("-" * 75)
    for method in methods:
        results = by_method[method]
        n = len(results)
        f1 = sum(r["f1"] for r in results) / n
        tc = sum(r["tc"] for r in results) / n
        rec = sum(r["recall"] for r in results) / n
        prec = sum(r["precision"] for r in results) / n
        t = sum(r["wall_time"] for r in results)
        print(f"{method:<25} {n:>4} {f1:>7.3f} {tc:>7.3f} {rec:>7.3f} {prec:>7.3f} {t:>7.1f}s")

    def bin_label(n):
        if n <= 5: return "2-5"
        if n <= 10: return "6-10"
        if n <= 20: return "11-20"
        if n <= 50: return "21-50"
        if n <= 100: return "51-100"
        return "101+"

    bin_order = ["2-5", "6-10", "11-20", "21-50", "51-100", "101+"]

    for metric in ["tc", "f1"]:
        print(f"\n=== {metric.upper()} BY SEQUENCE COUNT ===")
        header = f"{'Method':<25}"
        for b in bin_order:
            header += f" {b:>9}"
        header += f" {'ALL':>7}"
        print(header)
        print("-" * 90)

        for method in methods:
            bins = defaultdict(list)
            all_vals = []
            for r in by_method[method]:
                b = bin_label(r["nseq"])
                bins[b].append(r[metric])
                all_vals.append(r[metric])
            row = f"{method:<25}"
            for b in bin_order:
                if bins[b]:
                    mean = sum(bins[b]) / len(bins[b])
                    row += f" {mean:>.3f}({len(bins[b]):>2})"
                else:
                    row += f" {'—':>9}"
            if all_vals:
                row += f" {sum(all_vals)/len(all_vals):>7.3f}"
            print(row)

    rvs = ["RV11", "RV12", "RV20", "RV30", "RV40", "RV50"]
    for metric in ["tc", "f1"]:
        print(f"\n=== {metric.upper()} BY RV CATEGORY ===")
        header = f"{'Method':<25}"
        for rv in rvs:
            header += f" {rv:>9}"
        print(header)
        print("-" * 85)

        for method in methods:
            rv_bins = defaultdict(list)
            for r in by_method[method]:
                rv_bins[r["rv"]].append(r[metric])
            row = f"{method:<25}"
            for rv in rvs:
                if rv_bins[rv]:
                    mean = sum(rv_bins[rv]) / len(rv_bins[rv])
                    row += f" {mean:>.3f}({len(rv_bins[rv]):>2})"
                else:
                    row += f" {'—':>9}"
            print(row)

    # Key comparison
    print("\n=== KEY COMPARISON: MAFFT auto vs progressive-only ===")
    for rv in rvs:
        auto_tc = [r["tc"] for r in by_method.get("mafft_auto", []) if r["rv"] == rv]
        prog_tc = [r["tc"] for r in by_method.get("mafft_progressive", []) if r["rv"] == rv]
        if auto_tc and prog_tc:
            a = sum(auto_tc) / len(auto_tc)
            p = sum(prog_tc) / len(prog_tc)
            print(f"  {rv}: auto_TC={a:.3f}, progressive_TC={p:.3f}, delta={a-p:+.3f}")

    print("\n=== PROGRESSIVE-ONLY RANKING (TC) ===")
    prog_methods = ["kalign_baseline", "mafft_progressive", "muscle_default", "clustalo_default"]
    for rv in rvs + ["ALL"]:
        vals = {}
        for m in prog_methods:
            if rv == "ALL":
                tc_vals = [r["tc"] for r in by_method.get(m, [])]
            else:
                tc_vals = [r["tc"] for r in by_method.get(m, []) if r["rv"] == rv]
            if tc_vals:
                vals[m] = sum(tc_vals) / len(tc_vals)
        if vals:
            ranked = sorted(vals.items(), key=lambda x: -x[1])
            parts = [f"{m}={v:.3f}" for m, v in ranked]
            print(f"  {rv:>5}: {' > '.join(parts)}")


if __name__ == "__main__":
    main()
