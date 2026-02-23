"""Run external tools (mafft, muscle, clustalo) on BAliBASE inside container.

Scores against references using kalign's compare_detailed with XML masks.
Outputs JSON results compatible with probmsa_experiment.py format.

Usage (inside container):
    python -m benchmarks.external_balibase
"""

import json
import shutil
import statistics
import subprocess
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import kalign
from .datasets import get_cases
from .scoring import parse_balibase_xml


def _run_external(tool, unaligned, output):
    """Run an external alignment tool."""
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


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(case.reference), str(output_path), column_mask=mask)
    return kalign.compare_detailed(str(case.reference), str(output_path))


def main():
    tools = []
    for t in ["mafft", "muscle", "clustalo"]:
        if shutil.which(t):
            tools.append(t)
        else:
            print(f"  {t}: not found, skipping")

    if not tools:
        print("No external tools found. Run inside container.")
        return

    print(f"Tools: {tools}")

    cases = get_cases("balibase")
    print(f"{len(cases)} BAliBASE cases")

    results = []
    n_tasks = len(tools) * len(cases)
    done = 0
    t0 = time.perf_counter()

    for tool in tools:
        print(f"\n  Tool: {tool}", flush=True)
        for case in cases:
            with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
                tmp_path = tmp.name
            try:
                t1 = time.perf_counter()
                ok = _run_external(tool, case.unaligned, tmp_path)
                wall = time.perf_counter() - t1
                if ok:
                    scores = _score_case(case, tmp_path)
                    results.append({
                        "family": case.family, "dataset": case.dataset,
                        "method": tool,
                        "recall": scores["recall"], "precision": scores["precision"],
                        "f1": scores["f1"], "tc": scores["tc"],
                        "wall_time": wall,
                    })
                else:
                    results.append({
                        "family": case.family, "dataset": case.dataset,
                        "method": tool,
                        "recall": 0, "precision": 0, "f1": 0, "tc": 0,
                        "wall_time": 0, "error": f"{tool} failed",
                    })
            finally:
                Path(tmp_path).unlink(missing_ok=True)

            done += 1
            if done % 20 == 0:
                elapsed = time.perf_counter() - t0
                eta = elapsed / done * (n_tasks - done)
                print(f"    {done}/{n_tasks} ({elapsed:.0f}s, ETA {eta:.0f}s)", flush=True)

    elapsed = time.perf_counter() - t0
    print(f"\nAll done in {elapsed:.0f}s")

    # Summary
    groups = defaultdict(list)
    for r in results:
        groups[r["method"]].append(r)

    print(f"\n{'Method':>12} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8} {'Time':>7}")
    print("-" * 60)
    for tool in tools:
        entries = [r for r in groups[tool] if "error" not in r]
        if not entries:
            print(f"{tool:>12}  (no results)")
            continue
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        wt = sum(r["wall_time"] for r in entries)
        print(f"{tool:>12} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f} {wt:>6.1f}s")

    # Per-category
    all_cats = sorted({r["dataset"].replace("balibase_", "") for r in results})
    for cat in all_cats:
        cat_results = [r for r in results if cat in r["dataset"]]
        cat_groups = defaultdict(list)
        for r in cat_results:
            cat_groups[r["method"]].append(r)
        n = len(cat_groups.get(tools[0], []))
        print(f"\n=== {cat} ({n} cases) ===")
        print(f"{'Method':>12} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 50)
        for tool in tools:
            entries = [r for r in cat_groups.get(tool, []) if "error" not in r]
            if not entries:
                continue
            rec = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            print(f"{tool:>12} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")

    out = Path("benchmarks/data/external_balibase.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
