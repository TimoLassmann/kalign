"""Test PFASUM auto matrix switching at various thresholds.

Compares fixed PFASUM43, fixed PFASUM60, and auto switching
to find the best threshold for the distance-based matrix selection.

We can't easily change the threshold from Python (it's hardcoded in C),
so instead we test: PFASUM43, PFASUM60, and PFASUM_AUTO with the current
threshold. If we want to sweep thresholds, we modify the C code and rebuild.

Usage:
    uv run python -m benchmarks.pfasum_auto_sweep
"""
import statistics
import tempfile
import time
from pathlib import Path

import kalign
from .datasets import get_cases
from .scoring import parse_balibase_xml

CONFIGS = [
    {"name": "PFASUM43", "seq_type": "pfasum43"},
    {"name": "PFASUM60", "seq_type": "pfasum60"},
    {"name": "PFASUM_AUTO", "seq_type": "pfasum"},
]


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(case.reference), str(output_path), column_mask=mask)
    return kalign.compare_detailed(str(case.reference), str(output_path))


def main():
    cases = get_cases("balibase")
    print(f"{len(cases)} cases x {len(CONFIGS)} configs")

    results = []
    t0 = time.perf_counter()
    for ci, config in enumerate(CONFIGS):
        for case in cases:
            with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
                tmp_path = tmp.name
            try:
                kalign.align_file_to_file(
                    str(case.unaligned), tmp_path, format="fasta",
                    seq_type=config["seq_type"],
                    gap_open=7.0, gap_extend=1.5, terminal_gap_extend=1.0,
                    vsm_amax=2.0, refine="confident",
                )
                scores = _score_case(case, tmp_path)
                results.append({
                    "family": case.family, "dataset": case.dataset,
                    "config": config["name"],
                    "recall": scores["recall"], "precision": scores["precision"],
                    "f1": scores["f1"], "tc": scores["tc"],
                })
            except Exception as e:
                results.append({
                    "family": case.family, "dataset": case.dataset,
                    "config": config["name"],
                    "f1": 0, "recall": 0, "precision": 0, "tc": 0,
                    "error": str(e),
                })
            finally:
                Path(tmp_path).unlink(missing_ok=True)
        print(f"  {config['name']} done ({time.perf_counter()-t0:.0f}s)")

    print(f"\nAll done in {time.perf_counter()-t0:.0f}s\n")

    # Overall
    print("=" * 70)
    print("OVERALL (gpo=7.0, gpe=1.5, vsm_amax=2.0, refine=confident)")
    print("=" * 70)
    print(f"{'Config':<15} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 50)
    for config in CONFIGS:
        rows = [r for r in results if r["config"] == config["name"] and "error" not in r]
        if rows:
            rec = statistics.mean(r["recall"] for r in rows)
            prec = statistics.mean(r["precision"] for r in rows)
            f1 = statistics.mean(r["f1"] for r in rows)
            tc = statistics.mean(r["tc"] for r in rows)
            print(f"{config['name']:<15} {rec:8.4f} {prec:8.4f} {f1:8.4f} {tc:8.4f}")

    # Per category
    all_cats = sorted({r["dataset"] for r in results})
    for cat in all_cats:
        cat_label = cat.replace("balibase_", "")
        print(f"\n{cat_label}:")
        print(f"{'Config':<15} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 50)
        for config in CONFIGS:
            rows = [r for r in results if r["config"] == config["name"]
                    and r["dataset"] == cat and "error" not in r]
            if rows:
                rec = statistics.mean(r["recall"] for r in rows)
                prec = statistics.mean(r["precision"] for r in rows)
                f1 = statistics.mean(r["f1"] for r in rows)
                tc = statistics.mean(r["tc"] for r in rows)
                print(f"{config['name']:<15} {rec:8.4f} {prec:8.4f} {f1:8.4f} {tc:8.4f}")

    # Per-family comparison: how many families does AUTO match the best fixed?
    from collections import defaultdict
    by_family = defaultdict(dict)
    for r in results:
        if "error" not in r:
            by_family[r["family"]][r["config"]] = r

    auto_matches_best = 0
    auto_worse = 0
    for fam, configs in by_family.items():
        if "PFASUM_AUTO" in configs and "PFASUM43" in configs and "PFASUM60" in configs:
            best_fixed = max(configs["PFASUM43"]["f1"], configs["PFASUM60"]["f1"])
            auto_f1 = configs["PFASUM_AUTO"]["f1"]
            if abs(auto_f1 - best_fixed) < 0.001:
                auto_matches_best += 1
            elif auto_f1 < best_fixed:
                auto_worse += 1
    total = auto_matches_best + auto_worse + (len(by_family) - auto_matches_best - auto_worse)
    print(f"\nAUTO picks: matches best={auto_matches_best}, worse={auto_worse}, better={total - auto_matches_best - auto_worse}")


if __name__ == "__main__":
    main()
