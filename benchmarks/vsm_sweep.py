"""VSM amax sweep — sequential, low memory.

Usage:
    uv run python -m benchmarks.vsm_sweep
"""
import statistics
import tempfile
import time
from pathlib import Path

import kalign
from .datasets import get_cases
from .scoring import parse_balibase_xml

AMAX_VALUES = [0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(case.reference), str(output_path), column_mask=mask)
    return kalign.compare_detailed(str(case.reference), str(output_path))


def main():
    cases = get_cases("balibase")
    print(f"{len(cases)} cases x {len(AMAX_VALUES)} amax values")

    results = []
    t0 = time.perf_counter()
    for amax in AMAX_VALUES:
        for case in cases:
            with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
                tmp_path = tmp.name
            try:
                kalign.align_file_to_file(
                    str(case.unaligned), tmp_path, format="fasta",
                    seq_type="pfasum43",
                    gap_open=7.0, gap_extend=1.5, terminal_gap_extend=1.0,
                    vsm_amax=amax, refine="confident",
                )
                scores = _score_case(case, tmp_path)
                results.append({"dataset": case.dataset, "amax": amax,
                    "recall": scores["recall"], "precision": scores["precision"],
                    "f1": scores["f1"], "tc": scores["tc"]})
            except Exception:
                pass
            finally:
                Path(tmp_path).unlink(missing_ok=True)
        print(f"  amax={amax} done ({time.perf_counter()-t0:.0f}s)")

    print(f"\nAll done in {time.perf_counter()-t0:.0f}s\n")

    print(f"{'amax':>6} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
    print("-" * 40)
    for amax in AMAX_VALUES:
        rows = [r for r in results if r["amax"] == amax]
        rec = statistics.mean(r["recall"] for r in rows)
        prec = statistics.mean(r["precision"] for r in rows)
        f1 = statistics.mean(r["f1"] for r in rows)
        tc = statistics.mean(r["tc"] for r in rows)
        print(f"{amax:6.1f} {rec:8.4f} {prec:8.4f} {f1:8.4f} {tc:8.4f}")

    all_cats = sorted({r["dataset"] for r in results})
    for cat in all_cats:
        cat_label = cat.replace("balibase_", "")
        print(f"\n{cat_label}:")
        print(f"{'amax':>6} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 40)
        for amax in AMAX_VALUES:
            rows = [r for r in results if r["amax"] == amax and r["dataset"] == cat]
            if rows:
                rec = statistics.mean(r["recall"] for r in rows)
                prec = statistics.mean(r["precision"] for r in rows)
                f1 = statistics.mean(r["f1"] for r in rows)
                tc = statistics.mean(r["tc"] for r in rows)
                print(f"{amax:6.1f} {rec:8.4f} {prec:8.4f} {f1:8.4f} {tc:8.4f}")


if __name__ == "__main__":
    main()
