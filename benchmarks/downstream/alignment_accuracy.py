"""Pipeline: Alignment Accuracy (BAliBASE + BRAliBASE).

Scores alignment quality for all 6 methods against reference alignments.
BAliBASE cases use XML core-block masks; BRAliBASE uses all columns.

Usage::

    python -m benchmarks.downstream alignment_accuracy -j 8
    python -m benchmarks.downstream alignment_accuracy --quick -j 4
"""

from __future__ import annotations

import json
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

logger = logging.getLogger(__name__)

# Default methods for alignment accuracy (all 6 + true reference)
_DEFAULT_METHODS = [
    "kalign",
    "kalign_cons",
    "kalign_ens3",
    "mafft",
    "muscle",
    "clustalo",
]


# ---------------------------------------------------------------------------
# Single-case runner
# ---------------------------------------------------------------------------


def _run_case(
    case_dict: dict,
    method_name: str,
    work_dir: str,
    fingerprint: str | None = None,
) -> dict:
    """Align one case with one method and score against reference.

    Parameters
    ----------
    case_dict : dict
        Serialised BenchmarkCase with keys: family, dataset, unaligned,
        reference, seq_type.
    method_name : str
        Key into ``utils.METHODS``.
    work_dir : str
        Scratch directory.
    fingerprint : str or None
        Tool versions fingerprint for caching (None disables cache).

    Returns
    -------
    dict
        Per-case result with scoring metrics, timing, and metadata.
    """
    import kalign as _kalign

    from ..scoring import parse_balibase_xml
    from .utils import (
        cache_load,
        cache_save,
        clean_work_dir,
        parse_fasta,
        run_method,
        write_fasta,
    )

    wd = Path(work_dir)
    family = case_dict["family"]
    dataset = case_dict["dataset"]
    seq_type = case_dict["seq_type"]
    unaligned = Path(case_dict["unaligned"])
    reference = Path(case_dict["reference"])

    # Check cache
    if fingerprint:
        cached = cache_load(wd, fingerprint)
        if cached is not None:
            return cached

    clean_work_dir(wd)
    wd.mkdir(parents=True, exist_ok=True)

    try:
        # Align
        aln_result = run_method(
            method_name,
            unaligned,
            wd,
            seq_type=seq_type,
        )

        # Write alignment for scoring
        aln_path = wd / "aligned.fasta"
        write_fasta(aln_result.names, aln_result.sequences, aln_path)

        # Score: BAliBASE uses XML core-block mask, BRAliBASE uses all columns
        xml_path = reference.with_suffix(".xml")
        if xml_path.exists():
            mask = parse_balibase_xml(xml_path)
            scores = _kalign.compare_detailed(
                str(reference), str(aln_path), column_mask=mask
            )
        else:
            # Check for gapless reference (conserved RNA) — skip
            ref_names, ref_seqs = parse_fasta(reference)
            if ref_seqs and all("-" not in s for s in ref_seqs):
                return {
                    "family": family,
                    "dataset": dataset,
                    "benchmark": "bralibase" if "brali" in dataset else "balibase",
                    "method": method_name,
                    "skipped": True,
                    "reason": "gapless_reference",
                }

            scores = _kalign.compare_detailed(
                str(reference), str(aln_path), max_gap_frac=-1.0
            )

        result = {
            "family": family,
            "dataset": dataset,
            "benchmark": "bralibase" if "brali" in dataset else "balibase",
            "method": method_name,
            "recall": scores.get("recall", -1.0),
            "precision": scores.get("precision", -1.0),
            "f1": scores.get("f1", -1.0),
            "tc": scores.get("tc", -1.0),
            "wall_time": aln_result.wall_time,
            "peak_memory_mb": aln_result.peak_memory_mb,
        }

        if fingerprint:
            cache_save(wd, fingerprint, result)

        return result

    except Exception as exc:
        logger.warning("Case %s/%s failed: %s", family, method_name, exc)
        return {
            "family": family,
            "dataset": dataset,
            "benchmark": "bralibase" if "brali" in dataset else "balibase",
            "method": method_name,
            "error": str(exc),
        }


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------


def run_pipeline(params: dict) -> dict:
    """Run the alignment accuracy pipeline.

    Parameters
    ----------
    params : dict
        Configuration with keys:

        - ``data_dir`` : str or Path -- base data directory
        - ``results_dir`` : str or Path -- results output directory
        - ``methods`` : list[str] -- method names to evaluate
        - ``n_jobs`` : int -- parallel workers
        - ``quick`` : bool -- if True, limit to 5 cases per dataset
        - ``no_cache`` : bool -- disable caching

    Returns
    -------
    dict
        Result dict with provenance, cases, and summary.
    """
    from ..datasets import get_cases
    from .provenance import collect_provenance, result_path, update_latest_symlink
    from .utils import METHODS, tool_versions_fingerprint

    data_dir = Path(params.get("data_dir", "benchmarks/data"))
    results_dir = Path(params.get("results_dir", "benchmarks/results"))
    methods = params.get("methods", _DEFAULT_METHODS)
    n_jobs = params.get("n_jobs", 1)
    quick = params.get("quick", False)

    # Validate methods
    for m in methods:
        if m not in METHODS:
            raise ValueError(f"Unknown method {m!r}. Available: {sorted(METHODS.keys())}")

    # Load BAliBASE + BRAliBASE cases
    balibase_cases = get_cases("balibase")
    bralibase_cases = get_cases("bralibase")
    all_cases = balibase_cases + bralibase_cases

    if quick:
        # Take first 5 from each dataset
        balibase_subset = balibase_cases[:5]
        bralibase_subset = bralibase_cases[:5]
        all_cases = balibase_subset + bralibase_subset

    logger.info(
        "Alignment accuracy: %d cases x %d methods = %d total",
        len(all_cases),
        len(methods),
        len(all_cases) * len(methods),
    )

    # Compute fingerprint
    if params.get("no_cache"):
        fingerprint = None
        logger.info("Caching disabled (--no-cache)")
    else:
        fingerprint = tool_versions_fingerprint()
        logger.info("Tool versions fingerprint: %s", fingerprint)

    # Build work items
    work_items = []
    for case in all_cases:
        for method_name in methods:
            case_dict = {
                "family": case.family,
                "dataset": case.dataset,
                "unaligned": str(case.unaligned),
                "reference": str(case.reference),
                "seq_type": case.seq_type,
            }
            work_dir = str(
                data_dir / "alignment_accuracy_work" / case.dataset / case.family / method_name
            )
            work_items.append((case_dict, method_name, work_dir, fingerprint))

    # Execute
    from tqdm import tqdm

    results: list[dict] = []
    pbar = tqdm(total=len(work_items), desc="Alignment accuracy", unit="case")

    if n_jobs <= 1:
        for case_dict, method_name, work_dir, fp in work_items:
            r = _run_case(case_dict, method_name, work_dir, fp)
            results.append(r)
            pbar.update(1)
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {
                pool.submit(_run_case, cd, mn, wd, fp): (cd, mn)
                for cd, mn, wd, fp in work_items
            }
            for fut in as_completed(futures):
                try:
                    results.append(fut.result())
                except Exception as exc:
                    cd, mn = futures[fut]
                    logger.error("Worker failed: %s/%s: %s", cd["family"], mn, exc)
                    results.append({
                        "family": cd["family"],
                        "dataset": cd["dataset"],
                        "method": mn,
                        "error": str(exc),
                    })
                pbar.update(1)
    pbar.close()

    # Filter errors and skipped
    successful = [r for r in results if "error" not in r and not r.get("skipped")]
    skipped = [r for r in results if r.get("skipped")]
    errors = [r for r in results if "error" in r]
    logger.info(
        "%d successful, %d skipped, %d errors",
        len(successful), len(skipped), len(errors),
    )

    # Aggregate summary
    summary = _aggregate_summary(successful, methods)

    # Provenance
    provenance = collect_provenance(params)

    output = {
        "provenance": asdict(provenance),
        "cases": results,
        "summary": summary,
    }

    # Save
    out_path = result_path(results_dir, "alignment_accuracy")
    with open(out_path, "w") as fh:
        json.dump(output, fh, indent=2, default=str)
    update_latest_symlink(out_path)
    logger.info("Results saved to %s", out_path)

    # Print summary
    _print_summary(summary, methods)

    return output


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def _aggregate_summary(cases: list[dict], methods: list[str]) -> dict:
    """Compute per-method and per-category aggregated metrics."""
    from collections import defaultdict

    import numpy as np

    summary: dict = {}

    for method in methods:
        mc = [c for c in cases if c.get("method") == method]
        if not mc:
            summary[method] = {"n_cases": 0}
            continue

        recall_vals = [c["recall"] for c in mc if c.get("recall", -1) >= 0]
        prec_vals = [c["precision"] for c in mc if c.get("precision", -1) >= 0]
        f1_vals = [c["f1"] for c in mc if c.get("f1", -1) >= 0]
        tc_vals = [c["tc"] for c in mc if c.get("tc", -1) >= 0]
        time_vals = [c["wall_time"] for c in mc if "wall_time" in c]

        summary[method] = {
            "n_cases": len(mc),
            "mean_recall": float(np.mean(recall_vals)) if recall_vals else None,
            "mean_precision": float(np.mean(prec_vals)) if prec_vals else None,
            "mean_f1": float(np.mean(f1_vals)) if f1_vals else None,
            "mean_tc": float(np.mean(tc_vals)) if tc_vals else None,
            "mean_wall_time": float(np.mean(time_vals)) if time_vals else None,
        }

        # Per-category breakdown
        by_dataset: dict[str, list[dict]] = defaultdict(list)
        for c in mc:
            by_dataset[c["dataset"]].append(c)

        categories = {}
        for ds, ds_cases in sorted(by_dataset.items()):
            ds_recall = [c["recall"] for c in ds_cases if c.get("recall", -1) >= 0]
            ds_f1 = [c["f1"] for c in ds_cases if c.get("f1", -1) >= 0]
            ds_tc = [c["tc"] for c in ds_cases if c.get("tc", -1) >= 0]
            categories[ds] = {
                "n_cases": len(ds_cases),
                "mean_recall": float(np.mean(ds_recall)) if ds_recall else None,
                "mean_f1": float(np.mean(ds_f1)) if ds_f1 else None,
                "mean_tc": float(np.mean(ds_tc)) if ds_tc else None,
            }
        summary[method]["categories"] = categories

    return summary


def _print_summary(summary: dict, methods: list[str]) -> None:
    """Print a summary table to the logger."""
    logger.info("=" * 80)
    logger.info("Alignment Accuracy Summary")
    logger.info("=" * 80)
    logger.info(
        "%-18s %5s %8s %8s %8s %8s %8s",
        "Method", "N", "Recall", "Prec", "F1", "TC", "Time",
    )
    logger.info("-" * 80)
    for m in methods:
        s = summary.get(m, {})
        n = s.get("n_cases", 0)
        if n == 0:
            logger.info("%-18s %5d %8s", m, 0, "no data")
            continue
        logger.info(
            "%-18s %5d %8.3f %8.3f %8.3f %8.3f %7.2fs",
            m,
            n,
            s.get("mean_recall") or 0,
            s.get("mean_precision") or 0,
            s.get("mean_f1") or 0,
            s.get("mean_tc") or 0,
            s.get("mean_wall_time") or 0,
        )
    logger.info("=" * 80)
