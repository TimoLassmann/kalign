"""Pipeline 0: Confidence Calibration.

Validates that kalign ensemble confidence scores are meaningful by
comparing predicted per-column confidence to actual column correctness
against known true alignments from INDELible simulations.

Usage::

    python -m benchmarks.downstream.calibration -j 4
    python -m benchmarks.downstream.calibration -j 4 --quick
"""

from __future__ import annotations

import json
import logging
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-case result
# ---------------------------------------------------------------------------


@dataclass
class CaseResult:
    """Result from a single simulation x method calibration case."""

    sim_id: str
    method: str
    predicted_confidence: list[float]  # per-column confidence
    actual_correct: list[int]  # 1 if column matches true alignment, 0 otherwise
    sp_score: float
    tc_score: float
    wall_time: float
    peak_memory_mb: float


# ---------------------------------------------------------------------------
# Aggregated calibration result
# ---------------------------------------------------------------------------


@dataclass
class CalibrationResult:
    """Aggregated calibration result across all cases and methods."""

    provenance: dict
    cases: list[dict]  # serialised CaseResults
    summary: dict  # aggregated metrics per method


# ---------------------------------------------------------------------------
# Column correctness
# ---------------------------------------------------------------------------


def _column_correctness(
    test_seqs: list[str],
    true_seqs: list[str],
    test_names: list[str],
    true_names: list[str],
) -> list[int]:
    """Check per-column correctness of *test* alignment vs *true* alignment.

    For each column in the test alignment, determine whether the set of
    residue pairs in that column matches the corresponding column in the
    true alignment.  A column is correct (1) if every pair of non-gap
    characters in the test column appears at the same aligned positions in
    the true alignment; otherwise it is incorrect (0).

    Sequences are matched by name between test and true alignments.

    Parameters
    ----------
    test_seqs, true_seqs : list[str]
        Aligned sequences (all same length within each list).
    test_names, true_names : list[str]
        Corresponding sequence names.

    Returns
    -------
    list[int]
        One entry per column in the test alignment: 1 = correct, 0 = incorrect.
    """
    if not test_seqs:
        return []

    # Build name -> index mapping for the true alignment
    true_idx = {name: i for i, name in enumerate(true_names)}

    # Reorder true sequences to match test sequence order
    order = []
    for name in test_names:
        if name not in true_idx:
            raise ValueError(
                f"Sequence {name!r} in test alignment not found in true alignment"
            )
        order.append(true_idx[name])

    reordered_true = [true_seqs[i] for i in order]

    # Build a mapping: for each sequence, residue position -> true column index
    # "residue position" = count of non-gap characters seen so far
    n_seqs = len(test_seqs)
    true_ncols = len(reordered_true[0]) if reordered_true else 0

    # For each sequence, build: residue_index -> true_column
    true_res_to_col: list[dict[int, int]] = []
    for s in range(n_seqs):
        mapping: dict[int, int] = {}
        res_idx = 0
        for col in range(true_ncols):
            c = reordered_true[s][col]
            if c != "-" and c != ".":
                mapping[res_idx] = col
                res_idx += 1
        true_res_to_col.append(mapping)

    # For the test alignment, do the same: track residue positions
    test_ncols = len(test_seqs[0])
    test_res_counters = [0] * n_seqs

    result = []
    for col in range(test_ncols):
        # Collect the true column index for each non-gap character in this test column
        true_cols_for_residues: list[int] = []
        residue_indices_this_col: list[tuple[int, int]] = []  # (seq_idx, res_idx)

        for s in range(n_seqs):
            c = test_seqs[s][col]
            if c != "-" and c != ".":
                res_idx = test_res_counters[s]
                true_col = true_res_to_col[s].get(res_idx)
                if true_col is not None:
                    true_cols_for_residues.append(true_col)
                residue_indices_this_col.append((s, res_idx))
                test_res_counters[s] += 1

        # A column is correct if all non-gap residues in this column
        # are aligned to the same column in the true alignment
        if len(true_cols_for_residues) <= 1:
            # 0 or 1 non-gap characters: trivially correct
            result.append(1)
        elif all(tc == true_cols_for_residues[0] for tc in true_cols_for_residues):
            result.append(1)
        else:
            result.append(0)

    return result


# ---------------------------------------------------------------------------
# Calibration metrics
# ---------------------------------------------------------------------------


def brier_score(predicted: list[float], actual: list[int]) -> float:
    """Brier score: mean squared error between predicted confidence and binary correctness.

    Parameters
    ----------
    predicted : list[float]
        Predicted confidence values in [0, 1].
    actual : list[int]
        Binary correctness values (0 or 1).

    Returns
    -------
    float
        Brier score (lower is better). Returns 0.0 for empty inputs.
    """
    n = len(predicted)
    if n == 0:
        return 0.0
    total = 0.0
    for p, o in zip(predicted, actual):
        total += (p - o) ** 2
    return total / n


def expected_calibration_error(
    predicted: list[float],
    actual: list[int],
    n_bins: int = 10,
) -> float:
    """Expected Calibration Error (ECE).

    Bins predictions into *n_bins* equal-width bins over [0, 1],
    computes |avg_confidence - fraction_correct| per bin, weighted
    by the fraction of samples in each bin.

    Parameters
    ----------
    predicted : list[float]
        Predicted confidence values in [0, 1].
    actual : list[int]
        Binary correctness values (0 or 1).
    n_bins : int
        Number of equal-width bins.

    Returns
    -------
    float
        ECE value (lower is better). Returns 0.0 for empty inputs.
    """
    n = len(predicted)
    if n == 0:
        return 0.0

    bin_sums_conf = [0.0] * n_bins
    bin_sums_correct = [0.0] * n_bins
    bin_counts = [0] * n_bins

    for p, o in zip(predicted, actual):
        # Determine bin index
        b = int(p * n_bins)
        if b >= n_bins:
            b = n_bins - 1
        bin_sums_conf[b] += p
        bin_sums_correct[b] += o
        bin_counts[b] += 1

    ece = 0.0
    for b in range(n_bins):
        if bin_counts[b] > 0:
            avg_conf = bin_sums_conf[b] / bin_counts[b]
            frac_correct = bin_sums_correct[b] / bin_counts[b]
            ece += (bin_counts[b] / n) * abs(avg_conf - frac_correct)

    return ece


def calibration_curve(
    predicted: list[float],
    actual: list[int],
    n_bins: int = 10,
) -> tuple[list[float], list[float], list[int]]:
    """Compute calibration curve data for plotting.

    Parameters
    ----------
    predicted : list[float]
        Predicted confidence values in [0, 1].
    actual : list[int]
        Binary correctness values (0 or 1).
    n_bins : int
        Number of equal-width bins over [0, 1].

    Returns
    -------
    bin_centers : list[float]
        Midpoint of each bin.
    fraction_correct : list[float]
        Fraction of correct predictions in each bin (NaN for empty bins).
    bin_counts : list[int]
        Number of samples in each bin.
    """
    bin_width = 1.0 / n_bins
    bin_centers = [(b + 0.5) * bin_width for b in range(n_bins)]
    bin_sums_correct = [0.0] * n_bins
    bin_counts = [0] * n_bins

    for p, o in zip(predicted, actual):
        b = int(p * n_bins)
        if b >= n_bins:
            b = n_bins - 1
        bin_sums_correct[b] += o
        bin_counts[b] += 1

    fraction_correct = []
    for b in range(n_bins):
        if bin_counts[b] > 0:
            fraction_correct.append(bin_sums_correct[b] / bin_counts[b])
        else:
            fraction_correct.append(float("nan"))

    return bin_centers, fraction_correct, bin_counts


# ---------------------------------------------------------------------------
# Single-case runner
# ---------------------------------------------------------------------------


def run_calibration_case(
    sim_dataset,
    method_name: str,
    work_dir: Path,
    seq_type: str = "protein",
) -> CaseResult:
    """Run one calibration case: align with *method_name*, evaluate column correctness.

    Parameters
    ----------
    sim_dataset : SimulatedDataset
        A simulated dataset with known true alignment.
    method_name : str
        Key into ``utils.METHODS``.
    work_dir : Path
        Scratch directory for alignment output files.
    seq_type : str
        Sequence type hint for the aligner.

    Returns
    -------
    CaseResult
    """
    from .utils import alignment_accuracy, parse_fasta, run_method

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    sim_id = sim_dataset.params.get("sim_id", "unknown")
    if hasattr(sim_dataset, "params") and isinstance(sim_dataset.params, dict):
        sim_id = sim_dataset.params.get("sim_id", sim_id)

    # "true" method: use the true alignment directly as a performance ceiling
    if method_name == "true":
        true_names, true_seqs = parse_fasta(sim_dataset.true_alignment)
        aln_len = len(true_seqs[0]) if true_seqs else 0
        return CaseResult(
            sim_id=sim_id,
            method="true",
            predicted_confidence=[1.0] * aln_len,
            actual_correct=[1] * aln_len,
            sp_score=1.0,
            tc_score=1.0,
            wall_time=0.0,
            peak_memory_mb=0.0,
        )

    # Run the alignment method
    aln_result = run_method(method_name, sim_dataset.unaligned, work_dir, seq_type=seq_type)

    # Load the true alignment
    true_names, true_seqs = parse_fasta(sim_dataset.true_alignment)

    # Compute alignment accuracy (SP, TC)
    acc = alignment_accuracy(
        aln_result.sequences,
        true_seqs,
        aln_result.names,
        true_names,
    )

    # Compute per-column correctness
    correctness = _column_correctness(
        aln_result.sequences,
        true_seqs,
        aln_result.names,
        true_names,
    )

    # Extract confidence: use column_confidence if available, else uniform 1.0
    if aln_result.column_confidence is not None:
        confidence = list(aln_result.column_confidence)
        # After masking, the confidence list may be longer than the alignment;
        # truncate or pad to match the alignment length
        aln_len = len(aln_result.sequences[0]) if aln_result.sequences else 0
        if len(confidence) > aln_len:
            confidence = confidence[:aln_len]
        elif len(confidence) < aln_len:
            confidence.extend([1.0] * (aln_len - len(confidence)))
    else:
        # Methods without confidence get uniform 1.0
        aln_len = len(aln_result.sequences[0]) if aln_result.sequences else 0
        confidence = [1.0] * aln_len

    return CaseResult(
        sim_id=sim_id,
        method=method_name,
        predicted_confidence=confidence,
        actual_correct=correctness,
        sp_score=acc["sp_score"],
        tc_score=acc["tc_score"],
        wall_time=aln_result.wall_time,
        peak_memory_mb=aln_result.peak_memory_mb,
    )


def _run_one_calibration(sim_dataset, method_name, work_dir, fingerprint=None):
    """Worker for parallel calibration cases.  Returns a dict."""
    from dataclasses import asdict as _asdict

    from .utils import cache_load, cache_save, clean_work_dir

    work_dir = Path(work_dir)
    if fingerprint:
        cached = cache_load(work_dir, fingerprint)
        if cached is not None:
            return cached

    # Remove stale output files from previous runs
    clean_work_dir(work_dir)

    case = run_calibration_case(sim_dataset, method_name, work_dir, seq_type="protein")
    result_dict = _asdict(case)

    if fingerprint:
        cache_save(work_dir, fingerprint, result_dict)

    return result_dict


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------


def run_pipeline(params: dict) -> CalibrationResult:
    """Run the full confidence calibration pipeline.

    Parameters
    ----------
    params : dict
        Configuration with keys:

        - ``data_dir`` : str or Path -- base directory for simulation data
        - ``results_dir`` : str or Path -- base directory for result JSON files
        - ``methods`` : list[str] -- method names to evaluate (keys in METHODS)
        - ``n_simulations`` : int -- override number of simulation cases
          (default: use PROTEIN_GRID)
        - ``n_jobs`` : int -- parallelism level (default 1)
        - ``quick`` : bool -- if True, only 5 cases for smoke testing

    Returns
    -------
    CalibrationResult
    """
    from dataclasses import asdict as _asdict

    from .provenance import collect_provenance, result_path, update_latest_symlink
    from .simulation import (
        PROTEIN_GRID,
        PROTEIN_GRID_FULL,
        generate_indelible_dataset,
        iter_simulation_params,
        random_birth_death_tree,
    )
    from .utils import METHODS, cache_load, cache_save, clean_work_dir, tool_versions_fingerprint

    data_dir = Path(params.get("data_dir", "benchmarks/data/downstream/calibration"))
    results_dir = Path(params.get("results_dir", "benchmarks/results"))
    # Calibration evaluates raw confidence scores → masked methods are excluded
    _default_methods = ["kalign", "kalign_cons", "kalign_ens3", "true"]
    methods = params.get("methods", _default_methods)
    quick = params.get("quick", False)
    n_jobs = params.get("n_jobs", 1)

    data_dir.mkdir(parents=True, exist_ok=True)

    # Validate method names
    for m in methods:
        if m not in METHODS:
            raise ValueError(f"Unknown method {m!r}. Available: {sorted(METHODS.keys())}")

    # Generate simulation parameter sets
    grid = PROTEIN_GRID_FULL if params.get("full") else PROTEIN_GRID
    sim_params_list = list(iter_simulation_params(grid, "WAG"))

    # Allow override of the number of simulations
    n_simulations = params.get("n_simulations")
    if n_simulations is not None:
        sim_params_list = sim_params_list[:n_simulations]

    if quick:
        sim_params_list = sim_params_list[:5]

    logger.info(
        "Running calibration pipeline: %d simulations x %d methods = %d cases",
        len(sim_params_list),
        len(methods),
        len(sim_params_list) * len(methods),
    )

    # Step 1: Generate simulated datasets (sequential -- INDELible is fast)
    from tqdm import tqdm

    sim_datasets = []
    for sp in tqdm(sim_params_list, desc="Generating simulations", unit="sim"):
        sim_id = sp["sim_id"]
        sim_dir = data_dir / sim_id

        tree = random_birth_death_tree(
            n_taxa=sp["n_taxa"],
            target_depth=sp["target_depth"],
            seed=sp.get("seed", 42),
        )

        sim_dataset = generate_indelible_dataset(
            tree=tree,
            model=sp["model"],
            seq_length=sp["seq_length"],
            indel_rate=sp["indel_rate"],
            indel_length_mean=sp.get("indel_length_mean", 2.0),
            output_dir=sim_dir,
            seed=sp.get("seed", 42),
        )
        sim_dataset.params["sim_id"] = sim_id
        sim_datasets.append(sim_dataset)

    # Step 2: Build work items for alignment + scoring
    work_items = []
    for sim_dataset in sim_datasets:
        sim_id = sim_dataset.params["sim_id"]
        sim_dir = data_dir / sim_id
        for method_name in methods:
            work_dir = sim_dir / f"work_{method_name}"
            work_items.append((sim_dataset, method_name, work_dir))

    # Compute fingerprint once for all workers (None disables caching)
    if params.get("no_cache"):
        fingerprint = None
        logger.info("Caching disabled (--no-cache)")
    else:
        fingerprint = tool_versions_fingerprint()
        logger.info("Tool versions fingerprint: %s", fingerprint)

    # Step 3: Run alignment cases (parallel when n_jobs > 1)
    all_cases: list[dict] = []
    n_cached = 0
    pbar = tqdm(total=len(work_items), desc="Calibration", unit="case")
    if n_jobs <= 1:
        for sim_dataset, method_name, work_dir in work_items:
            try:
                cached = cache_load(work_dir, fingerprint) if fingerprint else None
                if cached is not None:
                    all_cases.append(cached)
                    n_cached += 1
                else:
                    clean_work_dir(work_dir)
                    case = run_calibration_case(
                        sim_dataset, method_name, work_dir, seq_type="protein",
                    )
                    result_dict = _asdict(case)
                    if fingerprint:
                        cache_save(work_dir, fingerprint, result_dict)
                    all_cases.append(result_dict)
            except Exception as exc:
                logger.error("Failed: %s / %s: %s",
                             sim_dataset.params["sim_id"], method_name, exc)
                all_cases.append({
                    "sim_id": sim_dataset.params["sim_id"],
                    "method": method_name,
                    "error": str(exc),
                })
            pbar.update(1)
            pbar.set_postfix(cached=n_cached, computed=pbar.n - n_cached)
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {}
            for sim_dataset, method_name, work_dir in work_items:
                fut = pool.submit(
                    _run_one_calibration,
                    sim_dataset, method_name, work_dir,
                    fingerprint=fingerprint,
                )
                futures[fut] = (sim_dataset.params["sim_id"], method_name)

            for fut in as_completed(futures):
                sim_id, method_name = futures[fut]
                try:
                    result_dict = fut.result()
                    all_cases.append(result_dict)
                except Exception as exc:
                    logger.error("Failed: %s / %s: %s", sim_id, method_name, exc)
                    all_cases.append({
                        "sim_id": sim_id,
                        "method": method_name,
                        "error": str(exc),
                    })
                pbar.update(1)
    pbar.close()

    # Aggregate results per method
    summary = _aggregate_summary(all_cases)

    # Collect provenance
    provenance = collect_provenance(params)

    result = CalibrationResult(
        provenance=_asdict(provenance),
        cases=all_cases,
        summary=summary,
    )

    # Save results
    out_path = result_path(results_dir, "calibration")
    with open(out_path, "w") as fh:
        json.dump(
            {
                "provenance": result.provenance,
                "cases": result.cases,
                "summary": result.summary,
            },
            fh,
            indent=2,
            default=str,
        )
    update_latest_symlink(out_path)
    logger.info("Results saved to %s", out_path)

    return result


def _aggregate_summary(cases: list[dict]) -> dict:
    """Compute per-method summary statistics from case result dicts.

    Returns a dict keyed by method name with aggregated metrics.
    """
    from collections import defaultdict

    method_preds: dict[str, list[float]] = defaultdict(list)
    method_actuals: dict[str, list[int]] = defaultdict(list)
    method_sp: dict[str, list[float]] = defaultdict(list)
    method_tc: dict[str, list[float]] = defaultdict(list)
    method_time: dict[str, list[float]] = defaultdict(list)
    method_mem: dict[str, list[float]] = defaultdict(list)

    for case in cases:
        if "error" in case:
            continue
        m = case["method"]
        method_preds[m].extend(case.get("predicted_confidence", []))
        method_actuals[m].extend(case.get("actual_correct", []))
        if case.get("sp_score", -1.0) >= 0:
            method_sp[m].append(case["sp_score"])
        if case.get("tc_score", -1.0) >= 0:
            method_tc[m].append(case["tc_score"])
        method_time[m].append(case.get("wall_time", 0.0))
        method_mem[m].append(case.get("peak_memory_mb", 0.0))

    summary: dict[str, dict] = {}
    for m in sorted(method_preds.keys()):
        preds = method_preds[m]
        actuals = method_actuals[m]
        bs = brier_score(preds, actuals)
        ece = expected_calibration_error(preds, actuals)
        bin_centers, frac_correct, bin_counts = calibration_curve(preds, actuals)

        sp_vals = method_sp[m]
        tc_vals = method_tc[m]
        time_vals = method_time[m]
        mem_vals = method_mem[m]

        summary[m] = {
            "brier_score": bs,
            "ece": ece,
            "calibration_curve": {
                "bin_centers": bin_centers,
                "fraction_correct": frac_correct,
                "bin_counts": bin_counts,
            },
            "n_columns_total": len(preds),
            "mean_sp": sum(sp_vals) / len(sp_vals) if sp_vals else None,
            "mean_tc": sum(tc_vals) / len(tc_vals) if tc_vals else None,
            "mean_time": sum(time_vals) / len(time_vals) if time_vals else None,
            "mean_peak_memory_mb": sum(mem_vals) / len(mem_vals) if mem_vals else None,
            "n_cases": len(time_vals),
        }

    return summary


# ---------------------------------------------------------------------------
# Result loading
# ---------------------------------------------------------------------------


def load_results(results_dir: Path | str) -> CalibrationResult:
    """Load calibration results from the latest.json symlink.

    Parameters
    ----------
    results_dir : Path or str
        Base results directory (the ``calibration/`` subdirectory is
        appended automatically).

    Returns
    -------
    CalibrationResult
    """
    results_dir = Path(results_dir)
    cal_dir = results_dir / "calibration"
    latest = cal_dir / "latest.json"
    with open(latest) as fh:
        data = json.load(fh)
    return CalibrationResult(
        provenance=data.get("provenance", {}),
        cases=data.get("cases", []),
        summary=data.get("summary", {}),
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import sys

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Pipeline 0: Confidence calibration benchmarks"
    )
    parser.add_argument(
        "-j", "--parallel", type=int, default=1, help="Number of parallel jobs"
    )
    parser.add_argument(
        "--quick", action="store_true", help="Quick mode: 5 cases only"
    )
    parser.add_argument(
        "--max-sims", type=int, default=None, help="Limit number of simulations"
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="benchmarks/data/downstream/calibration",
        help="Data directory for simulations",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="benchmarks/results",
        help="Results output directory",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=None,
        help="Methods to run (default: all)",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Ignore cached results and recompute everything.",
    )

    args = parser.parse_args()

    params = {
        "data_dir": args.data_dir,
        "results_dir": args.results_dir,
        "n_jobs": args.parallel,
        "quick": args.quick,
        "no_cache": args.no_cache,
    }
    if args.max_sims is not None:
        params["n_simulations"] = args.max_sims
    if args.methods is not None:
        params["methods"] = args.methods

    result = run_pipeline(params)

    # Print summary
    print("\n=== Calibration Summary ===\n")
    for method, stats in result.summary.items():
        print(f"  {method:25s}  Brier={stats['brier_score']:.4f}  "
              f"ECE={stats['ece']:.4f}  "
              f"SP={stats.get('mean_sp', 'N/A')!s:>6s}  "
              f"TC={stats.get('mean_tc', 'N/A')!s:>6s}  "
              f"n={stats['n_cases']}")
    print()
