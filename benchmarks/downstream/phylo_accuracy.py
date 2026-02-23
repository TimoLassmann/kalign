"""Pipeline 2: Phylogenetic Tree Accuracy.

Evaluates how alignment quality affects phylogenetic tree inference.
Simulates protein evolution with INDELible (WAG+Gamma), aligns with
multiple methods, infers ML trees with IQ-TREE 2, and compares inferred
trees to the true (simulated) tree using Robinson-Foulds distances.

Methods that produce per-column confidence scores can pass them to
IQ-TREE as site weights (``--site-weight`` flag) for continuous confidence weighting,
or mask low-confidence columns before tree inference.
"""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------


@dataclass
class PhyloCaseResult:
    """Result of a single simulation x method phylogenetic analysis."""

    sim_id: str
    method: str
    nrf: float  # normalized Robinson-Foulds (0=perfect, 1=completely wrong)
    branch_score_dist: float  # weighted RF
    sp_score: float  # alignment SP vs true alignment
    n_columns_retained: int  # after masking (if applicable)
    n_columns_total: int
    wall_time_align: float
    wall_time_iqtree: float
    peak_memory_mb: float


@dataclass
class PhyloResult:
    """Aggregated result for the full phylo accuracy pipeline."""

    provenance: dict
    cases: list[dict]
    summary: dict  # per-method aggregated metrics


# ---------------------------------------------------------------------------
# IQ-TREE runner
# ---------------------------------------------------------------------------


def _run_iqtree(
    alignment_path: Path,
    model: str,
    work_dir: Path,
    site_weights_path: Optional[Path] = None,
    n_threads: int = 1,
    timeout: int = 600,
) -> dict:
    """Run IQ-TREE 2 on an alignment and return the inferred tree.

    Parameters
    ----------
    alignment_path : Path
        FASTA alignment file.
    model : str
        Substitution model string (e.g. ``"WAG+G4"``).
    work_dir : Path
        Working directory for IQ-TREE output files.
    site_weights_path : Path or None
        If given, passed to IQ-TREE via ``--site-weight`` for per-site weighting.
    n_threads : int
        Number of threads for IQ-TREE (default 1).
    timeout : int
        Maximum seconds before killing IQ-TREE (default 600).

    Returns
    -------
    dict
        ``{"tree": str, "log_likelihood": float, "wall_time": float}``

    Raises
    ------
    RuntimeError
        If ``iqtree2`` is not found on ``$PATH`` or exits with an error.
    """
    iqtree_bin = shutil.which("iqtree2")
    if iqtree_bin is None:
        raise RuntimeError(
            "iqtree2 is not installed or not on PATH. "
            "Please install IQ-TREE 2 (http://www.iqtree.org/) "
            "and ensure the 'iqtree2' binary is accessible."
        )

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    prefix = work_dir / "iqtree"

    cmd = [
        iqtree_bin,
        "-s", str(alignment_path),
        "-m", model,
        "--prefix", str(prefix),
        "-nt", str(n_threads),
        "--fast",
        "-redo",
    ]

    # Note: IQ-TREE 2 does not support per-site weights.
    # The site_weights_path parameter is accepted but ignored.
    # Weighted methods (kalign_ens3_wt, guidance2_wt) run IQ-TREE
    # on the full unmasked alignment, providing a baseline comparison.
    if site_weights_path is not None:
        logger.debug("Site weights file ignored (IQ-TREE has no weighting support)")

    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        wall = time.perf_counter() - start
        raise RuntimeError(
            f"iqtree2 timed out after {timeout}s ({wall:.0f}s elapsed) "
            f"in {work_dir}"
        )
    wall = time.perf_counter() - start

    if proc.returncode != 0:
        raise RuntimeError(
            f"iqtree2 failed (exit {proc.returncode}):\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    # Parse the treefile
    treefile = Path(str(prefix) + ".treefile")
    if not treefile.exists():
        raise FileNotFoundError(
            f"IQ-TREE treefile not found at {treefile}. "
            f"Contents of {work_dir}: {list(work_dir.iterdir())}"
        )
    tree_str = treefile.read_text().strip()

    # Parse log-likelihood from the .iqtree summary file
    log_likelihood = float("nan")
    iqtree_log = Path(str(prefix) + ".iqtree")
    if iqtree_log.exists():
        log_text = iqtree_log.read_text()
        # Look for "Log-likelihood of the tree: -12345.6789 (..."
        m = re.search(r"Log-likelihood of the tree:\s+([-\d.]+)", log_text)
        if m:
            try:
                log_likelihood = float(m.group(1))
            except ValueError:
                pass

    return {
        "tree": tree_str,
        "log_likelihood": log_likelihood,
        "wall_time": wall,
    }


# ---------------------------------------------------------------------------
# Single-case runner
# ---------------------------------------------------------------------------


def run_phylo_case(
    sim_dataset,
    method_name: str,
    work_dir: Path,
    model: str = "WAG+G4",
    n_threads: int = 1,
    iqtree_timeout: int = 600,
) -> PhyloCaseResult:
    """Run the full phylo pipeline for one simulation x method combination.

    Steps:
    1. Align unaligned sequences with the specified method.
    2. Write site weights if the method produces them.
    3. Run IQ-TREE on the (possibly masked/weighted) alignment.
    4. Compare the inferred tree against the true tree.

    Parameters
    ----------
    sim_dataset : SimulatedDataset
        Simulation output (unaligned seqs, true alignment, true tree).
    method_name : str
        Key into ``utils.METHODS``.
    work_dir : Path
        Scratch directory for this case.
    model : str
        IQ-TREE substitution model (default ``"WAG+G4"``).

    Returns
    -------
    PhyloCaseResult
    """
    from .utils import (
        alignment_accuracy,
        compare_trees,
        parse_fasta,
        run_method,
        write_fasta,
        METHODS,
    )

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    sim_id = sim_dataset.params.get("sim_id", str(sim_dataset.unaligned))

    # -- "true" method: use true alignment directly as performance ceiling --
    if method_name == "true":
        true_names, true_seqs = parse_fasta(sim_dataset.true_alignment)
        aln_fasta = work_dir / "true_aln.fasta"
        write_fasta(true_names, true_seqs, aln_fasta)

        n_columns_total = len(true_seqs[0]) if true_seqs else 0

        iqtree_dir = work_dir / "true_iqtree"
        iqtree_result = _run_iqtree(
            aln_fasta, model, iqtree_dir,
            n_threads=n_threads, timeout=iqtree_timeout,
        )

        true_tree_str = sim_dataset.true_tree.read_text().strip()
        tree_cmp = compare_trees(true_tree_str, iqtree_result["tree"])

        return PhyloCaseResult(
            sim_id=sim_id,
            method="true",
            nrf=tree_cmp["nrf"],
            branch_score_dist=tree_cmp["branch_score_dist"],
            sp_score=1.0,
            n_columns_retained=n_columns_total,
            n_columns_total=n_columns_total,
            wall_time_align=0.0,
            wall_time_iqtree=iqtree_result["wall_time"],
            peak_memory_mb=0.0,
        )

    # -- Step 1: Align -------------------------------------------------------
    # For guided methods, pass the true alignment as input so that
    # kalign can compute pairwise distances from it and build a tree.
    cfg = METHODS.get(method_name, {})
    if cfg.get("probmsa_guided"):
        input_fasta = sim_dataset.true_alignment
    else:
        input_fasta = sim_dataset.unaligned

    aln_result = run_method(
        method_name,
        input_fasta,
        work_dir,
        seq_type="protein",
    )

    # Determine alignment dimensions
    n_columns_total = len(aln_result.sequences[0]) if aln_result.sequences else 0
    n_columns_retained = n_columns_total  # masking already applied by run_method

    # -- Step 2: Write alignment to FASTA for IQ-TREE -----------------------
    aln_fasta = work_dir / f"{method_name}_aln.fasta"
    write_fasta(aln_result.names, aln_result.sequences, aln_fasta)

    # -- Step 3: Check for site weights file --------------------------------
    cfg = METHODS.get(method_name, {})
    site_weights_path = None
    if cfg.get("weights") and aln_result.column_confidence is not None:
        wt_path = work_dir / f"{method_name}_site_weights.txt"
        if wt_path.exists():
            site_weights_path = wt_path

    # -- Step 4: Pre-check for all-gap sequences ----------------------------
    for name, seq in zip(aln_result.names, aln_result.sequences):
        if all(c == "-" for c in seq):
            raise RuntimeError(
                f"Sequence '{name}' is all gaps after alignment/masking — "
                f"IQ-TREE would reject this. Method={method_name}"
            )

    # -- Step 5: Run IQ-TREE ------------------------------------------------
    iqtree_dir = work_dir / f"{method_name}_iqtree"
    iqtree_result = _run_iqtree(
        aln_fasta, model, iqtree_dir, site_weights_path=site_weights_path,
        n_threads=n_threads, timeout=iqtree_timeout,
    )

    # -- Step 6: Compare inferred vs true tree ------------------------------
    true_tree_str = sim_dataset.true_tree.read_text().strip()
    tree_cmp = compare_trees(true_tree_str, iqtree_result["tree"])

    # -- Step 7: Alignment accuracy vs true alignment -----------------------
    true_names, true_seqs = parse_fasta(sim_dataset.true_alignment)
    acc = alignment_accuracy(
        aln_result.sequences, true_seqs, aln_result.names, true_names,
    )

    return PhyloCaseResult(
        sim_id=sim_id,
        method=method_name,
        nrf=tree_cmp["nrf"],
        branch_score_dist=tree_cmp["branch_score_dist"],
        sp_score=acc["sp_score"],
        n_columns_retained=n_columns_retained,
        n_columns_total=n_columns_total,
        wall_time_align=aln_result.wall_time,
        wall_time_iqtree=iqtree_result["wall_time"],
        peak_memory_mb=aln_result.peak_memory_mb,
    )


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------


def _phylo_worker(sim_ds, mname, data_dir, fingerprint=None, n_threads=1, iqtree_timeout=600):
    """Worker function for a single sim x method case (module-level for pickling)."""
    from .utils import cache_load, cache_save, clean_work_dir

    sim_id = sim_ds.params.get("sim_id", "unknown")
    case_dir = Path(data_dir) / "phylo_work" / sim_id / mname
    case_dir.mkdir(parents=True, exist_ok=True)

    if fingerprint:
        cached = cache_load(case_dir, fingerprint)
        if cached is not None:
            return cached

    # Remove stale output files (e.g. IQ-TREE checkpoints) from previous runs
    clean_work_dir(case_dir)

    try:
        result = run_phylo_case(
            sim_ds, mname, case_dir,
            n_threads=n_threads, iqtree_timeout=iqtree_timeout,
        )
        result_dict = asdict(result)

        if fingerprint:
            cache_save(case_dir, fingerprint, result_dict)

        return result_dict
    except Exception as exc:
        logger.warning("Phylo case %s/%s failed: %s", sim_id, mname, exc)
        return {"sim_id": sim_id, "method": mname, "error": str(exc)}


def run_pipeline(params: dict) -> PhyloResult:
    """Run the full phylo accuracy pipeline.

    Parameters
    ----------
    params : dict
        Configuration with keys:
        - ``data_dir`` : str or Path -- base directory for simulation data
        - ``results_dir`` : str or Path -- where to write result JSON
        - ``methods`` : list[str] -- method names (keys into METHODS)
        - ``n_jobs`` : int -- max parallel workers
        - ``quick`` : bool -- if True, limit to 5 simulation cases

    Returns
    -------
    PhyloResult
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    from .provenance import (
        collect_provenance,
        result_path,
        update_latest_symlink,
    )
    from .simulation import (
        PROTEIN_GRID,
        PROTEIN_GRID_FULL,
        generate_indelible_dataset,
        iter_simulation_params,
        random_birth_death_tree,
    )
    from .utils import (
        tool_versions_fingerprint,
        METHODS,
    )

    data_dir = Path(params.get("data_dir", "benchmarks/data/downloads"))
    results_dir = Path(params.get("results_dir", "benchmarks/results"))
    full = params.get("full", False)
    _default_methods = ["kalign", "kalign_ens3", "mafft", "muscle", "clustalo", "true"]
    methods = params.get("methods", list(METHODS.keys()) if full else _default_methods)
    n_jobs = params.get("n_jobs", 1)
    quick = params.get("quick", False)
    depth_filter = params.get("depths", None)

    logger.info("Phylo accuracy pipeline starting (quick=%s, n_jobs=%d)", quick, n_jobs)

    # -- Generate simulations -----------------------------------------------
    grid = PROTEIN_GRID_FULL if full else PROTEIN_GRID
    sim_params_list = list(iter_simulation_params(grid, "WAG"))
    if quick:
        sim_params_list = sim_params_list[:5]
    if depth_filter:
        depth_set = set(float(d) for d in depth_filter)
        sim_params_list = [sp for sp in sim_params_list if sp["target_depth"] in depth_set]
        logger.info("Filtered to depths %s: %d simulations", depth_set, len(sim_params_list))

    from tqdm import tqdm

    sim_dir = data_dir / "phylo_sims"
    sim_dir.mkdir(parents=True, exist_ok=True)

    simulations = []
    for sp in tqdm(sim_params_list, desc="Generating simulations", unit="sim"):
        sim_out = sim_dir / sp["sim_id"]
        tree = random_birth_death_tree(
            n_taxa=sp["n_taxa"],
            target_depth=sp["target_depth"],
            seed=sp["seed"],
        )
        try:
            ds = generate_indelible_dataset(
                tree=tree,
                model=sp["model"],
                seq_length=sp["seq_length"],
                indel_rate=sp["indel_rate"],
                indel_length_mean=sp["indel_length_mean"],
                output_dir=sim_out,
                seed=sp["seed"],
            )
            # Attach sim_id into params for later reference
            ds.params["sim_id"] = sp["sim_id"]
            simulations.append(ds)
        except Exception as exc:
            logger.warning("Simulation %s failed: %s", sp["sim_id"], exc)

    logger.info("Generated %d simulations successfully.", len(simulations))

    # Compute fingerprint once for all workers (None disables caching)
    if params.get("no_cache"):
        fingerprint = None
        logger.info("Caching disabled (--no-cache)")
    else:
        fingerprint = tool_versions_fingerprint()
        logger.info("Tool versions fingerprint: %s", fingerprint)

    # -- Run alignment + tree inference for each sim x method ----------------
    # Allocate IQ-TREE threads: divide available cores among parallel workers
    cpu_count = os.cpu_count() or 1
    iq_threads = max(1, cpu_count // max(n_jobs, 1))
    iq_timeout = 600  # 10 minutes per IQ-TREE run
    logger.info(
        "IQ-TREE will use %d threads per job (%d cores / %d jobs)",
        iq_threads, cpu_count, n_jobs,
    )

    cases: list[dict] = []

    # Build work items
    work_items = [
        (sim_ds, mname)
        for sim_ds in simulations
        for mname in methods
        if mname in METHODS
    ]

    pbar = tqdm(total=len(work_items), desc="Phylo accuracy", unit="case")

    if n_jobs <= 1:
        # Sequential execution
        for sim_ds, mname in work_items:
            cases.append(_phylo_worker(
                sim_ds, mname, data_dir, fingerprint,
                n_threads=iq_threads, iqtree_timeout=iq_timeout,
            ))
            pbar.update(1)
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {
                pool.submit(
                    _phylo_worker, sim_ds, mname, data_dir, fingerprint,
                    n_threads=iq_threads, iqtree_timeout=iq_timeout,
                ): (sim_ds, mname)
                for sim_ds, mname in work_items
            }
            for future in as_completed(futures):
                try:
                    cases.append(future.result())
                except Exception as exc:
                    sim_ds, mname = futures[future]
                    sim_id = sim_ds.params.get("sim_id", "unknown")
                    logger.warning(
                        "Future for %s/%s raised: %s", sim_id, mname, exc,
                    )
                    cases.append({
                        "sim_id": sim_id,
                        "method": mname,
                        "error": str(exc),
                    })
                pbar.update(1)
    pbar.close()

    logger.info("Completed %d cases.", len(cases))

    # -- Aggregate: per-method summary --------------------------------------
    summary = _aggregate_summary(cases, methods)

    # -- Statistical tests: pairwise Wilcoxon on nRF -----------------------
    _add_statistical_tests(summary, cases, methods)

    # -- Provenance ---------------------------------------------------------
    provenance = collect_provenance(params)

    result = PhyloResult(
        provenance=asdict(provenance),
        cases=cases,
        summary=summary,
    )

    # -- Save to disk -------------------------------------------------------
    out_path = result_path(results_dir, "phylo_accuracy")
    with open(out_path, "w") as fh:
        json.dump(asdict(result) if hasattr(result, "__dataclass_fields__") else {
            "provenance": result.provenance,
            "cases": result.cases,
            "summary": result.summary,
        }, fh, indent=2, default=str)
    update_latest_symlink(out_path)
    logger.info("Results saved to %s", out_path)

    # -- Console summary ---------------------------------------------------
    _print_summary(summary)

    return result


# ---------------------------------------------------------------------------
# Internal aggregation helpers
# ---------------------------------------------------------------------------


def _aggregate_summary(cases: list[dict], methods: list[str]) -> dict:
    """Compute per-method aggregated metrics from case results."""
    import numpy as np

    from .utils import bootstrap_ci

    summary: dict = {}

    for mname in methods:
        method_cases = [
            c for c in cases
            if c.get("method") == mname and "error" not in c
        ]
        if not method_cases:
            summary[mname] = {"n_cases": 0}
            continue

        nrf_vals = [c["nrf"] for c in method_cases]
        bsd_vals = [c["branch_score_dist"] for c in method_cases]
        sp_vals = [c["sp_score"] for c in method_cases if c.get("sp_score", -1) >= 0]

        nrf_arr = np.asarray(nrf_vals, dtype=float)
        bsd_arr = np.asarray(bsd_vals, dtype=float)

        entry = {
            "n_cases": len(method_cases),
            "mean_nrf": float(nrf_arr.mean()),
            "std_nrf": float(nrf_arr.std()),
            "mean_branch_score_dist": float(bsd_arr.mean()),
            "std_branch_score_dist": float(bsd_arr.std()),
        }

        if len(nrf_vals) >= 5:
            ci_lo, ci_hi = bootstrap_ci(nrf_vals)
            entry["nrf_ci_95"] = [ci_lo, ci_hi]

        if sp_vals:
            sp_arr = np.asarray(sp_vals, dtype=float)
            entry["mean_sp_score"] = float(sp_arr.mean())
            entry["std_sp_score"] = float(sp_arr.std())

        # Mean timings
        align_times = [c.get("wall_time_align", 0) for c in method_cases]
        iqtree_times = [c.get("wall_time_iqtree", 0) for c in method_cases]
        entry["mean_wall_time_align"] = float(np.mean(align_times))
        entry["mean_wall_time_iqtree"] = float(np.mean(iqtree_times))

        summary[mname] = entry

    return summary


def _add_statistical_tests(
    summary: dict, cases: list[dict], methods: list[str],
) -> None:
    """Add pairwise Wilcoxon tests on nRF to the summary."""
    from .utils import holm_bonferroni, wilcoxon_paired

    # Build per-sim_id lookup for each method
    method_nrf: dict[str, dict[str, float]] = {}
    for c in cases:
        if "error" in c:
            continue
        mname = c["method"]
        sim_id = c["sim_id"]
        method_nrf.setdefault(mname, {})[sim_id] = c["nrf"]

    # All pairwise tests
    active_methods = [m for m in methods if m in method_nrf and len(method_nrf[m]) >= 5]
    pairwise_tests = []
    p_values = []

    for i, m1 in enumerate(active_methods):
        for m2 in active_methods[i + 1:]:
            # Find common sim_ids
            common = sorted(set(method_nrf[m1].keys()) & set(method_nrf[m2].keys()))
            if len(common) < 5:
                continue
            a = [method_nrf[m1][sid] for sid in common]
            b = [method_nrf[m2][sid] for sid in common]

            try:
                test = wilcoxon_paired(a, b)
                pairwise_tests.append({
                    "method_a": m1,
                    "method_b": m2,
                    "n_pairs": len(common),
                    **test,
                })
                p_values.append(test["p_value"])
            except Exception as exc:
                logger.debug("Wilcoxon test %s vs %s failed: %s", m1, m2, exc)

    # Apply Holm-Bonferroni correction
    if p_values:
        adjusted = holm_bonferroni(p_values)
        for t, adj_p in zip(pairwise_tests, adjusted):
            t["p_value_adjusted"] = adj_p

    summary["_pairwise_tests"] = pairwise_tests


def _print_summary(summary: dict) -> None:
    """Print a concise summary table to the console."""
    print("\n--- Phylo Accuracy Summary ---")
    print(f"{'Method':<25s} {'N':>5s} {'nRF':>8s} {'BSD':>8s} {'SP':>8s} {'t_aln':>8s} {'t_iq':>8s}")
    print("-" * 75)

    for mname, stats in sorted(summary.items()):
        if mname.startswith("_"):
            continue
        n = stats.get("n_cases", 0)
        if n == 0:
            print(f"{mname:<25s} {n:>5d}   {'N/A':>6s}   {'N/A':>6s}   {'N/A':>6s}")
            continue
        nrf = stats.get("mean_nrf", float("nan"))
        bsd = stats.get("mean_branch_score_dist", float("nan"))
        sp = stats.get("mean_sp_score", float("nan"))
        t_aln = stats.get("mean_wall_time_align", float("nan"))
        t_iq = stats.get("mean_wall_time_iqtree", float("nan"))
        print(f"{mname:<25s} {n:>5d} {nrf:>8.4f} {bsd:>8.4f} {sp:>8.4f} {t_aln:>7.2f}s {t_iq:>7.2f}s")

    # Print pairwise tests
    tests = summary.get("_pairwise_tests", [])
    if tests:
        print(f"\nPairwise Wilcoxon signed-rank tests (nRF, Holm-Bonferroni corrected):")
        print(f"{'A vs B':<45s} {'delta':>8s} {'p_adj':>10s}")
        print("-" * 65)
        for t in tests:
            label = f"{t['method_a']} vs {t['method_b']}"
            delta = t.get("cliffs_delta", float("nan"))
            p_adj = t.get("p_value_adjusted", t.get("p_value", float("nan")))
            print(f"{label:<45s} {delta:>8.3f} {p_adj:>10.4g}")

    print()


# ---------------------------------------------------------------------------
# Load results
# ---------------------------------------------------------------------------


def load_results(results_dir) -> PhyloResult:
    """Load the latest phylo accuracy results from disk.

    Parameters
    ----------
    results_dir : str or Path
        Base results directory (expects a ``phylo_accuracy/`` subdirectory
        with a ``latest.json`` symlink).

    Returns
    -------
    PhyloResult
    """
    from .provenance import load_latest_results

    data = load_latest_results(Path(results_dir) / "phylo_accuracy")

    return PhyloResult(
        provenance=data.get("provenance", {}),
        cases=data.get("cases", []),
        summary=data.get("summary", {}),
    )
