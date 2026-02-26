"""Pipeline 1: Positive Selection Analysis.

Measures how alignment quality affects false-positive and false-negative
rates in dN/dS site-model tests.  Simulates codon evolution under M7/M8
models with INDELible, aligns with multiple methods, runs HyPhy FUBAR,
and compares detected positively-selected sites against ground truth.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Codon translation / back-translation helpers
# ---------------------------------------------------------------------------

# Standard genetic code (NCBI translation table 1)
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate_dna(dna_seq: str) -> str:
    """Translate a DNA sequence to protein using the standard genetic code.

    Stop codons are omitted from the output.  The DNA sequence length
    must be a multiple of 3.
    """
    dna_seq = dna_seq.upper()
    if len(dna_seq) % 3 != 0:
        raise ValueError(
            f"DNA sequence length {len(dna_seq)} is not a multiple of 3"
        )
    protein = []
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i : i + 3]
        aa = _CODON_TABLE.get(codon, "X")
        if aa != "*":
            protein.append(aa)
    return "".join(protein)


def _backtranslate(
    protein_aln: str,
    original_dna: str,
) -> str:
    """Back-translate an aligned protein sequence to a codon alignment.

    Each non-gap character in *protein_aln* is replaced by the
    corresponding original DNA codon (3 nucleotides).  Each gap
    character (``-``) is replaced by ``---``.

    Parameters
    ----------
    protein_aln : str
        Aligned protein sequence (may contain ``-`` gap characters).
    original_dna : str
        Original unaligned DNA sequence (ungapped, length must be a
        multiple of 3).

    Returns
    -------
    str
        Codon-aligned DNA sequence where every position is exactly 3
        characters (either a codon or ``---``).
    """
    original_dna = original_dna.upper()

    # Build list of codons from original DNA (excluding stop codons)
    codons: list[str] = []
    for i in range(0, len(original_dna), 3):
        codon = original_dna[i : i + 3]
        aa = _CODON_TABLE.get(codon, "X")
        if aa != "*":
            codons.append(codon)

    codon_idx = 0
    result: list[str] = []
    for ch in protein_aln:
        if ch == "-":
            result.append("---")
        else:
            if codon_idx >= len(codons):
                raise ValueError(
                    f"Protein alignment has more residues than available "
                    f"codons ({len(codons)})"
                )
            result.append(codons[codon_idx])
            codon_idx += 1

    if codon_idx != len(codons):
        raise ValueError(
            f"Protein alignment used {codon_idx} residues but original DNA "
            f"has {len(codons)} non-stop codons"
        )

    return "".join(result)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------


@dataclass
class SelectionCaseResult:
    """Result of a single (simulation x method) case."""

    sim_id: str
    method: str
    true_positives: int
    false_positives: int
    false_negatives: int
    true_negatives: int
    precision: float
    recall: float
    f1: float
    sp_score: float
    wall_time_align: float
    wall_time_hyphy: float
    peak_memory_mb: float
    n_sites_tested: int


@dataclass
class SelectionResult:
    """Aggregate result for the full positive-selection pipeline."""

    provenance: dict
    cases: list[dict]
    summary: dict  # per-method aggregated metrics


# ---------------------------------------------------------------------------
# HyPhy FUBAR runner
# ---------------------------------------------------------------------------


def _run_hyphy_fubar(
    alignment_path: Path,
    tree_path: Path,
    work_dir: Path,
) -> dict:
    """Run HyPhy FUBAR on a codon alignment.

    Parameters
    ----------
    alignment_path : Path
        Path to aligned codon sequences in FASTA format.
    tree_path : Path
        Path to Newick tree file.
    work_dir : Path
        Working directory for HyPhy output files.

    Returns
    -------
    dict
        ``{"site_posteriors": list[float], "n_sites": int}``
        where each element of *site_posteriors* is the posterior
        probability of positive selection (``Prob[alpha < beta]``)
        at each codon site.

    Raises
    ------
    RuntimeError
        If HyPhy is not found or exits with an error.
    """
    hyphy_bin = shutil.which("hyphy") or shutil.which("HYPHYMP")
    if hyphy_bin is None:
        raise RuntimeError(
            "HyPhy is not installed or not on PATH.  "
            "Please install HyPhy (https://www.hyphy.org/) "
            "and ensure the 'hyphy' binary is accessible."
        )

    alignment_path = Path(alignment_path).resolve()
    tree_path = Path(tree_path).resolve()
    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    # HyPhy FUBAR outputs a JSON file next to the alignment
    # with the suffix .FUBAR.json
    output_json = work_dir / (alignment_path.name + ".FUBAR.json")

    cmd = [
        hyphy_bin,
        "fubar",
        "--alignment",
        str(alignment_path),
        "--tree",
        str(tree_path),
        "--output",
        str(output_json),
    ]

    # Build environment with HYPHY_PATH so HyPhy can find its batch files.
    # HYPHY_LIB (set in Containerfile) is not the variable HyPhy checks;
    # it needs HYPHY_PATH.
    env = dict(os.environ)
    hyphy_lib = env.get("HYPHY_LIB") or env.get("HYPHY_PATH")
    if hyphy_lib:
        env["HYPHY_PATH"] = hyphy_lib

    logger.debug("Running HyPhy FUBAR: %s", " ".join(cmd))
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(work_dir),
        env=env,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"HyPhy FUBAR failed (exit {proc.returncode}).\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    # Parse results JSON
    if not output_json.exists():
        # HyPhy may place the JSON in the same directory as the alignment
        alt_json = alignment_path.parent / (alignment_path.name + ".FUBAR.json")
        if alt_json.exists():
            output_json = alt_json
        else:
            raise RuntimeError(
                f"HyPhy FUBAR output JSON not found at {output_json} "
                f"or {alt_json}"
            )

    with open(output_json) as fh:
        fubar_data = json.load(fh)

    # FUBAR JSON structure:
    #   "MLE" -> "content" -> {
    #     "0" -> {"headers": [...], "0": [row0], "1": [row1], ...}
    #   }
    # Headers typically include:
    #   "alpha", "beta", "Prob[alpha>beta]", "Prob[alpha<beta]", ...
    # We need "Prob[alpha<beta]" = posterior probability of positive selection.
    site_posteriors: list[float] = []

    mle = fubar_data.get("MLE", {})
    content = mle.get("content", {})

    # Headers are at MLE.headers as [[name, description], ...]
    headers = mle.get("headers", [])
    prob_col_idx = None
    for i, h in enumerate(headers):
        h_name = h[0] if isinstance(h, list) else h
        if isinstance(h_name, str) and "alpha<beta" in h_name.lower():
            prob_col_idx = i
            break

    if prob_col_idx is None:
        prob_col_idx = 4  # standard FUBAR column order
        logger.warning(
            "Could not find 'Prob[alpha<beta]' in FUBAR headers; "
            "using column index %d as fallback",
            prob_col_idx,
        )

    # Content is {partition_key: list_of_rows}
    partition_key = "0"
    if partition_key not in content:
        keys = sorted(content.keys())
        if not keys:
            raise RuntimeError("No partition data found in FUBAR JSON output")
        partition_key = keys[0]

    rows = content[partition_key]
    for row in rows:
        if isinstance(row, list) and len(row) > prob_col_idx:
            site_posteriors.append(float(row[prob_col_idx]))
        else:
            site_posteriors.append(0.0)

    return {
        "site_posteriors": site_posteriors,
        "n_sites": len(site_posteriors),
    }


# ---------------------------------------------------------------------------
# Site classification
# ---------------------------------------------------------------------------


def _classify_sites(
    posteriors: list[float],
    true_classes: list[int],
    threshold: float = 0.9,
) -> dict:
    """Classify sites as positively selected and compute confusion metrics.

    Parameters
    ----------
    posteriors : list[float]
        FUBAR posterior probability of positive selection per site.
    true_classes : list[int]
        Ground-truth site classes: 1 = positive selection, 0 = neutral.
    threshold : float
        A site is predicted positive if posterior >= threshold.

    Returns
    -------
    dict
        ``{"true_positives", "false_positives", "false_negatives",
          "true_negatives", "precision", "recall", "f1",
          "n_sites_tested"}``
    """
    n = min(len(posteriors), len(true_classes))
    tp = fp = fn = tn = 0

    for i in range(n):
        predicted_positive = posteriors[i] >= threshold
        actual_positive = true_classes[i] == 1

        if predicted_positive and actual_positive:
            tp += 1
        elif predicted_positive and not actual_positive:
            fp += 1
        elif not predicted_positive and actual_positive:
            fn += 1
        else:
            tn += 1

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    return {
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "true_negatives": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "n_sites_tested": n,
    }


# ---------------------------------------------------------------------------
# Single-case runner
# ---------------------------------------------------------------------------


def run_selection_case(
    sim_dataset,
    method_name: str,
    work_dir: Path,
    threshold: float = 0.9,
) -> SelectionCaseResult:
    """Run the full positive-selection pipeline for one (sim x method) case.

    Parameters
    ----------
    sim_dataset : SimulatedDataset
        A simulated dataset (from :mod:`.simulation`).
    method_name : str
        Key into ``utils.METHODS``.
    work_dir : Path
        Scratch directory for this case.
    threshold : float
        FUBAR posterior threshold for calling positive selection.

    Returns
    -------
    SelectionCaseResult
    """
    from .utils import (
        alignment_accuracy,
        mask_alignment_by_confidence,
        METHODS,
        parse_fasta,
        run_method,
        write_fasta,
    )

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    sim_id = sim_dataset.params.get("sim_id", "unknown")
    if hasattr(sim_dataset, "params") and isinstance(sim_dataset.params, dict):
        sim_id = sim_dataset.params.get("sim_id", sim_id)

    # "true" method: use the true codon alignment directly as performance ceiling
    if method_name == "true":
        true_names, true_seqs = parse_fasta(sim_dataset.true_alignment)
        aln_path = work_dir / "aligned.fasta"
        write_fasta(true_names, true_seqs, aln_path)

        t0_hyphy = time.perf_counter()
        fubar_result = _run_hyphy_fubar(
            alignment_path=aln_path,
            tree_path=sim_dataset.true_tree,
            work_dir=work_dir,
        )
        wall_time_hyphy = time.perf_counter() - t0_hyphy

        metrics = _classify_sites(
            posteriors=fubar_result["site_posteriors"],
            true_classes=sim_dataset.site_classes,
            threshold=threshold,
        )

        return SelectionCaseResult(
            sim_id=sim_id,
            method="true",
            true_positives=metrics["true_positives"],
            false_positives=metrics["false_positives"],
            false_negatives=metrics["false_negatives"],
            true_negatives=metrics["true_negatives"],
            precision=metrics["precision"],
            recall=metrics["recall"],
            f1=metrics["f1"],
            sp_score=1.0,
            wall_time_align=0.0,
            wall_time_hyphy=wall_time_hyphy,
            peak_memory_mb=0.0,
            n_sites_tested=metrics["n_sites_tested"],
        )

    # 1. Translate-align-backtranslate to preserve codon reading frame
    #    Standard DNA alignment can insert 1- or 2-nt gaps that break
    #    the codon frame, causing HyPhy FUBAR to fail.
    #    Masking is deferred to the codon level (skip_masking=True)
    #    so that backtranslation has the full protein alignment.

    # 1a. Read unaligned DNA sequences
    dna_names, dna_seqs = parse_fasta(sim_dataset.unaligned)

    # 1b. Translate DNA to protein
    protein_seqs = [_translate_dna(s) for s in dna_seqs]

    # 1c. Write temporary protein FASTA
    prot_fasta = work_dir / "protein_input.fasta"
    write_fasta(dna_names, protein_seqs, prot_fasta)

    # 1d. Align proteins (skip masking so we get the full alignment)
    aln_result = run_method(
        method_name,
        prot_fasta,
        work_dir,
        seq_type="protein",
        skip_masking=True,
    )

    # 1e. Back-translate protein alignment to codon alignment
    #     Map aligned protein names back to original DNA sequences
    dna_by_name = dict(zip(dna_names, dna_seqs))
    codon_seqs: list[str] = []
    for name, prot_aln_seq in zip(aln_result.names, aln_result.sequences):
        orig_dna = dna_by_name[name]
        codon_seqs.append(_backtranslate(prot_aln_seq, orig_dna))

    # Sanity check: all codon-aligned sequences must be the same length
    if codon_seqs:
        expected_len = len(codon_seqs[0])
        for i, cs in enumerate(codon_seqs):
            if len(cs) != expected_len:
                raise ValueError(
                    f"Codon alignment length mismatch: seq {i} has "
                    f"{len(cs)} nt, expected {expected_len}"
                )

    # 1f. Apply confidence masking at the codon level if configured
    #     Each protein column maps to one codon column (3 nt).
    mask_threshold = METHODS.get(method_name, {}).get("mask")
    if mask_threshold is not None and aln_result.column_confidence is not None:
        # Expand protein column_confidence to codon columns (3x)
        codon_confidence = [c for c in aln_result.column_confidence for _ in range(3)]
        codon_seqs, _ = mask_alignment_by_confidence(
            codon_seqs, codon_confidence, mask_threshold
        )

    # 2. Write codon-aligned output to work_dir
    aln_path = work_dir / "aligned.fasta"
    # Check if masking removed all columns (empty alignment)
    aln_len = len(codon_seqs[0]) if codon_seqs else 0
    if aln_len == 0:
        raise RuntimeError(
            f"Masking removed all columns for {method_name} — "
            "cannot run FUBAR on an empty alignment"
        )
    write_fasta(aln_result.names, codon_seqs, aln_path)

    # 3. Compute alignment accuracy vs true alignment (codon level)
    true_names, true_seqs = parse_fasta(sim_dataset.true_alignment)
    acc = alignment_accuracy(
        test_sequences=codon_seqs,
        true_sequences=true_seqs,
        test_names=aln_result.names,
        true_names=true_names,
    )
    sp_score = acc.get("sp_score", -1.0)

    # 4. Run HyPhy FUBAR
    t0_hyphy = time.perf_counter()
    fubar_result = _run_hyphy_fubar(
        alignment_path=aln_path,
        tree_path=sim_dataset.true_tree,
        work_dir=work_dir,
    )
    wall_time_hyphy = time.perf_counter() - t0_hyphy

    # 5. Classify sites vs true labels
    metrics = _classify_sites(
        posteriors=fubar_result["site_posteriors"],
        true_classes=sim_dataset.site_classes,
        threshold=threshold,
    )

    return SelectionCaseResult(
        sim_id=sim_id,
        method=method_name,
        true_positives=metrics["true_positives"],
        false_positives=metrics["false_positives"],
        false_negatives=metrics["false_negatives"],
        true_negatives=metrics["true_negatives"],
        precision=metrics["precision"],
        recall=metrics["recall"],
        f1=metrics["f1"],
        sp_score=sp_score,
        wall_time_align=aln_result.wall_time,
        wall_time_hyphy=wall_time_hyphy,
        peak_memory_mb=aln_result.peak_memory_mb,
        n_sites_tested=metrics["n_sites_tested"],
    )


# ---------------------------------------------------------------------------
# Worker function for parallel execution
# ---------------------------------------------------------------------------


def _run_one(args: dict) -> dict:
    """Worker: align + score one (simulation x method) case.

    Simulations must already be generated before calling this.
    Returns a dict (serialisable result or error).
    """
    from .utils import cache_load, cache_save, clean_work_dir

    sim_dataset = args["sim_dataset"]
    method_name = args["method_name"]
    work_dir = Path(args["work_dir"])
    threshold = args.get("threshold", 0.9)
    fingerprint = args.get("fingerprint")

    sim_id = sim_dataset.params.get("sim_id", "unknown")

    # Check cache
    if fingerprint:
        cached = cache_load(work_dir, fingerprint)
        if cached is not None:
            return cached

    # Remove stale output files from previous runs
    clean_work_dir(work_dir)

    try:
        result = run_selection_case(
            sim_dataset=sim_dataset,
            method_name=method_name,
            work_dir=work_dir,
            threshold=threshold,
        )
        result_dict = asdict(result)

        if fingerprint:
            cache_save(work_dir, fingerprint, result_dict)

        return result_dict

    except Exception as exc:
        logger.error(
            "Error in %s / %s: %s", sim_id, method_name, exc, exc_info=True
        )
        return {
            "sim_id": sim_id,
            "method": method_name,
            "error": str(exc),
        }


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------


def _aggregate_summary(
    cases: list[dict],
    methods: list[str],
) -> dict:
    """Compute per-method aggregated metrics from case results.

    Returns a dict keyed by method name, each with mean precision,
    recall, F1, SP score, and timing.
    """
    from .utils import bootstrap_ci

    summary: dict[str, dict] = {}

    for method in methods:
        method_cases = [c for c in cases if c.get("method") == method and "error" not in c]
        n = len(method_cases)
        if n == 0:
            summary[method] = {
                "n_cases": 0,
                "mean_precision": 0.0,
                "mean_recall": 0.0,
                "mean_f1": 0.0,
                "mean_sp_score": 0.0,
                "mean_wall_time_align": 0.0,
                "mean_wall_time_hyphy": 0.0,
            }
            continue

        precisions = [c["precision"] for c in method_cases]
        recalls = [c["recall"] for c in method_cases]
        f1s = [c["f1"] for c in method_cases]
        sp_scores = [c["sp_score"] for c in method_cases]
        align_times = [c["wall_time_align"] for c in method_cases]
        hyphy_times = [c["wall_time_hyphy"] for c in method_cases]

        mean_f1 = sum(f1s) / n
        mean_precision = sum(precisions) / n
        mean_recall = sum(recalls) / n
        mean_sp = sum(sp_scores) / n

        entry: dict[str, Any] = {
            "n_cases": n,
            "mean_precision": mean_precision,
            "mean_recall": mean_recall,
            "mean_f1": mean_f1,
            "mean_sp_score": mean_sp,
            "mean_wall_time_align": sum(align_times) / n,
            "mean_wall_time_hyphy": sum(hyphy_times) / n,
        }

        # Bootstrap CIs on F1 if we have enough cases
        if n >= 5:
            try:
                ci_low, ci_high = bootstrap_ci(f1s)
                entry["f1_ci_lower"] = ci_low
                entry["f1_ci_upper"] = ci_high
            except Exception:
                pass

        summary[method] = entry

    return summary


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

# Default methods for positive-selection pipeline (subset of all methods)
_DEFAULT_METHODS = [
    "kalign",
    "kalign_cons",
    "kalign_ens3",
    "mafft",
    "muscle",
    "clustalo",
    "true",
]

_FULL_METHODS = [
    "kalign",
    "kalign_cons",
    "kalign_ens3",
    "mafft",
    "muscle",
    "clustalo",
    "true",
]


def run_pipeline(params: dict) -> SelectionResult:
    """Run the positive-selection analysis pipeline.

    Parameters
    ----------
    params : dict
        Configuration with keys:

        - ``data_dir`` : str or Path -- directory for simulation data / scratch
        - ``results_dir`` : str or Path -- directory for output JSONs
        - ``methods`` : list[str] -- alignment methods to test
        - ``n_jobs`` : int -- parallel workers (default 1)
        - ``quick`` : bool -- if True, run only 5 simulations (default False)
        - ``threshold`` : float -- FUBAR posterior threshold (default 0.9)
        - ``max_sims`` : int or None -- max simulations to run

    Returns
    -------
    SelectionResult
    """
    from .provenance import collect_provenance, result_path, update_latest_symlink
    from .simulation import (
        CODON_GRID,
        CODON_GRID_FULL,
        SimulationGrid,
        generate_indelible_dataset,
        iter_simulation_params,
        random_birth_death_tree,
    )
    from .utils import tool_versions_fingerprint

    data_dir = Path(params.get("data_dir", "benchmarks/data"))
    results_dir = Path(params.get("results_dir", "benchmarks/results"))
    full = params.get("full", False)
    methods = params.get("methods", _FULL_METHODS if full else _DEFAULT_METHODS)
    n_jobs = params.get("n_jobs", 1)
    quick = params.get("quick", False)
    threshold = params.get("threshold", 0.9)
    max_sims = params.get("max_sims", None)

    data_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Positive selection pipeline starting")
    logger.info("Methods: %s", methods)
    logger.info("Quick mode: %s", quick)

    # Build the simulation parameter grid
    if quick:
        # Minimal grid for smoke-testing: fewer taxa, depths, reps
        grid = SimulationGrid(
            n_taxa=[8],
            tree_depths=[0.7],
            indel_rates=[0.05],
            psel_fractions=[0.0, 0.10],
            replicates=3,
            n_codons=200,
        )
    else:
        grid = CODON_GRID_FULL if full else CODON_GRID

    # Collect simulation parameter sets
    sim_params_list = list(iter_simulation_params(grid, model="M8"))

    if max_sims is not None:
        sim_params_list = sim_params_list[:max_sims]
    if quick and len(sim_params_list) > 5:
        sim_params_list = sim_params_list[:5]

    # Compute fingerprint once for all workers (None disables caching)
    if params.get("no_cache"):
        fingerprint = None
        logger.info("Caching disabled (--no-cache)")
    else:
        fingerprint = tool_versions_fingerprint()
        logger.info("Tool versions fingerprint: %s", fingerprint)

    # Step 1: Generate simulated datasets sequentially (avoids race conditions)
    from tqdm import tqdm

    sim_datasets = []
    for sp in tqdm(sim_params_list, desc="Generating codon simulations", unit="sim"):
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
            psel_fraction=sp.get("psel_fraction", 0.0),
            omega_positive=sp.get("omega_positive", 3.0),
            kappa=sp.get("kappa", 2.0),
        )
        sim_dataset.params["sim_id"] = sim_id
        sim_datasets.append(sim_dataset)

    logger.info(
        "Running %d simulations x %d methods = %d total cases",
        len(sim_datasets),
        len(methods),
        len(sim_datasets) * len(methods),
    )

    # Step 2: Build work items for alignment + scoring
    work_items = []
    for sim_dataset in sim_datasets:
        sim_id = sim_dataset.params["sim_id"]
        sim_dir = data_dir / sim_id
        for method_name in methods:
            work_dir = sim_dir / f"work_{method_name}"
            work_items.append(
                {
                    "sim_dataset": sim_dataset,
                    "method_name": method_name,
                    "work_dir": str(work_dir),
                    "threshold": threshold,
                    "fingerprint": fingerprint,
                }
            )

    # Step 3: Execute (parallel or serial)
    cases: list[dict] = []
    pbar = tqdm(total=len(work_items), desc="Positive selection", unit="case")

    if n_jobs > 1:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {pool.submit(_run_one, item): item for item in work_items}
            for future in as_completed(futures):
                try:
                    result = future.result()
                    cases.append(result)
                except Exception as exc:
                    item = futures[future]
                    logger.error("Worker failed: %s", exc)
                    cases.append(
                        {
                            "sim_id": item["sim_dataset"].params.get("sim_id", "unknown"),
                            "method": item["method_name"],
                            "error": str(exc),
                        }
                    )
                pbar.update(1)
    else:
        for item in work_items:
            result = _run_one(item)
            cases.append(result)
            pbar.update(1)
    pbar.close()

    # Filter out errors for summary
    successful = [c for c in cases if "error" not in c]
    error_count = len(cases) - len(successful)
    if error_count > 0:
        logger.warning("%d / %d cases had errors", error_count, len(cases))

    # Aggregate summary
    summary = _aggregate_summary(cases, methods)

    # Statistical tests: wilcoxon between kalign_ens3 and each other method
    stats_tests = _compute_statistical_tests(cases, methods)
    if stats_tests:
        summary["statistical_tests"] = stats_tests

    # Collect provenance
    provenance = collect_provenance(
        parameters={
            "pipeline": "positive_selection",
            "n_simulations": len(sim_datasets),
            "methods": methods,
            "n_jobs": n_jobs,
            "quick": quick,
            "threshold": threshold,
            "max_sims": max_sims,
        }
    )

    # Build result
    sel_result = SelectionResult(
        provenance=asdict(provenance),
        cases=cases,
        summary=summary,
    )

    # Save JSON
    out_path = result_path(results_dir, "positive_selection")
    with open(out_path, "w") as fh:
        json.dump(asdict(sel_result), fh, indent=2, default=str)
    update_latest_symlink(out_path)
    logger.info("Results saved to %s", out_path)

    # Print summary table
    _print_summary(summary, methods)

    return sel_result


def _compute_statistical_tests(
    cases: list[dict], methods: list[str]
) -> dict:
    """Run Wilcoxon signed-rank tests comparing kalign_ens3 to other methods."""
    from .utils import holm_bonferroni, wilcoxon_paired

    baseline = "kalign_ens3"
    if baseline not in methods:
        return {}

    # Group F1 scores by sim_id for paired comparisons
    baseline_f1: dict[str, float] = {}
    for c in cases:
        if c.get("method") == baseline and "error" not in c:
            baseline_f1[c["sim_id"]] = c["f1"]

    if len(baseline_f1) < 5:
        logger.warning(
            "Too few baseline cases (%d) for statistical tests", len(baseline_f1)
        )
        return {}

    comparisons = [m for m in methods if m != baseline]
    test_results: dict[str, dict] = {}
    raw_p_values: list[float] = []
    comp_names: list[str] = []

    for method in comparisons:
        method_f1: dict[str, float] = {}
        for c in cases:
            if c.get("method") == method and "error" not in c:
                method_f1[c["sim_id"]] = c["f1"]

        # Paired comparison: only use sim_ids present in both
        common_ids = sorted(set(baseline_f1.keys()) & set(method_f1.keys()))
        if len(common_ids) < 5:
            continue

        a = [baseline_f1[sid] for sid in common_ids]
        b = [method_f1[sid] for sid in common_ids]

        try:
            result = wilcoxon_paired(a, b)
            test_results[method] = {
                "n_pairs": len(common_ids),
                **result,
            }
            raw_p_values.append(result["p_value"])
            comp_names.append(method)
        except Exception as exc:
            logger.warning("Wilcoxon test failed for %s: %s", method, exc)

    # Apply Holm-Bonferroni correction
    if raw_p_values:
        adjusted = holm_bonferroni(raw_p_values)
        for name, adj_p in zip(comp_names, adjusted):
            test_results[name]["p_value_adjusted"] = adj_p

    return test_results


def _print_summary(summary: dict, methods: list[str]) -> None:
    """Print a summary table to the logger."""
    logger.info("=" * 72)
    logger.info("Positive Selection Pipeline — Summary")
    logger.info("=" * 72)
    logger.info(
        "%-22s %6s %8s %8s %8s %8s",
        "Method", "N", "Prec", "Recall", "F1", "SP",
    )
    logger.info("-" * 72)
    for m in methods:
        s = summary.get(m, {})
        if s.get("n_cases", 0) == 0:
            logger.info("%-22s %6d %8s", m, 0, "no data")
            continue
        logger.info(
            "%-22s %6d %8.3f %8.3f %8.3f %8.3f",
            m,
            s["n_cases"],
            s["mean_precision"],
            s["mean_recall"],
            s["mean_f1"],
            s["mean_sp_score"],
        )
    logger.info("=" * 72)


# ---------------------------------------------------------------------------
# Load results
# ---------------------------------------------------------------------------


def load_results(results_dir: Path) -> SelectionResult:
    """Load the latest positive-selection results from disk.

    Parameters
    ----------
    results_dir : Path
        Top-level results directory (e.g. ``benchmarks/results``).

    Returns
    -------
    SelectionResult
    """
    from .provenance import load_latest_results

    data = load_latest_results(Path(results_dir) / "positive_selection")
    return SelectionResult(
        provenance=data.get("provenance", {}),
        cases=data.get("cases", []),
        summary=data.get("summary", {}),
    )
