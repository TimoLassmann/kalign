"""Shared utilities for downstream benchmark pipelines.

Provides alignment runners, scoring helpers, statistical tests,
tree comparison wrappers, and consistent figure styling for the
kalign downstream benchmarks.
"""

from __future__ import annotations

import hashlib
import json as _json
import logging
import os
import resource
import subprocess
import tempfile
import time
from pathlib import Path
from typing import NamedTuple, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-case caching helpers
# ---------------------------------------------------------------------------


def tool_versions_fingerprint() -> str:
    """Return a short hex fingerprint of all external tool versions.

    Calls :func:`~.provenance.collect_tool_versions` once and hashes the
    sorted result.  The fingerprint changes whenever any tool is upgraded
    (or removed / added), causing cached results to be recomputed.
    """
    from .provenance import collect_tool_versions

    versions = collect_tool_versions()
    blob = _json.dumps(versions, sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def cache_load(work_dir: Path, fingerprint: str) -> Optional[dict]:
    """Load a cached result from *work_dir* if the fingerprint matches.

    Returns the cached result dict, or ``None`` on cache miss (file
    missing, fingerprint mismatch, or corrupt JSON).
    """
    cache_file = Path(work_dir) / "_cache.json"
    if not cache_file.exists():
        return None
    try:
        with open(cache_file) as fh:
            data = _json.load(fh)
        if data.get("fingerprint") == fingerprint:
            return data["result"]
    except (OSError, ValueError, KeyError):
        pass
    return None


def cache_save(work_dir: Path, fingerprint: str, result: dict) -> None:
    """Write a result dict to *work_dir/_cache.json* with the given fingerprint."""
    cache_file = Path(work_dir) / "_cache.json"
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    with open(cache_file, "w") as fh:
        _json.dump({"fingerprint": fingerprint, "result": result}, fh)


def clean_work_dir(work_dir: Path) -> None:
    """Remove stale output files from *work_dir*, preserving ``_cache.json``.

    Called before re-running a case to prevent leftover files (e.g. IQ-TREE
    checkpoints, partial FASTA files) from interfering with fresh runs.
    """
    work_dir = Path(work_dir)
    if not work_dir.is_dir():
        return
    for item in work_dir.iterdir():
        if item.name == "_cache.json":
            continue
        if item.is_dir():
            import shutil
            shutil.rmtree(item, ignore_errors=True)
        else:
            item.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# AlignResult
# ---------------------------------------------------------------------------


class AlignResult(NamedTuple):
    """Result returned by every alignment runner."""

    sequences: list[str]
    names: list[str]
    column_confidence: Optional[list[float]]
    residue_confidence: Optional[list[list[float]]]
    wall_time: float
    peak_memory_mb: float


# ---------------------------------------------------------------------------
# FASTA parser / writer helpers
# ---------------------------------------------------------------------------


def parse_fasta(path: Path) -> tuple[list[str], list[str]]:
    """Parse a FASTA file.

    Returns
    -------
    names : list[str]
        Sequence names (everything after '>' up to the first whitespace).
    sequences : list[str]
        Corresponding sequences with no whitespace.
    """
    names: list[str] = []
    sequences: list[str] = []
    current_seq: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
                # Name is the first whitespace-delimited token after '>'
                names.append(line[1:].split()[0] if line[1:].strip() else "")
            else:
                current_seq.append(line.strip())
        if current_seq:
            sequences.append("".join(current_seq))
    return names, sequences


def write_fasta(names: list[str], sequences: list[str], path: Path) -> None:
    """Write sequences to a FASTA file."""
    with open(path, "w") as fh:
        for name, seq in zip(names, sequences):
            fh.write(f">{name}\n")
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


# ---------------------------------------------------------------------------
# Confidence masking
# ---------------------------------------------------------------------------


def mask_alignment_by_confidence(
    sequences: list[str],
    column_confidence: list[float],
    threshold: float,
) -> tuple[list[str], int]:
    """Remove alignment columns where confidence is below *threshold*.

    Parameters
    ----------
    sequences : list[str]
        Aligned sequences (all same length).
    column_confidence : list[float]
        Per-column confidence values in [0, 1].
    threshold : float
        Columns with confidence < threshold are removed.

    Returns
    -------
    masked_sequences : list[str]
        Sequences with low-confidence columns stripped.
    n_columns_retained : int
        Number of columns that survived the filter.
    """
    if not sequences:
        return [], 0
    ncols = len(sequences[0])
    keep = [j for j in range(ncols) if column_confidence[j] >= threshold]
    masked = ["".join(seq[j] for j in keep) for seq in sequences]
    return masked, len(keep)


# ---------------------------------------------------------------------------
# IQ-TREE site weights
# ---------------------------------------------------------------------------


def write_site_weights(
    column_confidence: list[float],
    output_path: Path,
) -> None:
    """Write an IQ-TREE ``-a`` site-weight file.

    The file contains one weight per line.  Values are clamped to [0, 1].
    """
    with open(output_path, "w") as fh:
        for w in column_confidence:
            fh.write(f"{max(0.0, min(1.0, w)):.6f}\n")


# ---------------------------------------------------------------------------
# Alignment accuracy via kalign.compare
# ---------------------------------------------------------------------------


def alignment_accuracy(
    test_sequences: list[str],
    true_sequences: list[str],
    test_names: list[str],
    true_names: list[str],
) -> dict:
    """Compute SP and TC scores of *test* alignment vs *true* alignment.

    Writes both alignments to temporary FASTA files and calls
    ``kalign.compare_detailed``.  If ``kalign.compare_detailed`` is not
    available, falls back to ``kalign.compare`` for SP only.  If neither
    is available, returns sentinel values of -1.

    Returns
    -------
    dict
        ``{"sp_score": float, "tc_score": float}``
    """
    try:
        import kalign as _kalign
    except ImportError:
        return {"sp_score": -1.0, "tc_score": -1.0}

    # Quick validation: need at least 2 non-empty sequences
    if len(test_sequences) < 2 or len(true_sequences) < 2:
        return {"sp_score": -1.0, "tc_score": -1.0}
    if not test_sequences[0] or not true_sequences[0]:
        return {"sp_score": -1.0, "tc_score": -1.0}

    with tempfile.TemporaryDirectory() as tmpdir:
        ref_path = os.path.join(tmpdir, "ref.fasta")
        test_path = os.path.join(tmpdir, "test.fasta")
        write_fasta(true_names, true_sequences, Path(ref_path))
        write_fasta(test_names, test_sequences, Path(test_path))

        # Verify files were written (guard against race / disk issues)
        if os.path.getsize(ref_path) == 0 or os.path.getsize(test_path) == 0:
            logger.warning("alignment_accuracy: empty FASTA written, skipping")
            return {"sp_score": -1.0, "tc_score": -1.0}

        try:
            if hasattr(_kalign, "compare_detailed"):
                d = _kalign.compare_detailed(ref_path, test_path, max_gap_frac=-1.0)
                return {"sp_score": d.get("recall", -1.0), "tc_score": d.get("tc", -1.0)}
            elif hasattr(_kalign, "compare"):
                sp = _kalign.compare(ref_path, test_path)
                return {"sp_score": sp, "tc_score": -1.0}
            else:
                return {"sp_score": -1.0, "tc_score": -1.0}
        except Exception as exc:
            logger.warning("alignment_accuracy failed: %s", exc)
            return {"sp_score": -1.0, "tc_score": -1.0}


# ---------------------------------------------------------------------------
# Tree comparison (dendropy)
# ---------------------------------------------------------------------------


def compare_trees(true_tree: str, inferred_tree: str) -> dict:
    """Compare two Newick tree strings.

    Uses ``dendropy.calculate.treecompare`` (imported lazily).

    Returns
    -------
    dict
        ``{"nrf": float, "quartet_dist": float, "branch_score_dist": float}``
    """
    import dendropy
    from dendropy.calculate import treecompare

    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(data=true_tree, schema="newick", taxon_namespace=tns)
    t2 = dendropy.Tree.get(data=inferred_tree, schema="newick", taxon_namespace=tns)

    # Ensure the same leaf set
    t1.encode_bipartitions()
    t2.encode_bipartitions()

    nrf = treecompare.symmetric_difference(t1, t2) / (2 * (len(tns) - 3))
    quartet_dist = treecompare.euclidean_distance(t1, t2)
    bsd = treecompare.weighted_robinson_foulds_distance(t1, t2)

    return {
        "nrf": nrf,
        "quartet_dist": quartet_dist,
        "branch_score_dist": bsd,
    }


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------


def bootstrap_ci(
    values: list[float],
    n_bootstrap: int = 10000,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Compute a bootstrap confidence interval for the mean.

    Parameters
    ----------
    values : list[float]
        Observed values.
    n_bootstrap : int
        Number of bootstrap resamples.
    alpha : float
        Significance level (default 0.05 for a 95 % CI).

    Returns
    -------
    (lower, upper) : tuple[float, float]
    """
    import numpy as np

    arr = np.asarray(values, dtype=float)
    n = len(arr)
    rng = np.random.default_rng(seed=42)
    boot_means = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        sample = rng.choice(arr, size=n, replace=True)
        boot_means[i] = sample.mean()
    lower = float(np.percentile(boot_means, 100 * alpha / 2))
    upper = float(np.percentile(boot_means, 100 * (1 - alpha / 2)))
    return lower, upper


def wilcoxon_paired(a: list[float], b: list[float]) -> dict:
    """Wilcoxon signed-rank test with Cliff's delta effect size.

    Parameters
    ----------
    a, b : list[float]
        Paired observations (same length).

    Returns
    -------
    dict
        ``{"statistic": float, "p_value": float, "cliffs_delta": float}``
    """
    from scipy.stats import wilcoxon as _wilcoxon

    stat, pval = _wilcoxon(a, b)

    # Cliff's delta (non-parametric effect size)
    n = len(a)
    count = 0
    for ai, bi in zip(a, b):
        if ai > bi:
            count += 1
        elif ai < bi:
            count -= 1
    cliffs_d = count / n

    return {
        "statistic": float(stat),
        "p_value": float(pval),
        "cliffs_delta": cliffs_d,
    }


def holm_bonferroni(p_values: list[float]) -> list[float]:
    """Apply Holm-Bonferroni correction to a list of p-values.

    Returns adjusted p-values in the same order as the input.
    """
    m = len(p_values)
    if m == 0:
        return []

    # Sort indices by p-value
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * m
    cummax = 0.0
    for rank, (orig_idx, p) in enumerate(indexed):
        corrected = p * (m - rank)
        # Enforce monotonicity: adjusted value must be >= previous
        cummax = max(cummax, corrected)
        adjusted[orig_idx] = min(cummax, 1.0)
    return adjusted


# ---------------------------------------------------------------------------
# Method colours (for consistent figure styling)
# ---------------------------------------------------------------------------

METHOD_COLORS = {
    "kalign": "#1f77b4",       # blue
    "kalign_cons": "#2ca02c",  # green
    "kalign_ens3": "#ff7f0e",  # orange
    "mafft": "#d62728",        # red
    "muscle": "#9467bd",       # purple
    "clustalo": "#7f7f7f",     # grey
    "true": "#000000",         # black (reference)
}


# ---------------------------------------------------------------------------
# Internal helpers for runners
# ---------------------------------------------------------------------------


def _peak_mem_children_mb() -> float:
    """Return peak RSS (in MB) of child processes via getrusage.

    On macOS ``ru_maxrss`` is in bytes; on Linux it is in kilobytes.
    """
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    maxrss = ru.ru_maxrss
    import sys

    if sys.platform == "darwin":
        return maxrss / (1024 * 1024)
    return maxrss / 1024


def _parse_fasta_string(text: str) -> tuple[list[str], list[str]]:
    """Parse FASTA from a string.  Returns (names, sequences)."""
    names: list[str] = []
    sequences: list[str] = []
    current: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
            names.append(line[1:].split()[0] if line[1:].strip() else "")
        else:
            current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return names, sequences


# ---------------------------------------------------------------------------
# Alignment runners
# ---------------------------------------------------------------------------


def run_kalign(
    input_fasta: Path,
    ensemble: int = 0,
    seq_type: str = "auto",
    tgpo: float | None = None,
    terminal_dist_scale: float | None = None,
    refine: str | None = None,
    vsm_amax: float | None = None,
    realign: int | None = None,
    consistency: int | None = None,
    consistency_weight: float | None = None,
) -> AlignResult:
    """Run kalign via the Python API.

    Parameters
    ----------
    input_fasta : Path
        Path to unaligned FASTA.
    ensemble : int
        Number of ensemble runs (0 = off).
    seq_type : str
        Sequence type passed to ``kalign.align_from_file``.

    Returns
    -------
    AlignResult
    """
    import kalign as _kalign

    start = time.perf_counter()
    extra = {}
    if tgpo is not None:
        extra["tgpo"] = tgpo
    if terminal_dist_scale is not None:
        extra["terminal_dist_scale"] = terminal_dist_scale
    if refine is not None:
        extra["refine"] = refine
    if vsm_amax is not None:
        extra["vsm_amax"] = vsm_amax
    if realign is not None:
        extra["realign"] = realign
    if consistency is not None:
        extra["consistency"] = consistency
    if consistency_weight is not None:
        extra["consistency_weight"] = consistency_weight
    result = _kalign.align_from_file(
        str(input_fasta),
        seq_type=seq_type,
        ensemble=ensemble,
        **extra,
    )
    wall = time.perf_counter() - start

    # AlignedSequences supports unpacking and has .column_confidence etc.
    names = list(result.names)
    sequences = list(result.sequences)
    col_conf = result.column_confidence  # may be None
    res_conf = result.residue_confidence  # may be None

    ru = resource.getrusage(resource.RUSAGE_SELF)
    import sys

    if sys.platform == "darwin":
        mem_mb = ru.ru_maxrss / (1024 * 1024)
    else:
        mem_mb = ru.ru_maxrss / 1024

    return AlignResult(
        sequences=sequences,
        names=names,
        column_confidence=col_conf,
        residue_confidence=res_conf,
        wall_time=wall,
        peak_memory_mb=mem_mb,
    )


def run_mafft(input_fasta: Path, timeout: int | None = None) -> AlignResult:
    """Run MAFFT (``mafft --auto``) via subprocess."""
    start = time.perf_counter()
    proc = subprocess.run(
        ["mafft", "--auto", str(input_fasta)],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    wall = time.perf_counter() - start
    mem_mb = _peak_mem_children_mb()

    if proc.returncode != 0:
        raise RuntimeError(f"mafft failed (exit {proc.returncode}): {proc.stderr}")

    names, sequences = _parse_fasta_string(proc.stdout)
    return AlignResult(
        sequences=sequences,
        names=names,
        column_confidence=None,
        residue_confidence=None,
        wall_time=wall,
        peak_memory_mb=mem_mb,
    )


def run_muscle(input_fasta: Path, timeout: int | None = None) -> AlignResult:
    """Run MUSCLE v5 (``muscle -align ... -output ...``) via subprocess."""
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        out_path = tmp.name
    try:
        start = time.perf_counter()
        proc = subprocess.run(
            ["muscle", "-align", str(input_fasta), "-output", out_path],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        wall = time.perf_counter() - start
        mem_mb = _peak_mem_children_mb()

        if proc.returncode != 0:
            raise RuntimeError(
                f"muscle failed (exit {proc.returncode}): {proc.stderr}"
            )

        names, sequences = parse_fasta(Path(out_path))
    finally:
        if os.path.exists(out_path):
            os.unlink(out_path)

    return AlignResult(
        sequences=sequences,
        names=names,
        column_confidence=None,
        residue_confidence=None,
        wall_time=wall,
        peak_memory_mb=mem_mb,
    )


def run_clustalo(input_fasta: Path, timeout: int | None = None) -> AlignResult:
    """Run Clustal Omega (``clustalo -i ... -o ...``) via subprocess."""
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        out_path = tmp.name
    try:
        start = time.perf_counter()
        proc = subprocess.run(
            [
                "clustalo",
                "-i",
                str(input_fasta),
                "-o",
                out_path,
                "--outfmt=fasta",
                "--force",
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        wall = time.perf_counter() - start
        mem_mb = _peak_mem_children_mb()

        if proc.returncode != 0:
            raise RuntimeError(
                f"clustalo failed (exit {proc.returncode}): {proc.stderr}"
            )

        names, sequences = parse_fasta(Path(out_path))
    finally:
        if os.path.exists(out_path):
            os.unlink(out_path)

    return AlignResult(
        sequences=sequences,
        names=names,
        column_confidence=None,
        residue_confidence=None,
        wall_time=wall,
        peak_memory_mb=mem_mb,
    )


def run_guidance2(
    input_fasta: Path,
    output_dir: Path,
    seq_type: str = "aa",
    n_bootstrap: int = 100,
) -> AlignResult:
    """Run GUIDANCE2 via the perl script.

    Expects ``guidance.pl`` to be on ``$PATH`` or pointed to by the
    ``GUIDANCE2_PATH`` environment variable.

    Parameters
    ----------
    input_fasta : Path
        Unaligned FASTA input.
    output_dir : Path
        Working / output directory for GUIDANCE2 files.
    seq_type : str
        ``"aa"`` for protein, ``"nuc"`` for nucleotide (DNA/RNA).
    n_bootstrap : int
        Number of bootstrap guide-tree perturbations.

    Returns
    -------
    AlignResult
    """
    import shutil as _shutil

    guidance_bin = os.environ.get(
        "GUIDANCE2_PATH",
        _shutil.which("guidance2") or _shutil.which("guidance.pl") or "guidance.pl",
    )
    output_dir = Path(output_dir).resolve()
    # Clean any stale output from previous runs — GUIDANCE2 may skip steps
    # or produce inconsistent results if old output files are present.
    if output_dir.exists():
        _shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    input_fasta = Path(input_fasta).resolve()

    # GUIDANCE2 type flags
    g2_type = "aa" if seq_type in ("aa", "protein") else "nuc"

    cmd = [
        "perl",
        guidance_bin,
        "--seqFile",
        str(input_fasta),
        "--msaProgram",
        "MAFFT",
        "--seqType",
        g2_type,
        "--outDir",
        str(output_dir),
        "--bootstraps",
        str(n_bootstrap),
    ]

    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - start
    mem_mb = _peak_mem_children_mb()

    if proc.returncode != 0:
        raise RuntimeError(
            f"GUIDANCE2 failed (exit {proc.returncode}): {proc.stderr}"
        )

    # Parse the base MSA
    msa_path = output_dir / "MSA.MAFFT.aln.With_Names"
    if not msa_path.exists():
        # Fallback to the filtered alignment
        msa_path = output_dir / "MSA.MAFFT.Without_low_SP_Col.With_Names"
    if not msa_path.exists():
        raise FileNotFoundError(
            f"GUIDANCE2 output alignment not found in {output_dir}"
        )
    names, sequences = parse_fasta(msa_path)

    # Parse column confidence scores
    col_conf: Optional[list[float]] = None
    col_scores_path = output_dir / "MSA.MAFFT.Guidance2_col_col.scr"
    if col_scores_path.exists():
        col_conf = _parse_guidance2_col_scores(col_scores_path, len(sequences[0]))

    # Parse residue confidence scores
    res_conf: Optional[list[list[float]]] = None
    res_scores_path = output_dir / "MSA.MAFFT.Guidance2_res_pair_res.scr"
    if res_scores_path.exists():
        res_conf = _parse_guidance2_res_scores(
            res_scores_path, len(names), len(sequences[0])
        )

    return AlignResult(
        sequences=sequences,
        names=names,
        column_confidence=col_conf,
        residue_confidence=res_conf,
        wall_time=wall,
        peak_memory_mb=mem_mb,
    )


def _parse_guidance2_col_scores(
    path: Path, ncols: int
) -> list[float]:
    """Parse GUIDANCE2 per-column confidence file.

    The file has lines ``<col_index> <score>`` (1-based).
    """
    scores = [0.0] * ncols
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    col = int(parts[0]) - 1  # 1-based to 0-based
                    val = float(parts[1])
                    if 0 <= col < ncols:
                        scores[col] = val
                except ValueError:
                    continue
    return scores


def _parse_guidance2_res_scores(
    path: Path, nseqs: int, ncols: int
) -> list[list[float]]:
    """Parse GUIDANCE2 per-residue confidence file.

    The file has lines ``<seq_index> <col_index> <score>`` (1-based).
    """
    scores = [[0.0] * ncols for _ in range(nseqs)]
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    seq_i = int(parts[0]) - 1
                    col_j = int(parts[1]) - 1
                    val = float(parts[2])
                    if 0 <= seq_i < nseqs and 0 <= col_j < ncols:
                        scores[seq_i][col_j] = val
                except ValueError:
                    continue
    return scores


# ---------------------------------------------------------------------------
# METHODS registry
# ---------------------------------------------------------------------------

METHODS: dict[str, dict] = {
    "kalign": {
        "fn": run_kalign,
        "ensemble": 0,
        "mask": None,
        "vsm_amax": 2.0,
    },
    "kalign_cons": {
        "fn": run_kalign,
        "ensemble": 0,
        "mask": None,
        "vsm_amax": 2.0,
        "consistency": 5,
        "consistency_weight": 2.0,
    },
    "kalign_ens3": {
        "fn": run_kalign,
        "ensemble": 3,
        "mask": None,
        "vsm_amax": 2.0,
        "realign": 1,
    },
    "mafft": {"fn": run_mafft, "mask": None},
    "muscle": {"fn": run_muscle, "mask": None},
    "clustalo": {"fn": run_clustalo, "mask": None},
    "true": {"fn": None},
}


# ---------------------------------------------------------------------------
# High-level method dispatcher
# ---------------------------------------------------------------------------


def run_method(
    method_name: str,
    input_fasta: Path,
    work_dir: Path,
    seq_type: str = "auto",
    skip_masking: bool = False,
    timeout: int | None = None,
) -> AlignResult:
    """Run an alignment method by name.

    Looks up *method_name* in :data:`METHODS`, dispatches to the
    appropriate runner, and applies confidence masking or site-weight
    output when the method entry requests it.

    Parameters
    ----------
    method_name : str
        Key into :data:`METHODS`.
    input_fasta : Path
        Unaligned FASTA file.
    work_dir : Path
        Scratch / output directory (used by GUIDANCE2 and for weight
        files).
    seq_type : str
        Sequence type hint (``"auto"``, ``"protein"``, ``"dna"``, etc.).

    Returns
    -------
    AlignResult
    """
    if method_name not in METHODS:
        raise ValueError(
            f"Unknown method {method_name!r}.  "
            f"Available: {sorted(METHODS.keys())}"
        )

    cfg = METHODS[method_name]
    fn = cfg["fn"]

    if fn is None:
        raise ValueError(
            f"Method {method_name!r} has no runner function (fn=None). "
            f"It must be handled specially by each pipeline."
        )

    # -- dispatch to runner ------------------------------------------------
    if fn is run_kalign:
        result = fn(
            input_fasta,
            ensemble=cfg.get("ensemble", 0),
            seq_type=seq_type,
            tgpo=cfg.get("tgpo"),
            terminal_dist_scale=cfg.get("terminal_dist_scale"),
            refine=cfg.get("refine"),
            vsm_amax=cfg.get("vsm_amax"),
            realign=cfg.get("realign"),
            consistency=cfg.get("consistency"),
            consistency_weight=cfg.get("consistency_weight"),
        )
    elif fn is run_guidance2:
        g2_type = "nuc" if seq_type in ("dna", "rna") else "aa"
        g2_dir = Path(work_dir) / f"guidance2_{method_name}"
        result = fn(
            input_fasta,
            output_dir=g2_dir,
            seq_type=g2_type,
            n_bootstrap=cfg.get("n_bootstrap", 100),
        )
    else:
        # mafft, muscle, clustalo -- simple signature
        result = fn(input_fasta, timeout=timeout)

    # -- optional site-weight output ---------------------------------------
    if cfg.get("weights") and result.column_confidence is not None:
        wt_path = Path(work_dir) / f"{method_name}_site_weights.txt"
        write_site_weights(result.column_confidence, wt_path)

    # -- optional confidence masking ---------------------------------------
    mask_threshold = cfg.get("mask")
    if not skip_masking and mask_threshold is not None and result.column_confidence is not None:
        masked_seqs, _n_kept = mask_alignment_by_confidence(
            result.sequences, result.column_confidence, mask_threshold
        )
        result = AlignResult(
            sequences=masked_seqs,
            names=result.names,
            column_confidence=result.column_confidence,
            residue_confidence=result.residue_confidence,
            wall_time=result.wall_time,
            peak_memory_mb=result.peak_memory_mb,
        )

    return result
