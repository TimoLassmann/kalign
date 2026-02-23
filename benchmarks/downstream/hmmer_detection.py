"""Pipeline 3: HMMER Homology Detection benchmark.

Measures how alignment quality affects profile HMM sensitivity.  For each
Pfam seed family the pipeline strips gaps from the curated seed alignment,
re-aligns using each method, builds a profile HMM with ``hmmbuild``, and
searches a combined database with ``hmmsearch``.  True/false positive rates
are computed against known family memberships.
"""

from __future__ import annotations

import json
import logging
import os
import random
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pfam family list (hardcoded, well-known families)
# ---------------------------------------------------------------------------

PFAM_FAMILIES = [
    "PF00001",  # 7tm_1 (GPCR)
    "PF00002",  # 7tm_2
    "PF00003",  # 7tm_3
    "PF00004",  # AAA (ATPase)
    "PF00005",  # ABC_tran
    "PF00009",  # GTP_EFTU
    "PF00010",  # HLH
    "PF00012",  # HSP70
    "PF00013",  # KH domain
    "PF00014",  # Kunitz_BPTI
    "PF00018",  # SH3
    "PF00022",  # Actin
    "PF00023",  # Ank
    "PF00027",  # cNMP_binding
    "PF00028",  # Cadherin
    "PF00036",  # EF-hand
    "PF00041",  # Fibronectin_3
    "PF00042",  # Globin
    "PF00043",  # GST_C
    "PF00046",  # Homeodomain
    "PF00047",  # Ig
    "PF00048",  # IL8
    "PF00050",  # Kazal
    "PF00056",  # LDLA
    "PF00069",  # Pkinase
    "PF00071",  # Ras
    "PF00072",  # Response_reg
    "PF00076",  # RRM
    "PF00078",  # RVT
    "PF00079",  # Serpin
    "PF00081",  # Sod_Fe_N
    "PF00082",  # Subtilisin
    "PF00083",  # Sugar_tr
    "PF00084",  # Sushi
    "PF00085",  # Thioredoxin
    "PF00089",  # Trypsin
    "PF00096",  # Zinc_finger_C2H2
    "PF00100",  # Zona_pellucida
    "PF00102",  # Y_phosphatase
    "PF00104",  # Hormone_recep
    "PF00106",  # adh_short
    "PF00107",  # ADH_zinc_N
    "PF00109",  # ketoacyl-synt
    "PF00111",  # Fer2
    "PF00112",  # Pepsin
    "PF00115",  # COX1
    "PF00118",  # Cpn60_TCP1
    "PF00119",  # ATP-synt_A
    "PF00120",  # Glu_synthase
    "PF00121",  # TIM
]


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------


@dataclass
class HmmerCaseResult:
    """Result for a single family x method evaluation."""

    family_id: str
    method: str
    true_positives: int
    false_positives: int
    false_negatives: int
    sensitivity: float  # TP / (TP + FN)
    specificity: float  # precision = TP / (TP + FP)
    e_value_threshold: float
    hmm_length: int  # number of match states in the profile
    wall_time_align: float
    wall_time_hmmbuild: float
    wall_time_hmmsearch: float
    peak_memory_mb: float


@dataclass
class HmmerResult:
    """Aggregated result for the full HMMER detection pipeline."""

    provenance: dict
    cases: List[dict]
    summary: dict  # per-method aggregated metrics


# ---------------------------------------------------------------------------
# Pfam seed download / cache
# ---------------------------------------------------------------------------

_PFAM_SEED_URL_TEMPLATE = (
    "https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}"
    "?annotation=alignment:seed"
)


def _download_pfam_seeds(
    data_dir: Path,
    n_families: int = 50,
) -> List[dict]:
    """Download (or load from cache) Pfam seed alignments for benchmarking.

    Attempts to download seed alignments from the InterPro API for each
    family in :data:`PFAM_FAMILIES`.  Successfully downloaded alignments
    are cached as FASTA files under ``data_dir/pfam_seed/``.  If a
    download fails the cached version is used if available.

    Parameters
    ----------
    data_dir : Path
        Base data directory.  Seed files are stored in
        ``data_dir/pfam_seed/<accession>.fasta``.
    n_families : int
        Maximum number of families to return.

    Returns
    -------
    list[dict]
        Each dict has keys ``family_id``, ``seed_fasta_path``, ``n_members``.
    """
    import urllib.request
    import urllib.error

    seed_dir = Path(data_dir) / "pfam_seed"
    seed_dir.mkdir(parents=True, exist_ok=True)

    families_to_fetch = PFAM_FAMILIES[:n_families]
    results: List[dict] = []

    for accession in families_to_fetch:
        fasta_path = seed_dir / f"{accession}.fasta"

        # Try to use cache first
        if fasta_path.exists() and fasta_path.stat().st_size > 0:
            names, _seqs = _parse_fasta_simple(fasta_path)
            if names:
                results.append(
                    {
                        "family_id": accession,
                        "seed_fasta_path": str(fasta_path),
                        "n_members": len(names),
                    }
                )
                logger.debug("Using cached seed for %s (%d seqs)", accession, len(names))
                continue

        # Download from InterPro
        url = _PFAM_SEED_URL_TEMPLATE.format(accession=accession)
        logger.info("Downloading seed alignment for %s", accession)
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=30) as resp:
                raw = resp.read()

            # The InterPro API returns gzip-compressed Stockholm;
            # decompress if needed.
            import gzip as _gzip

            try:
                content = _gzip.decompress(raw)
            except (OSError, _gzip.BadGzipFile):
                content = raw

            # The API may return Stockholm or raw alignment text.
            # Convert to FASTA for uniform handling.
            text = content.decode("utf-8", errors="replace")
            names, seqs = _stockholm_or_text_to_fasta(text)

            if not names:
                logger.warning("No sequences parsed for %s", accession)
                continue

            # Strip gaps to get raw sequences, then write aligned FASTA
            # (we keep the aligned version since hmmbuild can use it)
            _write_fasta_simple(fasta_path, names, seqs)

            results.append(
                {
                    "family_id": accession,
                    "seed_fasta_path": str(fasta_path),
                    "n_members": len(names),
                }
            )
            logger.info("Downloaded seed for %s (%d seqs)", accession, len(names))

        except (urllib.error.URLError, OSError, ValueError) as exc:
            logger.warning("Failed to download %s: %s", accession, exc)
            # If there is no cache either, skip this family
            continue

    logger.info(
        "Loaded %d / %d Pfam seed families",
        len(results),
        len(families_to_fetch),
    )
    return results


# ---------------------------------------------------------------------------
# Stockholm / text parsing helpers
# ---------------------------------------------------------------------------


def _stockholm_or_text_to_fasta(text: str) -> tuple[list[str], list[str]]:
    """Parse Stockholm or aligned-text format into (names, sequences).

    Handles:
    - Stockholm format (lines with ``#=...`` annotations, ``//`` terminator)
    - Plain aligned FASTA
    - Clustal-like blocks
    """
    lines = text.splitlines()

    # Detect Stockholm
    if any(line.strip().startswith("# STOCKHOLM") for line in lines[:5]):
        return _parse_stockholm_text(lines)

    # Detect FASTA
    if any(line.strip().startswith(">") for line in lines[:10]):
        return _parse_fasta_lines(lines)

    # Fallback: treat as space-separated name+sequence blocks
    names: list[str] = []
    seqs: list[str] = []
    seq_dict: dict[str, list[str]] = {}
    name_order: list[str] = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            name, seq_part = parts
            if name not in seq_dict:
                seq_dict[name] = []
                name_order.append(name)
            seq_dict[name].append(seq_part.replace(" ", ""))

    for name in name_order:
        names.append(name)
        seqs.append("".join(seq_dict[name]))
    return names, seqs


def _parse_stockholm_text(lines: list[str]) -> tuple[list[str], list[str]]:
    """Parse Stockholm format lines into (names, sequences)."""
    seq_dict: dict[str, list[str]] = {}
    name_order: list[str] = []

    for line in lines:
        line = line.rstrip()
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            name, seq_part = parts
            if name not in seq_dict:
                seq_dict[name] = []
                name_order.append(name)
            seq_dict[name].append(seq_part.replace(" ", ""))

    names = []
    seqs = []
    for name in name_order:
        names.append(name)
        seqs.append("".join(seq_dict[name]))
    return names, seqs


def _parse_fasta_lines(lines: list[str]) -> tuple[list[str], list[str]]:
    """Parse FASTA from a list of lines."""
    names: list[str] = []
    seqs: list[str] = []
    current: list[str] = []
    for line in lines:
        line = line.rstrip()
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
                current = []
            name = line[1:].split()[0] if line[1:].strip() else ""
            names.append(name)
        else:
            current.append(line.strip())
    if current:
        seqs.append("".join(current))
    return names, seqs


def _parse_fasta_simple(path: Path) -> tuple[list[str], list[str]]:
    """Parse a FASTA file from disk.  Returns (names, sequences)."""
    names: list[str] = []
    seqs: list[str] = []
    current: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                    current = []
                names.append(line[1:].split()[0] if line[1:].strip() else "")
            else:
                current.append(line.strip())
    if current:
        seqs.append("".join(current))
    return names, seqs


def _write_fasta_simple(
    path: Path, names: list[str], sequences: list[str]
) -> None:
    """Write sequences in FASTA format."""
    with open(path, "w") as fh:
        for name, seq in zip(names, sequences):
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def _strip_gaps(sequences: list[str]) -> list[str]:
    """Remove gap characters from sequences."""
    return [s.replace("-", "").replace(".", "") for s in sequences]


# ---------------------------------------------------------------------------
# hmmbuild wrapper
# ---------------------------------------------------------------------------


def _run_hmmbuild(alignment_path: Path, hmm_path: Path) -> dict:
    """Run ``hmmbuild`` to construct a profile HMM from an alignment.

    Parameters
    ----------
    alignment_path : Path
        Input alignment (FASTA or Stockholm).
    hmm_path : Path
        Output HMM file path.

    Returns
    -------
    dict
        ``hmm_path``, ``hmm_length``, ``nseq``, ``wall_time``.

    Raises
    ------
    RuntimeError
        If ``hmmbuild`` is not found or exits with an error.
    """
    hmmbuild_bin = shutil.which("hmmbuild")
    if hmmbuild_bin is None:
        raise RuntimeError(
            "hmmbuild is not installed or not on PATH. "
            "Install HMMER3 (http://hmmer.org/)."
        )

    start = time.perf_counter()
    proc = subprocess.run(
        [hmmbuild_bin, "--amino", str(hmm_path), str(alignment_path)],
        capture_output=True,
        text=True,
    )
    wall = time.perf_counter() - start

    if proc.returncode != 0:
        raise RuntimeError(
            f"hmmbuild failed (exit {proc.returncode}): {proc.stderr}"
        )

    # Parse hmmbuild output for HMM length and nseq
    hmm_length = 0
    nseq = 0
    for line in proc.stdout.splitlines():
        # Example: "Number of sequences:  42"
        if "umber of sequences" in line:
            parts = line.split()
            for part in reversed(parts):
                try:
                    nseq = int(part)
                    break
                except ValueError:
                    continue
        # Example: "Model length:         285"
        # Also appears in HMM header: "LENG  285"
        if "odel length" in line or line.strip().startswith("LENG"):
            parts = line.split()
            for part in reversed(parts):
                try:
                    hmm_length = int(part)
                    break
                except ValueError:
                    continue

    # Fallback: parse the HMM file header for LENG
    if hmm_length == 0 and Path(hmm_path).exists():
        with open(hmm_path) as fh:
            for hline in fh:
                if hline.startswith("LENG"):
                    parts = hline.split()
                    if len(parts) >= 2:
                        try:
                            hmm_length = int(parts[1])
                        except ValueError:
                            pass
                    break

    return {
        "hmm_path": str(hmm_path),
        "hmm_length": hmm_length,
        "nseq": nseq,
        "wall_time": wall,
    }


# ---------------------------------------------------------------------------
# hmmsearch wrapper
# ---------------------------------------------------------------------------


def _run_hmmsearch(
    hmm_path: Path,
    target_db: Path,
    output_path: Path,
    e_value: float = 1e-5,
) -> dict:
    """Run ``hmmsearch`` to search a profile HMM against a sequence database.

    Parameters
    ----------
    hmm_path : Path
        Profile HMM file.
    target_db : Path
        FASTA database to search.
    output_path : Path
        Path for the ``--tblout`` tabular output.
    e_value : float
        E-value threshold for reporting hits.

    Returns
    -------
    dict
        ``hits`` (list of dicts with name, e_value, score),
        ``n_hits``, ``wall_time``.
    """
    hmmsearch_bin = shutil.which("hmmsearch")
    if hmmsearch_bin is None:
        raise RuntimeError(
            "hmmsearch is not installed or not on PATH. "
            "Install HMMER3 (http://hmmer.org/)."
        )

    start = time.perf_counter()
    proc = subprocess.run(
        [
            hmmsearch_bin,
            "--tblout",
            str(output_path),
            "--noali",
            "-E",
            str(e_value),
            str(hmm_path),
            str(target_db),
        ],
        capture_output=True,
        text=True,
    )
    wall = time.perf_counter() - start

    if proc.returncode != 0:
        raise RuntimeError(
            f"hmmsearch failed (exit {proc.returncode}): {proc.stderr}"
        )

    # Parse tblout
    hits = _parse_tblout(output_path)

    return {
        "hits": hits,
        "n_hits": len(hits),
        "wall_time": wall,
    }


def _parse_tblout(path: Path) -> List[dict]:
    """Parse HMMER3 ``--tblout`` tabular output.

    Columns (space-separated, first 6):
        target_name  accession  query_name  accession  E-value  score  bias ...

    Comment lines start with ``#``.

    Returns a list of dicts, each with ``name``, ``e_value``, ``score``.
    """
    hits: list[dict] = []
    if not Path(path).exists():
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                hit = {
                    "name": parts[0],
                    "e_value": float(parts[4]),
                    "score": float(parts[5]),
                }
                hits.append(hit)
            except (ValueError, IndexError):
                continue
    return hits


# ---------------------------------------------------------------------------
# Search database preparation
# ---------------------------------------------------------------------------


def _prepare_search_db(
    families: List[dict],
    data_dir: Path,
) -> Path:
    """Create a combined FASTA search database from family seed sequences.

    Pools seed sequences from all loaded families.  Each sequence is
    labelled ``<family_id>__<seq_name>`` so that true positives can be
    identified by family membership.

    Negative sequences are drawn from other families -- i.e. for any given
    query family, members of the remaining families act as negatives.

    Parameters
    ----------
    families : list[dict]
        Output of :func:`_download_pfam_seeds`.
    data_dir : Path
        Directory for writing the search DB.

    Returns
    -------
    Path
        Path to the combined FASTA database.
    """
    db_path = Path(data_dir) / "search_db.fasta"

    # If already exists and non-empty, check it covers all current families
    if db_path.exists() and db_path.stat().st_size > 0:
        needed = {f["family_id"] for f in families}
        present: set[str] = set()
        with open(db_path) as fh:
            for line in fh:
                if line.startswith(">"):
                    fam_id = line[1:].split("__")[0]
                    present.add(fam_id)
        if needed <= present:
            logger.info("Reusing existing search DB: %s (%d families)", db_path, len(present))
            return db_path
        logger.info(
            "Search DB has %d families but %d needed — rebuilding",
            len(present & needed), len(needed),
        )

    logger.info("Building combined search database from %d families", len(families))

    all_names: list[str] = []
    all_seqs: list[str] = []

    for fam in families:
        family_id = fam["family_id"]
        fasta_path = Path(fam["seed_fasta_path"])
        if not fasta_path.exists():
            continue
        names, seqs = _parse_fasta_simple(fasta_path)
        # Strip gaps to get unaligned sequences
        seqs = _strip_gaps(seqs)
        for name, seq in zip(names, seqs):
            # Label with family for ground truth tracking
            labelled_name = f"{family_id}__{name}"
            all_names.append(labelled_name)
            all_seqs.append(seq)

    if not all_names:
        raise RuntimeError("No sequences available to build search database")

    _write_fasta_simple(db_path, all_names, all_seqs)
    logger.info(
        "Search database written: %s (%d sequences)",
        db_path,
        len(all_names),
    )
    return db_path


def _get_family_members(search_db: Path, family_id: str) -> set[str]:
    """Return the set of sequence names in the search DB that belong to
    the given family (identified by the ``<family_id>__`` prefix).
    """
    members: set[str] = set()
    with open(search_db) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name.startswith(f"{family_id}__"):
                    members.add(name)
    return members


# ---------------------------------------------------------------------------
# Per-case pipeline
# ---------------------------------------------------------------------------


def run_hmmer_case(
    family: dict,
    method_name: str,
    search_db: Path,
    work_dir: Path,
    e_value_threshold: float = 1e-5,
    align_timeout: int | None = None,
) -> HmmerCaseResult:
    """Run the full HMMER pipeline for one family with one alignment method.

    Steps:
    1. Read Pfam seed sequences for the family
    2. Strip gaps to get unaligned sequences
    3. Re-align using the named method via ``utils.run_method()``
    4. Build profile HMM from the new alignment
    5. Search against the combined database
    6. Score: family members detected = TP, non-members detected = FP

    Parameters
    ----------
    family : dict
        Dict with ``family_id``, ``seed_fasta_path``, ``n_members``.
    method_name : str
        Alignment method name (key in ``METHODS``).
    search_db : Path
        Path to the combined search database FASTA.
    work_dir : Path
        Scratch directory for intermediate files.
    e_value_threshold : float
        E-value cutoff for hmmsearch.

    Returns
    -------
    HmmerCaseResult
    """
    from .utils import run_method, parse_fasta, write_fasta

    family_id = family["family_id"]
    seed_path = Path(family["seed_fasta_path"])
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Read seed alignment and strip gaps to get unaligned sequences
    names, seqs = _parse_fasta_simple(seed_path)
    unaligned_seqs = _strip_gaps(seqs)

    # Write unaligned FASTA for the aligner
    unaln_path = work_dir / f"{family_id}_{method_name}_unaligned.fasta"
    write_fasta(names, unaligned_seqs, unaln_path)

    # 2. Align (with timeout protection for slow external tools)
    try:
        aln_result = run_method(
            method_name=method_name,
            input_fasta=unaln_path,
            work_dir=work_dir,
            seq_type="protein",
            timeout=align_timeout,
        )
    except subprocess.TimeoutExpired:
        logger.warning(
            "Alignment timed out for %s/%s after %ss",
            family_id, method_name, align_timeout,
        )
        return HmmerCaseResult(
            family_id=family_id,
            method=method_name,
            true_positives=0,
            false_positives=0,
            false_negatives=family.get("n_members", 0),
            sensitivity=0.0,
            specificity=0.0,
            e_value_threshold=e_value_threshold,
            hmm_length=0,
            wall_time_align=float(align_timeout or 0),
            wall_time_hmmbuild=0.0,
            wall_time_hmmsearch=0.0,
            peak_memory_mb=0.0,
        )
    wall_time_align = aln_result.wall_time
    peak_memory_mb = aln_result.peak_memory_mb

    # Write alignment to FASTA for hmmbuild
    aln_path = work_dir / f"{family_id}_{method_name}_aligned.fasta"
    write_fasta(aln_result.names, aln_result.sequences, aln_path)

    # 3. Build HMM
    hmm_path = work_dir / f"{family_id}_{method_name}.hmm"
    try:
        hmm_info = _run_hmmbuild(aln_path, hmm_path)
    except RuntimeError as exc:
        logger.error("hmmbuild failed for %s/%s: %s", family_id, method_name, exc)
        return HmmerCaseResult(
            family_id=family_id,
            method=method_name,
            true_positives=0,
            false_positives=0,
            false_negatives=family.get("n_members", 0),
            sensitivity=0.0,
            specificity=0.0,
            e_value_threshold=e_value_threshold,
            hmm_length=0,
            wall_time_align=wall_time_align,
            wall_time_hmmbuild=0.0,
            wall_time_hmmsearch=0.0,
            peak_memory_mb=peak_memory_mb,
        )

    wall_time_hmmbuild = hmm_info["wall_time"]
    hmm_length = hmm_info["hmm_length"]

    # 4. Search
    tblout_path = work_dir / f"{family_id}_{method_name}_hits.tbl"
    try:
        search_info = _run_hmmsearch(
            hmm_path, search_db, tblout_path, e_value=e_value_threshold
        )
    except RuntimeError as exc:
        logger.error("hmmsearch failed for %s/%s: %s", family_id, method_name, exc)
        return HmmerCaseResult(
            family_id=family_id,
            method=method_name,
            true_positives=0,
            false_positives=0,
            false_negatives=family.get("n_members", 0),
            sensitivity=0.0,
            specificity=0.0,
            e_value_threshold=e_value_threshold,
            hmm_length=hmm_length,
            wall_time_align=wall_time_align,
            wall_time_hmmbuild=wall_time_hmmbuild,
            wall_time_hmmsearch=0.0,
            peak_memory_mb=peak_memory_mb,
        )

    wall_time_hmmsearch = search_info["wall_time"]

    # 5. Score against ground truth
    true_members = _get_family_members(search_db, family_id)
    detected_names = {hit["name"] for hit in search_info["hits"]}

    tp = len(detected_names & true_members)
    fp = len(detected_names - true_members)
    fn = len(true_members - detected_names)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0

    return HmmerCaseResult(
        family_id=family_id,
        method=method_name,
        true_positives=tp,
        false_positives=fp,
        false_negatives=fn,
        sensitivity=sensitivity,
        specificity=precision,  # using precision as the specificity metric
        e_value_threshold=e_value_threshold,
        hmm_length=hmm_length,
        wall_time_align=wall_time_align,
        wall_time_hmmbuild=wall_time_hmmbuild,
        wall_time_hmmsearch=wall_time_hmmsearch,
        peak_memory_mb=peak_memory_mb,
    )


# ---------------------------------------------------------------------------
# Worker function for parallel execution
# ---------------------------------------------------------------------------


def _run_one(args: dict) -> dict:
    """Worker function for a single family x method evaluation.

    Parameters
    ----------
    args : dict
        Must contain ``family``, ``method_name``, ``search_db``,
        ``work_dir``, ``e_value_threshold``.

    Returns
    -------
    dict
        Serialisable case result, or ``{"error": str}`` on failure.
    """
    from .utils import cache_load, cache_save, clean_work_dir

    fingerprint = args.get("fingerprint")
    work_dir = Path(args["work_dir"])

    if fingerprint:
        cached = cache_load(work_dir, fingerprint)
        if cached is not None:
            return cached

    # Remove stale output files from previous runs
    clean_work_dir(work_dir)

    try:
        result = run_hmmer_case(
            family=args["family"],
            method_name=args["method_name"],
            search_db=Path(args["search_db"]),
            work_dir=work_dir,
            e_value_threshold=args.get("e_value_threshold", 1e-5),
            align_timeout=args.get("align_timeout"),
        )
        result_dict = asdict(result)

        if fingerprint:
            cache_save(work_dir, fingerprint, result_dict)

        return result_dict
    except Exception as exc:
        logger.error(
            "Error in %s / %s: %s",
            args.get("family", {}).get("family_id", "?"),
            args.get("method_name", "?"),
            exc,
        )
        return {
            "error": str(exc),
            "family_id": args.get("family", {}).get("family_id", "unknown"),
            "method": args.get("method_name", "unknown"),
        }


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def _compute_summary(cases: List[dict], methods: List[str]) -> dict:
    """Compute per-method aggregated metrics from case results.

    Parameters
    ----------
    cases : list[dict]
        List of per-case result dicts (may include error entries).
    methods : list[str]
        Method names to summarise.

    Returns
    -------
    dict
        Keyed by method name, each value is a dict of aggregated metrics.
    """
    from .utils import bootstrap_ci, wilcoxon_paired, holm_bonferroni

    summary: dict = {}
    for method in methods:
        method_cases = [
            c for c in cases if c.get("method") == method and "error" not in c
        ]
        if not method_cases:
            summary[method] = {"n_cases": 0}
            continue

        sensitivities = [c["sensitivity"] for c in method_cases]
        specificities = [c["specificity"] for c in method_cases]

        n = len(method_cases)
        mean_sens = sum(sensitivities) / n if n > 0 else 0.0
        mean_spec = sum(specificities) / n if n > 0 else 0.0

        # Bootstrap CIs (only if we have enough data)
        sens_ci = bootstrap_ci(sensitivities) if n >= 3 else (mean_sens, mean_sens)
        spec_ci = bootstrap_ci(specificities) if n >= 3 else (mean_spec, mean_spec)

        total_tp = sum(c["true_positives"] for c in method_cases)
        total_fp = sum(c["false_positives"] for c in method_cases)
        total_fn = sum(c["false_negatives"] for c in method_cases)

        mean_align_time = (
            sum(c["wall_time_align"] for c in method_cases) / n if n > 0 else 0.0
        )
        mean_hmmbuild_time = (
            sum(c["wall_time_hmmbuild"] for c in method_cases) / n if n > 0 else 0.0
        )
        mean_hmmsearch_time = (
            sum(c["wall_time_hmmsearch"] for c in method_cases) / n if n > 0 else 0.0
        )

        summary[method] = {
            "n_cases": n,
            "mean_sensitivity": mean_sens,
            "sensitivity_ci_95": list(sens_ci),
            "mean_specificity": mean_spec,
            "specificity_ci_95": list(spec_ci),
            "total_tp": total_tp,
            "total_fp": total_fp,
            "total_fn": total_fn,
            "mean_wall_time_align": mean_align_time,
            "mean_wall_time_hmmbuild": mean_hmmbuild_time,
            "mean_wall_time_hmmsearch": mean_hmmsearch_time,
        }

    # Pairwise statistical tests: compare each method to "kalign"
    baseline = "kalign"
    baseline_cases = {
        c["family_id"]: c
        for c in cases
        if c.get("method") == baseline and "error" not in c
    }

    if baseline_cases:
        p_values = []
        comparisons = []
        for method in methods:
            if method == baseline:
                continue
            method_cases_map = {
                c["family_id"]: c
                for c in cases
                if c.get("method") == method and "error" not in c
            }
            # Find paired families
            shared_families = sorted(
                set(baseline_cases.keys()) & set(method_cases_map.keys())
            )
            if len(shared_families) >= 5:
                a_vals = [baseline_cases[f]["sensitivity"] for f in shared_families]
                b_vals = [method_cases_map[f]["sensitivity"] for f in shared_families]
                try:
                    test = wilcoxon_paired(a_vals, b_vals)
                    p_values.append(test["p_value"])
                    comparisons.append(
                        {
                            "comparison": f"{baseline}_vs_{method}",
                            "statistic": test["statistic"],
                            "p_value": test["p_value"],
                            "cliffs_delta": test["cliffs_delta"],
                            "n_paired": len(shared_families),
                        }
                    )
                except Exception:
                    pass

        if p_values:
            adjusted = holm_bonferroni(p_values)
            for comp, adj_p in zip(comparisons, adjusted):
                comp["p_value_adjusted"] = adj_p

        summary["_statistical_tests"] = comparisons

    return summary


# ---------------------------------------------------------------------------
# Main pipeline entry point
# ---------------------------------------------------------------------------


def run_pipeline(params: dict) -> HmmerResult:
    """Run the full HMMER homology detection benchmark.

    Parameters
    ----------
    params : dict
        Configuration with keys:

        - ``data_dir`` : str or Path -- base data directory
        - ``results_dir`` : str or Path -- where to save JSON output
        - ``methods`` : list[str] -- alignment methods to test
        - ``n_families`` : int (default 50) -- number of Pfam families
        - ``n_jobs`` : int (default 1) -- parallel workers
        - ``quick`` : bool (default False) -- use only 5 families
        - ``e_value_threshold`` : float (default 1e-5)

    Returns
    -------
    HmmerResult
    """
    from .provenance import collect_provenance, result_path, update_latest_symlink
    from .utils import tool_versions_fingerprint

    data_dir = Path(params.get("data_dir", "benchmarks/data/downloads"))
    results_dir = Path(params.get("results_dir", "benchmarks/results"))
    methods = params.get("methods", ["kalign", "kalign_ens3", "kalign_ens3_m70"])
    n_families = params.get("n_families", 50)
    n_jobs = params.get("n_jobs", 1)
    quick = params.get("quick", False)
    e_value_threshold = params.get("e_value_threshold", 1e-5)
    align_timeout = params.get("align_timeout", 300)  # 5 min per alignment

    if quick:
        n_families = min(n_families, 5)
        logger.info("Quick mode: using %d families", n_families)

    # Collect provenance
    provenance = collect_provenance(params)

    # Download / load Pfam seeds
    families = _download_pfam_seeds(data_dir, n_families=n_families)
    if not families:
        raise RuntimeError(
            "No Pfam seed families available. "
            "Check network connectivity or populate the cache manually."
        )

    # Prepare search database
    search_db = _prepare_search_db(families, data_dir)

    # Compute fingerprint once for all workers (None disables caching)
    if params.get("no_cache"):
        fingerprint = None
        logger.info("Caching disabled (--no-cache)")
    else:
        fingerprint = tool_versions_fingerprint()
        logger.info("Tool versions fingerprint: %s", fingerprint)

    # Build work items
    work_items: list[dict] = []
    for family in families:
        for method_name in methods:
            case_work_dir = (
                Path(data_dir)
                / "hmmer_work"
                / family["family_id"]
                / method_name
            )
            work_items.append(
                {
                    "family": family,
                    "method_name": method_name,
                    "search_db": str(search_db),
                    "work_dir": str(case_work_dir),
                    "e_value_threshold": e_value_threshold,
                    "fingerprint": fingerprint,
                    "align_timeout": align_timeout,
                }
            )

    logger.info(
        "Running %d cases (%d families x %d methods) with %d workers",
        len(work_items),
        len(families),
        len(methods),
        n_jobs,
    )

    # Execute
    from tqdm import tqdm

    all_cases: list[dict] = []
    pbar = tqdm(total=len(work_items), desc="HMMER detection", unit="case")
    if n_jobs <= 1:
        for item in work_items:
            result = _run_one(item)
            all_cases.append(result)
            pbar.update(1)
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {
                executor.submit(_run_one, item): item for item in work_items
            }
            for future in as_completed(futures):
                all_cases.append(future.result())
                pbar.update(1)
    pbar.close()

    # Aggregate
    summary = _compute_summary(all_cases, methods)

    # Log summary table
    logger.info("=== HMMER Detection Summary ===")
    for method in methods:
        s = summary.get(method, {})
        if s.get("n_cases", 0) > 0:
            logger.info(
                "  %-20s  sens=%.3f  prec=%.3f  (n=%d)",
                method,
                s["mean_sensitivity"],
                s["mean_specificity"],
                s["n_cases"],
            )

    # Save results
    provenance_dict = asdict(provenance)
    hmmer_result = HmmerResult(
        provenance=provenance_dict,
        cases=all_cases,
        summary=summary,
    )

    out_path = result_path(results_dir, "hmmer_detection")
    with open(out_path, "w") as fh:
        json.dump(asdict(hmmer_result), fh, indent=2, default=str)
    update_latest_symlink(out_path)
    logger.info("Results saved to %s", out_path)

    return hmmer_result


# ---------------------------------------------------------------------------
# Load results
# ---------------------------------------------------------------------------


def load_results(results_dir: Path) -> HmmerResult:
    """Load HMMER detection results from the latest result file.

    Parameters
    ----------
    results_dir : Path
        Base results directory (the ``hmmer_detection/`` subdirectory
        is appended automatically).

    Returns
    -------
    HmmerResult
    """
    from .provenance import load_latest_results

    data = load_latest_results(Path(results_dir) / "hmmer_detection")
    return HmmerResult(
        provenance=data.get("provenance", {}),
        cases=data.get("cases", []),
        summary=data.get("summary", {}),
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Pipeline 3: HMMER Homology Detection Benchmark"
    )
    parser.add_argument(
        "-j",
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)",
    )
    parser.add_argument(
        "--max-families",
        type=int,
        default=50,
        help="Maximum number of Pfam families to test (default: 50)",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="benchmarks/data/downloads",
        help="Data directory for downloads and caches",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="benchmarks/results",
        help="Results output directory",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick mode: only 5 families",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=None,
        help="Alignment methods to test (default: all)",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    from .utils import METHODS

    methods = args.methods or list(METHODS.keys())

    run_pipeline(
        {
            "data_dir": args.data_dir,
            "results_dir": args.results_dir,
            "methods": methods,
            "n_families": args.max_families,
            "n_jobs": args.parallel,
            "quick": args.quick,
        }
    )
