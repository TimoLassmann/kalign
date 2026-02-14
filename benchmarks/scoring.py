"""Alignment + scoring wrappers for benchmark runs."""

import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional

import kalign

from .datasets import BenchmarkCase


@dataclass
class AlignmentResult:
    """Result of aligning and scoring one benchmark case."""

    family: str
    dataset: str
    method: str
    sp_score: float
    wall_time: float
    seq_type: str
    refine: str = "none"
    ensemble: int = 0
    recall: float = 0.0
    precision: float = 0.0
    f1: float = 0.0
    tc: float = 0.0
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


def align_with_python_api(
    case: BenchmarkCase, output: Path, n_threads: int = 1, refine: str = "none",
    adaptive_budget: bool = False, ensemble: int = 0,
) -> float:
    """Align using kalign Python API. Returns wall time in seconds."""
    start = time.perf_counter()
    kalign.align_file_to_file(
        str(case.unaligned),
        str(output),
        format="fasta",
        seq_type=case.seq_type,
        n_threads=n_threads,
        refine=refine,
        adaptive_budget=adaptive_budget,
        ensemble=ensemble,
    )
    return time.perf_counter() - start


def align_with_cli(
    case: BenchmarkCase,
    output: Path,
    binary: str = "kalign",
    n_threads: int = 1,
    refine: str = "none",
    adaptive_budget: bool = False,
    ensemble: int = 0,
) -> float:
    """Align using kalign C binary via subprocess. Returns wall time in seconds."""
    cmd = [binary, "-i", str(case.unaligned), "-f", "fasta", "-o", str(output)]
    if n_threads > 1:
        cmd.extend(["--nthreads", str(n_threads)])
    if refine != "none":
        cmd.extend(["--refine", refine])
    if adaptive_budget:
        cmd.append("--adaptive-budget")
    if ensemble > 0:
        cmd.extend(["--ensemble", str(ensemble)])

    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start

    if result.returncode != 0:
        raise RuntimeError(
            f"kalign CLI failed (exit {result.returncode}): {result.stderr}"
        )

    return elapsed


EXTERNAL_TOOLS = {"clustalo", "mafft", "muscle"}


def align_with_external(
    case: BenchmarkCase, output: Path, tool: str, n_threads: int = 1,
) -> float:
    """Align using an external MSA tool. Returns wall time in seconds."""
    if tool == "clustalo":
        cmd = ["clustalo", "-i", str(case.unaligned), "-o", str(output),
               "--outfmt=fasta", "--force"]
        if n_threads > 1:
            cmd.extend(["--threads", str(n_threads)])
    elif tool == "mafft":
        cmd = ["mafft", "--auto"]
        if n_threads > 1:
            cmd.extend(["--thread", str(n_threads)])
        cmd.append(str(case.unaligned))
        # MAFFT writes to stdout
        start = time.perf_counter()
        with open(output, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
        elapsed = time.perf_counter() - start
        if result.returncode != 0:
            raise RuntimeError(
                f"mafft failed (exit {result.returncode}): {result.stderr}"
            )
        return elapsed
    elif tool == "muscle":
        cmd = ["muscle", "-align", str(case.unaligned), "-output", str(output)]
        if n_threads > 1:
            cmd.extend(["--threads", str(n_threads)])
    else:
        raise ValueError(f"Unknown external tool: {tool}")

    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start

    if result.returncode != 0:
        raise RuntimeError(
            f"{tool} failed (exit {result.returncode}): {result.stderr}"
        )

    return elapsed


def score_alignment(reference: Path, test_output: Path) -> float:
    """Score a test alignment against a reference. Returns SP score (0-100)."""
    return kalign.compare(str(reference), str(test_output))


def parse_balibase_xml(xml_path: Path) -> List[int]:
    """Parse BAliBASE XML annotation to extract core block column mask.

    The XML ``<column-score>/<colsco-data>`` element contains per-column
    integers: ``-1`` (flanking), ``0`` (non-core), ``1`` (core block).
    Returns a binary mask where only core columns (value == 1) are set.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()
    colsco = root.find(".//column-score/colsco-data")
    if colsco is None or colsco.text is None:
        raise ValueError(f"No <colsco-data> element found in {xml_path}")
    values = [int(v) for v in colsco.text.split()]
    return [1 if v == 1 else 0 for v in values]


def _fasta_ref_has_gaps(reference: Path) -> bool:
    """Check if a FASTA reference alignment contains any gap characters."""
    with open(reference) as f:
        for line in f:
            if not line.startswith('>') and '-' in line:
                return True
    return False


def score_alignment_detailed(reference: Path, test_output: Path) -> dict:
    """Score a test alignment with POAR metrics.

    If a BAliBASE XML annotation file exists alongside the reference
    (same name with .xml extension), the XML core block mask is used.
    Otherwise scores all columns (max_gap_frac=-1.0), appropriate for
    BRAliBASE and other benchmarks with hand-curated references.

    Gapless FASTA reference alignments (e.g. conserved RNA with no indels)
    are skipped — the C comparison code requires alnlen > 0.
    """
    xml_path = reference.with_suffix('.xml')
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(reference), str(test_output), column_mask=mask)

    # For non-BAliBASE references (FASTA format), check for gapless edge case
    if reference.suffix in ('.fa', '.fasta') and not _fasta_ref_has_gaps(reference):
        raise RuntimeError("Reference alignment has no gaps — skipping (trivially aligned)")

    return kalign.compare_detailed(str(reference), str(test_output), max_gap_frac=-1.0)


def run_case(
    case: BenchmarkCase,
    method: str = "python_api",
    binary: str = "kalign",
    n_threads: int = 1,
    refine: str = "none",
    adaptive_budget: bool = False,
    ensemble: int = 0,
) -> AlignmentResult:
    """Run alignment + scoring for a single benchmark case."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / f"{case.family}_aln.fasta"

        try:
            if method == "python_api":
                wall_time = align_with_python_api(case, output, n_threads, refine, adaptive_budget, ensemble)
            elif method == "cli":
                wall_time = align_with_cli(case, output, binary, n_threads, refine, adaptive_budget, ensemble)
            elif method in EXTERNAL_TOOLS:
                wall_time = align_with_external(case, output, method, n_threads)
            else:
                raise ValueError(f"Unknown method: {method}")

            sp_score = score_alignment(case.reference, output)
            detailed = score_alignment_detailed(case.reference, output)

            is_external = method in EXTERNAL_TOOLS
            return AlignmentResult(
                family=case.family,
                dataset=case.dataset,
                method=method,
                sp_score=sp_score,
                wall_time=wall_time,
                seq_type=case.seq_type,
                refine="n/a" if is_external else refine,
                ensemble=0 if is_external else ensemble,
                recall=detailed["recall"],
                precision=detailed["precision"],
                f1=detailed["f1"],
                tc=detailed["tc"],
            )
        except Exception as e:
            is_external = method in EXTERNAL_TOOLS
            return AlignmentResult(
                family=case.family,
                dataset=case.dataset,
                method=method,
                sp_score=0.0,
                wall_time=0.0,
                seq_type=case.seq_type,
                refine="n/a" if is_external else refine,
                ensemble=0 if is_external else ensemble,
                error=str(e),
            )
