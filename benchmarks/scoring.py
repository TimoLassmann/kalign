"""Alignment + scoring wrappers for benchmark runs."""

import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

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
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


def align_with_python_api(
    case: BenchmarkCase, output: Path, n_threads: int = 1, refine: str = "none",
    adaptive_budget: bool = False,
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
    )
    return time.perf_counter() - start


def align_with_cli(
    case: BenchmarkCase,
    output: Path,
    binary: str = "kalign",
    n_threads: int = 1,
    refine: str = "none",
    adaptive_budget: bool = False,
) -> float:
    """Align using kalign C binary via subprocess. Returns wall time in seconds."""
    cmd = [binary, "-i", str(case.unaligned), "-f", "fasta", "-o", str(output)]
    if n_threads > 1:
        cmd.extend(["--nthreads", str(n_threads)])
    if refine != "none":
        cmd.extend(["--refine", refine])
    if adaptive_budget:
        cmd.append("--adaptive-budget")

    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start

    if result.returncode != 0:
        raise RuntimeError(
            f"kalign CLI failed (exit {result.returncode}): {result.stderr}"
        )

    return elapsed


def score_alignment(reference: Path, test_output: Path) -> float:
    """Score a test alignment against a reference. Returns SP score (0-100)."""
    return kalign.compare(str(reference), str(test_output))


def run_case(
    case: BenchmarkCase,
    method: str = "python_api",
    binary: str = "kalign",
    n_threads: int = 1,
    refine: str = "none",
    adaptive_budget: bool = False,
) -> AlignmentResult:
    """Run alignment + scoring for a single benchmark case."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / f"{case.family}_aln.fasta"

        try:
            if method == "python_api":
                wall_time = align_with_python_api(case, output, n_threads, refine, adaptive_budget)
            elif method == "cli":
                wall_time = align_with_cli(case, output, binary, n_threads, refine, adaptive_budget)
            else:
                raise ValueError(f"Unknown method: {method}")

            sp_score = score_alignment(case.reference, output)

            return AlignmentResult(
                family=case.family,
                dataset=case.dataset,
                method=method,
                sp_score=sp_score,
                wall_time=wall_time,
                seq_type=case.seq_type,
                refine=refine,
            )
        except Exception as e:
            return AlignmentResult(
                family=case.family,
                dataset=case.dataset,
                method=method,
                sp_score=0.0,
                wall_time=0.0,
                seq_type=case.seq_type,
                refine=refine,
                error=str(e),
            )
