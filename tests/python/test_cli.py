import subprocess
import sys
from pathlib import Path


def test_cli_version() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "kalign.cli", "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip()


def test_cli_align_fasta_to_stdout() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    input_fasta = repo_root / "tests" / "data" / "tiny.fa"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "kalign.cli",
            "-i",
            str(input_fasta),
            "-o",
            "-",
            "--format",
            "fasta",
            "--type",
            "dna",
            "-n",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    # FASTA output should start with a header line.
    assert result.stdout.lstrip().startswith(">")
