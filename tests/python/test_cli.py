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


def test_cli_entry_point() -> None:
    """Verify kalign-py entry point is installed and runs."""
    result = subprocess.run(
        ["kalign-py", "--version"],
        capture_output=True,
        text=True,
    )
    # Entry point should exist and return version info
    assert result.returncode == 0
    assert result.stdout.strip()


_PROTEINS = {
    "a": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV",
    "b": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV",
    "c": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAPDAQFEVV",
    "d": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVA",
    "e": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKKALPDAQFEVV",
}


def _run_cli_stockholm(tmp_path: Path, output: str) -> str:
    fasta = tmp_path / "prot.fa"
    fasta.write_text("".join(f">{k}\n{v}\n" for k, v in _PROTEINS.items()))
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "kalign.cli",
            "-i",
            str(fasta),
            "-o",
            output,
            "--format",
            "stockholm",
            "--mode",
            "accurate",
            "-n",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def test_cli_stdout_keeps_confidence(tmp_path: Path) -> None:
    """Stockholm on stdout carries the same PP lines as a file."""
    out_file = tmp_path / "out.sto"
    _run_cli_stockholm(tmp_path, str(out_file))
    file_text = out_file.read_text()
    stdout_text = _run_cli_stockholm(tmp_path, "-")

    assert "#=GC PP_cons" in file_text
    assert "#=GC PP_cons" in stdout_text
    assert " PP " in stdout_text
    assert stdout_text.strip() == file_text.strip()


def test_cli_help_describes_modes_and_threads() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "kalign.cli", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    text = " ".join(result.stdout.split())
    assert "accurate" in text
    assert "5,000" in text
    assert "can only be used with fast" in text
    assert "default is 1" in text
