"""
Malformed input must raise a Python exception that says why, never crash the
interpreter. Each case runs in a subprocess so that a crash fails the test
instead of killing the test run.
"""

import subprocess
import sys
import textwrap

import pytest


def run_case(code):
    """Run `code` in a fresh interpreter; return (exit status, last output line)."""
    prog = textwrap.dedent(
        """
        import warnings
        warnings.simplefilter("ignore")
        import kalign
        try:
            result = CODE
            print("RETURNED", len(result))
        except Exception as e:
            print("RAISED", type(e).__name__, str(e))
        """
    ).replace("CODE", code)
    p = subprocess.run(
        [sys.executable, "-c", prog], capture_output=True, text=True, timeout=120
    )
    lines = p.stdout.strip().splitlines()
    return p.returncode, (lines[-1] if lines else "")


@pytest.fixture
def bad_files(tmp_path):
    files = {
        "empty": "",
        "header_only": ">a\n",
        "one": ">a\nACGTACGT\n",
        "no_header": "ACGTACGT\nACGTTCGT\n",
        "binary": None,
        "empty_seq": ">a\n\n>b\nACGTACGT\n>c\nACGTTCGT\n",
    }
    paths = {}
    for name, text in files.items():
        p = tmp_path / f"{name}.fa"
        if text is None:
            p.write_bytes(bytes(range(256)) * 8)
        else:
            p.write_text(text)
        paths[name] = str(p)
    return paths


@pytest.mark.parametrize(
    "name, message",
    [
        ("empty", "no sequences found"),
        ("header_only", "Only 1 sequence was found"),
        ("one", "Only 1 sequence was found"),
        ("no_header", "not FASTA, MSF or Clustal format"),
        ("binary", "not FASTA, MSF or Clustal format"),
    ],
)
@pytest.mark.parametrize("mode", ["fast", "accurate"])
def test_bad_file_raises_with_reason(bad_files, name, message, mode):
    rc, out = run_case(f"kalign.align_from_file({bad_files[name]!r}, mode={mode!r})")
    assert rc == 0, f"interpreter crashed (status {rc})"
    assert out.startswith("RAISED RuntimeError"), out
    assert message in out, out


@pytest.mark.parametrize("mode", ["fast", "accurate"])
def test_empty_sequence_in_file_is_dropped(bad_files, mode):
    rc, out = run_case(
        f"kalign.align_from_file({bad_files['empty_seq']!r}, mode={mode!r}).sequences"
    )
    assert rc == 0, f"interpreter crashed (status {rc})"
    assert out == "RETURNED 2", out


def test_file_to_file_bad_input_raises(bad_files, tmp_path):
    out_path = str(tmp_path / "out.fa")
    rc, out = run_case(f"kalign.align_file_to_file({bad_files['one']!r}, {out_path!r})")
    assert rc == 0, f"interpreter crashed (status {rc})"
    assert out.startswith("RAISED RuntimeError") and "Only 1 sequence" in out, out


def test_unwritable_output_raises(tmp_path):
    good = tmp_path / "ok.fa"
    good.write_text(">a\nACGTACGT\n>b\nACGTTCGT\n")
    bad_out = str(tmp_path / "no" / "such" / "dir" / "out.fa")
    rc, out = run_case(f"kalign.align_file_to_file({str(good)!r}, {bad_out!r})")
    assert rc == 0, f"interpreter crashed (status {rc})"
    assert out.startswith("RAISED") and "Could not open output file" in out, out


@pytest.mark.parametrize("mode", ["fast", "accurate"])
def test_non_ascii_sequence_raises(mode):
    rc, out = run_case(
        f"kalign.align(['ACGT\\u00c4CGT', 'ACGTACGT', 'ACGTTT'], mode={mode!r})"
    )
    assert rc == 0, f"interpreter crashed (status {rc})"
    assert out.startswith("RAISED RuntimeError") and "non-ASCII" in out, out
