"""Tests for kalign.compare() and kalign.align_file_to_file()."""

import os
import tempfile

import pytest

import kalign

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def _data(name):
    return os.path.join(DATA_DIR, name)


class TestCompare:
    def test_identical_alignment_scores_100(self):
        ref = _data("BB11001.msf")
        score = kalign.compare(ref, ref)
        assert score == pytest.approx(100.0)

    def test_compare_different_alignments(self):
        ref = _data("BB11001.msf")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "test_aln.fasta")
            kalign.align_file_to_file(_data("BB11001.tfa"), out)
            score = kalign.compare(ref, out)
            assert 0.0 <= score <= 100.0

    def test_compare_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            kalign.compare("/nonexistent/ref.msf", "/nonexistent/test.msf")

    def test_compare_ref_not_found(self):
        ref = _data("BB11001.msf")
        with pytest.raises(FileNotFoundError):
            kalign.compare("/nonexistent/ref.msf", ref)

    def test_compare_test_not_found(self):
        ref = _data("BB11001.msf")
        with pytest.raises(FileNotFoundError):
            kalign.compare(ref, "/nonexistent/test.msf")


class TestAlignFileToFile:
    def test_basic_align_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "out.fasta")
            kalign.align_file_to_file(_data("BB11001.tfa"), out)
            assert os.path.exists(out)
            assert os.path.getsize(out) > 0

    def test_align_msf_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "out.msf")
            kalign.align_file_to_file(_data("BB11001.tfa"), out, format="msf")
            assert os.path.exists(out)
            assert os.path.getsize(out) > 0

    def test_align_and_compare(self):
        """Align BB11001 and compare to reference - SP score should be reasonable."""
        ref = _data("BB11001.msf")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "aln.msf")
            kalign.align_file_to_file(_data("BB11001.tfa"), out, format="msf")
            score = kalign.compare(ref, out)
            assert score > 0.0, "SP score should be positive"

    def test_input_not_found(self):
        with pytest.raises(FileNotFoundError):
            kalign.align_file_to_file("/nonexistent/input.fa", "/tmp/out.fa")

    def test_align_preserves_sequence_names(self):
        """Output should contain FASTA headers from input."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "out.fasta")
            kalign.align_file_to_file(_data("BB11001.tfa"), out)
            with open(out) as f:
                content = f.read()
            assert content.startswith(">"), "Output should be FASTA format"
            headers = [l for l in content.splitlines() if l.startswith(">")]
            assert len(headers) >= 2, "Should have multiple sequences"
