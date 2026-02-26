"""Tests for the unified mode interface (default/fast/precise)."""

import os
import tempfile

import pytest
import kalign


TEST_SEQUENCES = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]

PROTEIN_SEQUENCES = [
    "MKFLILLFNILCLFPVLAADNHGVSLHCTTATAIP",
    "MKFLILLFNILCLFPVLAADNHGVSLHCTTATAIP",
    "MKFLILLFNILCLFPVLAADNHGVSLHCTTATAIP",
]

TEST_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "BB11001.tfa")


class TestModeConstants:
    """Test mode constant exports."""

    def test_mode_constants_exist(self):
        assert kalign.MODE_DEFAULT == "default"
        assert kalign.MODE_FAST == "fast"
        assert kalign.MODE_PRECISE == "precise"


class TestAlignModes:
    """Test mode parameter on align()."""

    def test_default_mode_no_kwarg(self):
        """Default mode (no mode kwarg) produces alignment."""
        result = kalign.align(TEST_SEQUENCES)
        assert isinstance(result, list)
        assert len(result) == len(TEST_SEQUENCES)
        assert all(len(s) == len(result[0]) for s in result)

    def test_default_mode_explicit(self):
        """mode='default' produces alignment."""
        result = kalign.align(TEST_SEQUENCES, mode="default")
        assert isinstance(result, list)
        assert len(result) == len(TEST_SEQUENCES)

    def test_fast_mode(self):
        """mode='fast' produces alignment."""
        result = kalign.align(TEST_SEQUENCES, mode="fast")
        assert isinstance(result, list)
        assert len(result) == len(TEST_SEQUENCES)
        assert all(len(s) == len(result[0]) for s in result)

    def test_precise_mode(self):
        """mode='precise' produces alignment (ensemble)."""
        result = kalign.align(TEST_SEQUENCES, mode="precise")
        # precise uses ensemble, which returns (seqs, confidence) tuple
        if isinstance(result, tuple):
            seqs = result[0]
        else:
            seqs = result
        assert len(seqs) == len(TEST_SEQUENCES)
        assert all(len(s) == len(seqs[0]) for s in seqs)

    def test_invalid_mode(self):
        """Invalid mode raises ValueError."""
        with pytest.raises(ValueError, match="Invalid mode"):
            kalign.align(TEST_SEQUENCES, mode="turbo")

    def test_explicit_param_overrides_mode(self):
        """Explicit consistency=10 overrides fast mode default (consistency=0)."""
        # This should not crash — fast base + explicit consistency
        result = kalign.align(TEST_SEQUENCES, mode="fast", consistency=10)
        assert isinstance(result, list)
        assert len(result) == len(TEST_SEQUENCES)

    def test_mode_case_insensitive(self):
        """Mode names are case-insensitive."""
        result = kalign.align(TEST_SEQUENCES, mode="FAST")
        assert isinstance(result, list)
        assert len(result) == len(TEST_SEQUENCES)


@pytest.mark.skipif(not os.path.exists(TEST_FILE), reason="Test data not found")
class TestAlignFromFileModes:
    """Test mode parameter on align_from_file()."""

    def test_default_mode(self):
        result = kalign.align_from_file(TEST_FILE)
        names, sequences = result
        assert len(names) > 0
        assert len(sequences) == len(names)

    def test_fast_mode(self):
        result = kalign.align_from_file(TEST_FILE, mode="fast")
        names, sequences = result
        assert len(names) > 0

    def test_precise_mode(self):
        result = kalign.align_from_file(TEST_FILE, mode="precise")
        names, sequences = result
        assert len(names) > 0


@pytest.mark.skipif(not os.path.exists(TEST_FILE), reason="Test data not found")
class TestAlignFileToFileModes:
    """Test mode parameter on align_file_to_file()."""

    def test_default_mode(self):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f:
            out = f.name
        try:
            kalign.align_file_to_file(TEST_FILE, out)
            assert os.path.getsize(out) > 0
        finally:
            os.unlink(out)

    def test_fast_mode(self):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f:
            out = f.name
        try:
            kalign.align_file_to_file(TEST_FILE, out, mode="fast")
            assert os.path.getsize(out) > 0
        finally:
            os.unlink(out)

    def test_precise_mode(self):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f:
            out = f.name
        try:
            kalign.align_file_to_file(TEST_FILE, out, mode="precise")
            assert os.path.getsize(out) > 0
        finally:
            os.unlink(out)
