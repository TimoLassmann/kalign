"""
Tests for bioinformatics ecosystem integration features.

This module tests the Biopython and scikit-bio integration capabilities,
as well as threading controls and I/O helpers.
"""

import pytest
import sys
from unittest.mock import Mock, patch
import kalign


# Test data
TEST_SEQUENCES = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]

TEST_IDS = ["seq1", "seq2", "seq3"]


class TestReturnFormats:
    """Test different return format options."""

    def test_plain_format_default(self):
        """Test default plain format returns list of strings."""
        result = kalign.align(TEST_SEQUENCES)
        assert isinstance(result, list)
        assert all(isinstance(seq, str) for seq in result)
        assert len(result) == len(TEST_SEQUENCES)

    def test_plain_format_explicit(self):
        """Test explicit plain format specification."""
        result = kalign.align(TEST_SEQUENCES, fmt="plain")
        assert isinstance(result, list)
        assert all(isinstance(seq, str) for seq in result)

    @pytest.mark.skipif(
        "biopython" not in sys.modules, reason="Biopython not available"
    )
    def test_biopython_format_with_biopython(self):
        """Test Biopython format when Biopython is available."""
        try:
            from Bio.Align import MultipleSeqAlignment

            result = kalign.align(TEST_SEQUENCES, fmt="biopython", ids=TEST_IDS)
            assert isinstance(result, MultipleSeqAlignment)
            assert len(result) == len(TEST_SEQUENCES)
            assert result[0].id == TEST_IDS[0]
        except ImportError:
            pytest.skip("Biopython not available")

    def test_biopython_format_without_biopython(self):
        """Test Biopython format fails gracefully without Biopython."""
        with patch.dict(
            "sys.modules", {"Bio.Align": None, "Bio.SeqRecord": None, "Bio.Seq": None}
        ):
            with pytest.raises(ImportError, match="Biopython not installed"):
                kalign.align(TEST_SEQUENCES, fmt="biopython")

    @pytest.mark.skipif("skbio" not in sys.modules, reason="scikit-bio not available")
    def test_skbio_format_with_skbio(self):
        """Test scikit-bio format when scikit-bio is available."""
        try:
            import skbio

            result = kalign.align(TEST_SEQUENCES, fmt="skbio")
            assert isinstance(result, skbio.TabularMSA)
            assert len(result) == len(TEST_SEQUENCES)
        except ImportError:
            pytest.skip("scikit-bio not available")

    def test_skbio_format_without_skbio(self):
        """Test scikit-bio format fails gracefully without scikit-bio."""
        with patch.dict("sys.modules", {"skbio": None, "skbio.sequence": None}):
            with pytest.raises(ImportError, match="scikit-bio not installed"):
                kalign.align(TEST_SEQUENCES, fmt="skbio")

    def test_invalid_format(self):
        """Test invalid format specification raises error."""
        with pytest.raises(ValueError, match="Unknown fmt"):
            kalign.align(TEST_SEQUENCES, fmt="invalid")

    def test_ids_validation(self):
        """Test ID validation for ecosystem formats."""
        # Too few IDs
        with pytest.raises(ValueError, match="Number of IDs"):
            kalign.align(TEST_SEQUENCES, fmt="plain", ids=["seq1"])

        # Too many IDs
        with pytest.raises(ValueError, match="Number of IDs"):
            kalign.align(
                TEST_SEQUENCES, fmt="plain", ids=["seq1", "seq2", "seq3", "seq4"]
            )

    def test_auto_generated_ids(self):
        """Test automatic ID generation."""
        # This test checks that IDs are generated even if not explicitly tested
        # by verifying the function doesn't crash
        result = kalign.align(TEST_SEQUENCES, fmt="plain")
        assert len(result) == len(TEST_SEQUENCES)


class TestThreadingControl:
    """Test global threading control functions."""

    def test_default_thread_count(self):
        """Test default thread count is 1."""
        assert kalign.get_num_threads() == 1

    def test_set_get_threads(self):
        """Test setting and getting thread count."""
        original = kalign.get_num_threads()
        try:
            kalign.set_num_threads(4)
            assert kalign.get_num_threads() == 4

            kalign.set_num_threads(8)
            assert kalign.get_num_threads() == 8
        finally:
            kalign.set_num_threads(original)

    def test_invalid_thread_count(self):
        """Test invalid thread counts raise errors."""
        with pytest.raises(ValueError, match="must be at least 1"):
            kalign.set_num_threads(0)

        with pytest.raises(ValueError, match="must be at least 1"):
            kalign.set_num_threads(-1)

    def test_thread_local_storage(self):
        """Test that thread settings are thread-local."""
        import threading
        import time

        results = {}

        def worker(thread_id, num_threads):
            kalign.set_num_threads(num_threads)
            time.sleep(0.1)  # Allow other threads to set different values
            results[thread_id] = kalign.get_num_threads()

        # Start multiple threads with different thread counts
        threads = []
        for i in range(3):
            t = threading.Thread(target=worker, args=(i, (i + 1) * 2))
            threads.append(t)
            t.start()

        # Wait for all threads to complete
        for t in threads:
            t.join()

        # Each thread should have maintained its own setting
        assert results[0] == 2
        assert results[1] == 4
        assert results[2] == 6

    def test_align_uses_default_threads(self):
        """Test that align() uses default thread count when n_threads=None."""
        original = kalign.get_num_threads()
        try:
            kalign.set_num_threads(4)

            # This should use the default (4 threads)
            # We can't easily test this directly, but we can ensure it doesn't crash
            result = kalign.align(TEST_SEQUENCES, n_threads=None)
            assert len(result) == len(TEST_SEQUENCES)

        finally:
            kalign.set_num_threads(original)

    def test_align_override_threads(self):
        """Test that align() can override default thread count."""
        original = kalign.get_num_threads()
        try:
            kalign.set_num_threads(4)

            # This should override the default
            result = kalign.align(TEST_SEQUENCES, n_threads=8)
            assert len(result) == len(TEST_SEQUENCES)

        finally:
            kalign.set_num_threads(original)


class TestIOModule:
    """Test I/O helper functions."""

    def test_io_module_imported(self):
        """Test that I/O module is properly imported."""
        assert hasattr(kalign, "io")
        assert hasattr(kalign.io, "read_fasta")
        assert hasattr(kalign.io, "write_fasta")
        assert hasattr(kalign.io, "write_clustal")
        assert hasattr(kalign.io, "write_stockholm")

    def test_io_functions_require_biopython(self):
        """Test that I/O functions require Biopython and fail gracefully."""
        with patch.dict("sys.modules", {"Bio": None, "Bio.SeqIO": None}):
            with pytest.raises(ImportError, match="Biopython required"):
                kalign.io.read_fasta("test.fasta")


class TestUtilsModule:
    """Test utility functions."""

    def test_utils_module_imported(self):
        """Test that utils module is properly imported."""
        assert hasattr(kalign, "utils")
        assert hasattr(kalign.utils, "to_array")
        assert hasattr(kalign.utils, "alignment_stats")
        assert hasattr(kalign.utils, "consensus_sequence")

    def test_to_array(self):
        """Test conversion to NumPy array."""
        aligned = kalign.align(TEST_SEQUENCES)
        arr = kalign.utils.to_array(aligned)

        assert arr.shape[0] == len(aligned)
        assert arr.shape[1] == len(aligned[0])
        assert arr.dtype.kind == "U"  # Unicode string

    def test_to_array_empty(self):
        """Test to_array with empty alignment."""
        with pytest.raises(ValueError, match="Empty alignment"):
            kalign.utils.to_array([])

    def test_alignment_stats(self):
        """Test alignment statistics calculation."""
        aligned = kalign.align(TEST_SEQUENCES)
        stats = kalign.utils.alignment_stats(aligned)

        required_keys = {
            "length",
            "n_sequences",
            "gap_fraction",
            "conservation",
            "identity",
        }
        assert set(stats.keys()) == required_keys

        assert stats["n_sequences"] == len(aligned)
        assert stats["length"] == len(aligned[0])
        assert 0 <= stats["gap_fraction"] <= 1
        assert 0 <= stats["conservation"] <= 1
        assert 0 <= stats["identity"] <= 1

    def test_consensus_sequence(self):
        """Test consensus sequence generation."""
        aligned = kalign.align(TEST_SEQUENCES)
        consensus = kalign.utils.consensus_sequence(aligned)

        assert isinstance(consensus, str)
        assert len(consensus) == len(aligned[0])

    def test_consensus_sequence_empty(self):
        """Test consensus with empty alignment."""
        with pytest.raises(ValueError, match="Empty alignment"):
            kalign.utils.consensus_sequence([])

    def test_consensus_sequence_invalid_threshold(self):
        """Test consensus with invalid threshold."""
        aligned = kalign.align(TEST_SEQUENCES)

        with pytest.raises(ValueError, match="Threshold must be between 0 and 1"):
            kalign.utils.consensus_sequence(aligned, threshold=-0.1)

        with pytest.raises(ValueError, match="Threshold must be between 0 and 1"):
            kalign.utils.consensus_sequence(aligned, threshold=1.1)


class TestBackwardCompatibility:
    """Test that all changes maintain backward compatibility."""

    def test_existing_api_unchanged(self):
        """Test that existing API still works exactly as before."""
        # All these calls should work exactly as they did before
        result1 = kalign.align(TEST_SEQUENCES)
        result2 = kalign.align(TEST_SEQUENCES, seq_type="dna")
        result3 = kalign.align(TEST_SEQUENCES, seq_type="dna", n_threads=2)
        result4 = kalign.align(TEST_SEQUENCES, gap_open=-5.0, gap_extend=-1.0)

        # All should return lists of strings
        for result in [result1, result2, result3, result4]:
            assert isinstance(result, list)
            assert all(isinstance(seq, str) for seq in result)
            assert len(result) == len(TEST_SEQUENCES)

    def test_module_exports(self):
        """Test that all expected exports are available."""
        expected_exports = {
            "align",
            "align_from_file",
            "write_alignment",
            "generate_test_sequences",
            "set_num_threads",
            "get_num_threads",
            "kalign",
            "io",
            "utils",
            "DNA",
            "DNA_INTERNAL",
            "RNA",
            "PROTEIN",
            "PROTEIN_DIVERGENT",
            "AUTO",
            "__version__",
            "__author__",
            "__email__",
        }

        actual_exports = set(kalign.__all__)
        assert expected_exports == actual_exports

    def test_convenience_alias(self):
        """Test that the kalign convenience alias still works."""
        result1 = kalign.align(TEST_SEQUENCES)
        result2 = kalign.kalign(TEST_SEQUENCES)  # Alias

        # Results should be identical
        assert result1 == result2
