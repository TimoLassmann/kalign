"""
Error handling and exception tests for Kalign Python package.
"""

import pytest
import kalign


class TestErrorHandling:
    """Test error conditions and exception handling."""

    def test_invalid_sequence_type_raises_error(self, dna_simple):
        """Test that invalid sequence type raises ValueError."""
        with pytest.raises(ValueError, match="Invalid seq_type"):
            kalign.align(dna_simple, seq_type="invalid_type")

    def test_zero_threads_raises_error(self, dna_simple):
        """Test that zero threads raises ValueError."""
        with pytest.raises(ValueError, match="at least 1"):
            kalign.align(dna_simple, n_threads=0)

    def test_negative_threads_raises_error(self, dna_simple):
        """Test that negative threads raises ValueError."""
        with pytest.raises(ValueError, match="at least 1"):
            kalign.align(dna_simple, n_threads=-1)

    def test_empty_sequence_in_list_raises_error(self):
        """Test that empty sequences raise ValueError."""
        with pytest.raises(
            ValueError, match="Sequence at index 1 is empty or contains only whitespace"
        ):
            kalign.align(["ATCG", "", "GCTA"])

    def test_non_string_sequence_raises_error(self):
        """Test that non-string sequences raise ValueError."""
        with pytest.raises(ValueError, match="must be strings"):
            kalign.align(["ATCG", 123, "GCTA"])

    def test_runtime_error_on_alignment_failure(self):
        """Test that alignment failures raise RuntimeError."""
        # Test with sequences containing invalid characters that should cause failure
        with pytest.raises((RuntimeError, ValueError)):
            kalign.align(["ATCGXYZ123", "ATCGXYZ456"], seq_type="dna")

    def test_file_not_found_error(self):
        """Test file not found error handling."""
        with pytest.raises(FileNotFoundError):
            kalign.align_from_file("/nonexistent/file.fasta")

    def test_invalid_file_format_error(self, temp_dir):
        """Test invalid file format handling."""
        invalid_file = temp_dir / "invalid.txt"
        invalid_file.write_text("This is not a sequence file")

        with pytest.raises(RuntimeError):
            kalign.align_from_file(str(invalid_file))
