"""
Input validation tests for Kalign Python package.
"""

import pytest
import kalign


class TestInputValidation:
    """Test input validation and error handling."""

    def test_empty_sequence_list(self):
        """Test empty sequence list raises ValueError."""
        with pytest.raises(ValueError, match="No sequences were found in the input"):
            kalign.align([])

    def test_empty_strings_in_list(self):
        """Test empty strings in sequence list."""
        with pytest.raises(
            ValueError, match="Sequence at index 1 is empty or contains only whitespace"
        ):
            kalign.align(["ATCG", "", "GCTA"])

    def test_whitespace_only_sequences(self):
        """Test sequences with only whitespace."""
        with pytest.raises(
            ValueError, match="Sequence at index 1 is empty or contains only whitespace"
        ):
            kalign.align(["ATCG", "   ", "GCTA"])

    def test_non_string_sequences(self):
        """Test non-string elements in sequence list."""
        with pytest.raises(ValueError, match="must be strings"):
            kalign.align(["ATCG", 123, "GCTA"])

    def test_none_in_sequence_list(self):
        """Test None values in sequence list."""
        with pytest.raises(ValueError, match="must be strings"):
            kalign.align(["ATCG", None, "GCTA"])

    def test_invalid_sequence_type_string(self):
        """Test invalid sequence type string."""
        with pytest.raises(ValueError, match="Invalid seq_type"):
            kalign.align(["ATCG", "GCTA"], seq_type="invalid")

    def test_invalid_thread_count_zero(self):
        """Test zero thread count."""
        with pytest.raises(ValueError, match="at least 1"):
            kalign.align(["ATCG", "GCTA"], n_threads=0)

    def test_invalid_thread_count_negative(self):
        """Test negative thread count."""
        with pytest.raises(ValueError, match="at least 1"):
            kalign.align(["ATCG", "GCTA"], n_threads=-1)

    def test_valid_thread_counts(self):
        """Test valid thread counts work."""
        sequences = ["ATCG", "GCTA"]

        # These should not raise errors
        for n_threads in [1, 2, 4]:
            aligned = kalign.align(sequences, n_threads=n_threads)
            assert len(aligned) == 2

    def test_case_insensitive_sequence_types(self):
        """Test case insensitive sequence type specification."""
        sequences = ["ATCG", "GCTA"]

        # These should all work
        for seq_type in ["DNA", "dna", "Dna", "DNA"]:
            aligned = kalign.align(sequences, seq_type=seq_type)
            assert len(aligned) == 2
