"""
Basic alignment functionality tests for Kalign Python package.
"""

import pytest
import sys
import kalign
from conftest import assert_valid_alignment, assert_alignment_preserves_characters


class TestBasicAlignment:
    """Test core alignment functionality."""

    def test_simple_dna_alignment(self, dna_simple):
        """Test basic DNA sequence alignment."""
        aligned = kalign.align(dna_simple, seq_type="dna")

        assert_valid_alignment(dna_simple, aligned)
        assert_alignment_preserves_characters(dna_simple, aligned)

        # Check that alignment adds gaps where needed
        assert len(aligned[0]) >= max(len(seq) for seq in dna_simple)

    def test_simple_protein_alignment(self, protein_simple):
        """Test basic protein sequence alignment."""
        aligned = kalign.align(protein_simple, seq_type="protein")

        assert_valid_alignment(protein_simple, aligned)
        assert_alignment_preserves_characters(protein_simple, aligned)

    def test_simple_rna_alignment(self, rna_simple):
        """Test basic RNA sequence alignment."""
        aligned = kalign.align(rna_simple, seq_type="rna")

        assert_valid_alignment(rna_simple, aligned)
        assert_alignment_preserves_characters(rna_simple, aligned)

    def test_auto_detection(self, dna_simple):
        """Test automatic sequence type detection."""
        aligned = kalign.align(dna_simple)  # No seq_type specified

        assert_valid_alignment(dna_simple, aligned)
        assert len(aligned) == len(dna_simple)

    def test_alignment_length_consistency(self, sequence_type_data):
        """Test that all aligned sequences have same length."""
        seq_type, sequences = sequence_type_data
        aligned = kalign.align(sequences, seq_type=seq_type)

        if aligned:
            expected_length = len(aligned[0])
            for seq in aligned:
                assert len(seq) == expected_length

    def test_empty_sequence_list(self):
        """Test behavior with empty sequence list."""
        with pytest.raises(ValueError, match="No sequences were found in the input"):
            kalign.align([])

    def test_single_sequence(self):
        """Test alignment with single sequence should raise error."""
        single_seq = ["ATCGATCG"]
        # Kalign requires at least 2 sequences for alignment
        with pytest.raises(
            ValueError,
            match="Only 1 sequence was found in the input - at least 2 sequences are required for alignment",
        ):
            kalign.align(single_seq, seq_type="dna")
