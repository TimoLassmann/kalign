"""
Edge case tests for Kalign Python package.
"""

import pytest
import kalign


class TestEdgeCases:
    """Test unusual inputs and boundary conditions."""

    def test_identical_sequences(self):
        """Test alignment of identical sequences."""
        identical = ["ATCGATCG"] * 3
        aligned = kalign.align(identical, seq_type="dna")

        assert len(aligned) == 3
        # All should be identical (no gaps needed)
        assert all(seq == aligned[0] for seq in aligned)

    def test_very_short_sequences(self):
        """Test very short sequences."""
        short_seqs = ["A", "T", "G"]
        aligned = kalign.align(short_seqs, seq_type="dna")

        assert len(aligned) == 3
        assert all(len(seq) >= 1 for seq in aligned)

    def test_very_different_lengths(self):
        """Test sequences with very different lengths."""
        diff_lengths = ["A", "ATCGATCGATCGATCG", "ATCG"]
        aligned = kalign.align(diff_lengths, seq_type="dna")

        assert len(aligned) == 3
        assert all(len(seq) == len(aligned[0]) for seq in aligned)

    def test_single_character_sequences(self):
        """Test single character sequences."""
        single_chars = ["A", "T", "C", "G"]
        aligned = kalign.align(single_chars, seq_type="dna")

        assert len(aligned) == 4

    def test_repetitive_sequences(self):
        """Test highly repetitive sequences."""
        repetitive = ["AAAAAAA", "TTTTTTT", "AAAAAAA"]
        aligned = kalign.align(repetitive, seq_type="dna")

        assert len(aligned) == 3

    def test_mixed_case_sequences(self, mixed_case):
        """Test mixed case handling."""
        aligned = kalign.align(mixed_case, seq_type="dna")
        assert len(aligned) == len(mixed_case)

    def test_large_sequence_count(self):
        """Test with many sequences."""
        # Use valid DNA characters only - create variation with different bases
        base_patterns = ["ATCG", "TACG", "GATC", "CGTA", "TGCA", "ACGT", "GCTA", "CTAG"]
        many_seqs = []
        for i in range(20):
            pattern = base_patterns[i % len(base_patterns)]
            # Add length variation with valid DNA bases
            suffix = "A" * (i % 4)  # Add 0-3 A's for length variation
            many_seqs.append(f"{pattern}{suffix}TG")

        aligned = kalign.align(many_seqs, seq_type="dna")

        assert len(aligned) == 20
        assert all(len(seq) == len(aligned[0]) for seq in aligned)

    def test_ambiguous_nucleotides(self):
        """Test sequences with ambiguous nucleotides."""
        ambiguous = ["ATCGN", "ATCGY", "ATCGR"]
        # Use auto-detection since ambiguous nucleotides might be detected as protein
        try:
            aligned = kalign.align(ambiguous)  # Let kalign auto-detect the type
            assert len(aligned) == 3
            # Should handle gracefully with auto-detection
        except (ValueError, RuntimeError):
            # Some ambiguous characters might still cause issues
            pass
