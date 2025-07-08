"""
Sequence type specific tests for Kalign Python package.
"""

import pytest
import kalign


class TestSequenceTypes:
    """Test different sequence types and their handling."""

    @pytest.mark.parametrize(
        "seq_type_str,seq_type_const",
        [
            ("auto", kalign.AUTO),
            ("dna", kalign.DNA),
            ("internal", kalign.DNA_INTERNAL),
        ],
    )
    def test_sequence_type_constants(self, seq_type_str, seq_type_const, dna_simple):
        """Test that string and constant sequence types give same results."""
        aligned_str = kalign.align(dna_simple, seq_type=seq_type_str)
        aligned_const = kalign.align(dna_simple, seq_type=seq_type_const)

        assert aligned_str == aligned_const

    @pytest.mark.parametrize(
        "seq_type_str,seq_type_const",
        [
            ("rna", kalign.RNA),
        ],
    )
    def test_rna_sequence_type_constants(
        self, seq_type_str, seq_type_const, rna_simple
    ):
        """Test RNA type constants with RNA sequences."""
        aligned_str = kalign.align(rna_simple, seq_type=seq_type_str)
        aligned_const = kalign.align(rna_simple, seq_type=seq_type_const)

        assert aligned_str == aligned_const

    @pytest.mark.parametrize(
        "seq_type_str,seq_type_const",
        [
            ("protein", kalign.PROTEIN),
            ("divergent", kalign.PROTEIN_DIVERGENT),
        ],
    )
    def test_protein_sequence_type_constants(
        self, seq_type_str, seq_type_const, protein_simple
    ):
        """Test protein type constants with protein sequences."""
        aligned_str = kalign.align(protein_simple, seq_type=seq_type_str)
        aligned_const = kalign.align(protein_simple, seq_type=seq_type_const)

        assert aligned_str == aligned_const

    def test_dna_alignment(self, dna_simple):
        """Test DNA-specific alignment."""
        aligned = kalign.align(dna_simple, seq_type="dna")

        assert len(aligned) == len(dna_simple)
        # Check that only valid DNA characters (including gaps) are present
        valid_chars = set("ATCG-")
        for seq in aligned:
            assert all(c.upper() in valid_chars for c in seq)

    def test_rna_alignment(self, rna_simple):
        """Test RNA-specific alignment."""
        aligned = kalign.align(rna_simple, seq_type="rna")

        assert len(aligned) == len(rna_simple)
        # Check that RNA characters are preserved
        for orig, align in zip(rna_simple, aligned):
            degapped = align.replace("-", "")
            assert degapped == orig

    def test_protein_alignment(self, protein_simple):
        """Test protein-specific alignment."""
        aligned = kalign.align(protein_simple, seq_type="protein")

        assert len(aligned) == len(protein_simple)
        # Check that amino acid characters are preserved
        for orig, align in zip(protein_simple, aligned):
            degapped = align.replace("-", "")
            assert degapped == orig

    def test_divergent_protein_type(self, protein_simple):
        """Test divergent protein sequence type."""
        aligned = kalign.align(protein_simple, seq_type="divergent")
        assert len(aligned) == len(protein_simple)

    def test_internal_dna_type(self, dna_simple):
        """Test internal DNA sequence type."""
        aligned = kalign.align(dna_simple, seq_type="internal")
        assert len(aligned) == len(dna_simple)
