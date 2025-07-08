"""
Integration tests for Kalign Python package.
"""

import pytest
import kalign
from pathlib import Path


class TestIntegration:
    """Test end-to-end workflows and integration scenarios."""

    @pytest.mark.integration
    def test_full_dna_workflow(self, temp_dir):
        """Test complete DNA alignment workflow."""
        # Create test sequences
        sequences = ["ATCGATCGATCGATCG", "ATCGATCGTCGATCG", "ATCGATCGATCATCG"]

        # Align sequences
        aligned = kalign.align(sequences, seq_type="dna", n_threads=2)

        # Validate results
        assert len(aligned) == 3
        assert all(len(seq) == len(aligned[0]) for seq in aligned)

        # Check that original sequences are preserved (without gaps)
        for orig, align in zip(sequences, aligned):
            degapped = align.replace("-", "")
            assert degapped == orig

    @pytest.mark.integration
    def test_full_protein_workflow(self, sample_protein_fasta_file):
        """Test complete protein alignment workflow."""
        # Read and align from file
        aligned = kalign.align_from_file(
            sample_protein_fasta_file,
            seq_type="protein",
            gap_open=-10.0,
            gap_extend=-1.0,
        )

        assert len(aligned) == 3
        assert all(isinstance(seq, str) for seq in aligned)

    @pytest.mark.integration
    def test_mixed_parameter_workflow(self, dna_simple):
        """Test workflow with various parameter combinations."""
        results = []

        # Test different parameter combinations
        param_sets = [
            {"seq_type": "dna", "n_threads": 1},
            {"seq_type": "dna", "n_threads": 2, "gap_open": -8.0},
            {"seq_type": "auto", "gap_extend": -2.0},
        ]

        for params in param_sets:
            aligned = kalign.align(dna_simple, **params)
            results.append(aligned)
            assert len(aligned) == len(dna_simple)

        # All results should be valid (though may differ)
        for result in results:
            assert all(len(seq) == len(result[0]) for seq in result)

    @pytest.mark.integration
    def test_error_recovery(self, dna_simple):
        """Test that system recovers from errors properly."""
        # First, cause an error
        try:
            kalign.align([], seq_type="dna")
        except ValueError:
            pass  # Expected

        # Then verify normal operation still works
        aligned = kalign.align(dna_simple, seq_type="dna")
        assert len(aligned) == len(dna_simple)

    @pytest.mark.integration
    def test_consistency_across_calls(self, dna_simple):
        """Test that results are consistent across multiple calls."""
        aligned1 = kalign.align(dna_simple, seq_type="dna", n_threads=1)
        aligned2 = kalign.align(dna_simple, seq_type="dna", n_threads=1)

        # Results should be identical for same input/parameters
        assert aligned1 == aligned2

    @pytest.mark.integration
    def test_real_biological_sequences(self, test_data_dir):
        """Test with real biological sequences if available."""
        # Use the test data files we created
        dna_file = test_data_dir / "dna_sequences.fasta"
        protein_file = test_data_dir / "protein_sequences.fasta"

        if dna_file.exists():
            dna_aligned = kalign.align_from_file(str(dna_file), seq_type="dna")
            assert len(dna_aligned) >= 1

        if protein_file.exists():
            protein_aligned = kalign.align_from_file(
                str(protein_file), seq_type="protein"
            )
            assert len(protein_aligned) >= 1
