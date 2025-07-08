"""
Parameter validation and behavior tests for Kalign Python package.
"""

import pytest
import kalign


class TestParameters:
    """Test parameter handling and validation."""

    def test_default_parameters(self, dna_simple):
        """Test that default parameters work."""
        aligned = kalign.align(dna_simple)
        assert len(aligned) == len(dna_simple)

    def test_custom_gap_penalties(self, dna_simple):
        """Test custom gap penalty parameters."""
        aligned = kalign.align(
            dna_simple,
            seq_type="dna",
            gap_open=-10.0,
            gap_extend=-1.0,
            terminal_gap_extend=0.0,
        )
        assert len(aligned) == len(dna_simple)

    def test_threading_parameters(self, dna_simple):
        """Test different thread counts."""
        for n_threads in [1, 2, 4]:
            aligned = kalign.align(dna_simple, n_threads=n_threads)
            assert len(aligned) == len(dna_simple)

    def test_none_parameters_use_defaults(self, dna_simple):
        """Test that None parameters use Kalign defaults."""
        aligned = kalign.align(
            dna_simple, gap_open=None, gap_extend=None, terminal_gap_extend=None
        )
        assert len(aligned) == len(dna_simple)

    @pytest.mark.parametrize("gap_penalty", [-20.0, -5.0, -1.0, 0.0])
    def test_gap_penalty_range(self, dna_simple, gap_penalty):
        """Test various gap penalty values."""
        aligned = kalign.align(
            dna_simple, gap_open=gap_penalty, gap_extend=gap_penalty / 2
        )
        assert len(aligned) == len(dna_simple)
