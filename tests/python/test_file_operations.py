"""
File I/O operation tests for Kalign Python package.
"""

import pytest
import kalign
from pathlib import Path


class TestFileOperations:
    """Test file-based alignment operations."""

    @pytest.mark.file_io
    def test_align_from_fasta_file(self, sample_fasta_file):
        """Test alignment from FASTA file."""
        result = kalign.align_from_file(sample_fasta_file, seq_type="dna")

        assert len(result.sequences) == 3
        assert len(result.names) == 3
        assert all(len(seq) == len(result.sequences[0]) for seq in result.sequences)

    @pytest.mark.file_io
    def test_align_protein_fasta_file(self, sample_protein_fasta_file):
        """Test protein alignment from FASTA file."""
        result = kalign.align_from_file(sample_protein_fasta_file, seq_type="protein")

        assert len(result.sequences) == 3
        assert len(result.names) == 3
        assert all(len(seq) == len(result.sequences[0]) for seq in result.sequences)

    @pytest.mark.file_io
    def test_auto_detect_from_file(self, sample_fasta_file):
        """Test auto-detection with file input."""
        result = kalign.align_from_file(sample_fasta_file)
        assert len(result.sequences) == 3
        assert len(result.names) == 3

    @pytest.mark.file_io
    def test_file_not_found(self, nonexistent_file):
        """Test handling of non-existent files."""
        with pytest.raises(FileNotFoundError):
            kalign.align_from_file(nonexistent_file)

    @pytest.mark.file_io
    def test_invalid_file_format(self, invalid_fasta_file):
        """Test handling of invalid file formats."""
        with pytest.raises(RuntimeError):
            kalign.align_from_file(invalid_fasta_file)

    @pytest.mark.file_io
    def test_file_with_threading(self, sample_fasta_file):
        """Test file operations with multiple threads."""
        result = kalign.align_from_file(sample_fasta_file, n_threads=2)
        assert len(result.sequences) == 3

    @pytest.mark.file_io
    def test_relative_file_path(self, test_data_dir):
        """Test relative file paths."""
        dna_file = test_data_dir / "dna_sequences.fasta"
        if dna_file.exists():
            result = kalign.align_from_file(str(dna_file))
            assert len(result.sequences) >= 1

    @pytest.mark.file_io
    def test_names_preserved(self, sample_fasta_file):
        """Test that sequence names from the file are preserved."""
        result = kalign.align_from_file(sample_fasta_file, seq_type="dna")
        assert result.names == ["seq1", "seq2", "seq3"]

    @pytest.mark.file_io
    def test_named_tuple_unpacking(self, sample_fasta_file):
        """Test that result can be unpacked as a tuple."""
        names, sequences = kalign.align_from_file(sample_fasta_file, seq_type="dna")
        assert len(names) == 3
        assert len(sequences) == 3
