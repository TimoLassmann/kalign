"""
Integration tests that exercise real ecosystem packages.

These tests require the actual packages installed (biopython, scikit-bio).
They are skipped when the packages are not available, but exercised in CI
via the test_ecosystem job which installs kalign[all].

Biopython and scikit-bio test groups are independently skippable.
"""

import os

import pytest
import kalign


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _skip_no_biopython():
    try:
        import Bio  # noqa: F401
        return False
    except ImportError:
        return True

def _skip_no_skbio():
    try:
        import skbio  # noqa: F401
        return False
    except ImportError:
        return True


def _assert_alignment_valid(aligned, input_seqs):
    """Verify aligned sequences are a plausible alignment of the inputs."""
    assert len(aligned) == len(input_seqs), (
        f"Aligned count ({len(aligned)}) != input count ({len(input_seqs)})"
    )

    # All aligned sequences must be the same length
    lengths = {len(s) for s in aligned}
    assert len(lengths) == 1, f"Inconsistent alignment lengths: {lengths}"
    aln_len = lengths.pop()

    # Alignment must be at least as long as the longest input
    max_input = max(len(s) for s in input_seqs)
    assert aln_len >= max_input

    # If inputs differ in length, at least one aligned seq must contain gaps
    if len({len(s) for s in input_seqs}) > 1:
        assert any("-" in s for s in aligned), "Inputs differ in length but no gaps in alignment"

    # Ungapped content must match the original sequences
    for orig, aln in zip(input_seqs, aligned):
        assert aln.replace("-", "") == orig


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

DNA_SEQS = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]
DNA_IDS = ["seq1", "seq2", "seq3"]

RNA_SEQS = ["AUCGAUCGAUCG", "AUCGUCGAUCG", "AUCGAUCAUCG"]

PROTEIN_SEQS = [
    "MKTAYIAKQRQISFVK",
    "MKTAYIAKQRQ",
    "MKTAYIAK",
]

# Path to test FASTA data shipped with the project
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
TEST_FASTA = os.path.join(TEST_DATA_DIR, "BB11001.tfa")


# ---------------------------------------------------------------------------
# Biopython tests  (skipped independently if Bio is not installed)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(_skip_no_biopython(), reason="Biopython not installed")
class TestBiopythonFormat:
    """Test fmt='biopython' with real Biopython."""

    def test_align_biopython_format_dna(self):
        from Bio.Align import MultipleSeqAlignment

        result = kalign.align(DNA_SEQS, fmt="biopython", ids=DNA_IDS)
        assert isinstance(result, MultipleSeqAlignment)
        assert len(result) == len(DNA_SEQS)
        assert result[0].id == "seq1"
        assert result[1].id == "seq2"
        assert result[2].id == "seq3"

        # Verify actual alignment content
        aligned_strs = [str(r.seq) for r in result]
        _assert_alignment_valid(aligned_strs, DNA_SEQS)

    def test_align_biopython_format_protein(self):
        from Bio.Align import MultipleSeqAlignment

        result = kalign.align(
            PROTEIN_SEQS,
            seq_type="protein",
            fmt="biopython",
            ids=["p1", "p2", "p3"],
        )
        assert isinstance(result, MultipleSeqAlignment)
        assert len(result) == len(PROTEIN_SEQS)

        aligned_strs = [str(r.seq) for r in result]
        _assert_alignment_valid(aligned_strs, PROTEIN_SEQS)

    def test_align_biopython_auto_ids(self):
        from Bio.Align import MultipleSeqAlignment

        result = kalign.align(DNA_SEQS, fmt="biopython")
        assert isinstance(result, MultipleSeqAlignment)
        assert result[0].id == "seq0"


@pytest.mark.skipif(_skip_no_biopython(), reason="Biopython not installed")
class TestIORead:
    """Test kalign.io read functions with real Biopython."""

    def test_read_fasta(self):
        sequences = kalign.io.read_fasta(TEST_FASTA)
        assert isinstance(sequences, list)
        assert len(sequences) > 0
        assert all(isinstance(s, str) for s in sequences)
        assert all(len(s) > 0 for s in sequences)

    def test_read_sequences(self):
        sequences, ids = kalign.io.read_sequences(TEST_FASTA)
        assert len(sequences) == len(ids)
        assert len(sequences) > 0
        assert all(len(s) > 0 for s in sequences)
        assert all(len(i) > 0 for i in ids)


@pytest.mark.skipif(_skip_no_biopython(), reason="Biopython not installed")
class TestIOWrite:
    """Test kalign.io write functions with real Biopython."""

    @pytest.fixture
    def aligned(self):
        return kalign.align(DNA_SEQS)

    def test_write_fasta_roundtrip(self, aligned, tmp_path):
        out = tmp_path / "out.fasta"
        kalign.io.write_fasta(aligned, str(out), ids=DNA_IDS)

        read_back = kalign.io.read_fasta(str(out))
        assert read_back == aligned

    def test_write_clustal(self, aligned, tmp_path):
        from Bio import AlignIO

        out = tmp_path / "out.aln"
        kalign.io.write_clustal(aligned, str(out), ids=DNA_IDS)

        aln = AlignIO.read(str(out), "clustal")
        assert len(aln) == len(aligned)
        assert aln.get_alignment_length() == len(aligned[0])
        # Verify content matches
        for orig, record in zip(aligned, aln):
            assert str(record.seq) == orig

    def test_write_stockholm(self, aligned, tmp_path):
        from Bio import AlignIO

        out = tmp_path / "out.sto"
        kalign.io.write_stockholm(aligned, str(out), ids=DNA_IDS)

        aln = AlignIO.read(str(out), "stockholm")
        assert len(aln) == len(aligned)
        assert aln.get_alignment_length() == len(aligned[0])
        for orig, record in zip(aligned, aln):
            assert str(record.seq) == orig

    def test_write_phylip(self, aligned, tmp_path):
        from Bio import AlignIO

        out = tmp_path / "out.phy"
        kalign.io.write_phylip(aligned, str(out), ids=DNA_IDS)

        aln = AlignIO.read(str(out), "phylip-sequential")
        assert len(aln) == len(aligned)
        for orig, record in zip(aligned, aln):
            assert str(record.seq) == orig


# ---------------------------------------------------------------------------
# scikit-bio tests  (skipped independently if skbio is not installed)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(_skip_no_skbio(), reason="scikit-bio not installed")
class TestSkbioFormat:
    """Test fmt='skbio' with real scikit-bio."""

    def test_align_skbio_format_dna(self):
        import skbio

        result = kalign.align(DNA_SEQS, seq_type="dna", fmt="skbio", ids=DNA_IDS)
        assert isinstance(result, skbio.TabularMSA)
        assert len(result) == len(DNA_SEQS)
        assert all(isinstance(seq, skbio.DNA) for seq in result)

        # Verify alignment content
        aligned_strs = [str(seq) for seq in result]
        _assert_alignment_valid(aligned_strs, DNA_SEQS)

    def test_align_skbio_format_rna(self):
        import skbio

        result = kalign.align(RNA_SEQS, seq_type="rna", fmt="skbio")
        assert isinstance(result, skbio.TabularMSA)
        assert len(result) == len(RNA_SEQS)
        assert all(isinstance(seq, skbio.RNA) for seq in result)

        aligned_strs = [str(seq) for seq in result]
        _assert_alignment_valid(aligned_strs, RNA_SEQS)

    def test_align_skbio_format_protein(self):
        import skbio

        result = kalign.align(
            PROTEIN_SEQS, seq_type="protein", fmt="skbio", ids=["p1", "p2", "p3"]
        )
        assert isinstance(result, skbio.TabularMSA)
        assert len(result) == len(PROTEIN_SEQS)
        assert all(isinstance(seq, skbio.Protein) for seq in result)

        aligned_strs = [str(seq) for seq in result]
        _assert_alignment_valid(aligned_strs, PROTEIN_SEQS)

    def test_align_skbio_auto_detect_dna(self):
        """AUTO seq_type should infer DNA and produce skbio.DNA objects."""
        import skbio

        result = kalign.align(DNA_SEQS, fmt="skbio")
        assert all(isinstance(seq, skbio.DNA) for seq in result)

    def test_align_skbio_auto_detect_protein(self):
        """AUTO seq_type should infer Protein and produce skbio.Protein objects."""
        import skbio

        result = kalign.align(PROTEIN_SEQS, fmt="skbio")
        assert all(isinstance(seq, skbio.Protein) for seq in result)

    def test_align_skbio_metadata(self):
        import skbio

        result = kalign.align(DNA_SEQS, seq_type="dna", fmt="skbio", ids=DNA_IDS)
        assert result[0].metadata["id"] == "seq1"
        assert result[1].metadata["id"] == "seq2"
        assert result[2].metadata["id"] == "seq3"
