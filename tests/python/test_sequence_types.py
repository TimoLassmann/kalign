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

    @pytest.mark.parametrize(
        "bad", ["internal", "divergent", "pfasum43", "pfasum60", "pfasum", 2, 3, 4, 6, 12, True]
    )
    def test_removed_seq_types_rejected(self, bad, dna_simple, tmp_path):
        """Only auto/protein/dna/rna (as the C --type) are accepted, everywhere."""
        fa = tmp_path / "in.fa"
        fa.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(dna_simple)))
        with pytest.raises(ValueError, match="Invalid seq_type"):
            kalign.align(dna_simple, seq_type=bad)
        with pytest.raises(ValueError, match="Invalid seq_type"):
            kalign.align_from_file(str(fa), seq_type=bad)
        with pytest.raises(ValueError, match="Invalid seq_type"):
            kalign.align_file_to_file(str(fa), str(tmp_path / "out.fa"), seq_type=bad)


def _nuc_family(n=30, length=80, seed=1):
    import random

    rng = random.Random(seed)
    root = [rng.choice("ACGT") for _ in range(length)]
    out = []
    for _ in range(n):
        s = []
        for c in root:
            r = rng.random()
            if r < 0.05:
                continue
            s.append(rng.choice("ACGT") if r < 0.2 else c)
        out.append("".join(s))
    return out


class TestSeqTypeOverride:
    """seq_type overrides detection and selects the protein or nucleotide presets."""

    @pytest.mark.parametrize("mode", ["fast", "accurate"])
    def test_override_warns_and_changes_presets(self, mode):
        seqs = _nuc_family()
        auto = kalign.align(seqs, mode=mode, n_threads=1)
        with pytest.warns(UserWarning, match="aligning as protein"):
            forced = kalign.align(seqs, seq_type="protein", mode=mode, n_threads=1)
        auto = auto[0] if isinstance(auto, tuple) else auto
        forced = forced[0] if isinstance(forced, tuple) else forced
        assert forced != auto
        for orig, a in zip(seqs, forced):
            assert a.replace("-", "") == orig

    def test_matching_type_is_silent_and_identical(self):
        import warnings

        seqs = _nuc_family()
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            same = kalign.align(seqs, seq_type="dna", n_threads=1)
            rna = kalign.align(seqs, seq_type="rna", n_threads=1)
        assert same == kalign.align(seqs, n_threads=1) == rna

    def test_add_to_alignment_seq_type(self, tmp_path):
        data = __import__("pathlib").Path(__file__).parents[1] / "data"
        out_auto, out_dna = tmp_path / "a.fa", tmp_path / "d.fa"
        kalign.add_to_alignment(str(data / "BB11001.msf"), str(data / "BB12006.tfa"),
                                str(out_auto), n_threads=1)
        with pytest.warns(UserWarning, match="aligning as DNA/RNA"):
            kalign.add_to_alignment(str(data / "BB11001.msf"), str(data / "BB12006.tfa"),
                                    str(out_dna), n_threads=1, seq_type="dna")
        assert out_auto.read_text() != out_dna.read_text()
        with pytest.raises(ValueError, match="Invalid seq_type"):
            kalign.add_to_alignment(str(data / "BB11001.msf"), str(data / "BB12006.tfa"),
                                    str(out_dna), seq_type="internal")


class TestAccurateSizeLimit:
    """accurate refuses more than 5,000 sequences for protein and nucleotide alike."""

    @pytest.mark.parametrize("seq,seq_type", [("ACGTACGTAC", "auto"), ("MKTAYIAKQR", "auto"),
                                              ("ACGTACGTAC", "protein"), ("MKTAYIAKQR", "dna")])
    def test_limit(self, seq, seq_type, tmp_path):
        import warnings

        fa = tmp_path / "big.fa"
        fa.write_text("".join(f">s{i}\n{seq}\n" for i in range(5001)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ValueError, match="not feasible above 5000"):
                kalign.align_from_file(str(fa), mode="accurate", seq_type=seq_type)
            with pytest.raises(ValueError, match="not feasible above 5000"):
                kalign.align([seq] * 5001, mode="accurate", seq_type=seq_type)
