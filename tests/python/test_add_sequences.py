"""Tests for adding sequences to an existing alignment."""
import pytest
import kalign
import os
import tempfile

TEST_FILE = os.path.join(
    os.path.dirname(__file__), "..", "data", "BB11001.tfa"
)


def _read_fasta(path):
    """Read a FASTA file, return list of (name, sequence) tuples."""
    entries = []
    name = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    entries.append((name, ''.join(seq_parts)))
                name = line.strip()[1:]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if name:
        entries.append((name, ''.join(seq_parts)))
    return entries


def _make_existing_and_new(test_file, n_holdout=1):
    """Align all sequences, then split into existing alignment + held-out sequences."""
    result = kalign.align_from_file(test_file, mode="fast")

    # Write existing alignment (all but last n_holdout sequences)
    existing_path = tempfile.mktemp(suffix=".fa")
    with open(existing_path, 'w') as f:
        for name, seq in zip(result.names[:-n_holdout], result.sequences[:-n_holdout]):
            f.write(f">{name}\n{seq}\n")

    # Write held-out sequences (unaligned, from original file)
    names_seqs = []
    with open(test_file) as f:
        name = None
        seq_parts = []
        for line in f:
            if line.startswith('>'):
                if name:
                    names_seqs.append((name, ''.join(seq_parts)))
                name = line.strip()[1:]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if name:
            names_seqs.append((name, ''.join(seq_parts)))

    new_path = tempfile.mktemp(suffix=".fa")
    with open(new_path, 'w') as f:
        for name, seq in names_seqs[-n_holdout:]:
            f.write(f">{name}\n{seq}\n")

    return existing_path, new_path, result


class TestAddSequences:

    def test_basic_add(self):
        """Basic add: 3 existing + 1 new = 4 sequences."""
        existing_path, new_path, full_result = _make_existing_and_new(TEST_FILE, 1)
        out_path = tempfile.mktemp(suffix=".fa")

        kalign.add_to_alignment(existing_path, new_path, out_path)

        result = _read_fasta(out_path)
        assert len(result) == len(full_result.sequences)

        # All sequences should have the same length
        lengths = set(len(seq) for _, seq in result)
        assert len(lengths) == 1, f"Unequal lengths: {lengths}"

        os.remove(existing_path)
        os.remove(new_path)
        os.remove(out_path)

    def test_existing_unchanged(self):
        """Existing sequences must not be modified (content, ignoring line wrapping)."""
        existing_path, new_path, _ = _make_existing_and_new(TEST_FILE, 1)
        out_path = tempfile.mktemp(suffix=".fa")

        existing_seqs = _read_fasta(existing_path)

        kalign.add_to_alignment(existing_path, new_path, out_path)

        output_seqs = _read_fasta(out_path)

        for i in range(len(existing_seqs)):
            assert existing_seqs[i][1] == output_seqs[i][1], \
                f"Existing sequence {i} was modified!\n  before: {existing_seqs[i][1][:60]}...\n  after:  {output_seqs[i][1][:60]}..."

        os.remove(existing_path)
        os.remove(new_path)
        os.remove(out_path)

    def test_residue_preservation(self):
        """Added sequences must preserve all residues (no residues lost)."""
        existing_path, new_path, _ = _make_existing_and_new(TEST_FILE, 1)
        out_path = tempfile.mktemp(suffix=".fa")

        # Read original new sequence
        with open(new_path) as f:
            lines = f.readlines()
        orig_seq = ''.join(l.strip() for l in lines if not l.startswith('>'))
        orig_residues = len(orig_seq)

        kalign.add_to_alignment(existing_path, new_path, out_path)

        # Read last sequence from output (the added one)
        result = _read_fasta(out_path)
        added_seq = result[-1][1]
        added_residues = sum(1 for c in added_seq if c != '-')

        assert added_residues == orig_residues, \
            f"Residue count changed: {orig_residues} -> {added_residues}"

        os.remove(existing_path)
        os.remove(new_path)
        os.remove(out_path)

    def test_alignment_length_matches(self):
        """Output alignment length must match existing alignment length."""
        existing_path, new_path, _ = _make_existing_and_new(TEST_FILE, 1)
        out_path = tempfile.mktemp(suffix=".fa")

        existing_seqs = _read_fasta(existing_path)
        existing_alnlen = len(existing_seqs[0][1])

        kalign.add_to_alignment(existing_path, new_path, out_path)

        output_seqs = _read_fasta(out_path)
        output_alnlen = len(output_seqs[0][1])

        assert output_alnlen == existing_alnlen, \
            f"Alignment length changed: {existing_alnlen} -> {output_alnlen}"

        os.remove(existing_path)
        os.remove(new_path)
        os.remove(out_path)

    def test_file_not_found(self):
        """Should raise FileNotFoundError for missing files."""
        with pytest.raises(FileNotFoundError):
            kalign.add_to_alignment("/nonexistent/file.fa", TEST_FILE, "/tmp/out.fa")
        with pytest.raises(FileNotFoundError):
            kalign.add_to_alignment(TEST_FILE, "/nonexistent/file.fa", "/tmp/out.fa")

    def test_larger_dataset(self):
        """Test with BB30014 (44 sequences, hold out 5)."""
        test_file = os.path.join(
            os.path.dirname(__file__), "..", "data", "BB30014.tfa"
        )
        if not os.path.exists(test_file):
            pytest.skip("BB30014.tfa not found")

        existing_path, new_path, full_result = _make_existing_and_new(test_file, 5)
        out_path = tempfile.mktemp(suffix=".fa")

        kalign.add_to_alignment(existing_path, new_path, out_path)

        result = _read_fasta(out_path)

        # Should have all 44 sequences
        assert len(result) == len(full_result.sequences)

        # All same length
        lengths = set(len(seq) for _, seq in result)
        assert len(lengths) == 1

        os.remove(existing_path)
        os.remove(new_path)
        os.remove(out_path)
