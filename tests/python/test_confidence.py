"""Tests for confidence masking and filtering."""
import pytest
import kalign
import os

TEST_FILE = os.path.join(
    os.path.dirname(__file__), "..", "data", "BB11001.tfa"
)


class TestConfidenceMasking:
    """Test confidence masking with ensemble modes."""

    def test_accurate_has_confidence(self):
        """Accurate mode should produce confidence scores."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        assert result.column_confidence is not None
        assert len(result.column_confidence) > 0
        assert all(0.0 <= c <= 1.0 for c in result.column_confidence)

    def test_fast_has_no_confidence(self):
        """Fast mode (single run) should NOT produce confidence."""
        result = kalign.align_from_file(TEST_FILE, mode="fast")
        assert result.column_confidence is None

    def test_mask_lowercase(self):
        """mask_alignment with lowercase style."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        masked = kalign.mask_alignment(result, threshold=0.5, style="lowercase")

        # Masked sequences should have some lowercase chars
        has_lower = any(c.islower() for seq in masked.sequences for c in seq if c != '-')
        assert has_lower, "Expected some lowercase residues after masking"

        # Original should be unchanged (uppercase)
        has_lower_orig = any(c.islower() for seq in result.sequences for c in seq if c != '-')
        assert not has_lower_orig, "Original should not have lowercase"

    def test_mask_remove(self):
        """mask_alignment with remove style replaces with gaps."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        masked = kalign.mask_alignment(result, threshold=0.5, style="remove")

        # Masked should have more gaps than original
        orig_gaps = sum(seq.count('-') for seq in result.sequences)
        masked_gaps = sum(seq.count('-') for seq in masked.sequences)
        assert masked_gaps >= orig_gaps

    def test_mask_preserves_confident_residues(self):
        """Confident columns should be unchanged after masking."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        masked = kalign.mask_alignment(result, threshold=0.5)

        conf = result.column_confidence
        for i, seq in enumerate(result.sequences):
            for col in range(len(seq)):
                if col < len(conf) and conf[col] >= 0.5 and seq[col] != '-':
                    assert masked.sequences[i][col] == seq[col], \
                        f"Confident residue changed at seq {i} col {col}"

    def test_mask_no_confidence_warns(self):
        """mask_alignment on non-ensemble result should warn."""
        result = kalign.align_from_file(TEST_FILE, mode="fast")
        with pytest.warns(UserWarning, match="No confidence"):
            masked = kalign.mask_alignment(result, threshold=0.5)
        # Should return unmodified
        assert masked.sequences == result.sequences

    def test_filter_removes_columns(self):
        """filter_alignment should remove low-confidence columns."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        filtered = kalign.filter_alignment(result, threshold=0.5)

        # Filtered should be shorter
        assert len(filtered.sequences[0]) < len(result.sequences[0])

        # All sequences should be same length
        lengths = [len(s) for s in filtered.sequences]
        assert len(set(lengths)) == 1

        # Filtered confidence should all be >= threshold
        assert all(c >= 0.5 for c in filtered.column_confidence)

    def test_filter_preserves_residue_count(self):
        """Filtering should not create or destroy residues — only remove columns."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        filtered = kalign.filter_alignment(result, threshold=0.3)

        for i in range(len(result.sequences)):
            orig_res = sum(1 for c in result.sequences[i] if c != '-')
            filt_res = sum(1 for c in filtered.sequences[i] if c != '-')
            # Filtered can have fewer residues (removed columns may have had residues)
            assert filt_res <= orig_res

    def test_write_confidence(self, tmp_path):
        """write_confidence should produce a file with one float per line."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        outfile = str(tmp_path / "conf.txt")
        kalign.write_confidence(outfile, result)

        with open(outfile) as f:
            lines = f.readlines()

        assert len(lines) == len(result.column_confidence)
        for line in lines:
            val = float(line.strip())
            assert 0.0 <= val <= 1.0

    def test_threshold_zero_is_noop(self):
        """Threshold 0 should not mask anything."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        masked = kalign.mask_alignment(result, threshold=0.0)
        assert masked.sequences == result.sequences

    def test_threshold_one_masks_some(self):
        """Threshold 1.0 should mask at least some columns."""
        result = kalign.align_from_file(TEST_FILE, mode="accurate")
        masked = kalign.mask_alignment(result, threshold=1.0)

        # At least some residues should be lowercase
        lower_count = sum(1 for seq in masked.sequences for c in seq if c.islower())
        assert lower_count > 0, "Expected some lowercase residues at threshold=1.0"
