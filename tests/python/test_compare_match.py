"""Sequence matching in kalign.compare / compare_detailed / compare_detailed_with_mask.

Reference sequences are matched to test sequences by name and residues (case and gaps ignored). The test may
hold extra sequences, which are ignored; a reference sequence missing from the test, or with different residues,
is an error. Expected numbers are worked out by hand in each test.
"""

import os
import tempfile

import pytest

import kalign

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def _write(d, name, recs):
    p = os.path.join(d, name)
    with open(p, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n{s}\n")
    return p


def _read(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    recs.append((name, "".join(buf)))
    return recs


def _cd(ref, test, **kw):
    kw.setdefault("max_gap_frac", -1.0)
    return kalign.compare_detailed(ref, test, **kw)


# Two sequences. Reference: s1 AC-, s2 A-C -> one aligned pair (A,A) = 2 directed pairs.
# Test: s1 AC, s2 AC -> pairs (A,A), (C,C) = 4 directed; common (A,A) = 2.
# recall 2/2 = 1, precision 2/4 = 0.5, TC: the one column with 2 residues (A) is correct -> 1.
REF2 = [("s1", "AC-"), ("s2", "A-C")]
TEST2 = [("s1", "AC"), ("s2", "AC")]

# Three sequences. Reference: s1 ACD, s2 ACD, s3 A-D. Test: s1 ACD, s2 ACD, s3 AD-.
# Reference pairs: col A 3 pairs, col C 1, col D 3 -> 7 pairs = 14 directed.
# Test pairs: col0 A,A,A 3; col1 C,C,D 3; col2 D,D 1 -> 7 = 14 directed.
# Common: A 3 + (s1C,s2C) 1 + (s1D,s2D) 1 = 5 -> 10 directed. recall = precision = 10/14.
# TC: columns with >= 2 residues: A, C, D; A and C correct, D split -> 2/3.
REF3 = [("s1", "ACD"), ("s2", "ACD"), ("s3", "A-D")]
TEST3 = [("s1", "ACD"), ("s2", "ACD"), ("s3", "AD-")]


def test_known_answer_pair():
    with tempfile.TemporaryDirectory() as d:
        r = _cd(_write(d, "r.fa", REF2), _write(d, "t.fa", TEST2))
    assert r["recall"] == pytest.approx(1.0)
    assert r["precision"] == pytest.approx(0.5)
    assert r["tc"] == pytest.approx(1.0)
    assert (r["ref_pairs"], r["test_pairs"], r["common_pairs"]) == (2, 4, 2)


def test_known_answer_three():
    with tempfile.TemporaryDirectory() as d:
        r = _cd(_write(d, "r.fa", REF3), _write(d, "t.fa", TEST3))
    assert r["recall"] == pytest.approx(10 / 14)
    assert r["precision"] == pytest.approx(10 / 14)
    assert r["tc"] == pytest.approx(2 / 3)
    assert (r["ref_pairs"], r["test_pairs"], r["common_pairs"]) == (14, 14, 10)


def test_extra_test_sequences_are_ignored():
    # extra sequences, placed first and last, with residues in every column: the answer stays the hand-computed one
    test = [("a0", "DDD")] + TEST3 + [("zz", "CCA")]   # "a0" sorts before the reference names
    with tempfile.TemporaryDirectory() as d:
        r = _cd(_write(d, "r.fa", REF3), _write(d, "t.fa", test))
        m = kalign.compare_detailed(_write(d, "r2.fa", REF3), _write(d, "t2.fa", test), column_mask=[1, 1, 1])
    for x in (r, m):
        assert x["recall"] == pytest.approx(10 / 14)
        assert x["precision"] == pytest.approx(10 / 14)
        assert x["tc"] == pytest.approx(2 / 3)
        assert (x["ref_pairs"], x["test_pairs"], x["common_pairs"]) == (14, 14, 10)


def test_extra_sequences_old_compare():
    # kalign.compare (the SP score in percent): identical rows plus extra test rows -> 100
    test = [("extra", "CCC")] + REF3
    with tempfile.TemporaryDirectory() as d:
        assert kalign.compare(_write(d, "r.fa", REF3), _write(d, "t.fa", test)) == pytest.approx(100.0)


def test_order_and_case_do_not_matter():
    test = [(n, s.lower()) for n, s in reversed(TEST3)]
    with tempfile.TemporaryDirectory() as d:
        r = _cd(_write(d, "r.fa", REF3), _write(d, "t.fa", test))
    assert r["recall"] == pytest.approx(10 / 14)
    assert r["tc"] == pytest.approx(2 / 3)


def test_duplicate_names_matched_by_residues():
    # two sequences named "d" with different residues; the test lists them in the other order
    ref = [("d", "AC-"), ("d", "A-W")]
    test = [("d", "A-W"), ("d", "AC-")]
    with tempfile.TemporaryDirectory() as d:
        r = _cd(_write(d, "r.fa", ref), _write(d, "t.fa", test))
    assert (r["recall"], r["precision"], r["tc"]) == pytest.approx((1.0, 1.0, 1.0))


def test_missing_reference_sequence_is_an_error():
    with tempfile.TemporaryDirectory() as d:
        ref = _write(d, "r.fa", REF3)
        renamed = _write(d, "t1.fa", [("s1", "ACD"), ("s2", "ACD"), ("other", "AD-")])
        fewer = _write(d, "t2.fa", TEST3[:2])
        for t in (renamed, fewer):
            with pytest.raises(RuntimeError):
                _cd(ref, t)
            with pytest.raises(RuntimeError):
                kalign.compare_detailed(ref, t, column_mask=[1, 1, 1])
            with pytest.raises(RuntimeError):
                kalign.compare(ref, t)


def test_different_residues_is_an_error():
    with tempfile.TemporaryDirectory() as d:
        ref = _write(d, "r.fa", REF3)
        changed = _write(d, "t1.fa", [("s1", "ACD"), ("s2", "ACE"), ("s3", "AD-")])
        longer = _write(d, "t2.fa", [("s1", "ACDW"), ("s2", "ACD-"), ("s3", "AD--")])
        for t in (changed, longer):
            with pytest.raises(RuntimeError):
                _cd(ref, t)
            with pytest.raises(RuntimeError):
                kalign.compare(ref, t)


def test_repeatable():
    # the same call gives the same answer every time (no reads of unset memory)
    with tempfile.TemporaryDirectory() as d:
        ref = _write(d, "r.fa", REF3)
        t = _write(d, "t.fa", [("x", "WWW")] + TEST3)
        first = _cd(ref, t)
        for _ in range(20):
            assert _cd(ref, t) == first


def test_real_alignment_subset_equals_projection():
    # BB11001: a "fast" alignment is the test; two rows of an "accurate" alignment (all-gap columns removed) are
    # the reference. Scoring the full test must give the same numbers as scoring the test projected onto the two
    # reference rows (all-gap columns removed), for every row pair tried.
    def drop_gap_cols(rows):
        cols = [j for j in range(len(rows[0])) if any(r[j].isalpha() for r in rows)]
        return ["".join(r[j] for j in cols) for r in rows]

    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(DATA_DIR, "BB11001.tfa")
        test_p, acc_p = os.path.join(d, "fast.fa"), os.path.join(d, "acc.fa")
        kalign.align_file_to_file(inp, test_p, mode="fast")
        kalign.align_file_to_file(inp, acc_p, mode="accurate")
        acc, test = _read(acc_p), dict(_read(test_p))
        assert len(acc) >= 4
        for keep in ([0, 1], [1, 3], [0, 2, 3]):
            names = [acc[i][0] for i in keep]
            sub_ref = _write(d, "sub_ref.fa", list(zip(names, drop_gap_cols([acc[i][1] for i in keep]))))
            proj = _write(d, "proj.fa", list(zip(names, drop_gap_cols([test[n] for n in names]))))
            for kw in ({"max_gap_frac": -1.0}, {"max_gap_frac": 0.2}):
                full = kalign.compare_detailed(sub_ref, test_p, **kw)
                assert full == kalign.compare_detailed(sub_ref, proj, **kw)
                assert 0.0 < full["recall"] <= 1.0
