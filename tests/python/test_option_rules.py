"""Option rules for the Python API: gap options need fast, fill options need accurate."""

import warnings

import pytest

import kalign

SEQS = [
    "MKVLAAGIVGLLLAQPSFA",
    "MKVLAGIVGLLLAQPSFA",
    "MKVLAAGIVGLLAQPSFA",
    "MKILAAGIVGLLLAQPSA",
    "MKVLAAGLVGLLLAQPSFA",
    "MKVLSAGIVGLLLAQPSFA",
]

FAST_ONLY = ["gap_open", "gap_extend", "terminal_gap_extend"]
ACCURATE_ONLY = {"fill_in": "realign", "fill_anchor_frac": 0.5, "fill_support": 2}

ACC_MSG = (
    "{name} cannot be used with mode='accurate': the accurate ensemble uses its "
    "fitted settings. Use mode='fast' to set alignment parameters."
)
FAST_MSG = (
    "{name} cannot be used with mode='fast': it only applies to the accurate "
    "ensemble. Use mode='accurate'."
)


@pytest.fixture
def fasta(tmp_path):
    p = tmp_path / "in.fa"
    p.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(SEQS)))
    return str(p)


def _call(kind, fasta, tmp_path, **kw):
    if kind == "align":
        return kalign.align(SEQS, **kw)
    if kind == "from_file":
        return kalign.align_from_file(fasta, **kw)
    return kalign.align_file_to_file(fasta, str(tmp_path / "out.fa"), **kw)


KINDS = ["align", "from_file", "file_to_file"]


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("mode", ["accurate", "recall", "ACCURATE"])
@pytest.mark.parametrize("name", FAST_ONLY)
def test_gap_option_with_accurate_raises(kind, mode, name, fasta, tmp_path):
    with pytest.raises(ValueError) as exc:
        _call(kind, fasta, tmp_path, mode=mode, **{name: 5.0})
    assert str(exc.value) == ACC_MSG.format(name=name)


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("mode", [None, "fast", "default"])
@pytest.mark.parametrize("name", list(ACCURATE_ONLY))
def test_fill_option_with_fast_raises(kind, mode, name, fasta, tmp_path):
    with pytest.raises(ValueError) as exc:
        _call(kind, fasta, tmp_path, mode=mode, **{name: ACCURATE_ONLY[name]})
    assert str(exc.value) == FAST_MSG.format(name=name)


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("name", FAST_ONLY)
@pytest.mark.parametrize("bad", [-1.0, -0.5, float("nan"), float("inf"), "5"])
def test_bad_gap_penalty_raises(kind, name, bad, fasta, tmp_path):
    with pytest.raises(ValueError, match=name):
        _call(kind, fasta, tmp_path, **{name: bad})


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("bad", [0, -1, 1.5, "2", True])
def test_bad_threads_raise(kind, bad, fasta, tmp_path):
    with pytest.raises(ValueError, match="n_threads"):
        _call(kind, fasta, tmp_path, n_threads=bad)


def test_add_to_alignment_bad_threads(fasta, tmp_path):
    with pytest.raises(ValueError, match="n_threads"):
        kalign.add_to_alignment(fasta, fasta, str(tmp_path / "o.fa"), n_threads=0)


def test_gap_penalties_with_fast_work():
    base = kalign.align(SEQS, mode="fast")
    res = kalign.align(SEQS, mode="default", gap_open=8.0, gap_extend=1.0)
    assert len(res) == len(SEQS)
    assert len({len(s) for s in res}) == 1
    assert len(base) == len(SEQS)
    # 0 is a valid penalty
    kalign.align(SEQS, gap_open=0, terminal_gap_extend=0.0)


def test_gap_penalties_do_not_switch_mode():
    # Before: a gap penalty silently turned accurate into fast.
    with pytest.raises(ValueError):
        kalign.align(SEQS, mode="accurate", gap_open=5.0)


def test_accurate_from_file_keeps_confidence(fasta):
    for mode in ("accurate", "recall"):
        res = kalign.align_from_file(fasta, mode=mode)
        assert res.column_confidence is not None
        assert len(res.column_confidence) == len(res.sequences[0])


def test_fast_has_no_confidence(fasta):
    for mode in (None, "fast", "default"):
        assert kalign.align_from_file(fasta, mode=mode).column_confidence is None


@pytest.mark.parametrize("mode", ["accurate", "recall"])
def test_fill_options_with_accurate_work(mode, fasta):
    res = kalign.align_from_file(
        fasta, mode=mode, fill_in="attach", fill_anchor_frac=0.5, fill_support=1
    )
    assert res.column_confidence is not None


def test_confidence_messages_name_accurate_mode(tmp_path):
    res = kalign.AlignedSequences(["a"], ["AC"])
    with pytest.warns(UserWarning, match="mode='accurate'"):
        kalign.mask_alignment(res)
    with pytest.warns(UserWarning, match="mode='accurate'"):
        kalign.filter_alignment(res)
    with pytest.raises(ValueError, match="mode='accurate'"):
        kalign.write_confidence(str(tmp_path / "c.txt"), res)


def test_mask_alignment_bad_style():
    res = kalign.AlignedSequences(["a", "b"], ["AC", "AC"], column_confidence=[1, 0])
    with pytest.raises(ValueError, match="style"):
        kalign.mask_alignment(res, style="upper")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert kalign.mask_alignment(res, style="remove").sequences == ["A-", "A-"]
        assert kalign.mask_alignment(res).sequences == ["Ac", "Ac"]


def test_no_high_thread_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        # Validation no longer warns above 1024; only the helper is checked here.
        assert kalign._resolve_threads(5000) == 5000


def test_default_threads_is_one():
    assert kalign.get_num_threads() == 1


def test_write_alignment_formats(tmp_path):
    pytest.importorskip("Bio")
    for fmt in ("fasta", "FA", "clustal", "aln", "stockholm", "sto", "phylip", "phy"):
        kalign.write_alignment(["AC-", "ACG"], str(tmp_path / f"o.{fmt}"), format=fmt)
    for fmt in ("msf", "clu"):
        with pytest.raises(ValueError):
            kalign.write_alignment(["AC-", "ACG"], str(tmp_path / "x"), format=fmt)


@pytest.mark.parametrize("fmt", ["fasta", "fa", "afa", "msf", "clu", "clustal"])
def test_file_to_file_formats(fmt, fasta, tmp_path):
    out = tmp_path / "o.out"
    kalign.align_file_to_file(fasta, str(out), format=fmt)
    assert out.stat().st_size > 0


def test_refine_constants_still_exported():
    for name in ("REFINE_NONE", "REFINE_ALL", "REFINE_CONFIDENT", "REFINE_INLINE"):
        assert hasattr(kalign, name)
    assert not hasattr(kalign, "_parse_refine_mode")


@pytest.mark.parametrize("kind", KINDS)
def test_numpy_scalars_accepted(kind, fasta, tmp_path):
    np = pytest.importorskip("numpy")
    _call(
        kind,
        fasta,
        tmp_path,
        mode="fast",
        gap_open=np.float32(5.0),
        gap_extend=np.int64(2),
        n_threads=np.int64(2),
    )


# These two depend on the exact, case-insensitive format matching in
# lib/src/msa_io.c (R4); they fail until that change is merged.
@pytest.mark.parametrize("fmt", ["FASTA", "MSF", "Clu"])
def test_file_to_file_formats_uppercase(fmt, fasta, tmp_path):
    out = tmp_path / "o.out"
    kalign.align_file_to_file(fasta, str(out), format=fmt)
    assert out.stat().st_size > 0


@pytest.mark.parametrize("fmt", ["clustalw", "xmsfx"])
def test_file_to_file_bad_format_raises(fmt, fasta, tmp_path):
    with pytest.raises(Exception):
        kalign.align_file_to_file(fasta, str(tmp_path / "o.out"), format=fmt)
