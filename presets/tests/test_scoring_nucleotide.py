"""Tests for recipe.cases_nucleotide (pure functions + a synthetic tree) and recipe.scoring_nucleotide."""
import random
import textwrap
from pathlib import Path

import pytest

from recipe import cases_nucleotide as cn
from recipe import scoring_nucleotide as sn


def test_ref_identity_pooled_and_tu_equal():
    assert cn.ref_identity(["ACGU", "ACGU"]) == 1.0
    assert cn.ref_identity(["ACGU", "ACGT"]) == 1.0          # T == U
    assert cn.ref_identity(["acgu", "ACGU"]) == 1.0          # case ignored
    # col1: A/A same; col2: C/G differ; col3 gap -> skipped; col4 U/U same -> 2 of 3 pairs
    assert cn.ref_identity(["ACGU", "AG-U"]) == pytest.approx(2 / 3)
    assert cn.ref_identity(["A-", "-A"]) == 0.0


def test_has_gap_and_normalise():
    assert cn.has_gap(["AC-G"]) and cn.has_gap(["AC.G"]) and not cn.has_gap(["ACGU"])
    assert cn.rfam_normalise("ac.g-U") == "AC-G-U"


STO = textwrap.dedent("""\
    # STOCKHOLM 1.0

    #=GF AC   RF00001
    #=GF ID   fam1
    #=GF TP   Gene; sRNA;
    s1/1-4        AC.GU-
    s2/1-4        ACUGUU
    s3/1-4        AC.GUU
    s4/1-4        ACAGUU
    #=GC SS_cons  ......
    //
    # STOCKHOLM 1.0

    #=GF AC   RF00002
    #=GF ID   fam2
    a/1-3         ACG
    b/1-3         ACG
    //
    """)


def test_parse_and_filter(tmp_path):
    p = tmp_path / "x.sto"
    p.write_text(STO)
    fams = list(cn.parse_stockholm(p))
    assert [f["ac"] for f in fams] == ["RF00001", "RF00002"]
    assert fams[0]["tp"] == "Gene; sRNA;" and len(fams[0]["rows"]) == 4
    refs, why = cn.rfam_filter(fams[0])
    assert why is None and refs[0] == ("s1/1-4", "AC-GU-")
    assert cn.rfam_filter(fams[1]) == (None, "n_seq<4")


def test_filter_reasons():
    mk = lambda rows: {"rows": rows}  # noqa: E731
    gapless = mk([(f"s{i}", "ACGU") for i in range(5)])
    assert cn.rfam_filter(gapless)[1] == "gapless"
    big = mk([(f"s{i}", "A" * 2001) for i in range(5)])
    assert cn.rfam_filter(big)[1] == "aln_len>2000"
    many = mk([(f"s{i}", "AC-") for i in range(201)])
    assert cn.rfam_filter(many)[1] == "n_seq>200"
    ragged = mk([("a", "AC-"), ("b", "AC"), ("c", "AC-"), ("d", "AC-")])
    assert cn.rfam_filter(ragged)[1] == "ragged"
    empty = mk([("a", "---"), ("b", "AC-"), ("c", "AC-"), ("d", "AC-")])
    assert cn.rfam_filter(empty)[1] == "empty_sequence"


def test_interleaved_blocks(tmp_path):
    p = tmp_path / "i.sto"
    p.write_text("# STOCKHOLM 1.0\n#=GF AC RF9\na AC\nb A-\n\na GU\nb GU\n//\n")
    (fam,) = list(cn.parse_stockholm(p))
    assert fam["rows"] == [("a", "ACGU"), ("b", "A-GU")]


def test_allocate_properties():
    sizes = {"a": 1000, "b": 100, "c": 5, "d": 2}
    al = cn.allocate(sizes, 100)
    assert sum(al.values()) == 100
    assert al["d"] == 2 and al["c"] == 3 and all(al[k] <= sizes[k] for k in sizes)
    assert al["a"] > al["b"] > al["c"]
    assert cn.allocate(sizes, 5000) == sizes


def test_allocate_exact_totals():
    rnd = random.Random(1)
    for _ in range(50):
        sizes = {i: rnd.randint(1, 300) for i in range(rnd.randint(1, 16))}
        tot = rnd.randint(1, sum(sizes.values()))
        al = cn.allocate(sizes, tot)
        assert sum(al.values()) == tot and all(0 <= al[k] <= sizes[k] for k in sizes)


def test_stratified_sample_deterministic():
    items = [f"RF{i:05d}" for i in range(500)]
    strata = [("x" if i % 3 else "y") for i in range(500)]
    a, al = cn.stratified_sample(items, strata, 60, seed=1234)
    b, _ = cn.stratified_sample(items, strata, 60, seed=1234)
    c, _ = cn.stratified_sample(items, strata, 60, seed=99)
    assert a == b and a != c and len(a) == 60
    assert sorted(a) == a                                    # input order kept


def test_bins():
    assert cn.stratum_label(4, 0.55) == "n4-5|id0.00-0.60"
    assert cn.stratum_label(200, 1.0) == "n25-200|id0.90-1.00"
    assert cn.nseq_bin(10) == 2 and cn.ident_bin(0.75) == 2


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def test_mdsa_filter_and_strip(tmp_path):
    ok = tmp_path / "ok.afa"
    _write(ok, ">a\nAC-GT\n>b\nACAGT\n>c\nAC-GA\n")
    assert cn.mdsa_old_filter(ok) == (3, 5, True)
    two = tmp_path / "two.afa"
    _write(two, ">a\nAC\n>b\nAC\n")
    assert cn.mdsa_old_filter(two)[2] is False
    long_ = tmp_path / "long.afa"
    _write(long_, ">a\n" + "A" * 401 + "\n>b\nA\n>c\nA\n")
    assert cn.mdsa_old_filter(long_)[2] is False
    out = tmp_path / "u" / "ok.afa"
    cn.strip_gaps_fasta(ok, out)
    assert out.read_text() == ">a\nACGT\n>b\nACAGT\n>c\nACGA\n"


def test_build_on_synthetic_tree(tmp_path):
    ws = tmp_path
    raw = ws / "data" / "raw"
    # bralibase: one gapped, one gapless case, plus a set2 pair
    d1 = raw / "bralibase/unpacked/data-set1/data-set1/tRNA"
    _write(d1 / "structural/aln1.fa", ">a\nAC-GU\n>b\nACAGU\n")
    _write(d1 / "unaligned/aln1.fa", ">a\nACGU\n>b\nACAGU\n")
    _write(d1 / "structural/aln2.fa", ">a\nACGU\n>b\nACGU\n")
    _write(d1 / "unaligned/aln2.fa", ">a\nACGU\n>b\nACGU\n")
    d2 = raw / "bralibase/unpacked/data-set2/data-set2"
    _write(d2 / "structural/aln9.fasta", ">a\nA-C\n>b\nAGC\n")
    _write(d2 / "unaligned/aln9.fasta", ">a\nAC\n>b\nAGC\n")
    # mdsa
    for db in cn.MDSA_DBS:
        (raw / f"mdsa/unpacked/{db}_mdsa_all/{db}_mdsa_all").mkdir(parents=True)
    _write(cn.mdsa_dir(raw / "mdsa/unpacked", "smart") / "X.afa", ">a\nAC-GT\n>b\nACAGT\n>c\nAC-GA\n")
    _write(cn.mdsa_dir(raw / "mdsa/unpacked", "smart") / "Y.afa", ">a\nACGTA\n>b\nACGTA\n>c\nACGTA\n")
    # rfam
    _write(raw / "rfam/Rfam.seed", STO)
    rows, st = cn.build_manifest(ws, sample=5)
    ids = [r["case_id"] for r in rows]
    assert ids == ["bralibase/tRNA/aln1", "bralibase/data-set2/aln9", "mdsa/smart/X", "rfam/RF00001"]
    assert st["bralibase"]["gapless_dropped"] == 1
    assert st["mdsa"]["per_db"]["smart"]["gapless_dropped"] == 1
    assert st["rfam"]["dropped"] == {"n_seq<4": 1}
    r = rows[-1]
    assert Path(r["reference"]).read_text().startswith(">s1/1-4\nAC-GU-\n")
    assert Path(r["unaligned"]).read_text().startswith(">s1/1-4\nACGU\n")
    assert r["score_rule"] == "all_columns" and r["n_seq"] == "4" and r["ref_len"] == "6"
    out = ws / "data" / "cases.tsv"
    cn.write_manifest(rows, out)
    from recipe.scoring import read_manifest
    assert len(read_manifest(out)) == 4


# ---------------------------------------------------------------- scorer (needs the kalign wheel)
kalign = pytest.importorskip("kalign")


def _case(tmp_path, ref_txt, unal_txt):
    ref, unal = tmp_path / "ref.fa", tmp_path / "unal.fa"
    ref.write_text(ref_txt)
    unal.write_text(unal_txt)
    return {"reference": str(ref), "unaligned": str(unal), "score_rule": "all_columns"}


def test_score_identity_and_wrong(tmp_path):
    case = _case(tmp_path, ">a\nACGU-A\n>b\nACGUUA\n>c\nAC-UUA\n", ">a\nACGUA\n>b\nACGUUA\n>c\nACUUA\n")
    same = tmp_path / "same.fa"
    same.write_text(Path(case["reference"]).read_text())
    s = sn.score_case(case, same)
    assert s["recall"] == pytest.approx(1.0) and s["precision"] == pytest.approx(1.0)
    assert s["tc"] == pytest.approx(1.0) and s["width_ratio"] == 1.0 and s["n_seq"] == 3
    assert s["integrity_bad"] == 0
    bad = tmp_path / "bad.fa"
    bad.write_text(">a\nACGUA--\n>b\nACGUUA-\n>c\nAC-UUA-\n".replace(">a\nACGUA--", ">a\n--ACGUA"))
    s2 = sn.score_case(case, bad)
    assert s2["recall"] < 1.0 and s2["width_ratio"] > 1.0


def test_integrity_detects_changed_and_missing(tmp_path):
    case = _case(tmp_path, ">a\nAC-G\n>b\nACAG\n", ">a\nACG\n>b\nACAG\n")
    aln = tmp_path / "o.fa"
    aln.write_text(">a\nACCG\n>b\nACAG\n")
    assert sn.integrity_bad(case["unaligned"], aln) == 1
    aln.write_text(">b\nACAG\n")
    assert sn.integrity_bad(case["unaligned"], aln) == 1


def test_unknown_rule_raises(tmp_path):
    case = _case(tmp_path, ">a\nA-\n>b\nAA\n", ">a\nA\n>b\nAA\n")
    case["score_rule"] = "core"
    with pytest.raises(ValueError):
        sn.score_case(case, case["reference"])


def test_score_raises_on_integrity_failure(tmp_path):
    case = _case(tmp_path, ">a\nACGU-A\n>b\nACGUUA\n", ">a\nACGUA\n>b\nACGUUA\n")
    aln = tmp_path / "o.fa"
    aln.write_text(">a\nACCUA-\n>b\nACGUUA\n")  # sequence a changed (G -> C)
    with pytest.raises(ValueError, match="differs from input"):
        sn.score_case(case, aln)
    aln.write_text(">b\nACGUUA\n")  # sequence a missing
    with pytest.raises(ValueError):
        sn.score_case(case, aln)
    assert sn.SCORER_VERSION == "nucleotide-v2"   # v2: score_rule "projected" added (6b2cc93)
