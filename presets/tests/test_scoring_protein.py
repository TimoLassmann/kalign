"""Tests for recipe.cases_protein and recipe.scoring_protein.

Synthetic cases need nothing. Tests that call kalign.compare_detailed are skipped without the wheel;
the whole-manifest self-test (reference scored against itself = 1.0 for every case) is skipped
unless <ws>/data/cases_protein.tsv exists (run it on the workstation).
"""
import importlib.util
import os
from pathlib import Path

import pytest

from recipe import cases_protein as cp
from recipe import scoring_protein as sp

HAVE_KALIGN = importlib.util.find_spec("kalign") is not None


def fa(path, recs):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(f">{n}\n{s}\n" for n, s in recs))


# ---------- pure helpers
def test_read_fasta_wrapped_and_header_token(tmp_path):
    p = tmp_path / "x"
    p.write_text(">a desc here\nAC\nGT\n>b\nTT\n")
    assert sp.read_fasta(p) == [("a", "ACGT"), ("b", "TT")]


def test_dash_gaps_normalises_dot_dash_bang():
    assert sp.dash_gaps([("a", "Ab.c-D!")]) == [("a", "Ab-c-D-")]


def test_core_columns_rule():
    ref = sp.dash_gaps([("a", "AAbB.-"), ("b", "A-bb.c")])
    # col0 upper/upper core; col1 upper+gap core; col2 lower/lower no; col3 upper+lower mixed no;
    # col4 gap/gap no; col5 gap+lower no
    mask, mixed = sp.core_mask(ref, "core_columns")
    assert mask == [1, 1, 0, 0, 0, 0]
    assert mixed == 1


def test_pair_core_rule_needs_both_upper():
    ref = sp.dash_gaps([("a", "AAb-"), ("b", "A-bC")])
    assert sp.core_mask(ref, "pair_core")[0] == [1, 0, 0, 0]
    with pytest.raises(ValueError):
        sp.core_mask(sp.dash_gaps([("a", "A"), ("b", "A"), ("c", "A")]), "pair_core")
    with pytest.raises(ValueError):
        sp.core_mask(ref, "nope")


def test_project_drops_allgap_columns_and_orders():
    test = [("c", "AC--"), ("a", "A-GT"), ("b", "A-G-")]
    assert sp.project(test, ["a", "b"]) == [("a", "AGT"), ("b", "AG-")]
    with pytest.raises(ValueError):
        sp.project(test, ["a", "zz"])
    with pytest.raises(ValueError):
        sp.project(test + [("a", "----")], ["a"])


# ---------- prepare(): consistency checks (no kalign)
def make_case(tmp_path, rule="core_columns"):
    ref = [("s1", "AC.DE"), ("s2", "ACgDE")]
    inp = [("s1", "ACDE"), ("s2", "ACGDE"), ("s3", "KKK")]
    fa(tmp_path / "ref", ref)
    fa(tmp_path / "in", inp)
    return {"case_id": "t.1", "reference": str(tmp_path / "ref"), "unaligned": str(tmp_path / "in"),
            "score_rule": rule}


def test_prepare_ok_projects_extra_sequences(tmp_path):
    row = make_case(tmp_path)
    fa(tmp_path / "aln", [("s1", "AC-DE----"), ("s2", "ACGDE----"), ("s3", "-----KKK-")][:2]
       + [("s3", "-----KKK-")])
    ref_dash, proj, mask, mixed = sp.prepare(row, tmp_path / "aln")
    assert [n for n, _ in proj] == ["s1", "s2"]
    assert proj == [("s1", "AC-DE"), ("s2", "ACGDE")]  # s3 and its all-gap columns removed
    assert mask == [1, 1, 0, 1, 1] and mixed == 0


@pytest.mark.parametrize("aln, msg", [
    ([("s1", "ACDE"), ("s2", "ACGDE")], "differs from input"),                       # s3 missing
    ([("s1", "ACDE"), ("s2", "ACGDE"), ("s3", "KKA")], "differs from input"),         # residues changed
    ([("s1", "ACDE"), ("s2", "ACGDE"), ("s3", "KKK"), ("s4", "A")], "differs from input"),
])
def test_prepare_rejects_bad_alignments(tmp_path, aln, msg):
    row = make_case(tmp_path)
    w = max(len(s) for _, s in aln)
    fa(tmp_path / "aln", [(n, s.ljust(w, "-")) for n, s in aln])
    with pytest.raises(ValueError, match=msg):
        sp.prepare(row, tmp_path / "aln")


def test_prepare_empty_mask(tmp_path):
    fa(tmp_path / "ref", [("a", "acd"), ("b", "acd")])
    fa(tmp_path / "in", [("a", "ACD"), ("b", "ACD")])
    fa(tmp_path / "aln", [("a", "ACD"), ("b", "ACD")])
    row = {"reference": str(tmp_path / "ref"), "unaligned": str(tmp_path / "in"), "score_rule": "core_columns"}
    with pytest.raises(ValueError, match="no core"):
        sp.prepare(row, tmp_path / "aln")


# ---------- scoring with the real compare_detailed
@pytest.mark.skipif(not HAVE_KALIGN, reason="kalign wheel not installed")
def test_score_perfect_and_shifted(tmp_path):
    ref = [("s1", "ACDEFGHIK"), ("s2", "ACDEFGHIK"), ("s3", "ACDEFGHIK")]
    fa(tmp_path / "ref", ref)
    fa(tmp_path / "in", ref)
    row = {"case_id": "t", "reference": str(tmp_path / "ref"), "unaligned": str(tmp_path / "in"),
           "score_rule": "core_columns", "n_seq": "3"}
    fa(tmp_path / "good", ref)
    r = sp.score_case(row, tmp_path / "good")
    assert (r["recall"], r["precision"], r["f1"], r["tc"]) == (1.0, 1.0, 1.0, 1.0)
    assert r["width_ratio"] == 1.0 and r["core_cols"] == 9 and r["n_seq"] == 3
    fa(tmp_path / "bad", [("s1", "ACDEFGHIK-"), ("s2", "ACDEFGHIK-"), ("s3", "-ACDEFGHIK")])
    r = sp.score_case(row, tmp_path / "bad")
    assert r["recall"] < 0.5 and r["tc"] < 0.5 and r["width_ratio"] > 1.0


@pytest.mark.skipif(not HAVE_KALIGN, reason="kalign wheel not installed")
def test_prefab_pair_scored_on_both_upper_columns(tmp_path):
    fa(tmp_path / "ref", [("p", "aCDe."), ("q", "aC-eF")])  # only the C/C column has both upper case
    fa(tmp_path / "in", [("p", "ACDE"), ("q", "ACEF"), ("z", "GGG")])
    row = {"case_id": "t", "reference": str(tmp_path / "ref"), "unaligned": str(tmp_path / "in"),
           "score_rule": "pair_core", "n_seq": "3"}
    fa(tmp_path / "aln", [("p", "ACDE-----"), ("q", "AC-EF----"), ("z", "-----GGG-")])
    r = sp.score_case(row, tmp_path / "aln")
    assert r["core_cols"] == 1  # only the C/C column has both residues upper case
    assert r["recall"] == 1.0 and r["tc"] == 1.0


# ---------- manifest builder
def test_pid_and_subset_helpers():
    assert [cp.pid_bin(x) for x in (0, 19.9, 20, 29.9, 30, 39.9, 40, 99)] == [
        "pid<20", "pid<20", "pid20-30", "pid20-30", "pid30-40", "pid30-40", "pid>=40", "pid>=40"]
    assert cp.bali_subset("BB11001") == "BB-RV11"
    assert cp.bali_subset("BBS50010") == "BBS-RV50"


def test_validate_case_reasons():
    ok_ref = [("a", "AC.D"), ("b", "ACgD")]
    ok_in = [("a", "ACD"), ("b", "ACGD")]
    assert cp.validate_case(ok_ref, ok_in, "core_columns") is None
    assert "fewer" in cp.validate_case(ok_ref[:1], ok_in, "core_columns")
    assert "not in input" in cp.validate_case(ok_ref, ok_in[:1], "core_columns")
    assert "differ from input" in cp.validate_case(ok_ref, [("a", "ACD"), ("b", "ACGG")], "core_columns")
    assert "duplicate names in reference" in cp.validate_case([("a", "AC"), ("a", "AC")], ok_in, "core_columns")
    assert "no core" in cp.validate_case([("a", "ac"), ("b", "ac")], [("a", "AC"), ("b", "AC")], "core_columns")
    # core columns exist but each holds a single residue: no pair to assess
    assert "no assessed" in cp.validate_case([("a", "A-"), ("b", "-C")], [("a", "A"), ("b", "C")], "core_columns")


def test_build_on_synthetic_workspace(tmp_path):
    ws = tmp_path / "ws"
    d5 = ws.joinpath(*cp.DRIVE5)
    for name, ids in (("bali3", ["BB11001", "BBS12002"]), ("prefab4", ["p1"]), ("ox", ["o1"]), ("sabre", ["s1"])):
        root = d5 / name
        (root / "info").mkdir(parents=True)
        (root / "info" / "ids.txt").write_text("\n".join(ids) + "\n")
        (root / "info" / "pctids.txt").write_text("".join(f"{i};25.5\n" for i in ids))
        for i in ids:
            if name == "prefab4":
                fa(root / "ref" / i, [("a", "AC.D"), ("b", "ACgD")])
                fa(root / "in" / i, [("a", "ACD"), ("b", "ACGD"), ("c", "kk")])
            else:
                fa(root / "ref" / i, [("a", "AC.D"), ("b", "ACgD")])
                fa(root / "in" / i, [("a", "acd"), ("b", "ACGD")])
    bf = ws.joinpath(*cp.BALIFAM)
    (bf / "info").mkdir(parents=True)
    (bf / "info" / "ids.txt").write_text("PF1.100\nPF2.100\n")
    fa(bf / "ref" / "PF1.100", [("a", "AC.D"), ("b", "ACgD")])
    fa(bf / "in" / "PF1.100", [("a", "ACD"), ("b", "ACGD")])
    fa(bf / "ref" / "PF2.100", [("a", "ac"), ("b", "ac")])  # no core -> dropped
    fa(bf / "in" / "PF2.100", [("a", "AC"), ("b", "AC")])
    rows, dropped, stats = cp.build(ws)
    assert len(rows) == 6 and [d[0] for d in dropped] == ["balifam.PF2.100"]
    man = (ws / "data" / "cases_protein.tsv").read_text().splitlines()
    assert man[0].split("\t") == list(cp.MANIFEST_COLUMNS)
    assert len(man) == 7
    r = {x["case_id"]: x for x in rows}
    assert r["bali3.BBS12002"]["subset"] == "BBS-RV12" and r["prefab4.p1"]["subset"] == "pid20-30"
    assert r["prefab4.p1"]["score_rule"] == "pair_core" and r["prefab4.p1"]["n_seq"] == "3"
    # normalised input: upper case, same sequences
    assert sp.read_fasta(r["bali3.BB11001"]["unaligned"]) == [("a", "ACD"), ("b", "ACGD")]
    assert r["bali3.BB11001"]["sha256_reference"] == cp.sha256_file(r["bali3.BB11001"]["reference"])
    assert stats["n_cases"] == 6 and stats["n_dropped"] == 1
    # deterministic
    before = (ws / "data" / "cases_protein.tsv").read_bytes()
    cp.build(ws)
    assert (ws / "data" / "cases_protein.tsv").read_bytes() == before


# ---------- whole manifest (workstation only)
MANIFEST = cp.default_ws() / "data" / "cases_protein.tsv"


@pytest.mark.skipif(not (HAVE_KALIGN and MANIFEST.exists()), reason="needs kalign wheel and the manifest")
def test_every_reference_scores_one_against_itself(tmp_path):
    from recipe.scoring import read_manifest

    rows = read_manifest(MANIFEST)
    assert len(rows) > 2000
    for r in rows:
        ref = sp.dash_gaps(sp.read_fasta(r["reference"]))
        names = {n for n, _ in ref}
        w = len(ref[0][1])
        test = list(ref)
        # extra input sequences (PREFAB, BaliFam): append as insertion columns of their own
        extra = [(n, s) for n, s in sp.read_fasta(r["unaligned"]) if n not in names]
        tot = w + sum(len(s) for _, s in extra)
        off = w
        full = [(n, s + "-" * (tot - w)) for n, s in test]
        for n, s in extra:
            full.append((n, "-" * off + s + "-" * (tot - off - len(s))))
            off += len(s)
        sp.write_fasta(tmp_path / "aln", [(n, s.upper()) for n, s in full])
        res = sp.score_case(r, tmp_path / "aln")
        assert (res["recall"], res["precision"], res["f1"], res["tc"]) == (1.0, 1.0, 1.0, 1.0), r["case_id"]
