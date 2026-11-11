import random

import numpy as np
import pytest

from recipe import grouping as g
from recipe.features import features_from_rows

AA = "ACDEFGHIKLMNPQRSTVWY"


def rnd(n, rng, alpha=AA):
    return "".join(rng.choice(alpha) for _ in range(n))


def mutate(s, frac, rng, alpha=AA):
    s = list(s)
    for i in rng.sample(range(len(s)), int(frac * len(s))):
        s[i] = rng.choice([c for c in alpha if c != s[i]])
    return "".join(s)


def make_rows(tmp_path, cases, pool="protein"):
    rows = []
    for name, seqs in cases.items():
        p = tmp_path / f"{name}.fa"
        p.write_text("".join(f">{name}_{i}\n{s}\n" for i, s in enumerate(seqs)))
        rows.append(dict(case_id=name, unaligned=str(p), reference=str(p), benchmark=name[0], pool=pool))
    return rows


def groups_of(rows, pool="protein", **kw):
    roots, edges, stats = g.group_cases(rows, pool, threads=2, log=lambda *_: None, **kw)
    out = {}
    for r, root in zip(rows, roots.tolist()):
        out.setdefault(root, set()).add(r["case_id"])
    return sorted(map(sorted, out.values())), edges, stats


def test_identity_kernel():
    a = g.encode("ACDEFGHIKLMNPQRSTVWY" * 3, "protein")
    assert g.sg_identity(a, a, 2, -1, -3) == 1.0
    sub = a[10:40]
    assert g.sg_identity(sub, a, 2, -1, -3) == 1.0           # containment
    assert g.sg_identity(a, sub, 2, -1, -3) == 1.0           # symmetric in argument order
    x = g.encode("A" * 30, "protein")
    y = g.encode("C" * 30, "protein")
    assert g.sg_identity(x, y, 2, -1, -3) == 0.0
    unk = g.encode("X" * 30, "protein")
    assert g.sg_identity(unk, unk, 2, -1, -3) == 0.0         # unknown residues never match


def test_normalise():
    assert g.normalise("ac-gu.n", "nucleotide") == "ACGTN"
    assert g.normalise("ac-gu", "protein") == "ACGU"


def test_rules_and_components(tmp_path):
    rng = random.Random(1)
    s1, s2, s3 = rnd(120, rng), rnd(150, rng), rnd(90, rng)
    cases = {
        "A1": [s1, rnd(100, rng)],
        "B1": [s1, rnd(100, rng)],                  # identical sequence to A1
        "C1": [s2[30:110], rnd(80, rng)],           # contained in D1's s2
        "D1": [s2],
        "E1": [mutate(s3, 0.30, rng)],              # 70 % identical to F1
        "F1": [s3],
        "G1": [mutate(s3[:60], 0.65, rng)],         # only 35 % identical: not linked
        "H1": [rnd(130, rng)],                      # unrelated
    }
    comps, edges, stats = groups_of(make_rows(tmp_path, cases))
    assert comps == [["A1", "B1"], ["C1", "D1"], ["E1", "F1"], ["G1"], ["H1"]]
    kinds = {k for *_, k, _ in edges}
    assert "exact" in kinds and ("contained/identical" in kinds or "identity" in kinds)


def test_transitive_chain(tmp_path):
    rng = random.Random(2)
    base = rnd(200, rng)
    a = base
    b = mutate(a, 0.3, rng)
    c = mutate(b, 0.3, rng)   # c is ~50 % from a, linked through b
    comps, _, _ = groups_of(make_rows(tmp_path, {"A": [a], "B": [b], "C": [c]}))
    assert comps == [["A", "B", "C"]]


def test_short_sequences_only_exact(tmp_path):
    rng = random.Random(3)
    short = rnd(12, rng)
    cases = {"A": [short, rnd(100, rng)], "B": [short + "A", rnd(100, rng)], "C": [short, rnd(100, rng)]}
    comps, _, _ = groups_of(make_rows(tmp_path, cases))
    assert comps == [["A", "C"], ["B"]]            # identical short seq links; containing it does not


def test_nucleotide_u_equals_t(tmp_path):
    rng = random.Random(4)
    dna = rnd(100, rng, "ACGT")
    rows = make_rows(tmp_path, {"A": [dna], "B": [dna.replace("T", "U")], "C": [rnd(100, rng, "ACGT")]}, "nucleotide")
    comps, _, _ = groups_of(rows, "nucleotide")
    assert comps == [["A", "B"], ["C"]]


def test_random_null_not_linked(tmp_path):
    rng = random.Random(5)
    cases = {f"R{i}": [rnd(rng.randint(60, 300), rng)] for i in range(40)}
    comps, _, _ = groups_of(make_rows(tmp_path, cases))
    assert len(comps) == 40


def test_relabel_order():
    roots = np.array([0, 0, 2, 2, 2, 5])
    ids = ["a", "b", "c", "d", "e", "f"]
    assert g.relabel(roots, ids) == ["G0002", "G0002", "G0001", "G0001", "G0001", "G0003"]


def test_prefilter_validation_runs(tmp_path):
    rng = random.Random(6)
    base = rnd(150, rng)
    cases = {f"M{i}": [mutate(base, 0.2, rng)] for i in range(6)}
    cases.update({f"U{i}": [rnd(150, rng)] for i in range(6)})
    _, _, stats = groups_of(make_rows(tmp_path, cases), validate=12)
    assert stats["prefilter_recall"]["pairs"] > 0


# ------------------------------------------------------------------ features
def test_features_basic():
    rows = ["AC--GT", "ACTTGT", "-C-TGA"]
    f = features_from_rows(rows, nucleotide=True)
    assert f["n_seq"] == 3 and f["ref_len"] == 6
    assert f["len_min"] == 4 and f["len_max"] == 6
    assert f["gap_frac"] == pytest.approx(4 / 18)
    assert f["n_gap_runs"] == 2 and f["gap_run_len_max"] == 2    # row0 "--" and row2 "-" (internal)
    assert f["lead_max"] == 1 and f["trail_max"] == 0
    assert f["term_ext_frac"] == pytest.approx(1 / 18)
    # pair identities: (0,1) 4/4, (0,2) 2/3 over C,G? cols with both: C(1),G(4)->T? see below
    assert 0.0 < f["mean_pid"] <= 1.0


def test_features_identity_exact():
    f = features_from_rows(["AAAA", "AAAC"], nucleotide=False)
    assert f["mean_pid"] == pytest.approx(0.75)
    f = features_from_rows(["AAAA", "aaac"], nucleotide=False, mask=[1, 1, 0, 0])
    assert f["mean_pid_core"] == 1.0 and f["core_cols"] == 2
    f = features_from_rows(["AUG", "ATG"], nucleotide=True)
    assert f["mean_pid"] == 1.0


def test_final_protein_groups_link_containment_with_ambiguous_residue(tmp_path):
    """Review 2026-10-05: a reference sequence contained in another but holding an X (never 'identical' in the
    identity kernel) must still link the two cases (rule: identical or contained)."""
    rng = random.Random(7)
    s = rnd(200, rng)
    s = s[:100] + "X" + s[101:]
    cases = {"A": [s[:150], rnd(90, rng)], "B": [s + "K"], "C": [rnd(150, rng)]}
    rows = make_rows(tmp_path, cases)
    roots, _, _ = g.group_cases(rows, "protein", threads=2, log=lambda *_: None, source="reference", identity=1.0)
    assert roots[0] != roots[1]                       # the identity kernel alone misses it (X never matches)
    roots, edges, stats = g.final_groups(rows, "protein", threads=2, log=lambda *_: None)
    assert roots[0] == roots[1] != roots[2]
    assert stats["containment_exact_merges"] == 1
