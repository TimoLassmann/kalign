"""Protein pool scorer: ``score_case(case_row, aln_path)`` for the drive5 benchmarks (BAliBASE 3,
PREFAB 4, OXBench, SABRE) and BaliFam100.

Every case is scored the same way: the test alignment is *projected* onto the reference sequences
(other sequences removed, columns that are gaps in all of them dropped), the reference is rewritten
with ``-`` gaps (the sources use ``.``), and ``kalign.compare_detailed(ref, projected,
column_mask=M)`` is called. Only the mask rule depends on the benchmark (``score_rule`` column of
the manifest):

  core_columns  M = reference columns with at least one upper-case residue and no lower-case
                residue. Upper case marks the core in drive5 and BaliFam references (drive5
                README; BaliFam README). Used for bali3, ox, sabre, balifam.
                (= ``score_balifam_case`` of experiments/redesign/tools/score_confirmation.py)
  pair_core     M = reference columns where BOTH residues of the reference pair are upper case.
                Used for prefab4, whose reference is the 2-sequence pair while the input holds
                the full set of ~50 sequences. (= ``score_prefab_case`` of the same file)

Returned dict: recall (SP), precision, f1, tc, width_ratio, n_seq (+ extras core_cols, mixed_cols,
ref_pairs, test_pairs). width_ratio = width of the projected test alignment / reference width.
A test alignment that does not contain every input sequence unchanged (residues compared case
insensitively) is an error, never silently scored.

The functions are pure Python except the call to ``kalign.compare_detailed``; everything else can be
tested without the kalign wheel via ``prepare()``.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

RULES = ("core_columns", "pair_core")
SCORER_VERSION = "protein-v1"  # part of the cache key: bump on ANY change to the scoring


def read_fasta(path) -> list[tuple[str, str]]:
    """List of (name, sequence) in file order. Name = first whitespace token of the header; sequence
    lines are concatenated (files may be wrapped)."""
    recs, name, buf = [], None, []
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                tok = line[1:].split()
                name, buf = (tok[0] if tok else ""), []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path, recs) -> None:
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n{s}\n")


def residues(s: str) -> str:
    return "".join(c for c in s.upper() if c.isalpha())


def dash_gaps(recs):
    """Reference rows with every non-letter (``.``, ``-``, ``!``...) written as ``-``."""
    return [(n, "".join(c if c.isalpha() else "-" for c in s)) for n, s in recs]


def check_rectangular(recs, what="alignment"):
    widths = {len(s) for _, s in recs}
    if len(widths) != 1:
        raise ValueError(f"{what}: rows differ in length {sorted(widths)}")
    return widths.pop()


def core_mask(ref_dash, rule: str):
    """(mask list, number of mixed-case columns) for a dash-gapped reference with case marks."""
    width = check_rectangular(ref_dash, "reference")
    if rule == "core_columns":
        mask, mixed = [], 0
        for j in range(width):
            col = [s[j] for _, s in ref_dash if s[j].isalpha()]
            up, lo = any(c.isupper() for c in col), any(c.islower() for c in col)
            mask.append(1 if (up and not lo) else 0)
            mixed += int(up and lo)
        return mask, mixed
    if rule == "pair_core":
        if len(ref_dash) != 2:
            raise ValueError(f"rule pair_core needs a reference pair, got {len(ref_dash)} sequences")
        a, b = ref_dash[0][1], ref_dash[1][1]
        return [1 if (x.isupper() and y.isupper()) else 0 for x, y in zip(a, b)], 0
    raise ValueError(f"unknown score_rule {rule!r}")


def project(test_recs, names):
    """Keep ``names`` (in that order), drop columns that are gaps in all of them."""
    by = {}
    for n, s in test_recs:
        if n in by:
            raise ValueError(f"duplicate sequence name {n} in test alignment")
        by[n] = s
    missing = [n for n in names if n not in by]
    if missing:
        raise ValueError(f"reference sequences missing from test alignment: {missing[:5]}")
    rows = [by[n] for n in names]
    check_rectangular([(n, r) for n, r in zip(names, rows)], "test alignment")
    keep = [i for i in range(len(rows[0])) if any(r[i].isalpha() for r in rows)]
    return [(n, "".join(r[i] for i in keep)) for n, r in zip(names, rows)]


def prepare(case_row, aln_path):
    """Everything except the kalign call: (ref_dash, projected test, mask, mixed_cols).
    Raises ValueError on any inconsistency between input, reference and test alignment."""
    ref = read_fasta(case_row["reference"])
    if len(ref) < 2:
        raise ValueError(f"{case_row['reference']}: fewer than 2 reference sequences")
    names = [n for n, _ in ref]
    if len(set(names)) != len(names):
        raise ValueError(f"{case_row['reference']}: duplicate sequence names")
    ref_dash = dash_gaps(ref)
    mask, mixed = core_mask(ref_dash, case_row["score_rule"])
    if sum(mask) == 0:
        raise ValueError(f"{case_row['reference']}: no core column (mask empty)")
    test = read_fasta(aln_path)
    # integrity: every input sequence present and unchanged
    inp = dict(read_fasta(case_row["unaligned"]))
    tdict = dict(test)
    bad = [k for k, v in inp.items() if k not in tdict or residues(tdict[k]) != residues(v)]
    extra = set(tdict) - set(inp)
    if bad or extra or len(test) != len(inp):
        raise ValueError(f"{aln_path}: output differs from input (changed/missing {bad[:3]}, "
                         f"extra {sorted(extra)[:3]}, {len(test)} vs {len(inp)} sequences)")
    proj = project(test, names)
    for (n, s), (_, r) in zip(proj, ref_dash):
        if residues(s) != residues(r):
            raise ValueError(f"{aln_path}: residues of {n} differ from the reference")
    return ref_dash, proj, mask, mixed


def score_case(case_row, aln_path) -> dict:
    import kalign

    ref_dash, proj, mask, mixed = prepare(case_row, aln_path)
    with tempfile.TemporaryDirectory() as td:
        rf, tf = Path(td) / "ref.fasta", Path(td) / "test.fasta"
        write_fasta(rf, ref_dash)
        write_fasta(tf, proj)
        r = kalign.compare_detailed(str(rf), str(tf), column_mask=mask)
    return {
        "recall": r["recall"], "precision": r["precision"], "f1": r["f1"], "tc": r["tc"],
        "width_ratio": len(proj[0][1]) / len(ref_dash[0][1]),
        "n_seq": len(read_fasta(case_row["unaligned"])),
        "core_cols": int(sum(mask)), "mixed_cols": mixed,
        "ref_pairs": r["ref_pairs"], "test_pairs": r["test_pairs"],
    }
