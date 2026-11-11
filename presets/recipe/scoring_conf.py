"""Scorer of the confidence re-fit (docs/PLAN_ensconf.md): the pool scorer (recipe.scoring) plus the pair Brier score
of the alignment's column confidence.

Pair Brier score B of one case. Population: every pair of scored residues from two different reference sequences;
scored residues are those in a core column of the reference (the protein scorer's core mask for the small pool, every
column for extHomFam, as the pool scorers). y = 1 when the reference aligns the pair; p = the column confidence of
the test column that holds both residues, 0 when the test alignment does not align them.

    B = sum (p - y)^2 / reference pairs
      = [sum over test columns (correct x (1 - c)^2 + wrong x c^2) + reference pairs the test leaves unaligned]
        / reference pairs

Pairs that neither alignment aligns contribute 0, so every alignment of a case is scored on the same population and
aligning less cannot lower B. Without confidence (no OUT.conf.json: a single alignment) c = 1 and
B = (missed + wrong pairs) / reference pairs. Also returned: the pair counts (correct_pairs / ref_pairs_b equals the
pool scorer's recall; tested in docs/PLAN_ensconf.md, identity checks).
"""
from __future__ import annotations

import json
from pathlib import Path

from .scoring import score_case as pool_score
from .scoring_protein import SCORER_VERSION as V_PROTEIN, core_mask, dash_gaps, read_fasta, residues
from .scoring_largefam import SCORER_VERSION as V_LARGE

SCORER_VERSION = f"conf-v1+{V_PROTEIN}+{V_LARGE}"   # part of the cache key: bump on ANY change to the scoring


def residue_cols(row):
    return [j for j, c in enumerate(row) if c.isalpha()]


def pair_counts(ref_recs, mask, test_recs):
    """-> (reference pairs, [(correct, wrong) per test column]). ref_recs dash-gapped, mask per reference column."""
    test = {n: s for n, s in test_recs}
    width = len(test_recs[0][1])
    cols = [[] for _ in range(width)]
    in_col = {}
    for n, r in ref_recs:
        if n not in test:
            raise ValueError(f"reference sequence {n} missing from the test alignment")
        rc, tc = residue_cols(r), residue_cols(test[n])
        if len(rc) != len(tc) or [r[j].upper() for j in rc] != [test[n][j].upper() for j in tc]:
            raise ValueError(f"residues of {n} differ between reference and test alignment")
        for a, b in zip(rc, tc):
            if mask[a]:
                cols[b].append(a)
                in_col[a] = in_col.get(a, 0) + 1
    n_ref = sum(m * (m - 1) // 2 for m in in_col.values())
    out = []
    for members in cols:
        k = len(members)
        cnt = {}
        for a in members:
            cnt[a] = cnt.get(a, 0) + 1
        correct = sum(v * (v - 1) // 2 for v in cnt.values())
        out.append((correct, k * (k - 1) // 2 - correct))
    return n_ref, out


def pair_brier(n_ref, per_col, conf):
    correct = sum(c for c, _ in per_col)
    num = sum(c * (1.0 - p) ** 2 + w * p * p for (c, w), p in zip(per_col, conf)) + (n_ref - correct)
    return num / n_ref


def score_case(case_row, aln_path) -> dict:
    base = pool_score(case_row, aln_path)
    ref = dash_gaps(read_fasta(case_row["reference"]))
    if case_row["pool"] == "protein":
        mask, _ = core_mask(ref, case_row["score_rule"])
    elif case_row["pool"] in ("exthomfam", "nucleotide"):
        mask = [1] * len(ref[0][1])          # every column (both pool scorers: max_gap_frac -1, all columns)
    else:
        raise ValueError(f"pool {case_row['pool']!r} not scored by the confidence scorer")
    test = [(n, "".join(c if c.isalpha() else "-" for c in s)) for n, s in read_fasta(aln_path)]
    width = len(test[0][1])
    cf = Path(f"{aln_path}.conf.json")
    if cf.exists():
        conf, source = json.loads(cf.read_text()), "ensemble"
        if len(conf) != width:
            raise ValueError(f"{cf}: {len(conf)} confidence values for {width} columns")
    else:
        conf, source = [1.0] * width, "hard"
    n_ref, per_col = pair_counts(ref, mask, test)
    if n_ref == 0:
        raise ValueError(f"{case_row['reference']}: no reference pair in the scored columns")
    return {**base, "brier_pair": pair_brier(n_ref, per_col, conf), "conf_source": source, "ref_pairs_b": n_ref,
            "correct_pairs": sum(c for c, _ in per_col), "wrong_pairs": sum(w for _, w in per_col),
            "mean_col_conf": sum(conf) / width}


__all__ = ["score_case", "pair_counts", "pair_brier", "residues", "SCORER_VERSION"]
