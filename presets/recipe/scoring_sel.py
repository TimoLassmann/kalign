"""Scorer of amendment 2 (docs/AMENDMENT_ensconf_2.md): recipe.scoring_conf plus
  - pair-level calibration counts: over the scored residue pairs the test alignment aligns, per equal-width bin of
    the column confidence (NB bins over [0, 1]): pairs, sum of confidence, correct pairs (aligned by the reference);
    a hard alignment (no OUT.conf.json) has confidence 1;
  - the selected member of a select_only ensemble (OUT.sel.json), -1 otherwise.
"""
from __future__ import annotations

import json
from pathlib import Path

from .scoring_conf import SCORER_VERSION as V_CONF, pair_counts, score_case as conf_score
from .scoring_protein import core_mask, dash_gaps, read_fasta

SCORER_VERSION = f"sel-v1+{V_CONF}"   # part of the cache key: bump on ANY change to the scoring
NB = 10


def calib_bins(per_col, conf):
    n, sp, ny = [0] * NB, [0.0] * NB, [0] * NB
    for (c, w), p in zip(per_col, conf):
        k = c + w
        if k == 0:
            continue
        b = min(int(p * NB), NB - 1)
        n[b] += k
        sp[b] += p * k
        ny[b] += c
    return {"n": n, "sum_p": sp, "correct": ny}


def score_case(case_row, aln_path) -> dict:
    base = conf_score(case_row, aln_path)
    ref = dash_gaps(read_fasta(case_row["reference"]))
    if case_row["pool"] == "protein":
        mask, _ = core_mask(ref, case_row["score_rule"])
    else:
        mask = [1] * len(ref[0][1])
    test = [(n, "".join(c if c.isalpha() else "-" for c in s)) for n, s in read_fasta(aln_path)]
    width = len(test[0][1])
    cf = Path(f"{aln_path}.conf.json")
    conf = json.loads(cf.read_text()) if cf.exists() else [1.0] * width
    _, per_col = pair_counts(ref, mask, test)
    sf = Path(f"{aln_path}.sel.json")
    sel = json.loads(sf.read_text())["selected_member"] if sf.exists() else -1
    return {**base, "calib": calib_bins(per_col, conf), "selected_member": sel}


__all__ = ["score_case", "calib_bins", "SCORER_VERSION"]
