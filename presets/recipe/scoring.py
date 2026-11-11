"""Scoring dispatcher. The pool scorers are separate modules written elsewhere:

  recipe.scoring_protein.score_case(case_row, aln_path)     -> dict
  recipe.scoring_nucleotide.score_case(case_row, aln_path)  -> dict

Contract: the returned dict has recall, precision, f1, tc, width_ratio, n_seq (floats / int; extra
keys are stored as JSON). Scorers use kalign.compare_detailed(reference, test, ...).
Case manifest columns (TSV): case_id, pool, benchmark, subset, unaligned, reference, score_rule,
n_seq, ref_len, sha256_unaligned, sha256_reference.
"""

from __future__ import annotations

import importlib
import json
import os
import subprocess
import sys
from pathlib import Path

MANIFEST_COLUMNS = ("case_id", "pool", "benchmark", "subset", "unaligned", "reference", "score_rule",
                    "n_seq", "ref_len", "sha256_unaligned", "sha256_reference")
POOL_MODULES = {"protein": "recipe.scoring_protein", "nucleotide": "recipe.scoring_nucleotide",
                "exthomfam": "recipe.scoring_largefam"}   # fast-fit plan (PLAN_fastfit.md)


# Cases with at least SUBPROCESS_MIN_SEQ reference rows are scored in a short-lived child process: the fits score
# from threads, and the compiled scorer holds the GIL, so in-process scoring of large references runs one at a
# time and leaves the process with the transient memory of every large case (2026-11-03). A child scores in
# parallel with the other threads and returns its memory on exit. The child runs the same scorer module
# (python = $KALIGN_SCORE_PYTHON, default this interpreter) and returns the same dict.
SUBPROCESS_MIN_SEQ = 500
SUBPROCESS_TIMEOUT_S = 3600


def _n_seq(case_row):
    try:
        return int(case_row.get("n_seq") or 0)
    except (TypeError, ValueError):
        return 0


def score_case(case_row, aln_path) -> dict:
    pool = case_row["pool"]
    if pool not in POOL_MODULES:
        raise ValueError(f"unknown pool {pool!r} for case {case_row.get('case_id')}")
    if _n_seq(case_row) >= SUBPROCESS_MIN_SEQ and not os.environ.get("KALIGN_SCORE_CHILD"):
        return _score_in_child(case_row, aln_path)
    mod = importlib.import_module(POOL_MODULES[pool])
    return mod.score_case(case_row, aln_path)


def _score_in_child(case_row, aln_path) -> dict:
    py = os.environ.get("KALIGN_SCORE_PYTHON") or sys.executable
    env = {**os.environ, "KALIGN_SCORE_CHILD": "1",
           "PYTHONPATH": os.pathsep.join([str(Path(__file__).resolve().parent.parent),
                                          os.environ.get("PYTHONPATH", "")]).rstrip(os.pathsep)}
    p = subprocess.run([py, "-m", "recipe.scoring", "--child"], input=json.dumps([case_row, str(aln_path)]),
                       capture_output=True, text=True, env=env, timeout=SUBPROCESS_TIMEOUT_S)
    if p.returncode != 0:
        tail = (p.stderr or "").strip().splitlines()[-1:] or [f"exit {p.returncode}"]
        raise RuntimeError(f"scoring child failed: {tail[0][:300]}")
    return json.loads(p.stdout)


def _child_main():
    case_row, aln_path = json.loads(sys.stdin.read())
    sys.stdout.write(json.dumps(score_case(case_row, aln_path)))


def read_manifest(path):
    """Read a case manifest TSV into a list of dicts (strings; n_seq/ref_len left as strings)."""
    import csv
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if rows:
        missing = [c for c in MANIFEST_COLUMNS if c not in rows[0]]
        if missing:
            raise ValueError(f"{path}: missing manifest columns {missing}")
    return rows


if __name__ == "__main__":
    if sys.argv[1:] == ["--child"]:
        _child_main()
    else:
        sys.exit("usage: python -m recipe.scoring --child  (case row and alignment path as JSON on stdin)")
