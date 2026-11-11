"""Fake scorer for the cache/evaluate tests (no kalign import)."""
from pathlib import Path

SCORER_VERSION = "fake-v1"


def score_case(case_row, aln_path):
    n = len(Path(aln_path).read_text())
    if case_row["case_id"].endswith("bad"):
        raise ValueError("synthetic scoring failure")
    return {"recall": 0.5, "precision": 0.25, "f1": 1 / 3, "tc": 0.1, "width_ratio": 1.0,
            "n_seq": int(case_row["n_seq"]), "aln_chars": n}
