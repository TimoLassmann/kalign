"""extHomFam scorer (fast-fit plan, PLAN_fastfit.md section 3): ``score_case(case_row, aln_path)`` for pool
"exthomfam".

The test alignment of a family subset (1,000 / 10,000 / ... sequences) is projected onto the reference sequences
(every other row removed, columns that are gaps in all reference rows dropped) and compared with the HOMSTRAD
reference alignment over all columns: ``kalign.compare_detailed(ref, projected, max_gap_frac=-1.0)``. This is the
scoring of the large-N explorations (largeN_mechanism / refitN ``score``), with the residue check of the protein
scorer added: a reference row whose residues differ from the test row is an error, never silently scored.

Returned dict: recall (SP), precision, f1, tc, width_ratio (projected test width / reference width), n_seq (rows of
the test alignment), n_ref, aln_width (full test alignment width).
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from .scoring_protein import check_rectangular, dash_gaps, read_fasta, residues, write_fasta

SCORER_VERSION = "exthomfam-v1"  # part of the cache key: bump on ANY change to the scoring


def project(ref_recs, test_recs):
    """-> (reference rows dash-gapped, projected test rows), both in reference order."""
    ref = dash_gaps(ref_recs)
    check_rectangular(ref, "reference")
    test = dict((n, s.upper().replace(".", "-")) for n, s in test_recs)
    check_rectangular(list(test.items()), "test alignment")
    rows = []
    for n, s in ref:
        if n not in test:
            raise ValueError(f"reference sequence {n} missing from the test alignment")
        if residues(s) != residues(test[n]):
            raise ValueError(f"reference sequence {n}: residues differ from the test alignment")
        rows.append((n, test[n]))
    keep = [j for j in range(len(rows[0][1])) if any(r[j] != "-" for _, r in rows)]
    return ref, [(n, "".join(r[j] for j in keep)) for n, r in rows]


def score_case(case_row, aln_path) -> dict:
    import kalign
    test_recs = read_fasta(aln_path)
    ref, proj = project(read_fasta(case_row["reference"]), test_recs)
    with tempfile.TemporaryDirectory(dir=None) as d:
        a, b = Path(d) / "ref.fa", Path(d) / "test.fa"
        write_fasta(a, ref)
        write_fasta(b, proj)
        s = kalign.compare_detailed(str(a), str(b), max_gap_frac=-1.0)
    return {"recall": s["recall"], "precision": s["precision"], "f1": s["f1"], "tc": s["tc"],
            "width_ratio": len(proj[0][1]) / len(ref[0][1]), "n_seq": len(test_recs), "n_ref": len(ref),
            "aln_width": len(test_recs[0][1])}
