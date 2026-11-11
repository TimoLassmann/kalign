"""Nucleotide-pool scorer: all reference columns, as the old BRAliBASE / MDSA scoring.

    score_case(case_row, aln_path) -> {recall, precision, f1, tc, width_ratio, n_seq, integrity_bad, ...}

recall = SP, precision, f1, tc come from kalign.compare_detailed(reference, aln, max_gap_frac=-1.0)
(-1.0: every reference column is scored, no gap-fraction filter, no core mask; the manuscript's BRAliBASE
rule; MDSA and Rfam use the same rule). width_ratio = alignment width / reference width
(experiments/redesign/tools/bench_bralibase.py). Extra keys (stored in the cache's `extra` JSON):
integrity_bad = number of sequences whose residues differ from the unaligned input (or that are missing /
extra); any non-zero value raises ValueError (stored as error), so stored rows always have 0; ref_pairs, test_pairs, common_pairs from compare_detailed.

score_rule "all_columns": the test alignment is compared with the reference as is. score_rule "projected" (the
nucfit simulations above 1,000 sequences): the reference holds a subset of the input sequences; the test alignment
is first projected onto those rows (other rows removed, columns that are gaps in every kept row dropped), as the
extHomFam scorer; width_ratio is then projected width / reference width. Anything else raises.
"""

from __future__ import annotations

MAX_GAP_FRAC = -1.0
SCORER_VERSION = "nucleotide-v2"   # v2: score_rule "projected" added; "all_columns" unchanged  # part of the cache key: bump on ANY change to the scoring


def read_fasta(path):
    """[(first word of header, sequence)] in file order (same as experiments/redesign/tools/krcommon)."""
    recs, name, buf = [], None, []
    with open(path, encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                words = line[1:].split()
                name, buf = (words[0] if words else ""), []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def _residues(s: str) -> str:
    return "".join(c for c in s.upper() if c.isalpha())


def integrity_bad(unaligned_path, aln_path) -> int:
    """Sequences of the unaligned input that are missing from the alignment or changed, plus extra
    sequences in the alignment (old krcommon.integrity)."""
    a = dict(read_fasta(unaligned_path))
    b = dict(read_fasta(aln_path))
    bad = [k for k in a if _residues(a[k]) != _residues(b.get(k, ""))]
    return len(bad) + len(set(b) - set(a))


def width(path) -> int:
    recs = read_fasta(path)
    return len(recs[0][1]) if recs else 0


def score_case(case_row, aln_path) -> dict:
    import kalign

    rule = case_row.get("score_rule", "all_columns")
    if rule not in ("all_columns", "projected"):
        raise ValueError(f"nucleotide scorer: unknown score_rule {rule!r}")
    ref = case_row["reference"]
    bad = integrity_bad(case_row["unaligned"], aln_path)
    if bad:  # never score an alignment that changed, dropped or added sequences (as the protein scorer)
        raise ValueError(f"{aln_path}: output differs from input ({bad} sequences changed/missing/extra)")
    if rule == "projected":
        import tempfile
        from pathlib import Path
        test = read_fasta(aln_path)
        refrecs = read_fasta(ref)
        tmap = dict(test)
        rows = []
        for n, s in refrecs:
            if n not in tmap or _residues(s) != _residues(tmap[n]):
                raise ValueError(f"{aln_path}: reference row {n} missing or changed")
            rows.append((n, tmap[n]))
        keep = [j for j in range(len(rows[0][1])) if any(r[j].isalpha() for _, r in rows)]
        with tempfile.TemporaryDirectory() as d:
            proj = Path(d) / "proj.fa"
            with open(proj, "w") as fh:
                for n, r in rows:
                    fh.write(f">{n}\n{''.join(r[j] for j in keep)}\n")
            r = kalign.compare_detailed(str(ref), str(proj), max_gap_frac=MAX_GAP_FRAC)
        w_test = len(keep)
    else:
        r = kalign.compare_detailed(str(ref), str(aln_path), max_gap_frac=MAX_GAP_FRAC)
        w_test = width(aln_path)
    wref = width(ref)
    out = {
        "recall": float(r["recall"]), "precision": float(r["precision"]), "f1": float(r["f1"]),
        "tc": float(r["tc"]),
        "width_ratio": (w_test / wref) if wref else float("nan"),
        "n_seq": len(read_fasta(aln_path)),
        "integrity_bad": bad,  # always 0 here (otherwise raised)
    }
    for k in ("ref_pairs", "test_pairs", "common_pairs"):
        if k in r:
            out[k] = r[k]
    return out
