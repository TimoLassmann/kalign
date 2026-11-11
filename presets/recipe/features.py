"""Case features computed from the REFERENCE alignment only (never from any aligner output).

    cd presets && uv run python -m recipe.features --pool protein [--ws ~/tmp/kalign-refit] [--jobs N]
    -> <ws>/data/features_<pool>.tsv

Per case (one row; columns are FEATURE_COLUMNS):
  n_seq_input      sequences in the unaligned input (manifest n_seq; the cost driver)
  n_seq            sequences in the reference alignment (prefab4: 2; the input holds ~50)
  ref_len          reference alignment width (columns)
  len_mean/sd/min/max   ungapped length of the reference sequences (sd = population sd)
  mean_pid         mean over all sequence pairs of (identical residues / columns where both have a
                   residue), case-insensitive, T=U for nucleotides; pairs without a shared column
                   are skipped (NaN if none)
  mean_pid_core    the same restricted to the scoring mask columns (protein: score_rule mask from
                   scoring_protein.core_mask; nucleotide: all columns, so equal to mean_pid)
  core_cols, core_frac   number / fraction of mask columns
  gap_frac         gap cells / all cells of the reference matrix (every non-letter is a gap)
  n_gap_runs       internal gap runs summed over sequences (a run = maximal stretch of gap cells in
                   one row strictly between its first and last residue)
  gap_runs_per_seq, gap_run_len_mean, gap_run_len_max   (len_mean/max NaN/0 if no run)
  internal_gap_frac   cells in internal gap runs / all cells
  lead_mean/lead_max, trail_mean/trail_max   gap cells before first / after last residue per
                   sequence (columns); term_ext_frac = (all lead + trail cells) / all cells
Because the reference of a drive5 case may hold fewer sequences than the input, all features describe
the reference only; n_seq_input is carried over from the manifest for convenience.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from .scoring import read_manifest
from .scoring_protein import core_mask, dash_gaps, read_fasta

FEATURE_COLUMNS = (
    "case_id", "benchmark", "subset", "n_seq_input", "n_seq", "ref_len", "len_mean", "len_sd",
    "len_min", "len_max", "mean_pid", "mean_pid_core", "core_cols", "core_frac", "gap_frac",
    "n_gap_runs", "gap_runs_per_seq", "gap_run_len_mean", "gap_run_len_max", "internal_gap_frac",
    "lead_mean", "lead_max", "trail_mean", "trail_max", "term_ext_frac",
)


def _matrix(rows, nucleotide: bool):
    """uint8 matrix (n, L) of upper-case letter codes, 0 = gap. T and U are merged for nucleotides."""
    n, L = len(rows), len(rows[0])
    m = np.zeros((n, L), dtype=np.uint8)
    for i, s in enumerate(rows):
        a = np.frombuffer(s.upper().encode("latin-1"), dtype=np.uint8)
        isl = ((a >= 65) & (a <= 90))
        if nucleotide:
            a = np.where(a == ord("U"), ord("T"), a)
        m[i] = np.where(isl, a, 0)
    return m


def mean_pairwise_identity(m: np.ndarray, mask=None):
    """Mean over pairs of identical/both-present; mask = optional bool vector over columns."""
    if mask is not None:
        m = m[:, mask]
    n = m.shape[0]
    if n < 2 or m.shape[1] == 0:
        return float("nan")
    pres = m > 0
    vals = []
    for i in range(n - 1):
        both = pres[i] & pres[i + 1:]
        nb = both.sum(axis=1)
        ident = (both & (m[i] == m[i + 1:])).sum(axis=1)
        ok = nb > 0
        vals.extend((ident[ok] / nb[ok]).tolist())
    return float(np.mean(vals)) if vals else float("nan")


def gap_structure(m: np.ndarray):
    """Per-row leading / trailing gap cells and internal gap runs (list of run lengths)."""
    pres = m > 0
    n, L = m.shape
    lead = np.zeros(n, dtype=np.int64)
    trail = np.zeros(n, dtype=np.int64)
    runs = []
    for i in range(n):
        idx = np.flatnonzero(pres[i])
        if idx.size == 0:  # row with no residue at all: all of it is "leading"
            lead[i] = L
            continue
        lead[i], trail[i] = idx[0], L - 1 - idx[-1]
        gaps = np.diff(idx) - 1  # gap cells between consecutive residues
        runs.extend(gaps[gaps > 0].tolist())
    return lead, trail, runs


def features_from_rows(rows, nucleotide: bool, mask=None) -> dict:
    """Features of a dash-gapped reference (list of equal-length strings; case allowed)."""
    m = _matrix(rows, nucleotide)
    n, L = m.shape
    cells = n * L
    pres = m > 0
    ulen = pres.sum(axis=1)
    lead, trail, runs = gap_structure(m)
    runs = np.asarray(runs, dtype=np.int64)
    mask_b = None if mask is None else np.asarray(mask, dtype=bool)
    core_cols = int(L if mask_b is None else mask_b.sum())
    return dict(
        n_seq=n, ref_len=L, len_mean=float(ulen.mean()), len_sd=float(ulen.std()),
        len_min=int(ulen.min()), len_max=int(ulen.max()),
        mean_pid=mean_pairwise_identity(m),
        mean_pid_core=mean_pairwise_identity(m, mask_b) if mask_b is not None else mean_pairwise_identity(m),
        core_cols=core_cols, core_frac=core_cols / L,
        gap_frac=float(1.0 - pres.sum() / cells),
        n_gap_runs=int(runs.size), gap_runs_per_seq=float(runs.size / n),
        gap_run_len_mean=float(runs.mean()) if runs.size else float("nan"),
        gap_run_len_max=int(runs.max()) if runs.size else 0,
        internal_gap_frac=float(runs.sum() / cells),
        lead_mean=float(lead.mean()), lead_max=int(lead.max()),
        trail_mean=float(trail.mean()), trail_max=int(trail.max()),
        term_ext_frac=float((lead.sum() + trail.sum()) / cells),
    )


def case_features(row: dict) -> dict:
    """Features of one manifest row (reads only the reference file)."""
    nucleotide = row["pool"] == "nucleotide"
    ref = read_fasta(row["reference"])
    if len(ref) < 2:
        raise ValueError(f"{row['reference']}: fewer than 2 reference sequences")
    rd = dash_gaps(ref)
    widths = {len(s) for _, s in rd}
    if len(widths) != 1:
        raise ValueError(f"{row['reference']}: ragged reference")
    mask = None
    if not nucleotide and row.get("score_rule") in ("core_columns", "pair_core"):
        mask, _ = core_mask(rd, row["score_rule"])
    f = features_from_rows([s for _, s in rd], nucleotide, mask)
    f.update(case_id=row["case_id"], benchmark=row["benchmark"], subset=row["subset"],
             n_seq_input=int(row["n_seq"]))
    return f


def _fmt(v):
    if isinstance(v, float):
        return "" if math.isnan(v) else f"{v:.6g}"
    return str(v)


def write_features(rows_out, path):
    tmp = str(path) + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(FEATURE_COLUMNS)
        for f in rows_out:
            w.writerow([_fmt(f[c]) for c in FEATURE_COLUMNS])
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--ws", default=os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit"))
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--final", action="store_true",
                    help="use cases_<pool>.final.tsv and write features_<pool>.final.tsv")
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    suffix = ".final" if a.final else ""
    rows = read_manifest(ws / "data" / f"cases_{a.pool}{suffix}.tsv")
    with ProcessPoolExecutor(a.jobs) as ex:
        out = list(ex.map(case_features, rows, chunksize=16))  # keeps manifest order
    dest = ws / "data" / f"features_{a.pool}{suffix}.tsv"
    write_features(out, dest)
    print(f"wrote {len(out)} rows to {dest}")


if __name__ == "__main__":
    sys.exit(main())
