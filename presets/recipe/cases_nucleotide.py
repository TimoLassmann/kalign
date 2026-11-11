"""Nucleotide case manifest: BRAliBASE II + MDSA + a stratified Rfam seed sample.

    python -m recipe.cases_nucleotide [--ws ~/tmp/kalign-refit] [--out WS/data/cases_nucleotide.tsv]
                                      [--rfam-sample 600] [--seed 1234]

Output: TSV with the shared contract columns (recipe.scoring.MANIFEST_COLUMNS) followed by extra columns
(seq_kind, family, ref_identity, stratum, note). Paths are absolute. Files derived here (MDSA unaligned
inputs, Rfam FASTA references and inputs) are written under WS/data/derived/nucleotide/.
Also writes WS/data/cases_nucleotide.stats.json (counts per benchmark, drops, distributions).

Selection rules (all documented in presets/docs/step0/cases_nucleotide.md):
  bralibase  data-set1 (5 families) + data-set2 (pairs); references without any '-' dropped (the old
             pipeline's rule: 'a "-" in any sequence line').
  mdsa       balibase/oxbench/smart "_mdsa_all"; 3..15 sequences, longest sequence line <= 400 columns
             (the old rule), gapless references dropped. Unaligned input = reference with '-' and '.'
             removed (the old datasets._strip_gaps_fasta), headers and case kept.
  rfam       Rfam.seed (Stockholm): 4..200 sequences, alignment length <= 2000, gapless dropped,
             duplicate sequence names dropped; stratified sample by (n_seq bin x identity bin), seed 1234.
Scoring is always on all reference columns (score_rule "all_columns" = max_gap_frac -1.0).

Stdlib only. Pure functions are separated from file I/O so that tests need no data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

MANIFEST_COLUMNS = ("case_id", "pool", "benchmark", "subset", "unaligned", "reference", "score_rule",
                    "n_seq", "ref_len", "sha256_unaligned", "sha256_reference")
EXTRA_COLUMNS = ("seq_kind", "family", "ref_identity", "stratum", "note")
SCORE_RULE = "all_columns"          # kalign.compare_detailed(..., max_gap_frac=-1.0)
POOL = "nucleotide"

MDSA_DBS = ("balibase", "oxbench", "smart")
MDSA_MIN_SEQS, MDSA_MAX_SEQS, MDSA_MAX_ALN_LEN = 3, 15, 400
RFAM_MIN_SEQS, RFAM_MAX_SEQS, RFAM_MAX_ALN_LEN = 4, 200, 2000
RFAM_SEED, RFAM_SAMPLE = 1234, 600
# Strata of the Rfam sample (bin edges chosen after looking at the distribution, see the docs):
NSEQ_BINS = ((4, 5), (6, 9), (10, 24), (25, 200))        # inclusive ranges
IDENT_BINS = ((0.0, 0.60), (0.60, 0.75), (0.75, 0.90), (0.90, 1.0001))   # [lo, hi)
MIN_PER_CELL = 3        # a non-empty cell contributes at least this many (capped by its size)


# ------------------------------------------------------------------ small helpers
def default_ws() -> Path:
    return Path(os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit")).expanduser()


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def read_fasta(path):
    """[(header_without_>, sequence)] in file order; wrapped lines joined; header = whole line."""
    recs, name, buf = [], None, []
    with open(path, encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def has_gap(seqs) -> bool:
    """True if any sequence has a gap character ('-' or '.')."""
    return any("-" in s or "." in s for s in seqs)


def ref_identity(seqs) -> float:
    """Column-pooled mean pairwise identity of an alignment: over all columns, the number of identical
    residue pairs divided by the number of residue pairs (both sequences non-gap). T and U are the same
    letter, case is ignored, '-' and '.' are gaps. Exact, O(columns x sequences) (no pair loop).
    Returns 0.0 if no pair of residues is ever aligned."""
    if not seqs:
        return 0.0
    ncol = max(len(s) for s in seqs)
    same = tot = 0
    for j in range(ncol):
        cnt = Counter()
        for s in seqs:
            if j < len(s):
                c = s[j].upper()
                if c in "-.":
                    continue
                cnt[("U" if c == "T" else c)] += 1
        m = sum(cnt.values())
        tot += m * (m - 1) // 2
        same += sum(k * (k - 1) // 2 for k in cnt.values())
    return same / tot if tot else 0.0


def write_fasta(path, records, width=0):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            if width:
                for i in range(0, len(seq), width):
                    f.write(seq[i:i + width] + "\n")
            else:
                f.write(seq + "\n")
    os.replace(tmp, path)


# ------------------------------------------------------------------ BRAliBASE II
def bralibase_cases(root: Path):
    """(family, case, unaligned_path, reference_path) for all BRAliBASE cases that have a reference.
    root = the 'unpacked' dir of the registry (contains data-set1/data-set1/<family>/ and data-set2/data-set2/)."""
    d1, d2 = root / "data-set1" / "data-set1", root / "data-set2" / "data-set2"
    out = []
    for u in sorted(d1.glob("*/unaligned/*.fa")):
        ref = u.parent.parent / "structural" / u.name
        if ref.exists():
            out.append((u.parent.parent.name, u.stem, u, ref))
    for u in sorted(d2.glob("unaligned/*.fasta")):
        ref = d2 / "structural" / u.name
        if ref.exists():
            out.append(("data-set2", u.stem, u, ref))
    return out


# ------------------------------------------------------------------ MDSA
def mdsa_dir(root: Path, db: str) -> Path:
    return root / f"{db}_mdsa_all" / f"{db}_mdsa_all"


def strip_gaps_fasta(ref_path, out_path):
    """Port of benchmarks/datasets._strip_gaps_fasta: one line per sequence, headers verbatim,
    '-' and '.' removed, case kept."""
    write_fasta(out_path, [(h, s.replace("-", "").replace(".", "")) for h, s in read_fasta(ref_path)])


def mdsa_old_filter(ref_path):
    """The old pipeline's test on one MDSA file. Returns (n_seq, max_line_len, passes_size_filter).
    Old rule (nested_nuc.discover): count '>' lines; longest non-header line (stripped) <= 400."""
    n = 0
    longest = 0
    with open(ref_path, encoding="latin-1") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
            else:
                longest = max(longest, len(line.strip()))
    return n, longest, MDSA_MIN_SEQS <= n <= MDSA_MAX_SEQS and longest <= MDSA_MAX_ALN_LEN


# ------------------------------------------------------------------ Rfam
def parse_stockholm(path):
    """Yield dicts {ac, id, de, tp, rows: [(name, aligned_seq)]} for each family of a Stockholm 1.0 file.
    Handles interleaved blocks (rows with the same name are concatenated, order of first appearance)."""
    cur = None
    with open(path, encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith("# STOCKHOLM"):
                cur = {"ac": "", "id": "", "de": "", "tp": "", "rows": {}}
            elif cur is None or not line.strip():
                continue
            elif line.startswith("//"):
                cur["rows"] = [(k, v) for k, v in cur["rows"].items()]
                yield cur
                cur = None
            elif line.startswith("#=GF"):
                parts = line.split(None, 2)
                if len(parts) == 3:
                    key = {"AC": "ac", "ID": "id", "DE": "de", "TP": "tp"}.get(parts[1])
                    if key and not cur[key]:
                        cur[key] = parts[2].strip()
            elif line.startswith("#"):
                continue
            else:
                name, _, seq = line.partition(" ")
                seq = seq.strip()
                cur["rows"][name] = cur["rows"].get(name, "") + seq


def rfam_clans(path) -> dict:
    """Rfam accession -> clan accession (`#=GF AC` / `#=GF CL` lines of Rfam.seed; families without a clan are absent)."""
    out, ac = {}, None
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#=GF AC"):
                ac = line.split()[2]
            elif line.startswith("#=GF CL") and ac:
                out[ac] = line.split()[2]
            elif line.startswith("//"):
                ac = None
    return out


def rfam_normalise(seq: str) -> str:
    """Rfam row -> reference row: '.' and '-' -> '-', residues upper-case, U/T kept as given."""
    return seq.replace(".", "-").upper()


def rfam_filter(fam):
    """Apply the Rfam family filters. Returns (None, reason) if dropped else (refs, None) where
    refs = [(name, normalised aligned seq)]. Reasons: n_seq<4, n_seq>200, aln_len>2000, ragged,
    duplicate_names, empty_sequence, gapless."""
    rows = fam["rows"]
    n = len(rows)
    if n < RFAM_MIN_SEQS:
        return None, f"n_seq<{RFAM_MIN_SEQS}"
    if n > RFAM_MAX_SEQS:
        return None, f"n_seq>{RFAM_MAX_SEQS}"
    lens = {len(s) for _, s in rows}
    if len(lens) != 1:
        return None, "ragged"
    if lens.pop() > RFAM_MAX_ALN_LEN:
        return None, f"aln_len>{RFAM_MAX_ALN_LEN}"
    refs = [(name, rfam_normalise(s)) for name, s in rows]
    if len({n_ for n_, _ in refs}) != len(refs):      # cannot happen with the dict above; kept as guard
        return None, "duplicate_names"
    if any(not s.replace("-", "") for _, s in refs):
        return None, "empty_sequence"
    if not has_gap([s for _, s in refs]):
        return None, "gapless"
    return refs, None


def nseq_bin(n: int) -> int:
    for i, (lo, hi) in enumerate(NSEQ_BINS):
        if lo <= n <= hi:
            return i
    raise ValueError(f"n_seq {n} outside the bins")


def ident_bin(x: float) -> int:
    for i, (lo, hi) in enumerate(IDENT_BINS):
        if lo <= x < hi:
            return i
    raise ValueError(f"identity {x} outside the bins")


def stratum_label(n: int, ident: float) -> str:
    lo, hi = NSEQ_BINS[nseq_bin(n)]
    ilo, ihi = IDENT_BINS[ident_bin(ident)]
    return f"n{lo}-{hi}|id{ilo:.2f}-{min(ihi, 1.0):.2f}"


def allocate(sizes: dict, total: int, min_per_cell: int = MIN_PER_CELL) -> dict:
    """Proportional allocation of `total` over cells with at least min(size, min_per_cell) each
    (largest remainder), never more than the cell size. If total >= sum(sizes) everything is taken."""
    pop = sum(sizes.values())
    if total >= pop:
        return dict(sizes)
    base = {k: min(v, min_per_cell) for k, v in sizes.items()}
    if sum(base.values()) > total:                       # tiny totals: proportional, no minimum
        base = {k: 0 for k in sizes}
    alloc = dict(base)
    # distribute the remaining quota proportionally to the *remaining* capacity, largest remainder
    while sum(alloc.values()) < total:
        left = total - sum(alloc.values())
        cap = {k: sizes[k] - alloc[k] for k in sizes}
        capsum = sum(cap.values())
        shares = {k: left * cap[k] / capsum for k in sizes if cap[k] > 0}
        given = {k: min(int(s), cap[k]) for k, s in shares.items()}
        for k, g in given.items():
            alloc[k] += g
        left = total - sum(alloc.values())
        if left <= 0:
            break
        rem = sorted(((shares[k] - int(shares[k]), k) for k in shares if alloc[k] < sizes[k]),
                     key=lambda t: (-t[0], str(t[1])))
        for _, k in rem[:left]:
            alloc[k] += 1
    return alloc


def stratified_sample(items, strata, total, seed=RFAM_SEED):
    """items: list of keys (sorted by caller for determinism), strata: parallel list of cell labels.
    Returns the chosen items (in input order) and the allocation dict."""
    cells = defaultdict(list)
    for it, s in zip(items, strata):
        cells[s].append(it)
    alloc = allocate({k: len(v) for k, v in cells.items()}, total)
    rng = random.Random(seed)
    chosen = set()
    for k in sorted(cells):                              # fixed cell order -> reproducible
        pick = rng.sample(cells[k], alloc[k]) if alloc[k] < len(cells[k]) else cells[k]
        chosen.update(pick)
    return [it for it in items if it in chosen], alloc


# ------------------------------------------------------------------ manifest builders
def _row(case_id, benchmark, subset, unal, ref, n_seq, ref_len, kind, family, ident, stratum="", note=""):
    return {
        "case_id": case_id, "pool": POOL, "benchmark": benchmark, "subset": subset,
        "unaligned": str(unal), "reference": str(ref), "score_rule": SCORE_RULE,
        "n_seq": str(n_seq), "ref_len": str(ref_len),
        "sha256_unaligned": sha256_file(unal), "sha256_reference": sha256_file(ref),
        "seq_kind": kind, "family": family, "ref_identity": f"{ident:.4f}", "stratum": stratum, "note": note,
    }


def build_bralibase(raw: Path, stats: dict):
    rows, dropped = [], []
    allc = bralibase_cases(raw / "bralibase" / "unpacked")
    for fam, case, unal, ref in allc:
        seqs = [s for _, s in read_fasta(ref)]
        if not has_gap(seqs):
            dropped.append(f"{fam}/{case}")
            continue
        sub = fam
        cid = f"bralibase/{fam}/{case}"
        rows.append(_row(cid, "bralibase", sub, unal, ref, len(seqs), max(len(s) for s in seqs), "rna", fam,
                         ref_identity(seqs)))
    stats["bralibase"] = {"found": len(allc), "gapless_dropped": len(dropped), "kept": len(rows),
                          "dropped_cases": dropped}
    return rows


def build_mdsa(raw: Path, derived: Path, stats: dict):
    rows = []
    st = {"found": 0, "per_db": {}}
    droplog = {"too_few_or_many_seqs": 0, "aln_len": 0, "gapless": 0}
    for db in MDSA_DBS:
        d = mdsa_dir(raw / "mdsa" / "unpacked", db)
        files = sorted(d.glob("*.afa"))
        pdb = {"found": len(files), "size_filter_pass": 0, "gapless_dropped": 0, "kept": 0}
        for ref in files:
            n, longest, ok = mdsa_old_filter(ref)
            if not ok:
                droplog["too_few_or_many_seqs" if not MDSA_MIN_SEQS <= n <= MDSA_MAX_SEQS else "aln_len"] += 1
                continue
            pdb["size_filter_pass"] += 1
            recs = read_fasta(ref)
            seqs = [s for _, s in recs]
            if not has_gap(seqs):
                pdb["gapless_dropped"] += 1
                droplog["gapless"] += 1
                continue
            unal = derived / "mdsa" / db / ref.name
            strip_gaps_fasta(ref, unal)
            rows.append(_row(f"mdsa/{db}/{ref.stem}", "mdsa", db, unal, ref, len(seqs), max(len(s) for s in seqs),
                             "dna", db, ref_identity(seqs)))
            pdb["kept"] += 1
        st["found"] += len(files)
        st["per_db"][db] = pdb
    st["dropped"] = droplog
    st["kept"] = len(rows)
    stats["mdsa"] = st
    return rows


def build_rfam(raw: Path, derived: Path, stats: dict, sample: int, seed: int):
    path = raw / "rfam" / "Rfam.seed"
    drops, passed, total = Counter(), [], 0
    for fam in parse_stockholm(path):
        total += 1
        refs, why = rfam_filter(fam)
        if why:
            drops[why] += 1
            continue
        ident = ref_identity([s for _, s in refs])
        passed.append((fam["ac"], fam, refs, ident))
    passed.sort(key=lambda t: t[0])
    labels = [stratum_label(len(r), i) for _, _, r, i in passed]
    cells_pop = Counter(labels)
    keys = [p[0] for p in passed]
    chosen, alloc = stratified_sample(keys, labels, sample, seed)
    chosen = set(chosen)
    rows, rng_unused = [], None
    for (ac, fam, refs, ident), lab in zip(passed, labels):
        if ac not in chosen:
            continue
        ref_path = derived / "rfam" / f"{ac}.ref.fa"
        unal_path = derived / "rfam" / f"{ac}.fa"
        write_fasta(ref_path, refs, width=0)
        write_fasta(unal_path, [(n, s.replace("-", "")) for n, s in refs], width=0)
        rows.append(_row(f"rfam/{ac}", "rfam", "seed", unal_path, ref_path, len(refs), len(refs[0][1]), "rna",
                         ac, ident, lab, f"{fam['id']}; {fam['tp'].rstrip(';')}"))
    stats["rfam"] = {
        "families_in_file": total, "dropped": dict(drops), "eligible": len(passed), "sample_target": sample,
        "sampled": len(rows), "seed": seed,
        "cells": {k: {"eligible": cells_pop[k], "sampled": alloc.get(k, 0)} for k in sorted(cells_pop)},
        "ident_eligible": _quantiles([p[3] for p in passed]),
        "ident_sampled": _quantiles([float(r["ref_identity"]) for r in rows]),
    }
    return rows


def _quantiles(xs):
    if not xs:
        return {}
    xs = sorted(xs)
    q = lambda p: xs[min(len(xs) - 1, int(p * (len(xs) - 1) + 0.5))]  # noqa: E731
    return {"min": xs[0], "q25": q(.25), "median": q(.5), "q75": q(.75), "max": xs[-1], "n": len(xs)}


def size_distribution(rows):
    """Per benchmark (and per subset for brevity only at benchmark level): n_seq / ref_len quantiles."""
    out = {}
    for b in sorted({r["benchmark"] for r in rows}):
        rs = [r for r in rows if r["benchmark"] == b]
        out[b] = {"cases": len(rs),
                  "n_seq": _quantiles([int(r["n_seq"]) for r in rs]),
                  "ref_len": _quantiles([int(r["ref_len"]) for r in rs]),
                  "ref_identity": _quantiles([float(r["ref_identity"]) for r in rs])}
    return out


def build_manifest(ws: Path, sample: int = RFAM_SAMPLE, seed: int = RFAM_SEED):
    raw = ws / "data" / "raw"
    derived = ws / "data" / "derived" / "nucleotide"
    stats = {}
    rows = build_bralibase(raw, stats) + build_mdsa(raw, derived, stats) + build_rfam(raw, derived, stats, sample, seed)
    ids = [r["case_id"] for r in rows]
    if len(set(ids)) != len(ids):
        dup = [k for k, v in Counter(ids).items() if v > 1]
        raise ValueError(f"duplicate case_ids: {dup[:5]}")
    stats["total"] = len(rows)
    stats["size_distribution"] = size_distribution(rows)
    return rows, stats


def write_manifest(rows, path):
    cols = list(MANIFEST_COLUMNS) + list(EXTRA_COLUMNS)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = f"{path}.tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ws", default=str(default_ws()))
    ap.add_argument("--out", default=None)
    ap.add_argument("--rfam-sample", type=int, default=RFAM_SAMPLE)
    ap.add_argument("--rfam-all", action="store_true",
                    help="use ALL eligible Rfam families (final pool; overrides --rfam-sample)")
    ap.add_argument("--seed", type=int, default=RFAM_SEED)
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    out = Path(a.out) if a.out else ws / "data" / "cases_nucleotide.tsv"
    rows, stats = build_manifest(ws, 10**9 if a.rfam_all else a.rfam_sample, a.seed)
    write_manifest(rows, out)
    sp = out.with_suffix(".stats.json")
    sp.write_text(json.dumps(stats, indent=1) + "\n")
    print(f"wrote {out} ({len(rows)} cases) and {sp}")
    for b, n in sorted(Counter(r["benchmark"] for r in rows).items()):
        print(f"  {b}: {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
