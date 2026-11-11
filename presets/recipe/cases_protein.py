"""Protein case manifest: WS/data/cases_protein.tsv (columns = recipe.scoring.MANIFEST_COLUMNS).

    cd presets && uv run python -m recipe.cases_protein [--ws ~/tmp/kalign-refit] [--check]

Sources (all under <ws>/data/raw/, see docs/step0/data.md):
  bali3    drive5 BAliBASE 3, BB (218) and BBS (168) cases; subset BB-RV11.. / BBS-RV11..
  prefab4  drive5 PREFAB 4, input = full set (~50 sequences), reference = the pair (prefab4ref is the
           same 1,681 references with a 2-sequence input; not used)
  ox       drive5 OXBench (oxx / MUSTANG variants / bali3pdb are not used: overlap or ill-defined)
  sabre    drive5 SABRE
  balifam  BaliFam100 (59 families, reference subset inside 104-242 Pfam sequences)
Subsets of ox, sabre, prefab4: percent identity of the core columns from info/pctids.txt
(pid<20, pid20-30, pid30-40, pid>=40).

Unaligned inputs are written as normalised copies (upper case letters only, names = first header
token) to <ws>/data/protein/in/<benchmark>/<id>.fa. References stay in data/raw untouched (the manifest
holds their sha256); the scorer normalises gaps on reading. Cases that cannot be scored are dropped
and listed with the reason in <ws>/data/cases_protein.dropped.tsv; statistics go to
<ws>/data/cases_protein.stats.json. Deterministic: same raw data -> same manifest bytes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import statistics
from pathlib import Path

from .scoring import MANIFEST_COLUMNS
from .scoring_protein import core_mask, dash_gaps, read_fasta, residues

DRIVE5 = ("data", "raw", "drive5_bench", "unpacked", "bench", "bench1.0")
BALIFAM = ("data", "raw", "balifam", "repo", "balifam100")
# benchmark -> (source dir name, score rule, subset kind)
DRIVE5_SETS = {
    "bali3": ("bali3", "core_columns", "bali"),
    "prefab4": ("prefab4", "pair_core", "pid"),
    "ox": ("ox", "core_columns", "pid"),
    "sabre": ("sabre", "core_columns", "pid"),
}


def default_ws() -> Path:
    return Path(os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit")).expanduser()


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def pid_bin(p: float) -> str:
    return "pid<20" if p < 20 else "pid20-30" if p < 30 else "pid30-40" if p < 40 else "pid>=40"


def bali_subset(case_id: str) -> str:
    kind = "BBS" if case_id.startswith("BBS") else "BB"
    digits = case_id[len(kind):]
    return f"{kind}-RV{digits[:2]}"


def read_pctids(info_dir: Path) -> dict:
    out = {}
    p = info_dir / "pctids.txt"
    if p.exists():
        for line in p.read_text().split():
            if ";" in line:
                k, v = line.split(";")
                out[k] = float(v)
    return out


def read_ids(info_dir: Path) -> list:
    return info_dir.joinpath("ids.txt").read_text().split()


def validate_case(ref_recs, in_recs, rule):
    """Return a drop reason (str) or None. Also checks the reference against the input."""
    if len(ref_recs) < 2:
        return "reference has fewer than 2 sequences"
    names = [n for n, _ in ref_recs]
    if len(set(names)) != len(names):
        return "duplicate names in reference"
    inn = [n for n, _ in in_recs]
    if len(set(inn)) != len(inn):
        return "duplicate names in input"
    if len({len(s) for _, s in ref_recs}) != 1:
        return "reference rows differ in length"
    ind = dict(in_recs)
    miss = [n for n in names if n not in ind]
    if miss:
        return f"reference sequences not in input ({miss[0]})"
    for n, s in ref_recs:
        if residues(s) != residues(ind[n]):
            return f"reference residues of {n} differ from input"
    rd = dash_gaps(ref_recs)
    mask, _ = core_mask(rd, rule)
    if sum(mask) == 0:
        return "no core (upper-case) column in reference"
    # at least one core column with >= 2 residues, else no pair is assessed
    if not any(m and sum(s[j].isalpha() for _, s in rd) >= 2 for j, m in enumerate(mask)):
        return "no assessed residue pair in core columns"
    return None


def write_input(path: Path, in_recs):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for n, s in in_recs:
            fh.write(f">{n}\n{residues(s)}\n")


def collect(ws: Path):
    """Yield (benchmark, subset, id, in_path, ref_path, rule, raw_in_path) in a deterministic order."""
    d5 = ws.joinpath(*DRIVE5)
    for bench, (dname, rule, kind) in DRIVE5_SETS.items():
        root = d5 / dname
        pid = read_pctids(root / "info") if kind == "pid" else {}
        for cid in sorted(read_ids(root / "info")):
            if kind == "bali":
                sub = bali_subset(cid)
            else:
                if cid not in pid:
                    raise ValueError(f"{bench}/{cid}: no percent identity in info/pctids.txt")
                sub = pid_bin(pid[cid])
            yield bench, sub, cid, root / "in" / cid, root / "ref" / cid, rule
    root = ws.joinpath(*BALIFAM)
    for cid in sorted(read_ids(root / "info")):
        yield "balifam", "balifam100", cid, root / "in" / cid, root / "ref" / cid, "core_columns"


def build(ws: Path, out: Path | None = None):
    out = out or ws / "data" / "cases_protein.tsv"
    rows, dropped = [], []
    for bench, sub, cid, in_p, ref_p, rule in collect(ws):
        in_recs, ref_recs = read_fasta(in_p), read_fasta(ref_p)
        why = validate_case(ref_recs, in_recs, rule)
        if why:
            dropped.append((f"{bench}.{cid}", bench, sub, why))
            continue
        norm = ws / "data" / "protein" / "in" / bench / f"{cid}.fa"
        write_input(norm, in_recs)
        rows.append({
            "case_id": f"{bench}.{cid}", "pool": "protein", "benchmark": bench, "subset": sub,
            "unaligned": str(norm), "reference": str(ref_p), "score_rule": rule,
            "n_seq": str(len(in_recs)), "ref_len": str(len(ref_recs[0][1])),
            "sha256_unaligned": sha256_file(norm), "sha256_reference": sha256_file(ref_p),
            # not manifest columns, used for the stats only
            "_n_ref": len(ref_recs), "_gapless": all("." not in s and "-" not in s for _, s in ref_recs),
            "_n_in_raw": len(in_recs),
        })
    ids = [r["case_id"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate case ids"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MANIFEST_COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in MANIFEST_COLUMNS])
    with open(out.with_name("cases_protein.dropped.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["case_id", "benchmark", "subset", "reason"])
        w.writerows(dropped)
    stats = summarise(rows, dropped)
    out.with_name("cases_protein.stats.json").write_text(json.dumps(stats, indent=1, sort_keys=True))
    return rows, dropped, stats


def _dup_key(row):
    """Content key of a case: ungapped input sequences + the reference alignment rows (names ignored)."""
    inp = tuple(sorted(sq.upper() for _, sq in read_fasta(row["unaligned"])))
    ref = tuple(sorted(sq for _, sq in read_fasta(row["reference"])))  # case = core marks
    return inp, ref


def find_exact_duplicates(rows):
    """-> {dropped case_id: kept case_id}. Cases with identical input sequences and identical reference
    alignment rows (names ignored) are duplicates; the first by case_id is kept (BB40016 over BB40017,
    sabre.sup_NNN over sabre.twi_MMM). Cases with equal reference file hashes are checked first (cheap)."""
    seen, drop = {}, {}
    for r in sorted(rows, key=lambda r: r["case_id"]):
        k = _dup_key(r)
        if k in seen:
            drop[r["case_id"]] = seen[k]
        else:
            seen[k] = r["case_id"]
    return drop


def make_final(ws: Path):
    """Final protein manifest: the step-0 manifest minus exact-duplicate cases."""
    from .scoring import read_manifest
    src = ws / "data" / "cases_protein.tsv"
    rows = read_manifest(src)
    drop = find_exact_duplicates(rows)
    keep = [r for r in rows if r["case_id"] not in drop]
    out = ws / "data" / "cases_protein.final.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MANIFEST_COLUMNS)
        for r in keep:
            w.writerow([r[c] for c in MANIFEST_COLUMNS])
    with open(ws / "data" / "cases_protein.final.dropped.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["case_id", "identical_to"])
        w.writerows(sorted(drop.items()))
    return keep, drop


def _dist(v):
    v = sorted(v)
    return {"min": v[0], "median": statistics.median(v), "max": v[-1]} if v else {}


def summarise(rows, dropped):
    stats = {"n_cases": len(rows), "n_dropped": len(dropped), "benchmarks": {}}
    for bench in sorted({r["benchmark"] for r in rows} | {d[1] for d in dropped}):
        rs = [r for r in rows if r["benchmark"] == bench]
        entry = {"n": len(rs), "dropped": sum(1 for d in dropped if d[1] == bench), "subsets": {}}
        entry["n_seq"] = _dist([int(r["n_seq"]) for r in rs])
        entry["n_ref_seq"] = _dist([r["_n_ref"] for r in rs])
        entry["ref_len"] = _dist([int(r["ref_len"]) for r in rs])
        entry["gapless_refs"] = sum(1 for r in rs if r["_gapless"])
        for sub in sorted({r["subset"] for r in rs}):
            ss = [r for r in rs if r["subset"] == sub]
            entry["subsets"][sub] = {"n": len(ss), "n_seq": _dist([int(r["n_seq"]) for r in ss]),
                                     "ref_len": _dist([int(r["ref_len"]) for r in ss]),
                                     "gapless_refs": sum(1 for r in ss if r["_gapless"])}
        stats["benchmarks"][bench] = entry
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ws", default=str(default_ws()))
    ap.add_argument("--out")
    ap.add_argument("--final", action="store_true",
                    help="write cases_protein.final.tsv = existing manifest minus exact-duplicate cases")
    a = ap.parse_args()
    if a.final:
        keep, drop = make_final(Path(a.ws).expanduser())
        print(f"final protein manifest: {len(keep)} cases, {len(drop)} exact duplicates dropped")
        for k, v in sorted(drop.items()):
            print("DROPPED", k, "identical to", v)
        return
    rows, dropped, stats = build(Path(a.ws).expanduser(), Path(a.out) if a.out else None)
    print(json.dumps(stats, indent=1, sort_keys=True))
    for d in dropped:
        print("DROPPED", *d)


if __name__ == "__main__":
    main()
