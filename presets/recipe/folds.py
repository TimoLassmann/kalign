"""Grouped, stratified 5-fold split of a pool and the homolog_in_train covariate.

    cd presets && python -m recipe.folds --pool protein      # needs numba for the covariate (WS/venv-grouping)
    -> <ws>/data/folds_<pool>.tsv   case_id, group_id, fold, benchmark, homolog_in_train, difficulty
       <ws>/data/folds_<pool>.stats.json

Rules (plan + data-final decisions):
  * a group (groups_<pool>.final.tsv) is never split: whole groups are assigned to folds;
  * 5 folds, stratified by benchmark x difficulty: greedy assignment, groups by size (large first, ties by a
    seeded hash), each group goes to the fold where it increases the (cell-count - target)^2 deviation least
    (target = cell total / 5); ties: fewer cases, lower fold index. Seed 1234 (only used for tie-break order);
  * nucleotide pool: a group holding >= FAMILY_MIN BRAliBASE cases is a "whole BRAliBASE family group"; these
    are placed first and one per fold (a fold is reused only after all folds hold one);
  * difficulty = identity bin of the reference (features mean_pid_core: protein core columns, nucleotide all
    columns): protein <20 / 20-30 / 30-40 / >=40 %, nucleotide <50 / 50-70 / 70-85 / >=85 %;
  * homolog_in_train: the case has a sequence linked by the step-0 identity rule (grouping.homolog_flags;
    protein >= 50 % identity over the shorter sequence, nucleotide the length-dependent rule) to a sequence of
    a case in another fold. It is a covariate for reporting, never used to build folds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

from .scoring import read_manifest

SEED = 1234
NFOLDS = 5
FAMILY_MIN = 20
BIN_EDGES = {"protein": (0.20, 0.30, 0.40), "nucleotide": (0.50, 0.70, 0.85)}
BIN_LABELS = {"protein": ("<20", "20-30", "30-40", ">=40"), "nucleotide": ("<50", "50-70", "70-85", ">=85")}


def default_ws() -> Path:
    return Path(os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit")).expanduser()


def hkey(seed, tag, case_id) -> str:
    """Deterministic pseudo-random sort key."""
    return hashlib.sha256(f"{seed}|{tag}|{case_id}".encode()).hexdigest()


def difficulty_bin(pool, pid) -> int:
    """Index of the identity bin; pid = fraction (nan/None -> lowest bin)."""
    if pid is None or pid != pid:
        return 0
    return sum(pid >= e for e in BIN_EDGES[pool])


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_pool(pool, ws=None):
    """Manifest rows joined with group id and features; deterministic order = manifest order."""
    ws = Path(ws) if ws else default_ws()
    rows = read_manifest(ws / "data" / f"cases_{pool}.final.tsv")
    groups = {r["case_id"]: r["group_id"] for r in read_tsv(ws / "data" / f"groups_{pool}.final.tsv")}
    feats = {r["case_id"]: r for r in read_tsv(ws / "data" / f"features_{pool}.final.tsv")}
    out = []
    for r in rows:
        f = feats[r["case_id"]]
        try:
            pid = float(f["mean_pid_core"])
        except ValueError:
            pid = float("nan")
        out.append(dict(case_id=r["case_id"], benchmark=r["benchmark"], subset=r["subset"], group_id=groups[r["case_id"]],
                        pid=pid, bin=difficulty_bin(pool, pid), n_seq=int(r["n_seq"]), ref_len=int(r["ref_len"]),
                        row=r))
    return out


def assign_folds(cases, pool, seed=SEED, nfolds=NFOLDS):
    """-> {group_id: fold}. See the module docstring."""
    members = defaultdict(list)
    for c in cases:
        members[c["group_id"]].append(c)
    cell_total = Counter((c["benchmark"], c["bin"]) for c in cases)
    target = {k: v / nfolds for k, v in cell_total.items()}
    counts = [Counter() for _ in range(nfolds)]
    sizes = [0] * nfolds
    fold_of = {}
    family_used = set()

    def is_family(ms):
        return pool == "nucleotide" and sum(m["benchmark"] == "bralibase" for m in ms) >= FAMILY_MIN

    order = sorted(members, key=lambda g: (not is_family(members[g]), -len(members[g]), hkey(seed, "fold", g)))
    for g in order:
        ms = members[g]
        gc = Counter((m["benchmark"], m["bin"]) for m in ms)
        allowed = list(range(nfolds))
        if is_family(ms):
            free = [f for f in allowed if f not in family_used]
            allowed = free or allowed
        best = None
        for f in allowed:
            inc = sum(((counts[f][k] + v - target[k]) ** 2 - (counts[f][k] - target[k]) ** 2) / (target[k] + 1.0)
                      for k, v in gc.items())
            cand = (round(inc, 9), sizes[f], f)
            if best is None or cand < best:
                best = cand
        f = best[2]
        fold_of[g] = f
        counts[f].update(gc)
        sizes[f] += len(ms)
        if is_family(ms):
            family_used.add(f)
            if len(family_used) == nfolds:
                family_used.clear()
    return fold_of


def build_folds(pool, ws=None, threads=16, log=print, seed=SEED):
    from . import grouping
    cases = load_pool(pool, ws)
    gfold = assign_folds(cases, pool, seed)
    folds = [gfold[c["group_id"]] for c in cases]
    rows = [c["row"] for c in cases]
    flags, info = grouping.homolog_flags(rows, pool, folds, threads=threads, log=log)
    out = []
    for c, f, h in zip(cases, folds, flags.tolist()):
        out.append(dict(case_id=c["case_id"], group_id=c["group_id"], fold=f, benchmark=c["benchmark"],
                        homolog_in_train=int(h), difficulty=BIN_LABELS[pool][c["bin"]]))
    return out, info


def read_folds(pool, ws=None):
    ws = Path(ws) if ws else default_ws()
    return {r["case_id"]: r for r in read_tsv(ws / "data" / f"folds_{pool}.tsv")}


def tables(out):
    """Fold balance tables for the report."""
    t = {}
    t["fold_sizes"] = dict(sorted(Counter(r["fold"] for r in out).items()))
    for key in ("benchmark", "difficulty", "homolog_in_train"):
        c = Counter((r["fold"], r[key]) for r in out)
        t[f"fold_x_{key}"] = {f"{f}|{k}": v for (f, k), v in sorted(c.items(), key=lambda kv: (kv[0][0], str(kv[0][1])))}
    t["homolog_by_benchmark"] = {
        b: dict(n=sum(1 for r in out if r["benchmark"] == b),
                homolog=sum(r["homolog_in_train"] for r in out if r["benchmark"] == b))
        for b in sorted({r["benchmark"] for r in out})}
    return t


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--ws", default=str(default_ws()))
    ap.add_argument("--threads", type=int, default=24)
    a = ap.parse_args(argv)
    try:
        import numba
        numba.set_num_threads(a.threads)
    except Exception:
        pass
    ws = Path(a.ws).expanduser()
    out, info = build_folds(a.pool, ws, threads=a.threads)
    cols = ["case_id", "group_id", "fold", "benchmark", "homolog_in_train", "difficulty"]
    dest = ws / "data" / f"folds_{a.pool}.tsv"
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in out:
            w.writerow([r[c] for c in cols])
    st = tables(out)
    st["homolog_search"] = info
    st["seed"] = SEED
    (ws / "data" / f"folds_{a.pool}.stats.json").write_text(json.dumps(st, indent=1))
    print(json.dumps(st, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
