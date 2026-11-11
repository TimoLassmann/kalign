"""Training cores and screening subsets for the nested-CV fits.

    cd presets && python -m recipe.core [--pool protein|nucleotide|both] [--ws WS]
    -> <ws>/data/cores/<pool>_fold<k>_core.txt     case ids of the training core of fold k (k = 0..4 or "all")
       <ws>/data/cores/<pool>_fold<k>_screen.txt   the screening subset (a subset of the core, search stage only)
       <ws>/data/cores/<pool>_fold<k>_guard.txt    protein only: the guard set (part of the core, not of screen)
       <ws>/data/cores/composition.json            counts per fold / benchmark / difficulty bin

Python API: training_core(pool, fold), screening_subset(pool, fold), guard_set(pool, fold); fold is an int or "all"
(= the final fit: the whole pool is the training population, same rules).

All draws are deterministic: a case is ranked by sha256(seed|purpose|case_id), so selections are stable and the
same case is preferred in every fold (nested where the populations overlap). Seed 1234.

Rules
  training population  = cases of the pool not in the test fold (all cases for "all");
  protein: cases whose 3.7.0 `accurate` run took > 10 s CPU (or failed) are excluded from cores (still tested);
           cap: at most 5 cases per group (survivors of a full group: balifam, bali3, ox, sabre, prefab4 in this
           priority order, then hash order); core = 800 cases + 40 guard cases:
           all BaliFam cases of the population, the rest (800 - n_balifam) split evenly over bali3, prefab4, ox,
           sabre (a benchmark with too few cases gives its share to the others); inside a benchmark the identity
           bins <20 / 20-30 / 30-40 / >=40 % get quota weights 1:1:1:0.5 (same redistribution); the guard set is
           40 further cases drawn at random (hash order) from the >=40 % bin of the capped population, drawn
           first and removed from the main draw.
           screening = 300 cases of the 800, stratified benchmark x bin (proportional, >= 3 per non-empty cell).
  nucleotide: cap: at most 25 cases per group; core = all remaining training cases;
           screening = 300 cases, equal share per benchmark (bralibase, mdsa, rfam; redistributed if a benchmark
           is small), inside a benchmark proportional over the difficulty bins (>= 3 per non-empty cell).
"""
from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

from .cases_nucleotide import allocate
from .folds import BIN_LABELS, SEED, NFOLDS, hkey, load_pool, read_folds

CAP = {"protein": 5, "nucleotide": 25}
CORE_SIZE = 800
GUARD_SIZE = 40
SCREEN_SIZE = 300
BIN_WEIGHTS = (1.0, 1.0, 1.0, 0.5)
CPU_LIMIT_S = 10.0
EVEN_BENCHMARKS = ("bali3", "prefab4", "ox", "sabre")
# Which cases of an over-full group survive the cap: the scarcer benchmarks first (the large ref-contained
# groups mix all benchmarks and a purely random cap would thin BaliFam / BAliBASE far below their quotas).
CAP_PRIORITY = {"balifam": 0, "bali3": 1, "ox": 2, "sabre": 3, "prefab4": 4}
FOLDS = [0, 1, 2, 3, 4, "all"]


def default_ws() -> Path:
    return Path(os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit")).expanduser()


def load_cpu_accurate(ws=None, cache_path=None, build_id=None):
    """{case_id: cpu_s of the 3.7.0 accurate run, None if it failed} from the step-0 cache."""
    ws = Path(ws) if ws else default_ws()
    cache_path = Path(cache_path) if cache_path else ws / "cache" / "cache.sqlite"
    key_json = '{"preset":"accurate"}'
    con = sqlite3.connect(f"file:{cache_path}?mode=ro", uri=True)
    try:
        key = con.execute("SELECT param_key FROM specs WHERE spec_json=?", (key_json,)).fetchone()[0]
        builds = [b for (b,) in con.execute("SELECT DISTINCT build_id FROM results WHERE param_key=?", (key,))]
        if build_id is None:
            if len(builds) != 1:
                raise ValueError(f"cache holds {len(builds)} builds for accurate; pass build_id")
            build_id = builds[0]
        cpu = {}
        for cid, c, err in con.execute("SELECT case_id, cpu_s, error FROM results WHERE build_id=? AND param_key=?",
                                       (build_id, key)):
            cpu[cid] = None if err else c
    finally:
        con.close()
    return cpu


def waterfill(total, avail, weights):
    """Integer allocation of `total` over keys proportional to weights, never above avail[key]; what a
    saturated key cannot take is redistributed over the others. Largest remainder, ties by key."""
    keys = sorted(avail)
    alloc = {k: 0 for k in keys}
    total = min(total, sum(avail.values()))
    open_ = [k for k in keys if avail[k] > 0 and weights[k] > 0]
    left = total
    while left > 0 and open_:
        wsum = sum(weights[k] for k in open_)
        share = {k: left * weights[k] / wsum for k in open_}
        sat = [k for k in open_ if alloc[k] + share[k] >= avail[k]]
        if sat:
            for k in sat:
                left -= avail[k] - alloc[k]
                alloc[k] = avail[k]
            open_ = [k for k in open_ if k not in sat]
            continue
        base = {k: int(share[k]) for k in open_}
        for k in open_:
            alloc[k] += base[k]
        left -= sum(base.values())
        rem = sorted(open_, key=lambda k: (-(share[k] - base[k]), str(k)))
        for k in rem[:left]:
            alloc[k] += 1
        left = 0
    return alloc


def training_population(cases, fold, folds):
    if fold == "all":
        return list(cases)
    return [c for c in cases if int(folds[c["case_id"]]["fold"]) != int(fold)]


def apply_cap(pop, pool, seed=SEED):
    by_group = defaultdict(list)
    for c in pop:
        by_group[c["group_id"]].append(c)
    kept = []
    for g, ms in by_group.items():
        ms.sort(key=lambda c: (CAP_PRIORITY.get(c["benchmark"], 0), hkey(seed, "cap", c["case_id"])))
        kept += ms[:CAP[pool]]
    return sorted(kept, key=lambda c: c["case_id"])


def protein_core(pop, seed=SEED):
    """-> dict(main=[cases], guard=[cases]) from the capped, CPU-filtered population. The guard set is drawn
    first (so it always has GUARD_SIZE cases when the >=40 % bin allows), the main core from the rest."""
    hi = sorted((c for c in pop if c["bin"] == 3), key=lambda c: hkey(seed, "guard", c["case_id"]))
    guard = hi[:GUARD_SIZE]
    gset = {c["case_id"] for c in guard}
    by_b = defaultdict(list)
    for c in pop:
        if c["case_id"] not in gset:
            by_b[c["benchmark"]].append(c)
    main = list(by_b.get("balifam", []))
    rest = CORE_SIZE - len(main)
    avail = {b: len(by_b.get(b, [])) for b in EVEN_BENCHMARKS}
    bq = waterfill(rest, avail, {b: 1.0 for b in EVEN_BENCHMARKS})
    for b in EVEN_BENCHMARKS:
        cells = defaultdict(list)
        for c in by_b.get(b, []):
            cells[c["bin"]].append(c)
        cq = waterfill(bq[b], {i: len(cells[i]) for i in range(4)}, {i: BIN_WEIGHTS[i] for i in range(4)})
        for i in range(4):
            cells[i].sort(key=lambda c: hkey(seed, "core", c["case_id"]))
            main += cells[i][:cq[i]]
    return dict(main=sorted(main, key=lambda c: c["case_id"]), guard=sorted(guard, key=lambda c: c["case_id"]))


def stratified_pick(cases, total, bin_key, seed, tag):
    cells = defaultdict(list)
    for c in cases:
        cells[(c["benchmark"], c["bin"])].append(c)
    alloc = allocate({k: len(v) for k, v in cells.items()}, total)
    out = []
    for k, v in cells.items():
        v = sorted(v, key=lambda c: hkey(seed, tag, c["case_id"]))
        out += v[:alloc[k]]
    return out


def screen_protein(main, seed=SEED):
    return sorted(stratified_pick(main, SCREEN_SIZE, "bin", seed, "screen"), key=lambda c: c["case_id"])


def screen_nucleotide(core, seed=SEED):
    by_b = defaultdict(list)
    for c in core:
        by_b[c["benchmark"]].append(c)
    bq = waterfill(SCREEN_SIZE, {b: len(v) for b, v in by_b.items()}, {b: 1.0 for b in by_b})
    out = []
    for b, v in by_b.items():
        out += stratified_pick(v, bq[b], "bin", seed, "screen")
    return sorted(out, key=lambda c: c["case_id"])


def build_core(pool, fold, ws=None, data=None, cpu=None, folds=None):
    """-> dict(core=[ids], screen=[ids], guard=[ids], info={...}) with all invariants asserted."""
    ws = Path(ws) if ws else default_ws()
    cases = data if data is not None else load_pool(pool, ws)
    folds = folds if folds is not None else read_folds(pool, ws)
    pop = training_population(cases, fold, folds)
    info = dict(pool=pool, fold=fold, population=len(pop))
    if pool == "protein":
        cpu = cpu if cpu is not None else load_cpu_accurate(ws)
        eligible = [c for c in pop if cpu.get(c["case_id"]) is not None and cpu[c["case_id"]] <= CPU_LIMIT_S]
        info["excluded_cpu"] = len(pop) - len(eligible)
        capped = apply_cap(eligible, pool)
        info["after_cap"] = len(capped)
        parts = protein_core(capped)
        main, guard = parts["main"], parts["guard"]
        core = sorted(main + guard, key=lambda c: c["case_id"])
        screen = screen_protein(main)
    else:
        capped = apply_cap(pop, pool)
        info["after_cap"] = len(capped)
        core, guard = capped, []
        screen = screen_nucleotide(core)
    # invariants
    test_ids = set() if fold == "all" else {k for k, r in folds.items() if int(r["fold"]) == int(fold)}
    ids = [c["case_id"] for c in core]
    assert len(set(ids)) == len(ids), "duplicate case in core"
    assert not (set(ids) & test_ids), f"test-fold case in the training core ({pool}, fold {fold})"
    assert {c["case_id"] for c in screen} <= set(ids), "screen not inside core"
    assert not ({c["case_id"] for c in guard} & {c["case_id"] for c in screen}), "guard overlaps screen"
    per_group = Counter(c["group_id"] for c in core)
    assert max(per_group.values()) <= CAP[pool], "cap violated"
    if pool == "protein":
        assert len(core) <= CORE_SIZE + GUARD_SIZE
        assert all(cpu[i] <= CPU_LIMIT_S for i in ids), "expensive case in core"
    info.update(core=len(core), screen=len(screen), guard=len(guard))
    return dict(core=ids, screen=[c["case_id"] for c in screen], guard=[c["case_id"] for c in guard],
                info=info, cases={c["case_id"]: c for c in core})


_CACHE = {}


def _get(pool, fold, ws):
    k = (pool, str(fold), str(ws))
    if k not in _CACHE:
        _CACHE[k] = build_core(pool, fold, ws)
    return _CACHE[k]


def training_core(pool, fold, ws=None):
    """Case ids of the training core (protein: 800 + guard)."""
    return list(_get(pool, fold, ws)["core"])


def screening_subset(pool, fold, ws=None):
    return list(_get(pool, fold, ws)["screen"])


def guard_set(pool, fold, ws=None):
    return list(_get(pool, fold, ws)["guard"])


def composition(res, pool):
    """benchmark x bin counts for core / screen / guard of one build_core result."""
    cmap = res["cases"]
    allc = {c["case_id"]: c for c in cmap.values()}
    out = {}
    for name in ("core", "screen", "guard"):
        cnt = Counter()
        for i in res[name]:
            c = allc[i]
            cnt[(c["benchmark"], BIN_LABELS[pool][c["bin"]])] += 1
        out[name] = {f"{b}|{lab}": n for (b, lab), n in sorted(cnt.items())}
        out[name + "_by_benchmark"] = dict(sorted(Counter(allc[i]["benchmark"] for i in res[name]).items()))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", default="both", choices=["protein", "nucleotide", "both"])
    ap.add_argument("--ws", default=str(default_ws()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    d = ws / "data" / "cores"
    d.mkdir(parents=True, exist_ok=True)
    comp = {}
    for pool in (["protein", "nucleotide"] if a.pool == "both" else [a.pool]):
        data = load_pool(pool, ws)
        folds = read_folds(pool, ws)
        cpu = load_cpu_accurate(ws) if pool == "protein" else None
        for fold in FOLDS:
            res = build_core(pool, fold, ws, data=data, cpu=cpu, folds=folds)
            for name in ("core", "screen") + (("guard",) if pool == "protein" else ()):
                (d / f"{pool}_fold{fold}_{name}.txt").write_text("".join(i + "\n" for i in res[name]))
            comp[f"{pool}_fold{fold}"] = dict(info=res["info"], **composition(res, pool))
            print(pool, fold, res["info"], flush=True)
    path = d / "composition.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(comp)
    path.write_text(json.dumps(old, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
