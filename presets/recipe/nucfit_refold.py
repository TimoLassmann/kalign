"""PLAN_nucfit.md amendment A4 (2026-11-04): folds without an oversized group.

Why: the grouping rule (shared sequence, length-dependent identity, Rfam clan; connected components) chained 1,161 of
the 2,988 fit cases (small + Rfam large) into one group, G0001: all 574 BRAliBASE cases, 474 Rfam seed cases, 72 of
the 92 Rfam large families and 41 MDSA cases. A group is never split, so all of it went to fold 0: fold 0's fit saw
no BRAliBASE and was tested on all of it, folds 1-4 were tested on none, and the cap of 25 cases per group left 15
BRAliBASE cases in the training cores while the objectives weight BRAliBASE as one third of Q_small.

Rule (deterministic):
  1. links: for every fit group over a limit (below), ALL case-level links inside it are computed (the grouping keeps
     only the links that merged two components): cases sharing a sequence; sequence pairs passing the grouping's
     length-dependent identity rule (input + reference sequences; k-mer prefilter, every candidate verified, none
     skipped); every pair of Rfam cases (seed or large) whose families are in one Rfam clan.
  2. split: while a component of such a group holds more than MAX_SHARE of the fit cases or more than BENCH_SHARE of
     the fit cases of one benchmark, a case of the first violating component (largest first) is removed. Candidates,
     in this order of benchmark tiers: Rfam seed + MDSA; Rfam large; BRAliBASE (a tier is used only when the earlier
     tiers have no case in the component); inside a tier the case with the most links inside the component, ties by
     case id. Rfam seed and MDSA are plentiful; Rfam large families are the scarce real large tier; BRAliBASE is the
     benchmark the split exists for.
  3. removed cases get role "excluded": in no core, no test fold, no lockbox (they would link folds).
     The pieces become groups <old id>.<n> (n by size, then first case id); the lockbox, every other group and the
     simulations are unchanged.
  4. folds: recipe.folds.assign_folds over the fit groups (as before: balanced benchmark x identity bin; groups with
     >= 20 BRAliBASE cases one per fold); simulations keep their folds. NFOLDS stays 5 (DECISIONS.md D14).
  5. cores and screens: the recipe.nucfit_data rules, unchanged.

    python -m recipe.nucfit_refold links  --ws WS --refit REFIT [--threads 60]    (numba: venv-grouping)
    python -m recipe.nucfit_refold refold --ws WS
    python -m recipe.nucfit_refold check  --ws WS --refit REFIT [--threads 60]    (numba)

links -> WS/data/refold/links_<group>.tsv; refold archives data/{cases_nucfit.tsv,cores,prepare.json} to
data/archive_folds_v1/ (once) and writes the new manifest, cores and data/refold/refold.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

MAX_SHARE = 0.10
BENCH_SHARE = 0.35
TIERS = (("rfam", "mdsa"), ("rfam_large",), ("bralibase",))
EXCLUDED = "excluded"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fit_rows(rows):
    """Fit cases of the grouped population (small + Rfam large; simulations are their own groups)."""
    return [r for r in rows if r["role"] == "fit" and r["benchmark"] != "sim_dna"]


def limits(rows):
    fit = fit_rows(rows)
    tot = Counter(r["benchmark"] for r in fit)
    return int(MAX_SHARE * len(fit)), {b: max(1, int(BENCH_SHARE * n)) for b, n in tot.items()}


def violates(members, bench, lim, blim):
    if len(members) > lim:
        return True
    c = Counter(bench[x] for x in members)
    return any(v > blim[b] for b, v in c.items())


def oversized_groups(rows):
    lim, blim = limits(rows)
    g = defaultdict(list)
    for r in fit_rows(rows):
        g[r["group_id"]].append(r["case_id"])
    bench = {r["case_id"]: r["benchmark"] for r in rows}
    return sorted(k for k, v in g.items() if violates(v, bench, lim, blim))


# ----------------------------------------------------------------------------------------------- 1. links
def full_links(rows, clans, threads=16, log=print):
    """All case-level links among rows (manifest rows with unaligned/reference paths). -> sorted [(a, b, kinds)]."""
    import numpy as np
    from . import grouping as G
    p = G.PARAMS["nucleotide"]
    ids = [r["case_id"] for r in rows]
    seqs, ucases = G.load_sequences(rows, "nucleotide", "both")
    nuid = len(seqs)
    kinds = defaultdict(set)
    for cs in ucases:
        cs = sorted(set(cs))
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                kinds[(cs[i], cs[j])].add("exact")
    usable = np.array([len(s) >= p["min_len"] for s in seqs])
    enc = [G.encode(s, "nucleotide") for s in seqs]
    lens = np.array([len(e) for e in enc], dtype=np.int64)
    cl = [e if usable[i] else e[:0] for i, e in enumerate(enc)]
    offs = np.zeros(nuid + 1, dtype=np.int64)
    offs[1:] = np.cumsum([len(e) for e in cl])
    index = G.build_index(np.concatenate(cl).astype(np.uint8), offs, p["k"], len(G.ALPHABETS["nucleotide"]))
    label = np.arange(nuid, dtype=np.int64)        # every sequence its own label: no candidate pair is skipped
    I, J, _ = G.candidates(index, label, p, nuid, threads)
    log(f"{len(rows)} cases, {nuid} unique sequences, {len(I)} candidate pairs")
    ac, ao = G.allcodes(enc), G.allofs(enc)
    step = 500000
    for st in range(0, len(I), step):
        ii, jj = I[st:st + step], J[st:st + step]
        idn = G.verify_pairs(ac, ao, ii, jj, p["match"], p["mismatch"], p["gap"])
        ok = idn >= G.required_identity("nucleotide", np.minimum(lens[ii], lens[jj]))
        for u, v in zip(ii[ok].tolist(), jj[ok].tolist()):
            for x in ucases[u]:
                for y in ucases[v]:
                    if x != y:
                        kinds[(min(x, y), max(x, y))].add("identity")
        if (st // step) % 20 == 0:
            log(f"  verified {min(st + step, len(I))} / {len(I)}")
    byclan = defaultdict(list)
    for i, r in enumerate(rows):
        if r["benchmark"] in ("rfam", "rfam_large") and r["family"] in clans:
            byclan[clans[r["family"]]].append(i)
    for m in byclan.values():
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                kinds[(min(m[i], m[j]), max(m[i], m[j]))].add("clan")
    return sorted((ids[x], ids[y], ",".join(sorted(k))) for (x, y), k in kinds.items())


# ----------------------------------------------------------------------------------------------- 2. split
def components(nodes, adj):
    seen, out = set(), []
    for s in sorted(nodes):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x]:
                if y in nodes and y not in seen:
                    seen.add(y)
                    stack.append(y)
        out.append(sorted(comp))
    return sorted(out, key=lambda c: (-len(c), c[0]))


def split_group(members, links, bench, lim, blim):
    """members: case ids of one group; links: [(a, b, kinds)]; bench: {case: benchmark}.
    -> (removed case ids in removal order, components of the rest)."""
    adj = defaultdict(set)
    for a, b, _ in links:
        adj[a].add(b)
        adj[b].add(a)
    alive, removed = set(members), []
    while True:
        comps = components(alive, adj)
        bad = [c for c in comps if violates(c, bench, lim, blim)]
        if not bad:
            return removed, comps
        cset = set(bad[0])
        order = sorted(bad[0], key=lambda x: (-len(adj[x] & cset), x))
        pick = next((x for tier in TIERS for x in order if bench[x] in tier), order[0])
        alive.discard(pick)
        removed.append(pick)


# ----------------------------------------------------------------------------------------------- 3-5. refold
def cmd_links(ws, refit, threads):
    from . import grouping as G
    from .cases_nucleotide import rfam_clans
    from .scoring import read_manifest
    if G.HAVE_NUMBA:
        G.numba.set_num_threads(threads)
    data = ws / "data"
    rows = read_tsv(data / "cases_nucfit.tsv")
    if any(r["role"] == EXCLUDED for r in rows):
        raise SystemExit("the manifest is already refolded; links are computed on the original groups")
    clans = rfam_clans(refit / "data" / "raw" / "rfam" / "Rfam.seed")
    ext = {r["case_id"]: r for r in rows}
    man = read_manifest(data / "cases_nucfit.tsv")
    (data / "refold").mkdir(exist_ok=True)
    for g in oversized_groups(rows):
        sub = [{**m, "benchmark": ext[m["case_id"]]["benchmark"], "family": ext[m["case_id"]]["family"]}
               for m in man if ext[m["case_id"]]["group_id"] == g and ext[m["case_id"]]["role"] == "fit"]
        links = full_links(sub, clans, threads)
        with open(data / "refold" / f"links_{g}.tsv", "w") as fh:
            for a, b, k in links:
                fh.write(f"{a}\t{b}\t{k}\n")
        print(f"{g}: {len(sub)} cases, {len(links)} links")
    return 0


def cmd_refold(ws):
    from . import core as CORE
    from .folds import assign_folds
    from .nucfit_data import EXTRA, NFOLDS, SEED, large_screen, sha256_file, sha256_lines
    from .scoring import MANIFEST_COLUMNS
    data = ws / "data"
    arch = data / "archive_folds_v1"
    if not arch.exists():
        arch.mkdir()
        for name in ("cases_nucfit.tsv", "prepare.json", "lockbox_nucfit.txt"):
            shutil.copy2(data / name, arch / name)
        shutil.copytree(data / "cores", arch / "cores")
    rows = read_tsv(arch / "cases_nucfit.tsv")          # always from the original manifest
    cols = list(rows[0].keys())
    assert cols == list(MANIFEST_COLUMNS) + list(EXTRA), "manifest columns changed"
    bench = {r["case_id"]: r["benchmark"] for r in rows}
    lim, blim = limits(rows)
    info = {"rule": "PLAN_nucfit.md amendment A4", "max_share": MAX_SHARE, "bench_share": BENCH_SHARE,
            "tiers": TIERS, "limit_cases": lim, "limit_per_benchmark": blim, "groups": {}}
    removed_all = set()
    newgid = {}
    for g in oversized_groups(rows):
        members = sorted(r["case_id"] for r in fit_rows(rows) if r["group_id"] == g)
        path = data / "refold" / f"links_{g}.tsv"
        links = [tuple(l.rstrip("\n").split("\t")) for l in open(path)]
        assert all(a in set(members) and b in set(members) for a, b, _ in links), f"{path}: foreign case"
        removed, comps = split_group(members, links, bench, lim, blim)
        removed_all |= set(removed)
        for n, comp in enumerate(comps, 1):
            for c in comp:
                newgid[c] = f"{g}.{n}"
        info["groups"][g] = {"cases": len(members), "links": len(links), "links_sha256": sha256_file(path),
                             "removed": len(removed), "removed_by_benchmark": dict(Counter(bench[c] for c in removed)),
                             "removed_ids": removed, "pieces": [len(c) for c in comps],
                             "pieces_by_benchmark": [dict(Counter(bench[x] for x in c)) for c in comps[:12]]}
    for r in rows:
        if r["case_id"] in removed_all:
            r["role"], r["fold"], r["group_id"] = EXCLUDED, "-1", r["group_id"]
        elif r["case_id"] in newgid:
            r["group_id"] = newgid[r["case_id"]]
    # folds: grouped fit cases re-assigned; simulations and lockbox unchanged
    grouped = fit_rows(rows)
    gfold = assign_folds([dict(case_id=r["case_id"], benchmark=r["benchmark"], bin=int(r["bin"]),
                               group_id=r["group_id"]) for r in grouped], "nucleotide", SEED)
    for r in grouped:
        r["fold"] = str(gfold[r["group_id"]])
    with open(data / "cases_nucfit.tsv.tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["case_id"]):
            w.writerow(r)
    (data / "cases_nucfit.tsv.tmp").replace(data / "cases_nucfit.tsv")
    lock = sorted(r["case_id"] for r in rows if r["role"] == "lockbox")
    old_lock = (arch / "lockbox_nucfit.txt").read_text().split()
    assert lock == sorted(old_lock), "lockbox changed"
    # cores (the recipe.nucfit_data rules)
    small_fit = [dict(case_id=r["case_id"], benchmark=r["benchmark"], bin=int(r["bin"]), group_id=r["group_id"])
                 for r in rows if r["role"] == "fit" and r["tier"] == "small"]
    small_folds = {r["case_id"]: {"fold": r["fold"]} for r in rows if r["role"] == "fit" and r["tier"] == "small"}
    large_fit = [dict(r, fold=int(r["fold"])) for r in rows if r["role"] == "fit" and r["tier"] == "large"]
    cores = {}
    for f in CORE.FOLDS:
        c = CORE.build_core("nucleotide", f, ws=ws, data=small_fit, folds=small_folds)
        train_l = [r for r in large_fit if f == "all" or r["fold"] != f]
        lists = {"nucsmall_core": c["core"], "nucsmall_screen": c["screen"],
                 "nuclarge_core": sorted(r["case_id"] for r in train_l), "nuclarge_screen": large_screen(train_l, f)}
        for name, ids in lists.items():
            pool, kind = name.split("_")
            (data / "cores" / f"{pool}_fold{f}_{kind}.txt").write_text("".join(f"{x}\n" for x in ids))
        cores[str(f)] = {k: {"n": len(v), "sha256": sha256_lines(v),
                             "by_benchmark": dict(Counter(bench[x] for x in v))} for k, v in lists.items()}
    fits = [r for r in rows if r["role"] == "fit"]
    info.update({
        "nfolds": NFOLDS, "excluded": len(removed_all),
        "excluded_by_benchmark": dict(Counter(bench[c] for c in removed_all)),
        "lockbox_sha256": sha256_lines(lock),
        "fold_x_benchmark": dict(sorted(Counter(f"{r['fold']}|{r['benchmark']}" for r in fits).items())),
        "largest_group_share": max(Counter(r["group_id"] for r in fit_rows(rows)).values()) / len(fit_rows(rows)),
        "cores": cores, "manifest_sha256": sha256_file(data / "cases_nucfit.tsv")})
    (data / "refold" / "refold.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("excluded", "excluded_by_benchmark", "fold_x_benchmark",
                                           "largest_group_share")}, indent=1))
    for f, c in cores.items():
        print(f"fold {f}: small core {c['nucsmall_core']['n']} {c['nucsmall_core']['by_benchmark']}, "
              f"large core {c['nuclarge_core']['n']}")
    return 0


def cmd_check(ws, refit, threads):
    """Identity checks of the refolded data: the recipe.nucfit_data checks, plus: no excluded case in any core;
    lockbox unchanged; every group (new ids) in one fold; no link (exact, identity, clan) between cases of
    different folds or between a fit case and a lockbox case (exhaustive: grouping.homolog_flags over all fit and
    lockbox cases, clans checked directly); BRAliBASE in every test fold."""
    from . import grouping as G
    from .cases_nucleotide import rfam_clans
    from .nucfit_data import cmd_check as data_check
    from .scoring import read_manifest
    if G.HAVE_NUMBA:
        G.numba.set_num_threads(threads)
    data = ws / "data"
    errs = []
    if data_check(ws) != 0:
        errs.append("recipe.nucfit_data check failed (see above)")
    rows = read_tsv(data / "cases_nucfit.tsv")
    excl = {r["case_id"] for r in rows if r["role"] == EXCLUDED}
    for p in sorted((data / "cores").glob("*.txt")):
        bad = set(p.read_text().split()) & excl
        if bad:
            errs.append(f"{p.name}: {len(bad)} excluded cases")
    old = (data / "archive_folds_v1" / "lockbox_nucfit.txt").read_text().split()
    if sorted(old) != sorted((data / "lockbox_nucfit.txt").read_text().split()):
        errs.append("lockbox file changed")
    for f in range(5):
        if not any(r["fold"] == str(f) and r["benchmark"] == "bralibase" for r in rows):
            errs.append(f"fold {f}: no BRAliBASE case")
    # exhaustive cross-fold links: fit folds 0-4, lockbox = fold 5; excluded cases left out
    ext = {r["case_id"]: r for r in rows}
    man = [m for m in read_manifest(data / "cases_nucfit.tsv")
           if ext[m["case_id"]]["role"] != EXCLUDED and ext[m["case_id"]]["benchmark"] != "sim_dna"]
    fold_of = [5 if ext[m["case_id"]]["role"] == "lockbox" else int(ext[m["case_id"]]["fold"]) for m in man]
    flags, finfo = G.homolog_flags(man, "nucleotide", fold_of, threads=threads)
    if flags.any():
        bad = [m["case_id"] for m, f in zip(man, flags.tolist()) if f]
        errs.append(f"{len(bad)} cases linked to another fold or the lockbox, e.g. {bad[:5]}")
    clans = rfam_clans(refit / "data" / "raw" / "rfam" / "Rfam.seed")
    cf = defaultdict(set)
    for m, f in zip(man, fold_of):
        r = ext[m["case_id"]]
        if r["benchmark"] in ("rfam", "rfam_large") and r["family"] in clans:
            cf[clans[r["family"]]].add(f)
    errs += [f"clan {c} in folds {sorted(v)}" for c, v in cf.items() if len(v) > 1]
    out = {"excluded": len(excl), "homolog_search": finfo, "errors": errs}
    (data / "refold" / "check.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))
    return 1 if errs else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["links", "refold", "check"])
    ap.add_argument("--ws", required=True)
    ap.add_argument("--refit", default="~/tmp/kalign-refit")
    ap.add_argument("--threads", type=int, default=32)
    a = ap.parse_args(argv)
    ws, refit = Path(a.ws).expanduser(), Path(a.refit).expanduser()
    if a.cmd == "links":
        return cmd_links(ws, refit, a.threads)
    if a.cmd == "refold":
        return cmd_refold(ws)
    return cmd_check(ws, refit, a.threads)


if __name__ == "__main__":
    sys.exit(main())
