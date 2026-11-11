"""Data of the nucfit plan (presets/docs/nucfit/PLAN_nucfit.md, section 2): small pool, large tier, groups, lockbox,
folds, training cores and screening subsets for the nucleotide fast and accurate fits.

    python -m recipe.nucfit_data simulate --ws WS [--jobs 32]     (needs dendropy + INDELible: container rootfs)
    python -m recipe.nucfit_data prepare  --ws WS --refit REFIT_WS [--threads 32]
    python -m recipe.nucfit_data check    --ws WS                 identity checks (no alignment program needed)

Inputs (read-only): REFIT_WS/data/cases_nucleotide.final.tsv + features_nucleotide.final.tsv (the gm1 small pool,
3,225 cases), REFIT_WS/data/raw/rfam/Rfam.seed (Rfam 15.1).

Outputs (WS/data):
  cases_nucfit.tsv            manifest: recipe.scoring columns + tier, kind, family, size_cat, n_seq, ident, bin,
                              group_id, role (fit | lockbox), fold (0-4; -1 for lockbox)
  lockbox_nucfit.txt          lockbox case ids (sha256 in prepare.json; committed at the freeze)
  cores/nucsmall_fold<f>_{core,screen}.txt   f = 0..4, all
  cores/nuclarge_fold<f>_{core,screen}.txt
  large/rfam/<AC>.{fa,ref.fa}, large/sim/<sim_id>/{unaligned.fasta,true.fasta,tree.nwk,control.txt,info.json}
  groups_nucfit.edges.tsv, prepare.json (counts, constants, sha256 of every input and output list)

Rules (deterministic; hash order = sha256("1234|<tag>|<id>")):
  * small pool: the gm1 nucleotide pool unchanged (benchmarks bralibase, mdsa, rfam).
  * large tier, Rfam: every Rfam seed family with more than 200 sequences (the families the small pool's
    n_seq <= 200 filter left out), other filters of the small pool unchanged (ragged, alignment > 2,000 columns,
    empty sequence, gapless); benchmark "rfam_large", full seed size.
  * large tier, simulated DNA (benchmark "sim_dna"): INDELible 1.03, HKY (kappa 2, base frequencies 0.3 0.2 0.2
    0.3), gamma 4 categories alpha 0.5, root length 500, indel model POW 1.7 max 10, birth-death trees (birth 1.0,
    death 0.5, seed edge removed, maximum root-to-tip = depth) as downstream_v2; grid n_taxa {1000, 10000} x depth
    {0.25, 0.5, 1.0} x indel rate {0.02, 0.05, 0.10}; replicates 5 at 1,000 and 3 at 10,000; the LAST replicate
    of every cell is lockbox (one new tree and sequence set per cell, never seen by a fit).
  * groups: recipe.grouping.final_groups (the gm1 nucleotide rule: shared or contained sequence, length-dependent
    identity, Rfam clan) over the small pool + Rfam large families together; each simulation is its own group.
  * lockbox (small + Rfam large): groups with at most CAP (25) cases are eligible. 15 % of the eligible groups
    without a large family, allocated over strata (dominant benchmark x identity bin of the group's first case,
    largest remainder, >= 1 per non-empty stratum), hash order (tag "lockbox") within a stratum; plus 20 % of the
    eligible groups with a large family, allocated over size categories (201-499, 500-999, >= 1000 sequences),
    hash order (tag "lockbox_large"). A lockbox group contributes all its cases.
  * folds: the remaining small + Rfam large cases by recipe.folds.assign_folds (balanced benchmark x bin, whole
    BRAliBASE family groups spread); simulations: within each grid cell, fit replicates in hash order (tag "fold")
    go to folds 0..4 in turn, starting where the previous cell stopped.
  * cores: small = recipe.core rule for nucleotide (cap 25 per group; core = all capped training cases; screen =
    300, equal share per benchmark, proportional over identity bins). Large core = all large training cases;
    large screen = SCREEN_LARGE cases: Rfam large and each simulation size proportional to their training counts
    (largest remainder, >= 2 each), hash order (tag "screen_large<f>").
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import random
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .cases_nucleotide import (allocate, has_gap, parse_stockholm, ref_identity, rfam_clans, rfam_normalise,
                               write_fasta, RFAM_MAX_ALN_LEN, RFAM_MAX_SEQS)
from .scoring import MANIFEST_COLUMNS

SEED = 1234
NFOLDS = 5
CAP = 25
LOCKBOX_SMALL_FRAC = 0.15
LOCKBOX_LARGE_FRAC = 0.20
SCREEN_LARGE = 24
RFAM_SIZE_CATS = ((201, 499), (500, 999), (1000, 10 ** 9))
SIM_GRID = {"n_taxa": (1000, 10000), "depth": (0.25, 0.5, 1.0), "indel_rate": (0.005, 0.01, 0.02)}
SIM_REPS = {1000: 5, 10000: 3}
SIM_LENGTH = 500
SIM_MODEL = {"kappa": 2.0, "freqs": (0.3, 0.2, 0.2, 0.3), "alpha": 0.5, "ncat": 4, "pow": 1.7, "indel_max": 10}
SCORE_SUBSET = 1000          # simulations above this size are scored on a fixed random subset (as extHomFam)
EXTRA = ("tier", "kind", "family", "size_cat", "ident", "bin", "group_id", "role", "fold")


def hkey(tag, x):
    return hashlib.sha256(f"{SEED}|{tag}|{x}".encode()).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sha256_lines(xs):
    return hashlib.sha256("".join(f"{x}\n" for x in xs).encode()).hexdigest()


def ident_bin(x):
    from .folds import difficulty_bin
    return difficulty_bin("nucleotide", x)


# ----------------------------------------------------------------------------------------------- simulations

def sims():
    out = []
    for n, d, ir in itertools.product(SIM_GRID["n_taxa"], SIM_GRID["depth"], SIM_GRID["indel_rate"]):
        for rep in range(SIM_REPS[n]):
            out.append({"sim_id": f"dna_t{n}_d{d}_ir{ir}_r{rep}", "n_taxa": n, "depth": d, "indel_rate": ir, "rep": rep,
                        "cell": f"t{n}_d{d}_ir{ir}", "lockbox": rep == SIM_REPS[n] - 1, "length": SIM_LENGTH})
    return out


def seed_of(sim_id, purpose):
    """32-bit seed per purpose and simulation (the downstream_v2 rule)."""
    return int(hashlib.sha256(f"{purpose}:{sim_id}".encode()).hexdigest()[:8], 16)


def birth_death_tree(n_taxa, depth, seed):
    """As downstream_v2 simulate.birth_death_tree: dendropy birth-death (1.0, 0.5), seed edge removed, maximum
    root-to-tip distance scaled to depth, fixed-point branch lengths, taxa T1..Tn."""
    from dendropy.simulate import treesim

    tree = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, num_extant_tips=n_taxa, rng=random.Random(seed))
    tree.seed_node.edge.length = None
    for i, leaf in enumerate(tree.leaf_node_iter(), 1):
        leaf.taxon.label = f"T{i}"
    for edge in tree.preorder_edge_iter():
        if edge.tail_node is not None and edge.length is None:
            edge.length = 0.0
    max_dist = max(leaf.distance_from_root() for leaf in tree.leaf_node_iter())
    for edge in tree.preorder_edge_iter():
        if edge.tail_node is not None:
            edge.length *= depth / max_dist
    nwk = tree.as_string(schema="newick", suppress_rooting=True).strip()
    nwk = re.sub(r"\[&[^\]]*\]", "", nwk).strip()
    return re.sub(r":([0-9eE.+-]+)", lambda m: f":{float(m.group(1)):.12f}", nwk)


def control_text(sim, nwk):
    m = SIM_MODEL
    return ("[TYPE] NUCLEOTIDE 1\n\n"
            f"[SETTINGS]\n  [randomseed] {seed_of(sim['sim_id'], 'indelible')}\n\n"
            "[MODEL] simmodel\n"
            f"  [submodel] HKY {m['kappa']}\n"
            f"  [statefreq] {' '.join(str(x) for x in m['freqs'])}\n"
            f"  [rates] 0 {m['alpha']} {m['ncat']}\n"
            f"  [indelmodel] POW {m['pow']} {m['indel_max']}\n"
            f"  [indelrate] {sim['indel_rate']}\n\n"
            f"[TREE] t1 {nwk}\n\n[PARTITIONS] p1\n  [t1 simmodel {sim['length']}]\n\n[EVOLVE] p1 1 output\n")


def read_phylip(path):
    names, seqs = [], []
    with open(path) as fh:
        n, width = (int(x) for x in fh.readline().split())
        for line in fh:
            parts = line.split(None, 1)
            if len(parts) == 2:
                names.append(parts[0])
                seqs.append(parts[1].replace(" ", "").strip())
    if len(names) != n or any(len(s) != width for s in seqs):
        raise ValueError(f"{path}: header says {n} x {width}, read {len(names)} rows")
    return names, seqs


def simulate_one(args):
    sim, root, indelible = args
    out = Path(root) / sim["sim_id"]
    if (out / "info.json").exists():
        return json.loads((out / "info.json").read_text())
    out.mkdir(parents=True, exist_ok=True)
    nwk = birth_death_tree(sim["n_taxa"], sim["depth"], seed_of(sim["sim_id"], "tree"))
    (out / "tree.nwk").write_text(nwk + "\n")
    (out / "control.txt").write_text(control_text(sim, nwk))
    p = subprocess.run([indelible, "control.txt"], cwd=out, capture_output=True, text=True, timeout=7200)
    if p.returncode != 0:
        raise RuntimeError(f"INDELible failed for {sim['sim_id']}: {p.stderr[-500:]}")
    true_path = next(x for x in (out / "output_TRUE_1.phy", out / "output_TRUE.phy") if x.exists())
    names, rows = read_phylip(true_path)
    if any(not s.replace("-", "") for s in rows):
        raise ValueError(f"{sim['sim_id']}: empty sequence")
    write_fasta(out / "true.fasta", list(zip(names, rows)))
    write_fasta(out / "unaligned.fasta", [(n, s.replace("-", "")) for n, s in zip(names, rows)])
    for f in out.glob("output*"):
        f.unlink()
    for f in ("trees.txt", "LOG.txt"):
        (out / f).unlink(missing_ok=True)
    lens = [len(s.replace("-", "")) for s in rows]
    info = {**sim, "true_width": len(rows[0]), "n_seq": len(rows), "len_min": min(lens), "len_max": max(lens),
            "gap_fraction": sum(s.count("-") for s in rows) / (len(rows) * len(rows[0])),
            "sha256_true": sha256_file(out / "true.fasta"), "sha256_unaligned": sha256_file(out / "unaligned.fasta")}
    (out / "info.json").write_text(json.dumps(info, indent=1))
    return info


def cmd_simulate(ws, jobs, indelible):
    root = ws / "data" / "large" / "sim"
    root.mkdir(parents=True, exist_ok=True)
    todo = sorted(sims(), key=lambda s: -s["n_taxa"])
    with ProcessPoolExecutor(jobs) as ex:
        infos = list(ex.map(simulate_one, [(s, str(root), indelible) for s in todo]))
    print(json.dumps({"n": len(infos), "max_width": max(i["true_width"] for i in infos),
                      "empty": sum(i["len_min"] == 0 for i in infos)}, indent=1))


# ---------------------------------------------------------------------------------------------- Rfam large

def rfam_large(seed_path, out_dir):
    """-> (rows, drops): Rfam seed families with > RFAM_MAX_SEQS sequences, small-pool filters otherwise."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rows, drops = [], Counter()
    for fam in parse_stockholm(seed_path):
        n = len(fam["rows"])
        if n <= RFAM_MAX_SEQS:
            continue
        lens = {len(s) for _, s in fam["rows"]}
        if len(lens) != 1:
            drops["ragged"] += 1
            continue
        if lens.pop() > RFAM_MAX_ALN_LEN:
            drops[f"aln_len>{RFAM_MAX_ALN_LEN}"] += 1
            continue
        refs = [(nm, rfam_normalise(s)) for nm, s in fam["rows"]]
        if any(not s.replace("-", "") for _, s in refs):
            drops["empty_sequence"] += 1
            continue
        if not has_gap([s for _, s in refs]):
            drops["gapless"] += 1
            continue
        ac = fam["ac"]
        ref_path, unal_path = out_dir / f"{ac}.ref.fa", out_dir / f"{ac}.fa"
        write_fasta(ref_path, refs)
        write_fasta(unal_path, [(nm, s.replace("-", "")) for nm, s in refs])
        ident = ref_identity([s for _, s in refs])
        cat = next(f"{lo}-{hi}" if hi < 10 ** 9 else f">={lo}" for lo, hi in RFAM_SIZE_CATS if lo <= n <= hi)
        rows.append({"case_id": f"rfam_large/{ac}", "pool": "nucleotide", "benchmark": "rfam_large", "subset": "seed",
                     "unaligned": str(unal_path), "reference": str(ref_path), "score_rule": "all_columns",
                     "n_seq": n, "ref_len": len(refs[0][1]), "sha256_unaligned": sha256_file(unal_path),
                     "sha256_reference": sha256_file(ref_path), "tier": "large", "kind": "rna", "family": ac,
                     "size_cat": cat, "ident": round(ident, 6)})
    return rows, drops


# ------------------------------------------------------------------------------------------------- prepare

def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def small_rows(refit):
    rows = read_tsv(refit / "data" / "cases_nucleotide.final.tsv")
    feats = {r["case_id"]: r for r in read_tsv(refit / "data" / "features_nucleotide.final.tsv")}
    out = []
    for r in rows:
        try:
            pid = float(feats[r["case_id"]]["mean_pid_core"])
        except ValueError:
            pid = float("nan")
        out.append({**{k: r[k] for k in MANIFEST_COLUMNS}, "tier": "small",
                    "kind": "dna" if r["seq_kind"] == "dna" else r["seq_kind"], "family": r["family"],
                    "size_cat": "small", "ident": pid})
    return out


def score_subset(d, sim_id, n_seq):
    """The scored reference of a simulation: true.fasta itself up to SCORE_SUBSET sequences; above, the true alignment
    of SCORE_SUBSET sequences drawn with a seed from the sim id (columns that are gaps in every drawn row removed),
    written once to ref_subset.fasta. The test alignment is projected onto these rows (score_rule "projected"), as the
    extHomFam reference rows inside a large family. -> (path, rule, sha256, width)."""
    if n_seq <= SCORE_SUBSET:
        return d / "true.fasta", "all_columns", None, None
    path = d / "ref_subset.fasta"
    recs = parse_fasta_simple(d / "true.fasta")
    pick = sorted(random.Random(seed_of(sim_id, "score_subset")).sample(range(len(recs)), SCORE_SUBSET))
    rows = [recs[i] for i in pick]
    keep = [j for j in range(len(rows[0][1])) if any(s[j] != "-" for _, s in rows)]
    sub = [(n, "".join(s[j] for j in keep)) for n, s in rows]
    if not path.exists():
        write_fasta(path, sub)
    elif parse_fasta_simple(path) != sub:
        raise ValueError(f"{path}: differs from the seeded subset")
    return path, "projected", sha256_file(path), len(keep)


def parse_fasta_simple(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def sim_rows(ws):
    out = []
    for s in sims():
        d = ws / "data" / "large" / "sim" / s["sim_id"]
        info = json.loads((d / "info.json").read_text())
        ref, rule, sha, w = score_subset(d, s["sim_id"], info["n_seq"])
        out.append({"case_id": f"sim_dna/{s['sim_id']}", "pool": "nucleotide", "benchmark": "sim_dna",
                    "subset": s["cell"], "unaligned": str(d / "unaligned.fasta"), "reference": str(ref),
                    "score_rule": rule, "n_seq": info["n_seq"], "ref_len": w or info["true_width"],
                    "sha256_unaligned": info["sha256_unaligned"], "sha256_reference": sha or info["sha256_true"],
                    "tier": "large", "kind": "dna", "family": s["sim_id"], "size_cat": str(s["n_taxa"]),
                    "ident": float("nan"), "_sim": s})
    return out


def draw_lockbox(groups_rows):
    """groups_rows: {group_id: [rows]} of small + Rfam large -> set of lockbox group ids, info."""
    elig = {g: ms for g, ms in groups_rows.items() if len(ms) <= CAP}
    has_large = {g for g, ms in elig.items() if any(m["benchmark"] == "rfam_large" for m in ms)}
    strata_s, strata_l = defaultdict(list), defaultdict(list)
    for g, ms in elig.items():
        if g in has_large:
            big = max(int(m["n_seq"]) for m in ms if m["benchmark"] == "rfam_large")
            strata_l[next(i for i, (lo, hi) in enumerate(RFAM_SIZE_CATS) if lo <= big <= hi)].append(g)
        else:
            dom = Counter(m["benchmark"] for m in ms).most_common(1)[0][0]
            first = min(ms, key=lambda m: m["case_id"])
            strata_s[(dom, first["bin"])].append(g)
    out, info = set(), {}
    for tag, strata, frac in (("lockbox", strata_s, LOCKBOX_SMALL_FRAC),
                              ("lockbox_large", strata_l, LOCKBOX_LARGE_FRAC)):
        total = sum(len(v) for v in strata.values())
        alloc = allocate({k: len(v) for k, v in strata.items()}, round(frac * total), 1)
        for k, v in strata.items():
            out |= set(sorted(v, key=lambda g: hkey(tag, g))[:alloc[k]])
        info[tag] = {"eligible_groups": total, "drawn": sum(alloc.values()),
                     "alloc": {str(k): alloc[k] for k in sorted(strata, key=str)}}
    info["ineligible_groups_over_cap"] = len(groups_rows) - len(elig)
    return out, info


def sim_folds(rows):
    by_cell = defaultdict(list)
    for r in rows:
        by_cell[r["_sim"]["cell"]].append(r)
    i, out = 0, {}
    for cell in sorted(by_cell):
        for r in sorted(by_cell[cell], key=lambda r: hkey("fold", r["case_id"])):
            out[r["case_id"]] = i % NFOLDS
            i += 1
    return out


def large_screen(train, fold):
    by = defaultdict(list)
    for r in train:
        by[r["size_cat"] if r["benchmark"] == "sim_dna" else "rfam_large"].append(r)
    alloc = allocate({k: len(v) for k, v in by.items()}, SCREEN_LARGE, 2)
    out = []
    for k, v in by.items():
        out += sorted(v, key=lambda r: hkey(f"screen_large{fold}", r["case_id"]))[:alloc[k]]
    return sorted(r["case_id"] for r in out)


def cmd_prepare(ws, refit, threads):
    from . import core as CORE
    from . import grouping
    from .folds import assign_folds

    data = ws / "data"
    seed_path = refit / "data" / "raw" / "rfam" / "Rfam.seed"
    small = small_rows(refit)
    rlarge, drops = rfam_large(seed_path, data / "large" / "rfam")
    srows = sim_rows(ws)
    for r in small + rlarge:
        r["bin"] = ident_bin(r["ident"])
    for r in srows:
        r["bin"] = 0
    # groups over small + Rfam large (the gm1 rule, with clans)
    grp_in = small + rlarge
    clans = rfam_clans(seed_path)
    # clan links apply to every Rfam case (grouping.clan_pairs looks at benchmark "rfam" only)
    view = [{**r, "benchmark": "rfam" if r["benchmark"] == "rfam_large" else r["benchmark"]} for r in grp_in]
    if grouping.HAVE_NUMBA:
        grouping.numba.set_num_threads(threads)
    roots, edges, gstats = grouping.final_groups(view, "nucleotide", threads=threads, clans=clans)
    ids = [r["case_id"] for r in grp_in]
    gids = grouping.relabel(roots, ids)
    edges = [(ids[x], ids[y], kind, idn) for x, y, kind, idn in edges]
    for r, g in zip(grp_in, gids):
        r["group_id"] = g
    for r in srows:
        r["group_id"] = "S:" + r["case_id"].split("/", 1)[1]
    members = defaultdict(list)
    for r in grp_in:
        members[r["group_id"]].append(r)
    lock_groups, linfo = draw_lockbox(members)
    for r in grp_in:
        r["role"] = "lockbox" if r["group_id"] in lock_groups else "fit"
    for r in srows:
        r["role"] = "lockbox" if r["_sim"]["lockbox"] else "fit"
    # folds
    fit_gs = [dict(case_id=r["case_id"], benchmark=r["benchmark"], bin=r["bin"], group_id=r["group_id"])
              for r in grp_in if r["role"] == "fit"]
    gfold = assign_folds(fit_gs, "nucleotide", SEED)
    sfold = sim_folds([r for r in srows if r["role"] == "fit"])
    allrows = grp_in + srows
    for r in allrows:
        r["fold"] = -1 if r["role"] == "lockbox" else (sfold[r["case_id"]] if r["benchmark"] == "sim_dna"
                                                        else gfold[r["group_id"]])
    # manifest
    cols = list(MANIFEST_COLUMNS) + list(EXTRA)
    with open(data / "cases_nucfit.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(allrows, key=lambda r: r["case_id"]):
            w.writerow({**r, "ident": "" if r["ident"] != r["ident"] else r["ident"]})
    lock_ids = sorted(r["case_id"] for r in allrows if r["role"] == "lockbox")
    (data / "lockbox_nucfit.txt").write_text("".join(f"{x}\n" for x in lock_ids))
    with open(data / "groups_nucfit.edges.tsv", "w") as fh:
        for e in edges:
            fh.write("\t".join(str(x) for x in e) + "\n")
    # cores
    (data / "cores").mkdir(exist_ok=True)
    small_fit = [dict(case_id=r["case_id"], benchmark=r["benchmark"], bin=r["bin"], group_id=r["group_id"])
                 for r in small if r["role"] == "fit"]
    small_folds = {r["case_id"]: {"fold": r["fold"]} for r in small if r["role"] == "fit"}
    large_fit = [r for r in rlarge + srows if r["role"] == "fit"]
    cores = {}
    for f in CORE.FOLDS:
        c = CORE.build_core("nucleotide", f, ws=ws, data=small_fit, folds=small_folds)
        train_l = [r for r in large_fit if f == "all" or r["fold"] != f]
        lscreen = large_screen(train_l, f)
        lists = {"nucsmall_core": c["core"], "nucsmall_screen": c["screen"],
                 "nuclarge_core": sorted(r["case_id"] for r in train_l), "nuclarge_screen": lscreen}
        for name, ids in lists.items():
            pool, kind = name.split("_")
            (data / "cores" / f"{pool}_fold{f}_{kind}.txt").write_text("".join(f"{x}\n" for x in ids))
        cores[str(f)] = {k: {"n": len(v), "sha256": sha256_lines(v)} for k, v in lists.items()}
    stats = {
        "constants": {"seed": SEED, "cap": CAP, "lockbox_small_frac": LOCKBOX_SMALL_FRAC,
                      "lockbox_large_frac": LOCKBOX_LARGE_FRAC, "screen_large": SCREEN_LARGE,
                      "sim_grid": SIM_GRID, "sim_reps": {str(k): v for k, v in SIM_REPS.items()},
                      "sim_length": SIM_LENGTH, "sim_model": SIM_MODEL},
        "inputs": {"cases_nucleotide.final.tsv": sha256_file(refit / "data" / "cases_nucleotide.final.tsv"),
                   "features_nucleotide.final.tsv": sha256_file(refit / "data" / "features_nucleotide.final.tsv"),
                   "Rfam.seed": sha256_file(seed_path)},
        "rfam_large": {"kept": len(rlarge), "dropped": dict(drops),
                       "size_cats": dict(Counter(r["size_cat"] for r in rlarge))},
        "groups": {"n_groups": len(set(gids)), "clan_links": gstats.get("clan_links"),
                   "clan_merges": gstats.get("clan_merges"), "links": gstats.get("links")},
        "lockbox": {**linfo, "n_cases": len(lock_ids), "sha256": sha256_lines(lock_ids),
                    "by_benchmark": dict(Counter(r["benchmark"] for r in allrows if r["role"] == "lockbox"))},
        "fold_x_benchmark": dict(sorted(Counter(f"{r['fold']}|{r['benchmark']}" for r in allrows).items())),
        "largest_group_share": max(Counter(r["group_id"] for r in grp_in if r["role"] == "fit").values())
        / sum(r["role"] == "fit" for r in grp_in),
        "cores": cores,
        "manifest_sha256": sha256_file(data / "cases_nucfit.tsv"),
    }
    (data / "prepare.json").write_text(json.dumps(stats, indent=1, default=str))
    print(json.dumps({k: stats[k] for k in ("rfam_large", "lockbox", "fold_x_benchmark", "largest_group_share")},
                     indent=1, default=str))


# --------------------------------------------------------------------------------------------------- check

def cmd_check(ws):
    """Identity checks: no lockbox case in any core; no group across folds or into the lockbox; screens inside
    cores; every training core excludes its test fold; files hash to the manifest; sims parse as alignments whose
    ungapped rows equal the unaligned input."""
    from .cases_nucleotide import read_fasta

    data = ws / "data"
    rows = read_tsv(data / "cases_nucfit.tsv")
    byid = {r["case_id"]: r for r in rows}
    lock = set((data / "lockbox_nucfit.txt").read_text().split())
    errs = []
    gf = defaultdict(set)
    for r in rows:
        gf[r["group_id"]].add(r["fold"])
    errs += [f"group {g} in folds {sorted(v)}" for g, v in gf.items() if len(v) > 1]
    for f in (0, 1, 2, 3, 4, "all"):
        for pool in ("nucsmall", "nuclarge"):
            core = set((data / "cores" / f"{pool}_fold{f}_core.txt").read_text().split())
            screen = set((data / "cores" / f"{pool}_fold{f}_screen.txt").read_text().split())
            if core & lock:
                errs.append(f"{pool} fold {f}: {len(core & lock)} lockbox cases in core")
            if not screen <= core:
                errs.append(f"{pool} fold {f}: screen not inside core")
            if f != "all" and any(byid[c]["fold"] == str(f) for c in core):
                errs.append(f"{pool} fold {f}: test-fold case in core")
    bad = 0
    for r in rows:
        if sha256_file(r["unaligned"]) != r["sha256_unaligned"] or sha256_file(r["reference"]) != r["sha256_reference"]:
            bad += 1
        if r["tier"] == "large":
            ref = read_fasta(r["reference"])
            unal = dict(read_fasta(r["unaligned"]))
            if len({len(s) for _, s in ref}) != 1 or any(unal.get(n) != s.replace("-", "") for n, s in ref):
                errs.append(f"{r['case_id']}: reference rows do not match the input")
    if bad:
        errs.append(f"{bad} files with a sha256 mismatch")
    print(json.dumps({"cases": len(rows), "lockbox": len(lock), "errors": errs}, indent=1))
    return 1 if errs else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["simulate", "prepare", "check"])
    ap.add_argument("--ws", required=True)
    ap.add_argument("--refit", default="~/tmp/kalign-refit")
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--indelible", default="indelible")
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    if a.cmd == "simulate":
        cmd_simulate(ws, a.jobs, a.indelible)
    elif a.cmd == "prepare":
        cmd_prepare(ws, Path(a.refit).expanduser(), a.threads)
    else:
        sys.exit(cmd_check(ws))


if __name__ == "__main__":
    main()
