"""extHomFam data of the fast-fit plan (PLAN_fastfit.md section 3): lockbox draw, subsets, folds, screening lists.

    python -m recipe.largefam prepare --ws WS --ehf EHF_DIR --xl SCALE_XL_DIR [--jobs 16]
    python -m recipe.largefam check   --ws WS                       identity I3 (needs the kalign wheel)

-> WS/data/cases_exthomfam.tsv               manifest (recipe.scoring columns), pool "exthomfam"
   WS/data/folds_exthomfam.tsv               case_id, family, k, size_cat, n_family, n_ref, role, fold
   WS/data/cores/exthomfam_fold<f>_core.txt  training case ids of fold f (0..4, all)
   WS/data/cores/exthomfam_fold<f>_screen.txt  the screening cases of fold f
   WS/data/large/in/<family>__k<k>.fa, WS/data/large/ref/<family>.ref.fa
   WS/data/large/prepare.json                counts, rule constants, sha256 of the metadata file

Rules (all deterministic; hash order = sha256("1234|<tag>|<family>")):
  * eligible: families with >= 1,000 sequences. The four 200,000-sequence inputs of the paper benchmark (cox,
    helicase_NC, oat, TyrKc) are lockbox only, as their existing files (SCALE_XL_DIR).
  * lockbox: 60 families drawn from the eligible families that no exploration used (EXPLORED) and that are not one
    of the four 200k families; allocated to the extHomFam size categories proportionally (largest remainder, at
    least 1 per category), hash order (tag "lockbox") within a category. Each is used at 1,000 and, with >= 10,000
    sequences, at 10,000; 10 of those with >= 50,000 (hash order, tag "lockbox50k") also at 50,000.
  * fit families: every other eligible family (EXPLORED included), at 1,000 and, with >= 10,000 sequences, at
    10,000. Folds: within each size category, families in hash order (tag "fold") go to folds 0..4 in turn,
    starting where the previous category stopped.
  * subsets: all reference sequences + random.Random(1).sample(other sequences, k - n_ref), written in family-file
    order (the exploration rule, with the total fixed at k).
  * screening of fold f: 40 training families at 1,000 and 12 at 10,000, proportional to size category (largest
    remainder), hash order (tags "screen1k<f>", "screen10k<f>"). The core of fold f = all training cases.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import shutil
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .cases_nucleotide import allocate
from .scoring import MANIFEST_COLUMNS

SEED = 1234
SUBSET_SEED = 1
NFOLDS = 5
N_LOCKBOX = 60
N_LOCKBOX_50K = 10
SCREEN_1K, SCREEN_10K = 40, 12
KS_FIT = (1000, 10000)
XL = ("cox", "helicase_NC", "oat", "TyrKc")
XL_K = 200000
# families used by any exploration (k4/explore/EXPLORATION_LOG.md, "Data used by exploration")
EXPLORED = ("EFTU_C", "EGF_Lam", "GATA", "RuBisCO_large_NC", "biotin_lipoyl", "ghf1", "kex", "lacI", "mdd",
            "oxidored_q6", "MHC_II_C", "WW", "cyt5", "cytb", "dCMP_cyt_deam", "igcon", "planttoxin", "rnr", "slectin",
            "GATase", "Glyco_hydro_18", "galpha", "thiored")
FOLD_COLS = ("case_id", "family", "k", "size_cat", "n_family", "n_ref", "role", "fold")


def hkey(tag, fam):
    return hashlib.sha256(f"{SEED}|{tag}|{fam}".encode()).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_meta(ehf):
    rows = list(csv.DictReader(open(Path(ehf) / "extHomFam-v37.0-metadata.csv")))
    return {r["Protein family"]: {"size_cat": r["Protein family size category"],
                                  "n": int(r["Total number of sequences"]),
                                  "n_ref": int(r["Number of reference sequences"])} for r in rows}


def strat_pick(fams, meta, total, tag, min_per_cell=1):
    cells = defaultdict(list)
    for f in fams:
        cells[meta[f]["size_cat"]].append(f)
    alloc = allocate({c: len(v) for c, v in cells.items()}, total, min_per_cell)
    out = []
    for c in sorted(cells):
        out += sorted(cells[c], key=lambda f: hkey(tag, f))[:alloc[c]]
    return sorted(out)


def design(meta):
    """-> dict(lockbox=[...], lockbox50k=[...], fit={family: fold})."""
    missing = sorted(set(EXPLORED + XL) - set(meta))
    if missing:
        raise SystemExit(f"families not in the metadata: {missing}")
    elig = sorted(f for f, m in meta.items() if m["n"] >= KS_FIT[0])
    pool = [f for f in elig if f not in EXPLORED and f not in XL]
    lock = strat_pick(pool, meta, N_LOCKBOX, "lockbox")
    big = [f for f in lock if meta[f]["n"] >= 50000]
    lock50 = sorted(sorted(big, key=lambda f: hkey("lockbox50k", f))[:N_LOCKBOX_50K])
    fit = [f for f in elig if f not in lock and f not in XL]
    folds, nxt = {}, 0
    for c in sorted({meta[f]["size_cat"] for f in fit}):
        for f in sorted((f for f in fit if meta[f]["size_cat"] == c), key=lambda f: hkey("fold", f)):
            folds[f] = nxt
            nxt = (nxt + 1) % NFOLDS
    return {"lockbox": lock, "lockbox50k": lock50, "fit": folds}


def read_fasta_names(path):
    """[(name, header line, sequence)] of an unaligned FASTA, in file order."""
    out, name, head, buf = [], None, None, []
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, head, "".join(buf)))
                tok = line[1:].split()
                name, head, buf = (tok[0] if tok else ""), line, []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        out.append((name, head, "".join(buf)))
    return out


def make_subsets(job):
    """(family, family file, reference file, [k...], out dir) -> {k: (path, n_seq)}. Idempotent."""
    fam, src, ref, ks, out = job
    out = Path(out)
    paths = {k: out / "in" / f"{fam}__k{k}.fa" for k in ks}
    if all(p.exists() for p in paths.values()):
        return fam, {k: str(p) for k, p in paths.items()}
    recs = read_fasta_names(src)
    names = [r[0] for r in recs]
    if len(set(names)) != len(names):
        raise SystemExit(f"{fam}: duplicate sequence names in the family file")
    refs = [l[1:].split()[0] for l in open(ref, encoding="latin-1") if l.startswith(">")]
    idx = {n: i for i, n in enumerate(names)}
    miss = [r for r in refs if r not in idx]
    if miss:
        raise SystemExit(f"{fam}: reference sequences not in the family file: {miss}")
    ref_i = {idx[r] for r in refs}
    others = [i for i in range(len(recs)) if i not in ref_i]
    for k in ks:
        if k > len(recs):
            raise SystemExit(f"{fam}: {len(recs)} sequences < k = {k}")
        keep = ref_i | set(random.Random(SUBSET_SEED).sample(others, k - len(ref_i)))
        tmp = paths[k].with_suffix(".tmp")
        with open(tmp, "w") as fh:
            for i in sorted(keep):
                fh.write(f"{recs[i][1]}\n{recs[i][2]}\n")
        tmp.replace(paths[k])
    return fam, {k: str(p) for k, p in paths.items()}


def ref_width(path):
    w, cur, seen = None, [], False
    for line in open(path, encoding="latin-1"):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if seen:
                return len("".join(cur))
            seen = True
        elif seen:
            cur.append(line.strip())
    return len("".join(cur)) if seen else w


def prepare(ws, ehf, xl, jobs):
    ws, ehf, xl = Path(ws), Path(ehf), Path(xl)
    meta = read_meta(ehf)
    d = design(meta)
    out = ws / "data" / "large"
    (out / "in").mkdir(parents=True, exist_ok=True)
    (out / "ref").mkdir(parents=True, exist_ok=True)
    want = {}
    for f in d["fit"]:
        want[f] = [k for k in KS_FIT if meta[f]["n"] >= k]
    for f in d["lockbox"]:
        want[f] = [k for k in KS_FIT if meta[f]["n"] >= k] + ([50000] if f in d["lockbox50k"] else [])
    for f in want:
        dst = out / "ref" / f"{f}.ref.fa"
        if not dst.exists():
            shutil.copyfile(ehf / "references" / f, dst)
    js = [(f, str(ehf / "families" / f), str(out / "ref" / f"{f}.ref.fa"), ks, str(out))
          for f, ks in sorted(want.items(), key=lambda kv: -meta[kv[0]]["n"])]
    got = {}
    with ProcessPoolExecutor(jobs) as ex:
        for fam, paths in ex.map(make_subsets, js):
            got[fam] = paths
            print(f"subsets {fam}: {sorted(paths)}", flush=True)
    rows, frows = [], []
    for f in sorted(want):
        role = "lockbox" if f in d["lockbox"] else "fit"
        ref = out / "ref" / f"{f}.ref.fa"
        for k in want[f]:
            cid = f"ehf.{f}.k{k}"
            inp = Path(got[f][k])
            rows.append({"case_id": cid, "pool": "exthomfam", "benchmark": f"exthomfam_k{k}",
                         "subset": meta[f]["size_cat"], "unaligned": str(inp), "reference": str(ref),
                         "score_rule": "ref_rows_all_columns", "n_seq": k, "ref_len": ref_width(ref),
                         "sha256_unaligned": sha256_file(inp), "sha256_reference": sha256_file(ref)})
            frows.append({"case_id": cid, "family": f, "k": k, "size_cat": meta[f]["size_cat"],
                          "n_family": meta[f]["n"], "n_ref": meta[f]["n_ref"], "role": role,
                          "fold": "lockbox" if role == "lockbox" else d["fit"][f]})
    for f in XL:
        inp, ref = xl / f"{f}__{XL_K}.fa", xl / f"{f}.ref.fa"
        cid = f"ehf.{f}.k{XL_K}"
        rows.append({"case_id": cid, "pool": "exthomfam", "benchmark": f"exthomfam_k{XL_K}",
                     "subset": meta[f]["size_cat"], "unaligned": str(inp), "reference": str(ref),
                     "score_rule": "ref_rows_all_columns", "n_seq": XL_K, "ref_len": ref_width(ref),
                     "sha256_unaligned": sha256_file(inp), "sha256_reference": sha256_file(ref)})
        frows.append({"case_id": cid, "family": f, "k": XL_K, "size_cat": meta[f]["size_cat"],
                      "n_family": meta[f]["n"], "n_ref": meta[f]["n_ref"], "role": "lockbox", "fold": "lockbox"})
    _write(ws / "data" / "cases_exthomfam.tsv", MANIFEST_COLUMNS, rows)
    _write(ws / "data" / "folds_exthomfam.tsv", FOLD_COLS, frows)
    cores = ws / "data" / "cores"
    cores.mkdir(parents=True, exist_ok=True)
    info = {}
    for fold in list(range(NFOLDS)) + ["all"]:
        train = [r for r in frows if r["role"] == "fit" and (fold == "all" or r["fold"] != fold)]
        fams = sorted({r["family"] for r in train})
        f10 = sorted({r["family"] for r in train if r["k"] == 10000})
        s1 = strat_pick(fams, meta, SCREEN_1K, f"screen1k{fold}")
        s10 = strat_pick(f10, meta, SCREEN_10K, f"screen10k{fold}")
        core = sorted(r["case_id"] for r in train)
        screen = sorted([f"ehf.{f}.k1000" for f in s1] + [f"ehf.{f}.k10000" for f in s10])
        assert set(screen) <= set(core)
        (cores / f"exthomfam_fold{fold}_core.txt").write_text("".join(c + "\n" for c in core))
        (cores / f"exthomfam_fold{fold}_screen.txt").write_text("".join(c + "\n" for c in screen))
        info[str(fold)] = {"core_cases": len(core), "core_families": len(fams), "core_10k": len(f10),
                           "screen_1k": len(s1), "screen_10k": len(s10)}
    summary = {
        "metadata_sha256": sha256_file(ehf / "extHomFam-v37.0-metadata.csv"),
        "constants": {"seed": SEED, "subset_seed": SUBSET_SEED, "n_lockbox": N_LOCKBOX, "n_lockbox_50k": N_LOCKBOX_50K,
                      "screen_1k": SCREEN_1K, "screen_10k": SCREEN_10K, "explored": list(EXPLORED), "xl": list(XL)},
        "n_eligible": sum(1 for m in meta.values() if m["n"] >= KS_FIT[0]),
        "lockbox": d["lockbox"], "lockbox50k": d["lockbox50k"],
        "lockbox_by_size": dict(Counter(meta[f]["size_cat"] for f in d["lockbox"])),
        "fit_families": len(d["fit"]), "fit_by_fold": dict(Counter(d["fit"].values())),
        "fit_10k": sum(1 for f in d["fit"] if meta[f]["n"] >= 10000),
        "explored_in_fit": sorted(f for f in EXPLORED if f in d["fit"]),
        "cores": info, "n_cases": len(rows),
    }
    (out / "prepare.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("lockbox", "cores")}, indent=1))
    return summary


def _write(path, cols, rows):
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def check(ws):
    """I3: every subset holds all reference sequences, has n_seq sequences and unique names; each reference scores
    SP = TC = 1 against itself."""
    from .scoring import read_manifest
    from .scoring_largefam import score_case
    rows = read_manifest(Path(ws) / "data" / "cases_exthomfam.tsv")
    bad = []
    for r in rows:
        names = [l[1:].split()[0] for l in open(r["unaligned"], encoding="latin-1") if l.startswith(">")]
        refs = [l[1:].split()[0] for l in open(r["reference"], encoding="latin-1") if l.startswith(">")]
        if len(names) != int(r["n_seq"]) or len(set(names)) != len(names) or set(refs) - set(names):
            bad.append((r["case_id"], "subset"))
    seen = set()
    for r in rows:
        if r["reference"] in seen:
            continue
        seen.add(r["reference"])
        s = score_case(r, r["reference"])
        if abs(s["recall"] - 1.0) > 1e-9 or abs(s["tc"] - 1.0) > 1e-9:
            bad.append((r["case_id"], f"self-score SP={s['recall']} TC={s['tc']}"))
    print(f"I3: {len(rows)} cases, {len(seen)} references, failures {len(bad)}", flush=True)
    for b in bad[:20]:
        print("  ", b)
    return not bad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--ws", required=True)
    p.add_argument("--ehf", required=True)
    p.add_argument("--xl", required=True)
    p.add_argument("--jobs", type=int, default=16)
    c = sub.add_parser("check")
    c.add_argument("--ws", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "prepare":
        prepare(a.ws, a.ehf, a.xl, a.jobs)
        return 0
    return 0 if check(a.ws) else 1


if __name__ == "__main__":
    sys.exit(main())
