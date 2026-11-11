"""Before launch (PLAN_ensemble.md sections 3-4): per fold, the 24 selection families at 10,000 and the 10,000
timing-guard family.

    python -m recipe.ensfit_prep select10k [--ws WS]

Rules:
  * selection families of fold f (0..4, all): 24 families drawn from the fold's training families that have 10,000
    sequences, seeded 1234 and stratified by extHomFam size category, with the fast fit's draw rule
    (recipe.largefam.strat_pick: largest-remainder allocation, at least 1 per category, hash order
    sha256("1234|select10k<f>|<family>") within a category);
  * guard of fold f: the fast setting (single thread, the fit's evaluator and cache) aligns the 24 inputs; the family
    with the median CPU (24 values: the lower median, the 12th smallest; ties by family name) is the guard.
Training data only: every input is a training case of its fold (load_fold core); no lockbox input is touched.

Output: WS/runs/ensfit/select10k.json (written once; refuses to overwrite a different selection).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

from . import ensfit as EF
from .fastfit import load_fold
from .largefam import strat_pick

FOLDS = [0, 1, 2, 3, 4, "all"]


def size_cats(ws):
    out = {}
    for r in csv.DictReader(open(Path(ws) / "data" / "folds_exthomfam.tsv"), delimiter="\t"):
        out[r["family"]] = {"size_cat": r["size_cat"]}
    return out


def draw(ws):
    cats = size_cats(ws)
    sel, allowed, rows_all = {}, set(), None
    for f in FOLDS:
        rows, meta, core, screen, test, lockbox = load_fold(ws, f)
        rows_all = rows
        f10 = sorted({meta[c]["family"] for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 10000})
        fams = strat_pick(f10, cats, EF.N_SELECT_10K, f"select10k{f}")
        ids = [f"ehf.{x}.k10000" for x in fams]
        if not set(ids) <= set(core) or set(ids) & (set(test) | set(lockbox)):
            raise EF.LeakError(f"fold {f}: selection outside the training data")
        sel[str(f)] = {"families": fams, "n_training_10k": len(f10),
                       "by_size_cat": {c: sum(1 for x in fams if cats[x]["size_cat"] == c)
                                       for c in sorted({cats[x]["size_cat"] for x in fams})}}
        allowed |= set(ids)
    return sel, allowed, rows_all


def select10k(ws):
    ws = Path(ws)
    out = EF.select10k_path(ws)
    out.parent.mkdir(parents=True, exist_ok=True)
    sel, allowed, rows = draw(ws)
    ev, b, cache = EF.make_evaluator(ws, rows, allowed, threads=24)
    fast = EF.fast_setting(ws)
    kf = EF.pkey(fast)
    t0 = time.time()
    res = ev.evaluate({kf: fast}, sorted(allowed))[kf]
    bad = {c: r.get("error") for c, r in res.items() if r.get("error")}
    if bad:
        raise SystemExit(f"the fast setting failed on {len(bad)} inputs: {bad}")
    cpu = {c.split(".")[1]: float(r["cpu_s"]) for c, r in res.items()}
    for f, d in sel.items():
        order = sorted(d["families"], key=lambda x: (cpu[x], x))
        d["fast_cpu_s"] = {x: cpu[x] for x in order}
        d["guard"] = order[(len(order) - 1) // 2]
        d["guard_fast_cpu_s"] = cpu[d["guard"]]
    doc = {"plan": EF.PLAN, "plan_sha256": EF.PLAN_SHA256, "rule": __doc__.split("Rules:")[1].split("Output:")[0],
           "build": b.build_id, "key_fast": kf, "cache": str(cache), "n_inputs": len(allowed),
           "measured": time.strftime("%Y-%m-%dT%H:%M:%S"), "wall_s": time.time() - t0, "folds": sel}
    if out.exists():
        old = json.loads(out.read_text())
        if {f: (d["families"], d["guard"]) for f, d in old["folds"].items()} != \
                {f: (d["families"], d["guard"]) for f, d in sel.items()}:
            raise SystemExit(f"{out} exists with a different selection; refusing to overwrite")
        print(f"{out} exists with the same selection")
    else:
        out.write_text(json.dumps(doc, indent=1) + "\n")
    for f, d in sel.items():
        print(f"fold {f}: {len(d['families'])} of {d['n_training_10k']} training families at 10,000 "
              f"{d['by_size_cat']}; guard {d['guard']} (fast CPU {d['guard_fast_cpu_s']:.1f} s; range "
              f"{min(d['fast_cpu_s'].values()):.1f}-{max(d['fast_cpu_s'].values()):.1f} s)")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["select10k"])
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    return select10k(Path(a.ws).expanduser())


if __name__ == "__main__":
    sys.exit(main())
