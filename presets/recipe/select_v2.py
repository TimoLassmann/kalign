"""Amendment v2 (PLAN.md): two modes per pool, chosen at the knees of the (SP+TC)/2 Pareto curve.

    python -m recipe.select_v2 --pool protein --src WS/runs/v1_1/refit --dst WS/runs/v2/refit \
        --externals WS/runs/v1_1/externals_protein.tsv

For every fit of the pool (folds 0-4 and the final fit "all"), reads the fit's own core-scored Pareto front
(selection.json "front": stage-1 and stage-2 points scored on that fit's training core) and writes a v2 fit
directory: the same run.json / cases_core.txt / refs.json, the stage directories linked (read only), and a new
selection.json whose modes come from the rule below. `recipe.nested test` / `final` then align the selected
presets exactly as for v1.

Rule (fixed in PLAN.md Amendment v2 before any v2 number was computed). Quality q = (SP + TC) / 2 on the
training core, CPU = total single-thread CPU on the training core, both as stored in the front.
  marginal(p, limit): p is acceptable when no front point q with cpu(p) < cpu(q) <= limit gains
                      q.q - p.q >= GAIN_PER_DOUBLING * log2(cpu(q) / cpu(p)).
  accurate = the cheapest acceptable point with limit = infinity.
  fast     = the cheapest acceptable point among points with cpu <= limit, limit = the total CPU on the same
             core of the fastest of MAFFT L-INS-i, MUSCLE 5 and Clustal Omega, divided by FAST_VS_OTHERS.
             If no point is within that limit, fast = the cheapest front point (flagged).
  default = recall = accurate (the command line keeps all four names).
Ties: lower CPU, then key.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
from pathlib import Path

GAIN_PER_DOUBLING = 0.002
FAST_VS_OTHERS = 5.0
OTHER_TOOLS = ("mafft_linsi", "muscle5", "clustalo")
MODES = ("fast", "default", "recall", "accurate")
FITS = ("0", "1", "2", "3", "4", "all")


def q(e):
    return (e["sp"] + e["tc"]) / 2.0


def acceptable(p, pts, limit):
    for o in pts:
        if p["cpu"] < o["cpu"] <= limit and q(o) - q(p) >= GAIN_PER_DOUBLING * math.log2(o["cpu"] / p["cpu"]):
            return False
    return True


def cheapest_acceptable(pts, limit):
    c = sorted((p for p in pts if p["cpu"] <= limit), key=lambda p: (p["cpu"], -q(p), p["key"]))
    for p in c:
        if acceptable(p, c, limit):
            return p
    return None


def knee_modes(front, fastest_other_cpu):
    """-> ({"fast": entry, "accurate": entry}, info)."""
    acc = cheapest_acceptable(front, math.inf)
    limit = fastest_other_cpu / FAST_VS_OTHERS
    fast = cheapest_acceptable(front, limit)
    flag = None
    if fast is None:
        fast = min(front, key=lambda p: (p["cpu"], -q(p), p["key"]))
        flag = "no front point within the fast limit: cheapest point used"
    if fast["cpu"] > acc["cpu"]:
        fast = acc
    return {"fast": fast, "accurate": acc}, {"fast_limit_cpu": limit, "fastest_other_cpu": fastest_other_cpu,
                                             "flag": flag}


def fastest_other(externals_tsv, core_ids):
    tot = {t: 0.0 for t in OTHER_TOOLS}
    seen = {t: set() for t in OTHER_TOOLS}
    with open(externals_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["tool"] in tot and r["case_id"] in core_ids:
                tot[r["tool"]] += float(r["cpu_s"])
                seen[r["tool"]].add(r["case_id"])
    for t in OTHER_TOOLS:
        if seen[t] != core_ids:
            raise SystemExit(f"{t}: {len(core_ids - seen[t])} core cases without an external row")
    name = min(tot, key=lambda t: tot[t])
    return name, tot[name], tot


def params_by_key(src_dir):
    out = {}
    for st in ("stage1", "stage2"):
        f = src_dir / st / "front_core.json"
        for e in json.loads(f.read_text())["front"]:
            out[e["key"]] = e["params"]
    return out


def fit_dir(root, pool, fold):
    return Path(root) / pool / ("final" if fold == "all" else f"fold{fold}")


def select_pool(pool, src, dst, externals, log=print):
    summary = {}
    for fold in FITS:
        s, d = fit_dir(src, pool, fold), fit_dir(dst, pool, fold)
        d.mkdir(parents=True, exist_ok=True)
        for f in ("run.json", "cases_core.txt", "cases_screen.txt", "refs.json"):
            shutil.copy2(s / f, d / f)
        for st in ("stage1", "stage2"):
            link = d / st
            if not link.exists():
                os.symlink(os.path.relpath(s / st, d), link)
        old = json.loads((s / "selection.json").read_text())
        core = set((s / "cases_core.txt").read_text().split())
        oname, ocpu, otot = fastest_other(externals, core)
        picks, info = knee_modes(old["front"], ocpu)
        params = params_by_key(s)
        sel = {}
        for m in MODES:
            e = picks["fast" if m == "fast" else "accurate"]
            sel[m] = {"key": e["key"], "params": params[e["key"]], "source": e["source"], "is_370": False,
                      "fallback_370": False, "core": {k: e[k] for k in ("sp", "tc", "precision", "f1", "cpu")},
                      "q": q(e)}
        out = {"pool": pool, "fold": old["fold"], "registered": old["registered"], "rule": "amendment-v2",
               "gain_per_doubling": GAIN_PER_DOUBLING, "fast_vs_others": FAST_VS_OTHERS,
               "fastest_other": {"tool": oname, "cpu": ocpu, "all": otot}, **info,
               "refs_core": old["refs_core"], "selection": sel, "front": old["front"],
               "v1_1_selection": {m: old["selection"][m].get("key") for m in old["selection"]}}
        (d / "selection.json").write_text(json.dumps(out, indent=1) + "\n")
        summary[fold] = {m: (round(sel[m]["q"], 4), round(sel[m]["core"]["cpu"], 1), sel[m]["source"])
                         for m in ("fast", "accurate")}
        log(f"{pool} fit {fold}: fast q={summary[fold]['fast'][0]} cpu={summary[fold]['fast'][1]}s "
            f"({summary[fold]['fast'][2]}); accurate q={summary[fold]['accurate'][0]} "
            f"cpu={summary[fold]['accurate'][1]}s ({summary[fold]['accurate'][2]}); fast limit "
            f"{info['fast_limit_cpu']:.1f}s = {oname}/{FAST_VS_OTHERS:g}" + (f"; {info['flag']}" if info["flag"] else ""))
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--src", required=True, help="v1.1 refit directory (WS/runs/v1_1/refit)")
    ap.add_argument("--dst", required=True, help="v2 refit directory (WS/runs/v2/refit)")
    ap.add_argument("--externals", required=True, help="recipe.externals rows of the pool")
    a = ap.parse_args(argv)
    select_pool(a.pool, Path(a.src), Path(a.dst), a.externals)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
