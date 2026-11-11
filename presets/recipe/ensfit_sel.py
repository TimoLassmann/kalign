"""Amendment 2 of the confidence re-fit (docs/AMENDMENT_ensconf_2.md): selection mode instead of the consensus.

  check  S2-S4 on 10 small cases and 2 families at 1,000 of fold 0's test cases (final pick, select_only = 1)
         -> runs/enssel/check.json
  oof    for each fold k: fast, the fold's amendment-1 pick ("pick") and the same with select_only = 1 ("select") on
         fold k's test cases (small, 1,000), build 5d56876 -> runs/enssel/oof/rows.jsonl
  gate   S1 (fast and pick reproduce the ensconf out-of-fold rows of build ea02ca9) and the decision rule
         -> runs/enssel/oof/gate.json

    ~/tmp/kalign-refit/venv-fit/bin/python -m recipe.ensfit_sel check|oof|gate [--ws WS]
"""
from __future__ import annotations

import argparse
import copy
import json
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from . import ensfit as EF
from . import ensfit_conf as EC
from . import ensfit_oof as OOF
from . import evaluate as E
from .fit import SlotPool, _atomic_json, _py, log, pkey

AMENDMENT = "AMENDMENT_ensconf_2.md"
BUILD = "5d56876"
SCORER = "recipe.scoring_sel:score_case"
WORKER = Path(__file__).with_name("_align_worker_sel.py")
CHECK = Path(__file__).with_name("_sel_check.py")
ARMS = ("fast", "pick", "select")
EPS = 0.005
EXTRA = ("brier_pair", "conf_source", "ref_pairs_b", "correct_pairs", "wrong_pairs", "mean_col_conf", "calib",
         "selected_member")


def out_dir(ws, *parts):
    d = Path(ws) / "runs" / "enssel" / Path(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d


def select(p):
    q = copy.deepcopy(p)
    q["select_only"] = 1
    return _py(q)


def amend1_pick(ws, fit):
    sel = json.loads((Path(ws) / "runs" / "ensconf" / EC.fit_dir_name(fit) / "selection.json").read_text())
    if sel["pick"] is None:
        return None
    if pkey(sel["pick"]["params"]) != sel["pick"]["key"]:
        raise SystemExit(f"fit {fit}: selection params do not hash to their key")
    return sel["pick"]["params"]


def make_evaluator(ws, rows, allowed, logfh=None):
    E.WORKER_SCRIPT = WORKER
    EF.use_ws_tmp(ws)
    b = E.Build.load(Path(ws) / "kalign" / BUILD)
    slots = SlotPool(Path(ws) / "runs" / "slots", EF.SLOTS)
    ev = EF.EnsEvaluator(b, rows.values(), allowed, cache_path=Path(ws) / "cache" / "enssel.sqlite", slots=slots,
                         threads=EF.SLOTS, timeout_s=EF.TIMEOUT_S, scorer=SCORER, logfh=logfh)
    return ev, b


# ------------------------------------------------------------------------------------ S2-S4
def check(ws):
    rows, meta, test = OOF.fold_test(ws, 0)
    small = sorted(c for c in test if meta[c]["kind"] == "small")
    k1 = sorted(c for c in test if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)
    cases = small[:: max(1, len(small) // 10)][:10] + k1[:: max(1, len(k1) // 2)][:2]
    p = select(amend1_pick(ws, "all"))
    od = out_dir(ws, "check")
    job = {"params": p, "cases": [[c, rows[c]["unaligned"], str(od / c)] for c in cases]}
    (od / "job.json").write_text(json.dumps(job, indent=1))
    b = E.Build.load(Path(ws) / "kalign" / BUILD)
    subprocess.run([b.python, str(CHECK), str(od / "job.json"), str(od / "result.json")], check=True,
                   cwd=str(CHECK.parent))
    res = json.loads((od / "result.json").read_text())
    summary = {"amendment": AMENDMENT, "build": b.build_id, "n_cases": len(res),
               **{s: sum(r[s] for r in res) for s in ("S2", "S3", "S4")},
               "selected": {str(k): sum(r["selected"] == k for r in res) for k in range(5)},
               "cases": res, "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    summary["passed"] = all(summary[s] == len(res) for s in ("S2", "S3", "S4"))
    _atomic_json(out_dir(ws) / "check.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k != "cases"}, indent=1))
    return 0 if summary["passed"] else 1


# ------------------------------------------------------------------------------------ out of fold
def oof(ws, folds):
    od = out_dir(ws, "oof")
    logfh = open(od / "progress.log", "a")
    rows_path = od / "rows.jsonl"
    done = set()
    if rows_path.exists():
        done = {(r["fold"], r["arm"], r["case_id"]) for r in map(json.loads, open(rows_path))}
    fast = EF.fast_setting(ws)
    for k in folds:
        rows, meta, test = OOF.fold_test(ws, k)
        test = [c for c in test if meta[c]["kind"] == "small" or meta[c]["k"] == 1000]
        arms = {"fast": fast}
        pick = amend1_pick(ws, str(k))
        if pick is not None:
            arms["pick"], arms["select"] = pick, select(pick)
        keys = {a: pkey(p) for a, p in arms.items()}
        log(f"[enssel oof fold {k}] {len(test)} test cases; " + ", ".join(f"{a}={v[:10]}" for a, v in keys.items()),
            logfh)
        ev, b = make_evaluator(ws, rows, test, logfh)
        res = ev.evaluate({keys[a]: p for a, p in arms.items()}, test)
        with open(rows_path, "a") as fh:
            for a in arms:
                for cid in test:
                    if (k, a, cid) in done:
                        continue
                    r = res[keys[a]][cid]
                    fh.write(json.dumps({"fold": k, "arm": a, "key": keys[a], "case_id": cid, "build": b.build_id,
                                         **{m: r.get(m) for m in ("precision", "recall", "tc", "width_ratio",
                                                                   "cpu_s", "aln_sha256", "error")},
                                         **{m: (r.get("extra") or {}).get(m) for m in EXTRA}}) + "\n")
    return 0


# ------------------------------------------------------------------------------------ S1 + decision rule
def s1(ws, rows):
    """fast and pick on build 5d56876 == the ensconf out-of-fold rows (build ea02ca9) and their cached alignments."""
    old = {(r["fold"], r["arm"], r["case_id"]): r for r in map(json.loads, open(Path(ws) / "runs" / "ensconf" /
                                                                             "oof" / "rows.jsonl"))}
    con = sqlite3.connect(Path(ws) / "cache" / "ensconf.sqlite")
    old_build = E.Build.load(Path(ws) / "kalign" / EC.BUILD).build_id
    out = {"compared": 0, "metric_mismatch": [], "sha_mismatch": [], "sha_missing": 0}
    for r in rows:
        if r["arm"] not in ("fast", "pick"):
            continue
        o = old.get((r["fold"], r["arm"], r["case_id"]))
        if o is None:
            out["metric_mismatch"].append([r["fold"], r["arm"], r["case_id"], "missing in ensconf"])
            continue
        out["compared"] += 1
        for m in ("precision", "recall", "tc", "width_ratio", "brier_pair"):
            if (o[m] is None) != (r[m] is None) or (o[m] is not None and abs(o[m] - r[m]) > 1e-12):
                out["metric_mismatch"].append([r["fold"], r["arm"], r["case_id"], m, o[m], r[m]])
        q = con.execute("SELECT aln_sha256 FROM results WHERE build_id=? AND param_key=? AND case_id=?",
                        (old_build, o["key"], r["case_id"])).fetchall()
        shas = {x[0] for x in q if x[0]}
        if not shas:
            out["sha_missing"] += 1
        elif r["aln_sha256"] not in shas:
            out["sha_mismatch"].append([r["fold"], r["arm"], r["case_id"]])
    out["passed"] = not out["metric_mismatch"] and not out["sha_mismatch"]
    return out


def ece(rows, cases):
    n, sp, ny = np.zeros(10), np.zeros(10), np.zeros(10)
    for r in rows:
        if r["case_id"] in cases and r.get("calib"):
            n += r["calib"]["n"]
            sp += r["calib"]["sum_p"]
            ny += r["calib"]["correct"]
    tot = n.sum()
    return {"ece": float(np.abs(sp - ny).sum() / tot) if tot else None, "pairs": int(tot),
            "bins": {"n": n.astype(int).tolist(), "mean_p": (sp / np.maximum(n, 1)).tolist(),
                     "frac_correct": (ny / np.maximum(n, 1)).tolist()}}


def gate(ws, folds):
    od = out_dir(ws, "oof")
    rows = [r for r in map(json.loads, open(od / "rows.jsonl")) if r["fold"] in folds]
    meta = EF.load_ens_fold(ws, 0)[1]
    check_s1 = s1(ws, rows)
    cases = {(r["fold"], r["case_id"]) for r in rows if r["arm"] == "fast"}
    has = {f: any(r["arm"] == "pick" and r["fold"] == f for r in rows) for f in folds}
    for a in ("pick", "select"):
        if {(r["fold"], r["case_id"]) for r in rows if r["arm"] == a} != {x for x in cases if has[x[0]]}:
            raise SystemExit(f"arm {a}: test rows missing")
    # procedure: a fold without a pick counts as the fast setting for pick and select
    proc = rows + [dict(r, arm=a) for r in rows if r["arm"] == "fast" and not has[r["fold"]] for a in ("pick",
                                                                                                         "select")]
    with_pick = [r for r in rows if has[r["fold"]]]

    def boot_pool(rr):
        P = EC.Pool(rr, meta, list(ARMS))
        G, F = len(P.groups), len(P.fams)
        rng = np.random.default_rng(EC.SEED)
        wg = np.stack([np.bincount(rng.integers(0, G, G), minlength=G) for _ in range(EC.N_BOOT)]).astype(float)
        wf = np.stack([np.bincount(rng.integers(0, F, F), minlength=F) for _ in range(EC.N_BOOT)]).astype(float)
        ok = np.all(wg @ P.n_small > 0, axis=1)
        wg, wf = wg[ok], wf[ok]
        point = {a: {k: float(v[0]) for k, v in P.metrics(a, np.ones((1, G)), np.ones((1, F))).items()}
                 for a in ARMS}
        boot = {a: P.metrics(a, wg, wf) for a in ARMS}

        def diff(a, b):
            return {m: {"diff": point[a][m] - point[b][m],
                        "lo": float(np.percentile(boot[a][m] - boot[b][m], 2.5)),
                        "hi": float(np.percentile(boot[a][m] - boot[b][m], 97.5))} for m in point[a]}
        return P, point, diff, int(ok.sum())

    P, point, diff, nb = boot_pool(proc)
    d = diff("select", "fast")
    g = {
        "g1_sp_small": d["recall_small"]["diff"] >= 0 and d["recall_small"]["lo"] >= -EC.GATE_EPS,
        "g1_tc_small": d["tc_small"]["diff"] >= 0 and d["tc_small"]["lo"] >= -EC.GATE_EPS,
        "g2_sp_1k": d["recall_1k"]["diff"] >= 0 and d["recall_1k"]["lo"] >= -EC.GATE_EPS,
        "g2_tc_1k": d["tc_1k"]["diff"] >= 0 and d["tc_1k"]["lo"] >= -EC.GATE_EPS,
        "g3_drop_families": P.drops["select"]["below_-0.2"] <= EC.MAX_DROP_FAMILIES,
        "g4_brier_small_below_fast": d["brier_pair_small"]["hi"] < 0,
        "g4_brier_1k_below_fast": d["brier_pair_1k"]["hi"] < 0,
        "g5_every_input_completes": P.err["select"]["small"] == 0 and P.err["select"]["1k"] == 0,
    }
    P4, point4, diff4, nb4 = boot_pool(with_pick)
    sp = diff4("select", "pick")
    rule = {
        "r1_select_passes_gate": all(g.values()),
        "r2_sp_tc_noninferior": all(sp[m]["lo"] >= -EPS for m in ("recall_small", "tc_small", "recall_1k", "tc_1k")),
        "r3_brier_noninferior": all(sp[m]["hi"] <= EPS for m in ("brier_pair_small", "brier_pair_1k")),
        "r4_drop_families": P.drops["select"]["below_-0.2"] <= EC.MAX_DROP_FAMILIES,
    }
    small = {c for f, c in cases if meta[c]["kind"] == "small"}
    k1 = {c for f, c in cases if meta[c]["kind"] == "large"}
    calib = {a: {"small": ece([r for r in with_pick if r["arm"] == a], small),
                 "1k": ece([r for r in with_pick if r["arm"] == a], k1)} for a in ARMS}
    selm = {"small": {}, "1k": {}}
    for r in with_pick:
        if r["arm"] == "select":
            pool = "small" if r["case_id"] in small else "1k"
            selm[pool][str(r["selected_member"])] = selm[pool].get(str(r["selected_member"]), 0) + 1
    res = {"amendment": AMENDMENT, "folds": folds, "folds_with_pick": [f for f in folds if has[f]],
           "S1": check_s1, "n_boot": nb, "n_boot_4fold": nb4,
           "procedure": {a: {**point[a], "errors": P.err[a], "cpu_small": P.cpu[a][0], "cpu_1k": P.cpu[a][1],
                             "drops_1k": P.drops[a]} for a in ARMS},
           "procedure_vs_fast": {a: diff(a, "fast") for a in ("pick", "select")},
           "four_fold": {a: {**point4[a], "cpu_small": P4.cpu[a][0], "cpu_1k": P4.cpu[a][1], "drops_1k": P4.drops[a]}
                         for a in ARMS},
           "select_minus_pick_4fold": sp, "gate_select": g, "rule": rule,
           "select_replaces_consensus": all(rule.values()) and check_s1["passed"],
           "calibration_4fold": calib, "selected_member": selm, "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _atomic_json(od / "gate.json", res)
    print(json.dumps({"S1": {k: (v if not isinstance(v, list) else len(v)) for k, v in check_s1.items()},
                      "gate_select": g, "rule": rule, "select_replaces_consensus": res["select_replaces_consensus"],
                      "selected_member": selm}, indent=1))
    for a in ARMS:
        s = point4[a]
        print(f"4-fold {a:6s} SP {s['recall_small']:.4f}/{s['recall_1k']:.4f} TC {s['tc_small']:.4f}/{s['tc_1k']:.4f} "
              f"B {s['brier_pair_small']:.4f}/{s['brier_pair_1k']:.4f} CPU {P4.cpu[a][0]:.0f}/{P4.cpu[a][1]:.0f} "
              f"ECE {calib[a]['small']['ece']:.4f}/{calib[a]['1k']['ece']:.4f}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["check", "oof", "gate"])
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    EF.use_ws_tmp(a.ws)
    if a.cmd == "check":
        return check(a.ws)
    folds = [int(f) for f in a.folds.split(",")]
    return oof(a.ws, folds) if a.cmd == "oof" else gate(a.ws, folds)


if __name__ == "__main__":
    sys.exit(main())
