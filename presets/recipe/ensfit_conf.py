"""Confidence re-fit of protein accurate (docs/PLAN_ensconf.md, decision D5 of the Kalign 4 paper).

Per fit (folds 0-4 and the final fit), the members of that fit's ensemble pick (runs/ensfit/<fit>/selection_a4.json)
are kept; every member takes member 1's guide tree (the fast setting's); a grid over the single-linkage tree noise of
members 2-5 (SIGMAS) and the consensus support (SUPPORTS) is scored on the fit's training core (small core + every
training family at 1,000) through the fit's evaluator with build ea02ca9 and the confidence scorer
(recipe.scoring_conf: pair Brier score). Selection: eligible = no failed case, SP and TC not below the fast setting's
(small and 1,000, point estimates on the core) and no family at 1,000 with SP more than DROP below the fast setting's;
pick = lowest J = (B_small + B_1k) / 2; ties: smaller noise, then smaller support. No eligible point: the fit has no
pick (counts as the fast setting out of fold).

  fit    [--fits 0,1,2,3,4,all]  grid + selection -> runs/ensconf/<fit>/{core.json, selection.json}
  oof    for each fold k: the pick, the fast setting, the fold's ensemble pick as fitted ("asfit") and the same with
         member 1's tree in every member ("tree1") on fold k's test cases (small, 1,000) -> runs/ensconf/oof/rows.jsonl
  gate   pooled out-of-fold rows -> paired differences with 95 % bootstrap intervals (small: groups, 1,000: families;
         10,000 resamples, seed 1234) and the gate of docs/PLAN_ensconf.md -> runs/ensconf/oof/gate.json
  check  identity checks of docs/PLAN_ensconf.md (before the fit) -> runs/ensconf/identity.json

    ~/tmp/kalign-refit/venv-fit/bin/python -m recipe.ensfit_conf fit|oof|gate|check [--ws WS]
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from . import ensfit as EF
from . import ensfit_oof as OOF
from . import evaluate as E
from .fit import LeakError, SlotPool, _atomic_json, _py, log, pkey

PLAN = "PLAN_ensconf.md"
BUILD = "ea02ca9"
SCORER = "recipe.scoring_conf:score_case"
WORKER = Path(__file__).with_name("_align_worker_conf.py")
SIGMAS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5)
SUPPORTS = (1, 2)
DROP = 0.2
FITS = ("0", "1", "2", "3", "4", "all")
FINAL_KEY = "09413e6c"
TREE = ("run_tree_method", "run_tree_dist", "run_tree_alpha", "run_tree_maxlen", "run_tree_leaf",
        "run_tree_kmeans_dist", "run_realign_tree", "run_guide_dist")
N_BOOT, SEED = 10000, 1234
OOF_ARMS = ("pick", "fast", "asfit", "tree1")
MET = ("recall", "tc", "brier_pair")
GATE_EPS = 0.005          # SP / TC lower bounds (vs the fast setting) >= -0.005
MAX_DROP_FAMILIES = 3     # out-of-fold families at 1,000 with SP < fast - DROP (the as-fitted ensemble: 9 of 308)


def out_dir(ws, *parts):
    d = Path(ws) / "runs" / "ensconf" / Path(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d


def fold_arg(fit):
    return "all" if fit == "all" else int(fit)


def fit_dir_name(fit):
    return "final" if fit == "all" else f"fold{fit}"


def ens_pick(ws, fit):
    """The fit's ensemble pick (as fitted, PLAN_ensemble.md) -> (key, params)."""
    sel = json.loads((Path(ws) / "runs" / "ensfit" / fit_dir_name(fit) / "selection_a4.json").read_text())
    p = sel["pick"]["params"]
    if pkey(p) != sel["pick"]["key"]:
        raise SystemExit(f"fit {fit}: pick params do not hash to their key")
    if fit == "all" and not sel["pick"]["key"].startswith(FINAL_KEY):
        raise SystemExit(f"final fit pick {sel['pick']['key'][:8]} is not the shipped accurate {FINAL_KEY}")
    return sel["pick"]["key"], _py(copy.deepcopy(p))


def tree1(p):
    q = copy.deepcopy(p)
    for k in TREE:
        if k in q:
            q[k] = [p[k][0]] * p["n_runs"]
    return q


def config(pick, sigma, support):
    q = tree1(pick)
    q["min_support"] = int(support)
    q["run_tree_sl_noise"] = [0.0] + [float(sigma)] * (q["n_runs"] - 1)
    return _py(q)


def grid(pick):
    return {(s, m): config(pick, s, m) for s in SIGMAS for m in SUPPORTS}


def make_evaluator(ws, rows, allowed, logfh=None):
    """recipe.ensfit's evaluator with build ea02ca9, the confidence worker and the confidence scorer."""
    E.WORKER_SCRIPT = WORKER        # ParamWorker reads the module global when it starts a worker
    EF.use_ws_tmp(ws)
    b = E.Build.load(Path(ws) / "kalign" / BUILD)
    slots = SlotPool(Path(ws) / "runs" / "slots", EF.SLOTS)
    ev = EF.EnsEvaluator(b, rows.values(), allowed, cache_path=Path(ws) / "cache" / "ensconf.sqlite", slots=slots,
                         threads=EF.SLOTS, timeout_s=EF.TIMEOUT_S, scorer=SCORER, logfh=logfh)
    return ev, b


def val(r, m):
    if m == "brier_pair":
        return float((r.get("extra") or {})[m])
    return float(r[m])


def summarise(res, rows, small, k1, fast_res=None):
    """Benchmark-weighted means over the small cases, family means at 1,000; per-family SP drops vs fast."""
    errs = [c for c in small + k1 if res[c].get("error")]
    s = {"errors": len(errs), "error_cases": errs[:10],
         "cpu_small": sum(float(res[c].get("cpu_s") or 0) for c in small),
         "cpu_1k": sum(float(res[c].get("cpu_s") or 0) for c in k1)}
    if errs:
        return s
    bench = sorted({rows[c]["benchmark"] for c in small})
    for m in MET:
        s[f"{m}_small"] = float(np.mean([np.mean([val(res[c], m) for c in small if rows[c]["benchmark"] == b])
                                         for b in bench]))
        s[f"{m}_1k"] = float(np.mean([val(res[c], m) for c in k1]))
    if fast_res is not None:
        d = [val(res[c], "recall") - val(fast_res[c], "recall") for c in k1]
        s["drops_1k"] = {"below_-0.2": sum(x < -DROP for x in d), "below_-0.1": sum(x < -0.1 for x in d),
                         "above_+0.1": sum(x > 0.1 for x in d)}
    s["J"] = (s["brier_pair_small"] + s["brier_pair_1k"]) / 2 if "brier_pair_small" in s else None
    return s


def eligible(s, f):
    return (s["errors"] == 0 and s["recall_small"] >= f["recall_small"] and s["tc_small"] >= f["tc_small"]
            and s["recall_1k"] >= f["recall_1k"] and s["tc_1k"] >= f["tc_1k"] and s["drops_1k"]["below_-0.2"] == 0)


# ------------------------------------------------------------------------------------ fit
def fit_one(ws, fit):
    od = out_dir(ws, fit_dir_name(fit))
    logfh = open(od / "progress.log", "a")
    rows, meta, core, _, _, lockbox, _ = EF.load_ens_fold(ws, fold_arg(fit))
    small = [c for c in core if meta[c]["kind"] == "small"]
    k1 = [c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 1000]
    if set(small + k1) & set(lockbox):
        raise LeakError(f"fit {fit}: lockbox case in the core")
    akey, pick = ens_pick(ws, fit)
    g = grid(pick)
    fast = EF.fast_setting(ws)
    specs = {pkey(fast): fast, **{pkey(p): p for p in g.values()}}
    log(f"[ensconf {fit}] {len(small)} small + {len(k1)} families at 1,000; ensemble pick {akey[:8]}; "
        f"{len(g)} grid points", logfh)
    ev, b = make_evaluator(ws, rows, small + k1, logfh)
    res = ev.evaluate(specs, small + k1)
    fs = summarise(res[pkey(fast)], rows, small, k1)
    if fs["errors"]:
        raise SystemExit(f"fit {fit}: the fast setting failed on {fs['error_cases']}")
    pts = []
    for (sig, sup), p in g.items():
        k = pkey(p)
        s = summarise(res[k], rows, small, k1, res[pkey(fast)])
        s.update({"sigma": sig, "min_support": sup, "key": k})
        s["eligible"] = bool(s["errors"] == 0 and eligible(s, fs))
        pts.append(s)
    el = sorted((s for s in pts if s["eligible"]), key=lambda s: (s["J"], s["sigma"], s["min_support"]))
    sel = el[0] if el else None
    _atomic_json(od / "core.json", {"fit": fit, "build": b.build_id, "ensemble_pick": akey, "fast": fs,
                                    "points": pts, "n_small": len(small), "n_1k": len(k1)})
    out = {"fit": fit, "plan": PLAN, "build": b.build_id, "ensemble_pick": akey,
           "pick": None if sel is None else {"key": sel["key"], "sigma": sel["sigma"],
                                             "min_support": sel["min_support"], "params": g[(sel["sigma"],
                                                                                            sel["min_support"])],
                                             "sha256": hashlib.sha256(json.dumps(
                                                 g[(sel["sigma"], sel["min_support"])], sort_keys=True)
                                                 .encode()).hexdigest(), "core": sel},
           "n_eligible": len(el), "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _atomic_json(od / "selection.json", out)
    log(f"[ensconf {fit}] eligible {len(el)}/{len(pts)}; pick "
        + ("none" if sel is None else f"sigma {sel['sigma']} support {sel['min_support']} J {sel['J']:.4f}"), logfh)
    for s in sorted(pts, key=lambda s: (s["sigma"], s["min_support"])):
        if s["errors"]:
            log(f"   sigma {s['sigma']:.2f} ms {s['min_support']}: {s['errors']} errors", logfh)
            continue
        log(f"   sigma {s['sigma']:.2f} ms {s['min_support']}: SP {s['recall_small']:.4f}/{s['recall_1k']:.4f} "
            f"TC {s['tc_small']:.4f}/{s['tc_1k']:.4f} B {s['brier_pair_small']:.4f}/{s['brier_pair_1k']:.4f} "
            f"J {s['J']:.4f} drops {s['drops_1k']['below_-0.2']} {'*' if s['eligible'] else ''}", logfh)
    log(f"   fast: SP {fs['recall_small']:.4f}/{fs['recall_1k']:.4f} TC {fs['tc_small']:.4f}/{fs['tc_1k']:.4f} "
        f"B {fs['brier_pair_small']:.4f}/{fs['brier_pair_1k']:.4f}", logfh)
    return out


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
        sel = json.loads((out_dir(ws, fit_dir_name(str(k))) / "selection.json").read_text())
        _, asfit = ens_pick(ws, str(k))
        arms = {"fast": fast, "asfit": asfit, "tree1": tree1(asfit)}
        if sel["pick"] is not None:
            if pkey(sel["pick"]["params"]) != sel["pick"]["key"]:
                raise SystemExit(f"fold {k}: selection params do not hash to their key")
            arms["pick"] = sel["pick"]["params"]
        keys = {a: pkey(p) for a, p in arms.items()}
        log(f"[ensconf oof fold {k}] {len(test)} test cases; " + ", ".join(f"{a}={v[:10]}" for a, v in keys.items()),
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
                                                                   "cpu_s", "error")},
                                         **{m: (r.get("extra") or {}).get(m) for m in
                                            ("brier_pair", "conf_source", "ref_pairs_b", "correct_pairs",
                                             "wrong_pairs", "mean_col_conf")}}) + "\n")
    return 0


class Pool:
    """Small: per (group, benchmark) sums; 1,000: per family sums. A failed case scores SP = TC = 0 and B = 1 (every
    reference pair missed) and is counted."""

    def __init__(self, rows, meta, arms):
        small = sorted({r["case_id"] for r in rows if meta[r["case_id"]]["kind"] == "small"})
        k1 = sorted({r["case_id"] for r in rows if meta[r["case_id"]]["kind"] == "large"})
        self.groups = sorted({meta[c]["group"] for c in small})
        self.bench = sorted({meta[c]["benchmark"] for c in small})
        self.fams = sorted({meta[c]["family"] for c in k1})
        gi = {g: i for i, g in enumerate(self.groups)}
        bi = {b: i for i, b in enumerate(self.bench)}
        fi = {f: i for i, f in enumerate(self.fams)}
        by = {(r["arm"], r["case_id"]): r for r in rows}
        G, B, F = len(self.groups), len(self.bench), len(self.fams)
        self.n_small, self.n_1k = np.zeros((G, B)), np.zeros(F)
        for c in small:
            self.n_small[gi[meta[c]["group"]], bi[meta[c]["benchmark"]]] += 1
        for c in k1:
            self.n_1k[fi[meta[c]["family"]]] += 1
        self.S, self.L, self.err, self.cpu, self.drops = {}, {}, {}, {}, {}
        for a in arms:
            S, L = np.zeros((len(MET), G, B)), np.zeros((len(MET), F))
            ne, cpu = {"small": 0, "1k": 0}, [0.0, 0.0]
            for c in small + k1:
                r = by[(a, c)]
                e = bool(r["error"])
                v = [(1.0 if m == "brier_pair" else 0.0) if e else float(r[m]) for m in MET]
                m = meta[c]
                if m["kind"] == "small":
                    S[:, gi[m["group"]], bi[m["benchmark"]]] += v
                    ne["small"] += e
                    cpu[0] += float(r["cpu_s"] or 0)
                else:
                    L[:, fi[m["family"]]] += v
                    ne["1k"] += e
                    cpu[1] += float(r["cpu_s"] or 0)
            self.S[a], self.L[a], self.err[a], self.cpu[a] = S, L, ne, cpu
            d = [(0.0 if by[(a, c)]["error"] else float(by[(a, c)]["recall"])) - float(by[("fast", c)]["recall"])
                 for c in k1]
            self.drops[a] = {"below_-0.2": sum(x < -DROP for x in d), "below_-0.1": sum(x < -0.1 for x in d),
                             "above_+0.1": sum(x > 0.1 for x in d), "n": len(d)}

    def metrics(self, a, wg, wf):
        ns = wg @ self.n_small
        out = {}
        for j, m in enumerate(MET):
            out[f"{m}_small"] = np.mean((wg @ self.S[a][j]) / ns, axis=1)
            out[f"{m}_1k"] = (wf @ self.L[a][j]) / (wf @ self.n_1k)
        return out


def gate(ws, folds):
    od = out_dir(ws, "oof")
    rows = [r for r in map(json.loads, open(od / "rows.jsonl")) if r["fold"] in folds]
    meta = EF.load_ens_fold(ws, 0)[1]
    cases = {(r["fold"], r["case_id"]) for r in rows if r["arm"] == "fast"}
    has_pick = {f: any(r["arm"] == "pick" and r["fold"] == f for r in rows) for f in folds}
    for a in ("fast", "asfit", "tree1"):
        if {(r["fold"], r["case_id"]) for r in rows if r["arm"] == a} != cases:
            raise SystemExit(f"arm {a}: test rows missing")
    # a fold without a pick counts as the fast setting (the procedure as it would ship)
    rows = [r for r in rows if r["arm"] != "pick" or has_pick[r["fold"]]] + \
        [dict(r, arm="pick") for r in rows if r["arm"] == "fast" and not has_pick[r["fold"]]]
    arms = list(OOF_ARMS)
    P = Pool(rows, meta, arms)
    G, F = len(P.groups), len(P.fams)
    rng = np.random.default_rng(SEED)
    wg = np.stack([np.bincount(rng.integers(0, G, G), minlength=G) for _ in range(N_BOOT)]).astype(float)
    wf = np.stack([np.bincount(rng.integers(0, F, F), minlength=F) for _ in range(N_BOOT)]).astype(float)
    ok = np.all(wg @ P.n_small > 0, axis=1)
    wg, wf = wg[ok], wf[ok]
    point = {a: {k: float(v[0]) for k, v in P.metrics(a, np.ones((1, G)), np.ones((1, F))).items()} for a in arms}
    boot = {a: P.metrics(a, wg, wf) for a in arms}

    def diff(a, b):
        return {m: {"diff": point[a][m] - point[b][m], "lo": float(np.percentile(boot[a][m] - boot[b][m], 2.5)),
                    "hi": float(np.percentile(boot[a][m] - boot[b][m], 97.5))} for m in point[a]}

    summary = {a: {**point[a], "errors": P.err[a], "cpu_small": P.cpu[a][0], "cpu_1k": P.cpu[a][1],
                   "drops_1k": P.drops[a]} for a in arms}
    for a in arms:
        if a != "fast":
            summary[a]["vs_fast"] = diff(a, "fast")
    summary["pick"]["vs_asfit"] = diff("pick", "asfit")
    d = summary["pick"]["vs_fast"]
    crit = {
        "g1_sp_small": d["recall_small"]["diff"] >= 0 and d["recall_small"]["lo"] >= -GATE_EPS,
        "g1_tc_small": d["tc_small"]["diff"] >= 0 and d["tc_small"]["lo"] >= -GATE_EPS,
        "g2_sp_1k": d["recall_1k"]["diff"] >= 0 and d["recall_1k"]["lo"] >= -GATE_EPS,
        "g2_tc_1k": d["tc_1k"]["diff"] >= 0 and d["tc_1k"]["lo"] >= -GATE_EPS,
        "g3_drop_families": P.drops["pick"]["below_-0.2"] <= MAX_DROP_FAMILIES,
        "g4_brier_small_below_fast": d["brier_pair_small"]["hi"] < 0,
        "g4_brier_1k_below_fast": d["brier_pair_1k"]["hi"] < 0,
        "g5_every_input_completes": P.err["pick"]["small"] == 0 and P.err["pick"]["1k"] == 0,
    }
    res = {"plan": PLAN, "folds": folds, "folds_without_pick": [f for f, v in has_pick.items() if not v],
           "n_boot": int(ok.sum()), "seed": SEED, "n_groups": G, "n_families_1k": F, "arms": summary,
           "gate": crit, "passed": all(crit.values()), "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _atomic_json(od / "gate.json", res)
    print(json.dumps({"gate": crit, "passed": res["passed"]}, indent=1))
    for a in arms:
        s = summary[a]
        print(f"{a:6s} SP {s['recall_small']:.4f}/{s['recall_1k']:.4f} TC {s['tc_small']:.4f}/{s['tc_1k']:.4f} "
              f"B {s['brier_pair_small']:.4f}/{s['brier_pair_1k']:.4f} drops {s['drops_1k']} err {s['errors']}")
    return 0


# ------------------------------------------------------------------------------------ identity checks
def check(ws, n_small=20, n_1k=2):
    """I1 the confidence worker + ea02ca9 reproduce the fit's cached rows of the shipped evaluator (build 5b93950,
    recipe/_align_worker.py) for the fast setting and the final ensemble pick: same alignment sha256 and scores;
    I2 sigma 0 with run_tree_sl_noise == the same dict without it (same alignment sha256);
    I3 correct_pairs / ref_pairs_b == the pool scorer's recall (1e-12) on every scored row;
    I4 the fast setting is scored as a hard alignment (conf_source "hard") and an ensemble with confidence
       (conf_source "ensemble"), and B of a hard alignment == (missed + wrong pairs) / reference pairs."""
    from . import cache as C
    rows, meta, core, *_ = EF.load_ens_fold(ws, "all")
    small = sorted(c for c in core if meta[c]["kind"] == "small")[:: max(1, 840 // n_small)][:n_small]
    k1 = sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)[:n_1k]
    ids = small + k1
    _, acc = ens_pick(ws, "all")
    fast = EF.fast_setting(ws)
    z = config(acc, 0.0, 2)
    z_none = {k: v for k, v in z.items() if k != "run_tree_sl_noise"}
    ev, b = make_evaluator(ws, rows, ids)
    res = ev.evaluate({pkey(p): p for p in (fast, acc, z, z_none)}, ids)
    old = C.Cache(Path(ws) / "cache" / "ensfit.sqlite")
    b0 = E.Build.load(Path(ws) / "kalign" / EF.BUILD)
    out = {"I1": [], "I2": [], "I3": [], "I4": []}
    for name, p in (("fast", fast), ("accurate", acc)):
        for c in ids:
            o = old.get(b0.build_id, pkey(p), C.ident(rows[c], E.scorer_version("recipe.scoring:score_case",
                                                                                 rows[c])))
            n = res[pkey(p)][c]
            out["I1"].append({"arm": name, "case": c, "old_cached": o is not None,
                              "same": o is not None and o["aln_sha256"] == n["aln_sha256"]
                              and all(o[m] == n[m] for m in ("recall", "precision", "tc"))})
    for c in ids:
        out["I2"].append({"case": c, "same": res[pkey(z)][c]["aln_sha256"] == res[pkey(z_none)][c]["aln_sha256"]})
    for p in (fast, acc, z):
        for c in ids:
            r = res[pkey(p)][c]
            x = r["extra"]
            out["I3"].append({"case": c, "ok": abs(x["correct_pairs"] / x["ref_pairs_b"] - r["recall"]) < 1e-12})
    for c in ids:
        xf, xa = res[pkey(fast)][c]["extra"], res[pkey(acc)][c]["extra"]
        hard = (xf["ref_pairs_b"] - xf["correct_pairs"] + xf["wrong_pairs"]) / xf["ref_pairs_b"]
        out["I4"].append({"case": c, "ok": xf["conf_source"] == "hard" and xa["conf_source"] == "ensemble"
                          and abs(xf["brier_pair"] - hard) < 1e-12})
    summ = {"I1_same_as_shipped_evaluator": f"{sum(r['same'] for r in out['I1'])}/{len(out['I1'])}",
            "I1_old_rows_found": f"{sum(r['old_cached'] for r in out['I1'])}/{len(out['I1'])}",
            "I2_sigma0_eq_none": f"{sum(r['same'] for r in out['I2'])}/{len(out['I2'])}",
            "I3_pairs_eq_recall": f"{sum(r['ok'] for r in out['I3'])}/{len(out['I3'])}",
            "I4_hard_and_ensemble": f"{sum(r['ok'] for r in out['I4'])}/{len(out['I4'])}"}
    _atomic_json(out_dir(ws) / "identity.json", {"summary": summ, **out, "build": b.build_id})
    print(json.dumps(summ, indent=1))
    passed = all(r.get("same", r.get("ok")) for v in out.values() for r in v)
    return 0 if passed else 1


# ------------------------------------------------------------------------------------ amendment 1
def amend1(ws, src=Path("~/tmp/ensconf_timeouts/rows.jsonl").expanduser()):
    """AMENDMENT_ensconf_1.md: replace each grid "timeout" row whose idle re-run finished within the limit by that
    re-run (scored with the plan's scorer); log every replacement to runs/ensconf/amend1_rows.jsonl."""
    import sqlite3
    from . import cache as C
    from .scoring_conf import score_case
    rows = EF.load_ens_fold(ws, "all")[0]
    b = E.Build.load(Path(ws) / "kalign" / BUILD)
    cpath = Path(ws) / "cache" / "ensconf.sqlite"
    cache = C.Cache(cpath)
    spec = {k: json.loads(v) for k, v in sqlite3.connect(cpath).execute("SELECT param_key, spec_json FROM specs")}
    log_path = out_dir(ws) / "amend1_rows.jsonl"
    if log_path.exists():
        raise SystemExit(f"{log_path} exists: amendment 1 was applied already")
    n = {"replaced": 0, "kept_failure": 0}
    with open(log_path, "w") as fh:
        for r in map(json.loads, open(src)):
            key, cid = r["key"], r["case_id"]
            idt = C.ident(rows[cid], E.scorer_version(SCORER, rows[cid]))
            old = cache.get(b.build_id, key, idt)
            if old is None or old["error"] != "timeout":
                raise SystemExit(f"{key[:12]} {cid}: cached row is not a timeout ({old and old['error']})")
            m = json.loads(r["out"]) if r["rc"] == 0 else {"error": f"rc {r['rc']}", "wall": None, "cpu": None}
            entry = {"key": key, "case_id": cid, "fit": r["fit"], "sigma": r["sigma"], "min_support": r["min_support"],
                     "old_cpu_s": old["cpu_s"], "old_wall_s": old["wall_s"], "idle_cpu_s": m["cpu"],
                     "idle_wall_s": m["wall"]}
            if m["error"] is None and m["wall"] < EF.TIMEOUT_S:
                fa = src.parent / key[:12] / f"{cid}.fa"
                metrics = score_case(rows[cid], str(fa))
                if not cache.put(b.build_id, key, idt, metrics, cpu_s=m["cpu"], wall_s=m["wall"],
                                 aln_sha256=C.sha256_file(fa), error=None, spec=spec[key],
                                 harness=E.HARNESS + "+amend1", replace_error=True):
                    raise SystemExit(f"{key[:12]} {cid}: replacement not written")
                entry.update(action="replaced", recall=metrics["recall"], brier_pair=metrics["brier_pair"])
                n["replaced"] += 1
            else:
                entry["action"] = "kept_failure"
                n["kept_failure"] += 1
            fh.write(json.dumps(entry) + "\n")
    print(json.dumps(n))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["fit", "oof", "gate", "check", "amend1"])
    ap.add_argument("--fits", default=",".join(FITS))
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    EF.use_ws_tmp(a.ws)
    if a.cmd == "check":
        return check(a.ws)
    if a.cmd == "amend1":
        return amend1(a.ws)
    if a.cmd == "fit":
        for f in a.fits.split(","):
            fit_one(a.ws, f)
        return 0
    folds = [int(f) for f in a.folds.split(",")]
    return oof(a.ws, folds) if a.cmd == "oof" else gate(a.ws, folds)


if __name__ == "__main__":
    sys.exit(main())
