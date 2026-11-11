"""Evaluation of the nucfit plan (presets/docs/nucfit/PLAN_nucfit.md sections 1, 4, 5): out-of-fold (OOF),
in-sample (IS) and lockbox (LB) scores of the fast and accurate settings, the overfitting criteria O1 / O2 and the
gates.

  run  --setting fast|accurate --part oof|is|lockbox|infold
       oof      fold k's pick and the comparison arms on fold k's test cases (k = 0..4)
       infold   fold k's pick and the comparison arms on fold k's training cases: every fit case outside fold k (the
                cases fold k's fit could train on; reported beside OOF, not used by the gate) (k = 0..4)
       is       the final fit's pick and the comparison arms on all fit data (= the pooled OOF cases)
       lockbox  the final fit's pick and the comparison arms on the lockbox (refused until both final picks exist;
                their sha256 are recorded with the rows)
  gate --setting fast|accurate   -> runs/nuceval/<setting>/gate.json

Arms. fast: pick (the fit's nucfast pick), ref (pre-release nucleotide fast 4ff7e4d1, tree written out).
accurate: pick (the fit's nucconf pick; a fit without one counts as its fast pick), fast (the fit's nucfast pick),
ref (pre-release nucleotide accurate bd03cdb7). Accurate cases: no input above 5,000 sequences (Rfam large and
simulations at 1,000 in the large tier). Accurate arms are scored with the confidence scorer (pair Brier score;
single runs are scored as hard alignments).

Statistics: a failed case scores SP = TC = precision = 0 and B = 1, and is counted. 95 % bootstrap intervals over
groups (small) and families (large; every simulation is its own family), 10,000 resamples, seed 1234. Small values
are benchmark-weighted; large values are the means of the category means. IS and OOF are on the same cases and use
the same resamples (paired). O1: IS minus OOF of the pick (OOF minus IS for B) <= 0.01 or its interval contains 0.
O2: the lockbox difference pick - comparison lies inside the OOF interval of the same difference.

    python -m recipe.nucfit_eval run|gate --setting fast|accurate [--part oof|is|lockbox] [--ws WS]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from . import fastfit as FF
from . import fit as FIT
from . import nucconf as NC
from . import nucens as NE
from . import nucfast as NF
from .fit import Evaluator, LeakError, SlotPool, _atomic_json, _py, log, pkey
from .ensfit import use_ws_tmp

TIMEOUT_S = 1800.0
N_BOOT, SEED = 10000, 1234
O1_TOL = 0.01
PRIMARY = {"fast": ("sp_small", "tc_small", "q_small", "q_large"),
           "accurate": ("sp_small", "tc_small", "sp_large", "tc_large", "brier_small", "brier_large")}
COMPARE = {"fast": ("ref",), "accurate": ("fast", "ref")}


def ev_dir(ws, *parts):
    d = Path(ws) / "runs" / "nuceval" / Path(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d


def fit_name(fold):
    return "final" if str(fold) == "all" else f"fold{fold}"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def arms_of(ws, setting, fold):
    """-> {arm: params} for one fit."""
    fkey, fast = NC.fast_pick(ws, fold)
    if setting == "fast":
        return {"pick": fast, "ref": _py(NF.with_gap(FF.with_tree(NF.v2_params(NF.NUC_FAST_KEY), FF.TREE_400), 0))}
    sel = json.loads((Path(ws) / "runs" / "nucconf" / fit_name(fold) / "selection.json").read_text())
    pick = sel["pick"]["params"] if sel.get("pick") else fast
    if sel.get("pick") and pkey(pick) != sel["pick"]["key"]:
        raise SystemExit(f"nucconf {fold}: pick params do not hash to their key")
    return {"pick": pick, "fast": fast, "ref": NF.v2_params(NF.NUC_ACC_KEY)}


def cases_of(setting, meta, rows, ids):
    if setting == "fast":
        return sorted(ids)
    return sorted(c for c in ids if meta[c]["kind"] == "small"
                  or (meta[c]["k"] in NE.LARGE_K and int(rows[c]["n_seq"]) <= NE.MAX_SEQS))


def evaluator(ws, setting, rows, allowed, logfh):
    if setting == "accurate":
        return NC.make_evaluator(ws, rows, allowed, logfh)[0]
    from .evaluate import Build
    b = Build.load(Path(ws) / "kalign" / NF.BUILD)
    return Evaluator(b, rows.values(), allowed, cache_path=Path(ws) / "cache" / "cache.sqlite",
                     slots=SlotPool(Path(ws) / "runs" / "slots", FF.SLOTS), threads=FF.SLOTS, timeout_s=TIMEOUT_S,
                     logfh=logfh, rerun_short_timeouts=True)


def final_hashes(ws):
    paths = {"nucfast": Path(ws) / "runs" / "nucfast" / "final" / "selection.json",
             "nucconf": Path(ws) / "runs" / "nucconf" / "final" / "selection.json"}
    missing = [k for k, p in paths.items() if not p.exists()]
    if missing:
        raise SystemExit(f"lockbox refused: final selection missing for {missing}")
    return {k: sha(p) for k, p in paths.items()}


def write_rows(path, fold, res, keys, cases, extra=None):
    done = set()
    if path.exists():
        done = {(r["fold"], r["arm"], r["case_id"]) for r in map(json.loads, open(path))}
    with open(path, "a") as fh:
        for arm, k in keys.items():
            for c in cases:
                if (fold, arm, c) in done:
                    continue
                r = res[k][c]
                x = r.get("extra") or {}
                fh.write(json.dumps({"fold": fold, "arm": arm, "key": k, "case_id": c, "recall": r.get("recall"),
                                     "precision": r.get("precision"), "tc": r.get("tc"), "cpu_s": r.get("cpu_s"),
                                     "error": r.get("error"), "brier": x.get("brier_pair"), **(extra or {})}) + "\n")


def run(ws, setting, part):
    od = ev_dir(ws, setting, part)
    logfh = open(od / "progress.log", "a")
    jobs = []
    if part == "oof":
        for k in range(5):
            rows, meta, core, _, test, lockbox = NF.load_fold(ws, k)
            for name, bad in (("training core", set(test) & set(core)), ("lockbox", set(test) & set(lockbox))):
                if bad:
                    raise LeakError(f"fold {k}: {len(bad)} test cases in the {name}")
            jobs.append((k, rows, meta, cases_of(setting, meta, rows, test), None))
    elif part == "infold":
        for k in range(5):
            rows, meta, core, _, test, lockbox = NF.load_fold(ws, k)
            train = [c for c, m in meta.items() if m["fold"] not in ("lockbox", "excluded", str(k))]
            if set(train) & (set(test) | set(lockbox)):
                raise LeakError(f"fold {k}: in-fold cases overlap the test fold or the lockbox")
            if not set(core) <= set(train):
                raise LeakError(f"fold {k}: training core outside the in-fold cases")
            jobs.append((k, rows, meta, cases_of(setting, meta, rows, train), None))
    else:
        rows, meta, _, _, _, lockbox = NF.load_fold(ws, "all")
        if part == "is":
            ids = [c for c, m in meta.items() if m["fold"] not in ("lockbox", "excluded")]
            jobs.append(("all", rows, meta, cases_of(setting, meta, rows, ids), None))
        else:
            jobs.append(("all", rows, meta, cases_of(setting, meta, rows, lockbox), final_hashes(ws)))
    for fold, rows, meta, cases, extra in jobs:
        arms = arms_of(ws, setting, fold)
        keys = {a: pkey(p) for a, p in arms.items()}
        log(f"[nuceval {setting} {part} {fold}] {len(cases)} cases; " +
            ", ".join(f"{a}={k[:10]}" for a, k in keys.items()), logfh)
        ev = evaluator(ws, setting, rows, set(cases), logfh)
        res = ev.evaluate({keys[a]: p for a, p in arms.items()}, cases)
        write_rows(od / "rows.jsonl", fold, res, keys, cases, extra)
    log(f"[nuceval {setting} {part}] done", logfh)
    return 0


# ------------------------------------------------------------------------------------ statistics
MET = ("sp", "tc", "q", "p", "brier")
SIM_KS = ("1000", "10000")   # simulated size categories; reporting rule A2: <m>_sim = mean of their category means


class Pool:
    """Rows of several (part, arm) labels on one case set, as matrices over resampling units."""

    def __init__(self, rows, meta, labels):
        cases = sorted({r["case_id"] for r in rows})
        small = [c for c in cases if meta[c]["kind"] == "small"]
        large = [c for c in cases if meta[c]["kind"] == "large"]
        self.groups = sorted({meta[c]["group"] for c in small})
        self.fams = sorted({meta[c]["family"] for c in large})
        self.bench = sorted({meta[c]["benchmark"] for c in small})
        self.ks = sorted({meta[c]["k"] for c in large})
        gi, fi = {g: i for i, g in enumerate(self.groups)}, {f: i for i, f in enumerate(self.fams)}
        bi, ki = {b: i for i, b in enumerate(self.bench)}, {k: i for i, k in enumerate(self.ks)}
        G, B, F, K = len(self.groups), len(self.bench), len(self.fams), len(self.ks)
        self.n_small, self.n_large = np.zeros((G, B)), np.zeros((F, K))
        for c in small:
            self.n_small[gi[meta[c]["group"]], bi[meta[c]["benchmark"]]] += 1
        for c in large:
            self.n_large[fi[meta[c]["family"]], ki[meta[c]["k"]]] += 1
        by = {(r["label"], r["case_id"]): r for r in rows}
        self.S, self.L, self.err, self.cpu = {}, {}, {}, {}
        for lab in labels:
            S, L = np.zeros((len(MET), G, B)), np.zeros((len(MET), F, K))
            ne, cpu = {"small": 0, "large": 0}, [0.0, 0.0]
            for c in cases:
                r = by.get((lab, c))
                if r is None:
                    raise SystemExit(f"{lab}: no row for {c}")
                e = bool(r["error"])
                sp, tc = (0.0, 0.0) if e else (float(r["recall"]), float(r["tc"]))
                p = 0.0 if e else float(r["precision"])
                br = 1.0 if e else (float(r["brier"]) if r.get("brier") is not None else np.nan)
                v = (sp, tc, 0.5 * (sp + tc), p, br)
                m = meta[c]
                if m["kind"] == "small":
                    S[:, gi[m["group"]], bi[m["benchmark"]]] += v
                    ne["small"] += e
                    cpu[0] += float(r["cpu_s"] or 0)
                else:
                    L[:, fi[m["family"]], ki[m["k"]]] += v
                    ne["large"] += e
                    cpu[1] += float(r["cpu_s"] or 0)
            self.S[lab], self.L[lab], self.err[lab], self.cpu[lab] = S, L, ne, cpu

    def weights(self, n_boot, seed):
        G, F = len(self.groups), len(self.fams)
        rng = np.random.default_rng(seed)
        wg = np.stack([np.bincount(rng.integers(0, G, G), minlength=G) for _ in range(n_boot)]).astype(float)
        wf = np.stack([np.bincount(rng.integers(0, F, F), minlength=F) for _ in range(n_boot)]).astype(float)
        ok = np.all(wg @ self.n_small > 0, axis=1) & np.all(wf @ self.n_large > 0, axis=1)
        return wg[ok], wf[ok], int((~ok).sum())

    def metrics(self, lab, wg, wf):
        ns, nl = wg @ self.n_small, wf @ self.n_large
        out = {}
        for j, m in enumerate(MET):
            out[f"{m}_small"] = np.mean((wg @ self.S[lab][j]) / ns, axis=1)
            per_k = (wf @ self.L[lab][j]) / nl
            for t, k in enumerate(self.ks):
                out[f"{m}_{k}"] = per_k[:, t]
            out[f"{m}_large"] = per_k.mean(axis=1)
            sim = [t for t, k in enumerate(self.ks) if k in SIM_KS]
            if sim:
                out[f"{m}_sim"] = per_k[:, sim].mean(axis=1)
        return out


def summarise(P, labels, pairs):
    """-> point summaries per label and paired differences (with 95 % intervals) for (a, b) in pairs."""
    one_g, one_f = np.ones((1, len(P.groups))), np.ones((1, len(P.fams)))
    wg, wf, dropped = P.weights(N_BOOT, SEED)
    point = {lab: {k: float(v[0]) for k, v in P.metrics(lab, one_g, one_f).items()} for lab in labels}
    boot = {lab: P.metrics(lab, wg, wf) for lab in labels}
    out = {"n_small": int(P.n_small.sum()), "n_large": int(P.n_large.sum()), "n_groups": len(P.groups),
           "n_families": len(P.fams), "categories": P.ks, "n_boot": len(wg), "n_boot_dropped": dropped,
           "point": {lab: {**point[lab], "errors": P.err[lab], "cpu_small": P.cpu[lab][0],
                           "cpu_large": P.cpu[lab][1]} for lab in labels}, "diff": {}}
    for a, b in pairs:
        d = {}
        for m in point[a]:
            x = boot[a][m] - boot[b][m]
            if np.all(np.isnan(x)):
                continue
            d[m] = {"diff": point[a][m] - point[b][m], "lo": float(np.nanpercentile(x, 2.5)),
                    "hi": float(np.nanpercentile(x, 97.5))}
        out["diff"][f"{a} - {b}"] = d
    return out


def load_rows(ws, setting, part):
    p = ev_dir(ws, setting, part) / "rows.jsonl"
    return [dict(r, label=f"{part}:{r['arm']}") for r in map(json.loads, open(p))] if p.exists() else []


def gate(ws, setting):
    _, meta, *_ = NF.load_fold(ws, "all")
    oof, is_, lb = (load_rows(ws, setting, p) for p in ("oof", "is", "lockbox"))
    arms = ["pick", *COMPARE[setting]]
    oof_cases = {r["case_id"] for r in oof if r["arm"] == "pick"}
    if len({r["fold"] for r in oof}) != 5:
        raise SystemExit("OOF rows of all five folds are needed")
    res = {"setting": setting, "written": time.strftime("%Y-%m-%dT%H:%M:%S"), "primary": PRIMARY[setting]}
    # OOF and IS on the same cases, paired
    rows = oof + is_
    labels = [f"oof:{a}" for a in arms] + ([f"is:{a}" for a in arms] if is_ else [])
    if is_ and {r["case_id"] for r in is_ if r["arm"] == "pick"} != oof_cases:
        raise SystemExit("IS and OOF case sets differ")
    pairs = [(f"oof:pick", f"oof:{b}") for b in COMPARE[setting]]
    if is_:
        pairs += [(f"is:pick", f"is:{b}") for b in COMPARE[setting]] + [("is:pick", "oof:pick")]
    res["fit_data"] = summarise(Pool(rows, meta, labels), labels, pairs)
    if lb:
        lab = [f"lockbox:{a}" for a in arms]
        res["lockbox"] = summarise(Pool(lb, meta, lab), lab, [("lockbox:pick", f"lockbox:{b}")
                                                             for b in COMPARE[setting]])
        res["lockbox_final_sha256"] = lb[0].get("nucfast"), lb[0].get("nucconf")
    D = res["fit_data"]["diff"]
    crit = {}
    if is_:
        o = D["is:pick - oof:pick"]
        for m in PRIMARY[setting]:
            opt = -o[m]["diff"] if m.startswith("brier") else o[m]["diff"]
            crit[f"O1_{m}"] = bool(opt <= O1_TOL or (o[m]["lo"] <= 0 <= o[m]["hi"]))
    if lb:
        for b in COMPARE[setting]:
            L, O = res["lockbox"]["diff"][f"lockbox:pick - lockbox:{b}"], D[f"oof:pick - oof:{b}"]
            for m in PRIMARY[setting]:
                crit[f"O2_{m}_vs_{b}"] = bool(O[m]["lo"] <= L[m]["diff"] <= O[m]["hi"])
    if setting == "fast":
        d = D["oof:pick - oof:ref"]
        crit.update({"g1_sp_small_lo_ge_-0.005": d["sp_small"]["lo"] >= -0.005,
                     "g1_tc_small_ge_0": d["tc_small"]["diff"] >= 0, "g2_q_large_lo_gt_0": d["q_large"]["lo"] > 0})
        if lb:
            L = res["lockbox"]["diff"]["lockbox:pick - lockbox:ref"]
            e = res["lockbox"]["point"]["lockbox:pick"]["errors"]
            crit.update({"g4_lb_q_small_above_ref": L["q_small"]["diff"] > 0,
                         "g4_lb_q_large_above_ref": L["q_large"]["diff"] > 0, "g4_lb_large_complete": e["large"] == 0})
    else:
        d = D["oof:pick - oof:fast"]
        for m in ("sp_small", "tc_small", "sp_large", "tc_large"):
            crit[f"g_{m}"] = d[m]["diff"] >= 0 and d[m]["lo"] >= -0.005
        crit["g_brier_small_hi_lt_0"] = d["brier_small"]["hi"] < 0
        crit["g_brier_large_hi_lt_0"] = d["brier_large"]["hi"] < 0
        e = res["fit_data"]["point"]["oof:pick"]["errors"]
        crit["g_every_input_completes"] = e["small"] == 0 and e["large"] == 0
        if lb:
            L = res["lockbox"]["diff"]["lockbox:pick - lockbox:fast"]
            for m in ("sp_small", "tc_small", "sp_large", "tc_large"):
                crit[f"g_lb_{m}_ge_fast_-0.005"] = L[m]["diff"] >= -0.005
            crit["g_lb_brier_below_fast"] = L["brier_small"]["diff"] < 0 and L["brier_large"]["diff"] < 0
    res["criteria"] = crit
    res["complete"] = bool(is_ and lb)
    res["passed"] = bool(res["complete"] and all(crit.values()))
    _atomic_json(ev_dir(ws, setting) / "gate.json", res)
    print(json.dumps({"criteria": crit, "complete": res["complete"], "passed": res["passed"]}, indent=1))
    for lab, s in res["fit_data"]["point"].items():
        print(f"{lab:14s} SP {s['sp_small']:.4f}/{s['sp_large']:.4f} TC {s['tc_small']:.4f}/{s['tc_large']:.4f} "
              f"Q {s['q_small']:.4f}/{s['q_large']:.4f} err {s['errors']}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["run", "gate"])
    ap.add_argument("--setting", required=True, choices=["fast", "accurate"])
    ap.add_argument("--part", choices=["oof", "is", "lockbox", "infold"])
    ap.add_argument("--ws", default=str(NF.default_ws()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    use_ws_tmp(ws)   # temporary files under WS/tmp, never /tmp (as nucens / nucconf)
    FIT.MAX_SLOTS = FF.SLOTS
    if a.cmd == "run":
        if not a.part:
            raise SystemExit("--part is required for run")
        return run(ws, a.setting, a.part)
    return gate(ws, a.setting)


if __name__ == "__main__":
    sys.exit(main())
