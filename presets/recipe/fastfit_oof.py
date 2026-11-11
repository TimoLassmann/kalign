"""Out-of-fold evaluation and gate of the fast-fit (PLAN_fastfit.md section 7, DECISION_fastfit_selection.md).

  run    for each fold k: the fold's picks (knee rule, its fixed-delta sensitivity picks, the registered fast and
         ultra-fast picks; from <fold>/knee_selection.json), the gm1 fast preset (new build) and 4.0.0 fast (4.0.0
         build) on every test case of fold k (its small cases; its families at 1,000 and 10,000). Rows ->
         runs/fastfit/oof/rows.jsonl. Refuses any case outside fold k's test set (and checks none is in fold k's
         training core or the lockbox).
  gate   pooled out-of-fold rows -> per-arm summaries, paired differences to 4.0.0 fast with 95 % bootstrap
         intervals (small: groups, large: families; 10,000 resamples, seed 1234) and the gate criteria 1-2 (fast) and
         the out-of-fold ultra-fast criteria -> runs/fastfit/oof/gate.json.

A failed case (error or timeout, 1,800 s) scores SP = TC = 0 and is counted.

    python -m recipe.fastfit_oof run|gate [--folds 0,1,2,3,4] [--ws WS]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from . import fit as FIT
from .fastfit import P400_FAST, REF_BUILD, NEW_BUILD, SLOTS, load_fold
from .fit import Evaluator, LeakError, SlotPool, _atomic_json, log, pkey

FIT.MAX_SLOTS = SLOTS
OOF_TIMEOUT_S = 1800.0
N_BOOT = 10000
SEED = 1234
ARMS = ["knee", "knee_d0.005", "knee_d0.01", "registered_fast", "registered_ultra", "gm1"]
REF_ARM = "p400"


def oof_dir(ws):
    d = ws / "runs" / "fastfit" / "oof"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------------------------------ run
def run(ws, folds):
    from .evaluate import Build
    out = oof_dir(ws)
    logfh = open(out / "progress.log", "a")
    b_new, b_ref = Build.load(ws / "kalign" / NEW_BUILD), Build.load(ws / "kalign" / REF_BUILD)
    slots = SlotPool(ws / "runs" / "slots", SLOTS)
    cache = ws / "cache" / "cache.sqlite"
    rows_path = out / "rows.jsonl"
    done = set()
    if rows_path.exists():
        for l in open(rows_path):
            r = json.loads(l)
            done.add((r["fold"], r["arm"]))
    for k in folds:
        rows, meta, core, _, test, lockbox = load_fold(ws, k)
        test = sorted(test)
        for name, bad in (("training core", set(test) & set(core)), ("lockbox", set(test) & set(lockbox))):
            if bad:
                raise LeakError(f"fold {k}: {len(bad)} test cases in the {name}, e.g. {sorted(bad)[:3]}")
        if {meta[c]["fold"] for c in test} != {str(k)}:
            raise LeakError(f"fold {k}: test set holds cases of another fold")
        sel = json.load(open(ws / "runs" / "fastfit" / f"fold{k}" / "knee_selection.json"))
        arms = {a: sel[a] for a in ARMS if sel.get(a)}
        ev_new = Evaluator(b_new, rows.values(), test, cache_path=cache, slots=slots, threads=SLOTS,
                           timeout_s=OOF_TIMEOUT_S, logfh=logfh, rerun_short_timeouts=True)
        ev_ref = Evaluator(b_ref, rows.values(), test, cache_path=cache, slots=slots, threads=SLOTS,
                           timeout_s=OOF_TIMEOUT_S, logfh=logfh, rerun_short_timeouts=True)
        n_small = sum(meta[c]["kind"] == "small" for c in test)
        log(f"[oof fold {k}] test: {n_small} small, {len(test) - n_small} large; arms "
            f"{', '.join(f'{a}={v['key'][:10]}' for a, v in arms.items())}", logfh)
        specs = {v["key"]: v["params"] for v in arms.values()}
        for kk, p in specs.items():
            if pkey(p) != kk:
                raise SystemExit(f"fold {k}: params do not hash to their key {kk[:12]}")
        res_new = ev_new.evaluate(specs, test)
        k400 = pkey(P400_FAST)
        res_ref = ev_ref.evaluate({k400: P400_FAST}, test)
        with open(rows_path, "a") as fh:
            for arm, v in [*arms.items(), (REF_ARM, {"key": k400})]:
                if (k, arm) in done:
                    continue
                res = res_ref[k400] if arm == REF_ARM else res_new[v["key"]]
                for cid in test:
                    r = res[cid]
                    fh.write(json.dumps({"fold": k, "arm": arm, "key": v["key"], "case_id": cid,
                                         "recall": r.get("recall"), "tc": r.get("tc"), "cpu_s": r.get("cpu_s"),
                                         "error": r.get("error")}) + "\n")
        log(f"[oof fold {k}] done", logfh)
    return 0


# ------------------------------------------------------------------------------------ gate
class Pool:
    """Pooled out-of-fold rows of all arms, as matrices over resampling units, so every bootstrap replicate is a
    matrix product. Small: per (group, benchmark) sums; large: per (family, k) sums."""

    def __init__(self, rows, meta, arms):
        self.arms = arms
        small = sorted({r["case_id"] for r in rows if meta[r["case_id"]]["kind"] == "small"})
        large = sorted({r["case_id"] for r in rows if meta[r["case_id"]]["kind"] == "large"})
        self.groups = sorted({meta[c]["group"] for c in small})
        self.fams = sorted({meta[c]["family"] for c in large})
        self.bench = sorted({meta[c]["benchmark"] for c in small})
        self.ks = sorted({meta[c]["k"] for c in large})
        gi = {g: i for i, g in enumerate(self.groups)}
        fi = {f: i for i, f in enumerate(self.fams)}
        bi = {b: i for i, b in enumerate(self.bench)}
        ki = {k: i for i, k in enumerate(self.ks)}
        by = {(r["arm"], r["case_id"]): r for r in rows}
        G, B, F, K = len(self.groups), len(self.bench), len(self.fams), len(self.ks)
        self.n_small = np.zeros((G, B))
        self.n_large = np.zeros((F, K))
        for c in small:
            self.n_small[gi[meta[c]["group"]], bi[meta[c]["benchmark"]]] += 1
        for c in large:
            self.n_large[fi[meta[c]["family"]], ki[meta[c]["k"]]] += 1
        self.S, self.L, self.cpu, self.n_err = {}, {}, {}, {}
        for a in arms:
            S = np.zeros((3, G, B))             # sp, tc, q
            L = np.zeros((3, F, K))
            cpu = np.zeros((2, max(G, F)))      # small per group, large per family
            ne = 0
            for c in small + large:
                r = by.get((a, c))
                if r is None:
                    raise SystemExit(f"arm {a}: no row for {c}")
                err = bool(r["error"])
                ne += err
                sp = 0.0 if err else float(r["recall"])
                tc = 0.0 if err else float(r["tc"])
                m = meta[c]
                if m["kind"] == "small":
                    S[:, gi[m["group"]], bi[m["benchmark"]]] += (sp, tc, 0.5 * (sp + tc))
                    cpu[0, gi[m["group"]]] += float(r["cpu_s"] or 0.0)
                else:
                    L[:, fi[m["family"]], ki[m["k"]]] += (sp, tc, 0.5 * (sp + tc))
                    cpu[1, fi[m["family"]]] += float(r["cpu_s"] or 0.0)
            self.S[a], self.L[a], self.cpu[a], self.n_err[a] = S, L, cpu, ne

    def metrics(self, a, wg, wf):
        """wg: (R, G) group multiplicities, wf: (R, F) family multiplicities -> dict of (R,) arrays."""
        G, F = len(self.groups), len(self.fams)
        ns = wg @ self.n_small                                  # (R, B)
        nl = wf @ self.n_large                                  # (R, K)
        out = {}
        for j, m in enumerate(("sp", "tc", "q")):
            out[f"{m}_small"] = np.mean((wg @ self.S[a][j]) / ns, axis=1)      # benchmark-weighted
            per_k = (wf @ self.L[a][j]) / nl
            for t, k in enumerate(self.ks):
                out[f"{m}_{k}"] = per_k[:, t]
            out[f"{m}_large"] = per_k.mean(axis=1)
        out["cpu_small"] = wg @ self.cpu[a][0, :G]
        out["cpu_large"] = wf @ self.cpu[a][1, :F]
        return out


def gate(ws, folds):
    out = oof_dir(ws)
    rows = [json.loads(l) for l in open(out / "rows.jsonl")]
    rows = [r for r in rows if r["fold"] in folds]
    _, meta, _, _, _, _ = load_fold(ws, "all")
    arms = [a for a in [REF_ARM, *ARMS] if any(r["arm"] == a for r in rows)]
    # an arm is complete only if it has a row for every test case of every fold
    cases = {(r["fold"], r["case_id"]) for r in rows if r["arm"] == REF_ARM}
    for a in arms:
        got = {(r["fold"], r["case_id"]) for r in rows if r["arm"] == a}
        if got != cases:
            raise SystemExit(f"arm {a}: {len(cases - got)} missing test rows")
    P = Pool(rows, meta, arms)
    G, F = len(P.groups), len(P.fams)
    one_g, one_f = np.ones((1, G)), np.ones((1, F))
    rng = np.random.default_rng(SEED)
    wg = np.stack([np.bincount(rng.integers(0, G, G), minlength=G) for _ in range(N_BOOT)]).astype(float)
    wf = np.stack([np.bincount(rng.integers(0, F, F), minlength=F) for _ in range(N_BOOT)]).astype(float)
    # a replicate that drops every group of a benchmark has no benchmark-weighted Q_small; such replicates are
    # left out of every interval and counted (negligible when all folds are pooled)
    ok = np.all(wg @ P.n_small > 0, axis=1)
    wg, wf = wg[ok], wf[ok]
    point = {a: {k: float(v[0]) for k, v in P.metrics(a, one_g, one_f).items()} for a in arms}
    boot = {a: P.metrics(a, wg, wf) for a in arms}
    ref, refb = point[REF_ARM], boot[REF_ARM]
    summary = {}
    for a in arms:
        s = dict(point[a])
        s["n_err"] = P.n_err[a]
        s["cpu_ratio_small"] = s["cpu_small"] / ref["cpu_small"]
        s["cpu_ratio_large"] = s["cpu_large"] / ref["cpu_large"]
        s["R"] = math.sqrt(s["cpu_ratio_small"] * s["cpu_ratio_large"])
        diffs = {}
        if a != REF_ARM:
            for m in ("sp_small", "tc_small", "q_small", "sp_large", "tc_large", "q_large",
                      *[f"{x}_{k}" for x in ("sp", "tc", "q") for k in P.ks]):
                d = boot[a][m] - refb[m]
                diffs[m] = {"diff": s[m] - ref[m], "lo": float(np.percentile(d, 2.5)),
                            "hi": float(np.percentile(d, 97.5))}
            rb = np.sqrt((boot[a]["cpu_small"] / refb["cpu_small"]) * (boot[a]["cpu_large"] / refb["cpu_large"]))
            diffs["R"] = {"value": s["R"], "lo": float(np.percentile(rb, 2.5)), "hi": float(np.percentile(rb, 97.5))}
        s["vs_p400"] = diffs
        summary[a] = s
    crit = {}
    for a in ("knee", "registered_fast", "knee_d0.005", "knee_d0.01"):
        if a in summary:
            d = summary[a]["vs_p400"]
            crit[a] = {"c1_small_sp_lo_ge_-0.005": d["sp_small"]["lo"] >= -0.005,
                       "c2_q_large_lo_gt_0": d["q_large"]["lo"] > 0.0}
    if "registered_ultra" in summary:
        d = summary["registered_ultra"]["vs_p400"]
        crit["registered_ultra"] = {"small_sp_lo_ge_-0.015": d["sp_small"]["lo"] >= -0.015,
                                    "q_large_diff_ge_0": d["q_large"]["diff"] >= 0.0,
                                    "R_le_0.8": summary["registered_ultra"]["R"] <= 0.8}
    res = {"folds": folds, "n_boot": N_BOOT, "n_boot_dropped": int((~ok).sum()), "seed": SEED, "n_groups": G, "n_families": F,
           "n_small_cases": int(P.n_small.sum()), "n_large_cases": int(P.n_large.sum()), "arms": summary,
           "criteria_oof": crit, "note": "criteria 3-4 (lockbox) are not part of the out-of-fold gate",
           "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _atomic_json(out / "gate.json", res)
    print(f"pooled out-of-fold: {res['n_small_cases']} small cases ({G} groups), {res['n_large_cases']} large "
          f"({F} families); bootstrap replicates dropped (a benchmark absent): {res['n_boot_dropped']}")
    hdr = f"{'arm':17s} {'SPs':>6s} {'TCs':>6s} {'SP1k':>6s} {'TC1k':>6s} {'SP10k':>6s} {'TC10k':>6s} {'R':>5s} err"
    print(hdr)
    for a, s in summary.items():
        print(f"{a:17s} {s['sp_small']:6.3f} {s['tc_small']:6.3f} {s['sp_1000']:6.3f} {s['tc_1000']:6.3f} "
              f"{s['sp_10000']:6.3f} {s['tc_10000']:6.3f} {s['R']:5.2f} {s['n_err']}")
    for a, s in summary.items():
        if a == REF_ARM:
            continue
        d = s["vs_p400"]
        print(f"{a:17s} vs 4.0.0: small SP {d['sp_small']['diff']:+.4f} [{d['sp_small']['lo']:+.4f}, "
              f"{d['sp_small']['hi']:+.4f}]  Q_large {d['q_large']['diff']:+.4f} [{d['q_large']['lo']:+.4f}, "
              f"{d['q_large']['hi']:+.4f}]  R {s['R']:.2f} [{d['R']['lo']:.2f}, {d['R']['hi']:.2f}]")
    for a, c in crit.items():
        print(f"{a:17s} {c}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["run", "gate"])
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--ws", default=str(Path("~/tmp/kalign-fastfit").expanduser()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    folds = [int(f) for f in a.folds.split(",")]
    return run(ws, folds) if a.cmd == "run" else gate(ws, folds)


if __name__ == "__main__":
    sys.exit(main())
