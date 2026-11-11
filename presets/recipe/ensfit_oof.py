"""Out-of-fold evaluation and gate of the ensemble fit (PLAN_ensemble.md section 7 + AMENDMENT_ensemble_4.md).

  run    for each fold k: the fold's pick (<fold>/selection_a4.json), the fast setting, the pre-release design with
         the fast setting as member 1 (the declared fallback) and the pre-release accurate preset on every test case
         of fold k (small cases; families at 1,000 and 10,000), 1 thread, 1,800 s per case; FAMSA 2 (default options,
         1 thread) on the same cases. Rows -> runs/ensfit/oof/rows.jsonl. Refuses any case outside fold k's test set.
  gate   pooled out-of-fold rows -> per-arm summaries and paired differences to the fast setting with 95 % bootstrap
         intervals (small: groups, large: families; 10,000 resamples, seed 1234) and the gate (criteria 1-3 and the
         1,000 part of 5; 10,000 reported only, amendment 4) for the fold picks and for the fallback
         -> runs/ensfit/oof/gate.json.

Small and 1,000: a failed case scores 0 (precision, SP, TC) and is counted. 10,000 (reported): completion counts and
means over the families every compared arm completed.

    python -m recipe.ensfit_oof run|gate [--folds 0,1,2,3,4] [--ws WS]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import ensfit as EF
from .fit import LeakError, SlotPool, _atomic_json, log, pkey
from .scoring import score_case

OOF_TIMEOUT_S = 1800.0
N_BOOT, SEED = 10000, 1234
FAMSA = Path.home() / "tmp/k4bench/tools/famsa/famsa"
KALIGN_ARMS = ["pick", "fast", "design_fast_m1", "prerelease_accurate"]
ARMS = KALIGN_ARMS + ["famsa"]
REF = "fast"


def oof_dir(ws):
    d = Path(ws) / "runs" / "ensfit" / "oof"
    d.mkdir(parents=True, exist_ok=True)
    return d


def fold_test(ws, k):
    rows, meta, core, screen, test, lockbox, guard = EF.load_ens_fold(ws, k)
    test = sorted(test)
    for name, bad in (("training core", set(test) & set(core)), ("lockbox", set(test) & set(lockbox))):
        if bad:
            raise LeakError(f"fold {k}: {len(bad)} test cases in the {name}, e.g. {sorted(bad)[:3]}")
    if {meta[c]["fold"] for c in test} != {str(k)}:
        raise LeakError(f"fold {k}: test set holds cases of another fold")
    test = [c for c in test if meta[c]["kind"] == "small" or meta[c]["k"] in (1000, 10000)]
    return rows, meta, test


def famsa_one(row, tmpdir, timeout_s):
    with tempfile.TemporaryDirectory(dir=tmpdir) as d:
        out = Path(d) / "aln.fa"
        errf = Path(d) / "stderr.txt"
        t0 = time.time()
        with open(errf, "wb") as eh:
            p = subprocess.Popen([str(FAMSA), "-t", "1", row["unaligned"], str(out)],
                                 stdout=subprocess.DEVNULL, stderr=eh)
            status, ru = None, None
            while True:
                pid, st, r = os.wait4(p.pid, os.WNOHANG)
                if pid:
                    status, ru = st, r
                    break
                if time.time() - t0 > timeout_s:
                    p.kill()
                    _, status, ru = os.wait4(p.pid, 0)
                    return {"error": "timeout", "cpu_s": ru.ru_utime + ru.ru_stime, "wall_s": time.time() - t0}
                time.sleep(0.05)
        wall, cpu = time.time() - t0, ru.ru_utime + ru.ru_stime
        rc = os.waitstatus_to_exitcode(status)
        p.returncode = rc
        if rc != 0 or not out.exists():
            return {"error": f"famsa exit {rc}: {errf.read_text(errors='replace')[-300:]}", "cpu_s": cpu,
                    "wall_s": wall}
        try:
            m = score_case(row, str(out))
        except Exception as e:  # noqa: BLE001
            return {"error": f"score {type(e).__name__}: {e}"[:300], "cpu_s": cpu, "wall_s": wall}
        return {"precision": m.get("precision"), "recall": m.get("recall"), "tc": m.get("tc"),
                "width_ratio": m.get("width_ratio"), "error": None, "cpu_s": cpu, "wall_s": wall}


def run(ws, folds):
    ws = Path(ws)
    out = oof_dir(ws)
    logfh = open(out / "progress.log", "a")
    rows_path = out / "rows.jsonl"
    done = set()
    if rows_path.exists():
        for line in open(rows_path):
            r = json.loads(line)
            done.add((r["fold"], r["arm"], r["case_id"]))
    lock = threading.Lock()
    fast, acc = EF.fast_setting(ws), EF.prerelease_accurate()
    design = EF.Space(fast, acc).seed_params()[0]
    for k in folds:
        rows, meta, test = fold_test(ws, k)
        sel = json.loads((ws / "runs" / "ensfit" / f"fold{k}" / "selection_a4.json").read_text())
        arms = {"fast": fast, "design_fast_m1": design, "prerelease_accurate": acc}
        if sel.get("pick"):
            arms["pick"] = sel["pick"]["params"]
        keys = {a: pkey(p) for a, p in arms.items()}
        if sel.get("pick") and keys["pick"] != sel["pick"]["key"]:
            raise SystemExit(f"fold {k}: pick params do not hash to their key")
        n_small = sum(meta[c]["kind"] == "small" for c in test)
        log(f"[oof fold {k}] test: {n_small} small, {len(test) - n_small} large; arms "
            f"{', '.join(f'{a}={v[:10]}' for a, v in keys.items())}", logfh)
        ev, b, cache = EF.make_evaluator(ws, rows, test, threads=EF.SLOTS, timeout_s=OOF_TIMEOUT_S)
        ev.logfh = logfh
        res = ev.evaluate({keys[a]: p for a, p in arms.items()}, test)
        with lock, open(rows_path, "a") as fh:
            for a in arms:
                for cid in test:
                    if (k, a, cid) in done:
                        continue
                    r = res[keys[a]][cid]
                    fh.write(json.dumps({"fold": k, "arm": a, "key": keys[a], "case_id": cid,
                                         **{m: r.get(m) for m in ("precision", "recall", "tc", "width_ratio",
                                                                   "cpu_s", "error")}}) + "\n")
        # FAMSA (comparison row; not cached in the fit's cache), under the machine-wide slot pool
        todo = [c for c in test if (k, "famsa", c) not in done]
        slots = SlotPool(ws / "runs" / "slots", EF.SLOTS)
        tmpdir = str(ws / "tmp")

        def job(cid):
            with slots.slot():
                r = famsa_one(rows[cid], tmpdir, OOF_TIMEOUT_S)
            with lock, open(rows_path, "a") as fh:
                fh.write(json.dumps({"fold": k, "arm": "famsa", "key": "famsa-2.5.2-default", "case_id": cid,
                                     **{m: r.get(m) for m in ("precision", "recall", "tc", "width_ratio", "cpu_s",
                                                               "error")}}) + "\n")
        with ThreadPoolExecutor(16) as ex:
            list(ex.map(job, sorted(todo, key=lambda c: -int(rows[c].get("n_seq") or 0))))
        log(f"[oof fold {k}] done", logfh)
    return 0


# ------------------------------------------------------------------------------------ gate
class Pool:
    """Small: per (group, benchmark) sums; 1,000: per family sums (failed case = 0)."""

    MET = ("precision", "recall", "tc")

    def __init__(self, rows, meta, arms):
        self.arms = arms
        small = sorted({r["case_id"] for r in rows if meta[r["case_id"]]["kind"] == "small"})
        k1 = sorted({r["case_id"] for r in rows if meta[r["case_id"]]["kind"] == "large"
                     and meta[r["case_id"]]["k"] == 1000})
        self.groups = sorted({meta[c]["group"] for c in small})
        self.bench = sorted({meta[c]["benchmark"] for c in small})
        self.fams = sorted({meta[c]["family"] for c in k1})
        gi = {g: i for i, g in enumerate(self.groups)}
        bi = {b: i for i, b in enumerate(self.bench)}
        fi = {f: i for i, f in enumerate(self.fams)}
        by = {(r["arm"], r["case_id"]): r for r in rows}
        G, B, F = len(self.groups), len(self.bench), len(self.fams)
        self.n_small = np.zeros((G, B))
        self.n_1k = np.zeros(F)
        for c in small:
            self.n_small[gi[meta[c]["group"]], bi[meta[c]["benchmark"]]] += 1
        for c in k1:
            self.n_1k[fi[meta[c]["family"]]] += 1
        self.S, self.L, self.cpu, self.err = {}, {}, {}, {}
        for a in arms:
            S, L = np.zeros((3, G, B)), np.zeros((3, F))
            cpu = [0.0, 0.0]
            ne = {"small": 0, "1k": 0}
            for c in small + k1:
                r = by.get((a, c))
                if r is None:
                    raise SystemExit(f"arm {a}: no row for {c}")
                e = bool(r["error"])
                v = [0.0 if e else float(r[m]) for m in self.MET]
                m = meta[c]
                if m["kind"] == "small":
                    S[:, gi[m["group"]], bi[m["benchmark"]]] += v
                    cpu[0] += float(r["cpu_s"] or 0.0)
                    ne["small"] += e
                else:
                    L[:, fi[m["family"]]] += v
                    cpu[1] += float(r["cpu_s"] or 0.0)
                    ne["1k"] += e
            self.S[a], self.L[a], self.cpu[a], self.err[a] = S, L, cpu, ne

    def metrics(self, a, wg, wf):
        ns = wg @ self.n_small
        out = {}
        for j, m in enumerate(("p", "sp", "tc")):
            out[f"{m}_small"] = np.mean((wg @ self.S[a][j]) / ns, axis=1)
            out[f"{m}_1k"] = (wf @ self.L[a][j]) / (wf @ self.n_1k)
        return out


def report_10k(rows, meta, arms):
    k10 = sorted({r["case_id"] for r in rows if meta[r["case_id"]]["kind"] == "large"
                  and meta[r["case_id"]]["k"] == 10000})
    by = {(r["arm"], r["case_id"]): r for r in rows}
    out = {}
    for a in arms:
        rs = [by.get((a, c)) for c in k10]
        out[a] = {"n": len(k10), "n_done": sum(1 for r in rs if r and not r["error"]),
                  "cpu_s": float(sum(float(r["cpu_s"] or 0) for r in rs if r))}
    for a in arms:
        if a == REF:
            continue
        both = [c for c in k10 if by.get((a, c)) and not by[(a, c)]["error"] and not by[(REF, c)]["error"]]
        out[a]["paired_with_fast"] = {
            "n": len(both),
            **{f"{m}_{w}": (float(np.mean([float(by[(w2, c)][m]) for c in both])) if both else None)
               for m in ("precision", "recall", "tc") for w, w2 in (("arm", a), ("fast", REF))}}
    return out


def gate(ws, folds):
    ws = Path(ws)
    out = oof_dir(ws)
    rows = [json.loads(line) for line in open(out / "rows.jsonl")]
    rows = [r for r in rows if r["fold"] in folds]
    meta = EF.load_ens_fold(ws, 0)[1]                  # meta covers every case of the data
    cases = {(r["fold"], r["case_id"]) for r in rows if r["arm"] == REF}
    arms = [a for a in ARMS if a in {r["arm"] for r in rows}]
    picks_complete = {f: any(r["arm"] == "pick" and r["fold"] == f for r in rows) for f in folds}
    for a in arms:
        got = {(r["fold"], r["case_id"]) for r in rows if r["arm"] == a}
        if a == "pick":
            need = {x for x in cases if picks_complete[x[0]]}
        else:
            need = cases
        if got != need:
            raise SystemExit(f"arm {a}: {len(need - got)} missing test rows")
    small_1k = [r for r in rows if meta[r["case_id"]]["kind"] == "small" or meta[r["case_id"]]["k"] == 1000]
    # a fold without a pick (no eligible point) counts as the fast setting for the procedure (plan section 7: the
    # procedure is evaluated as it would ship; no eligible point = no accurate setting from that fold)
    proc = []
    for r in small_1k:
        if r["arm"] == "pick":
            proc.append(r)
        elif r["arm"] == REF and not picks_complete[r["fold"]]:
            proc.append(dict(r, arm="pick"))
    rows_g = [r for r in small_1k if r["arm"] != "pick"] + proc
    P = Pool(rows_g, meta, arms)
    G, F = len(P.groups), len(P.fams)
    rng = np.random.default_rng(SEED)
    wg = np.stack([np.bincount(rng.integers(0, G, G), minlength=G) for _ in range(N_BOOT)]).astype(float)
    wf = np.stack([np.bincount(rng.integers(0, F, F), minlength=F) for _ in range(N_BOOT)]).astype(float)
    ok = np.all(wg @ P.n_small > 0, axis=1)
    wg, wf = wg[ok], wf[ok]
    one_g, one_f = np.ones((1, G)), np.ones((1, F))
    point = {a: {k: float(v[0]) for k, v in P.metrics(a, one_g, one_f).items()} for a in arms}
    boot = {a: P.metrics(a, wg, wf) for a in arms}
    summary = {}
    for a in arms:
        s = dict(point[a])
        s["errors"] = P.err[a]
        s["cpu_small"], s["cpu_1k"] = P.cpu[a]
        if a != REF:
            s["cpu_ratio_small"] = s["cpu_small"] / point_cpu(P, REF, 0)
            s["cpu_ratio_1k"] = s["cpu_1k"] / point_cpu(P, REF, 1)
            s["vs_fast"] = {m: {"diff": s[m] - point[REF][m],
                                "lo": float(np.percentile(boot[a][m] - boot[REF][m], 2.5)),
                                "hi": float(np.percentile(boot[a][m] - boot[REF][m], 97.5))}
                            for m in ("p_small", "sp_small", "tc_small", "p_1k", "sp_1k", "tc_1k")}
        summary[a] = s

    def crit(a):
        d, e = summary[a]["vs_fast"], summary[a]["errors"]
        return {"c1_small_precision_lo_gt_0": d["p_small"]["lo"] > 0,
                "c2_small_sp_diff_ge_-0.005": d["sp_small"]["diff"] >= -0.005,
                "c2_small_sp_lo_ge_-0.015": d["sp_small"]["lo"] >= -0.015,
                "c3_1k_precision_lo_gt_0": d["p_1k"]["lo"] > 0,
                "c3_1k_sp_diff_ge_-0.005": d["sp_1k"]["diff"] >= -0.005,
                "c5_every_1k_input_completes": e["1k"] == 0}
    gates = {a: crit(a) for a in ("pick", "design_fast_m1") if a in summary}
    passed = {a: all(v.values()) for a, v in gates.items()}
    decision = ("final pick" if passed.get("pick") else
                "fallback: pre-release design with the fast setting as member 1" if passed.get("design_fast_m1")
                else "no accurate setting from this experiment")
    res = {"folds": folds, "folds_without_pick": [f for f, v in picks_complete.items() if not v],
           "n_boot": int(ok.sum()), "n_boot_dropped": int((~ok).sum()), "seed": SEED, "n_groups": G,
           "n_families_1k": F, "arms": summary, "gate": gates, "passed": passed, "decision": decision,
           "reported_10k": report_10k(rows, meta, arms),
           "note": "AMENDMENT_ensemble_4.md: 10,000 reported only; lockbox not used (amendment 1)",
           "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _atomic_json(out / "gate.json", res)
    print(json.dumps({k: res[k] for k in ("folds_without_pick", "passed", "decision")}, indent=1))
    for a, s in summary.items():
        print(f"{a:22s} P_small {s['p_small']:.4f} SP_small {s['sp_small']:.4f} P_1k {s['p_1k']:.4f} "
              f"SP_1k {s['sp_1k']:.4f} errors {s['errors']}")
    return 0


def point_cpu(P, a, i):
    return max(P.cpu[a][i], 1e-9)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["run", "gate"])
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    EF.use_ws_tmp(a.ws)
    folds = [int(f) for f in a.folds.split(",")]
    return run(a.ws, folds) if a.cmd == "run" else gate(a.ws, folds)


if __name__ == "__main__":
    sys.exit(main())
