"""Nested cross-validation of the preset fit (5 grouped folds per pool) and the final fit.

    python -m recipe.nested cv    --pool protein  [--folds 0,1,2,3,4] [--out WS/runs/refit]
    python -m recipe.nested run   --pool protein  [--folds ...]      fits only (recipe.fit per fold)
    python -m recipe.nested test  --pool protein  [--folds ...]      out-of-fold evaluation only
    python -m recipe.nested final --pool protein                     final fit ("all") + in-sample scores
    python -m recipe.nested insample --pool protein                  in-sample scores of an existing final
                                                                     selection only (Amendment v2: no new fit)

run    starts one `python -m recipe.fit` process per fold, all at once. They share the machine-wide slot
       pool (WS/runs/slots, 30 lock files), so together they never run more than 30 aligner processes,
       and the machine stays at 30 concurrent aligners while any fold has work.
test   for fold k: the fold's selected presets (selection.json; a mode without a choice falls back to
       3.7.0, flagged) AND the 3.7.0 presets (parameter dicts, same worker path) are aligned on ALL test
       cases of fold k (never reduced). Guard: the test cases are exactly fold k of folds_<pool>.tsv and
       none is in the fit's training core (cases_core.txt, hash in run.json); the evaluator refuses any
       other case. -> <out>/<pool>/fold<k>/oof_rows.tsv
final  the identical fit on fold "all" (whole pool as training population), then the final selection and
       the 3.7.0 presets on every case of the pool (in-sample; for the optimism estimate)
       -> <out>/<pool>/final/insample_rows.tsv
A fit with non-registered settings (smoke test) is refused unless --allow-smoke; rows are then marked.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from . import fit as FIT

ROW_COLS = ("pool", "fold", "scope", "mode", "arm", "key", "is_370", "fallback_370", "smoke", "case_id",
            "benchmark", "difficulty", "homolog_in_train", "recall", "tc", "precision", "f1", "width_ratio",
            "cpu_s", "wall_s", "error")
TEST_TIMEOUT_S = 1800.0
PRESETS_DIR = Path(__file__).resolve().parents[1]


def default_out(ws):
    return Path(ws) / "runs" / "refit"


def fold_dir(out, pool, fold):
    return Path(out) / pool / ("final" if str(fold) == "all" else f"fold{fold}")


# ------------------------------------------------------------------------------------ run fits
def run_fits(pool, folds, out, ws, fit_args=(), log=FIT.log):
    procs = {}
    for k in folds:
        d = fold_dir(out, pool, k)
        d.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, "-m", "recipe.fit", "--pool", pool, "--fold", str(k), "--out", str(d),
               "--ws", str(ws), *fit_args]
        fh = open(d / "fit.out", "a")
        procs[k] = (subprocess.Popen(cmd, cwd=PRESETS_DIR, stdout=fh, stderr=subprocess.STDOUT), fh)
        log(f"{pool} fold {k}: fit started (pid {procs[k][0].pid}) -> {d}")
    bad = []
    for k, (p, fh) in procs.items():
        rc = p.wait()
        fh.close()
        log(f"{pool} fold {k}: fit exited {rc}")
        if rc != 0:
            bad.append(k)
    if bad:
        raise SystemExit(f"{pool}: fits of folds {bad} failed (see fit.out)")


# ------------------------------------------------------------------------------------ evaluation
def _make_evaluator(ws, rows, allowed, a, logfh=None):
    from .evaluate import Build
    build = Build.load(a.build_dir or Path(ws) / "kalign" / "960a4f6")
    slots = FIT.SlotPool(a.slots_dir or Path(ws) / "runs" / "slots", FIT.MAX_SLOTS)
    return FIT.Evaluator(build, rows.values(), allowed, cache_path=a.cache or Path(ws) / "cache" / "cache.sqlite",
                         slots=slots, threads=a.threads, timeout_s=TEST_TIMEOUT_S, rerun_short_timeouts=True, logfh=logfh)


def arms_of(selection, pool):
    """-> [(mode, arm, params, key, is_370, fallback)] for the selected presets and the 3.7.0 presets."""
    from . import presets_370 as P370
    out = []
    for m in FIT.MODES[pool]:
        s = selection["selection"][m]
        ref = FIT._py(P370.params_370(pool, m))
        if s.get("key") is None:
            out.append((m, "new", ref, FIT.pkey(ref), True, True))
        else:
            if FIT.pkey(s["params"]) != s["key"]:
                raise SystemExit(f"selection {m}: stored key does not match its parameters")
            out.append((m, "new", s["params"], s["key"], s["key"] == FIT.pkey(ref), False))
        out.append((m, "v370", ref, FIT.pkey(ref), True, False))
    return out


def write_rows(path, rows):
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ROW_COLS)
        for r in rows:
            w.writerow(["" if r.get(c) is None else r.get(c) for c in ROW_COLS])
    tmp.replace(path)


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def evaluate_arms(ev, arms, case_ids, meta, pool, fold, scope, smoke):
    specs = {k: p for _, _, p, k, _, _ in arms}
    res = ev.evaluate(specs, case_ids)
    rows = []
    for m, arm, p, k, is370, fb in arms:
        for cid in case_ids:
            r = res[k][cid]
            mm = meta[cid]
            rows.append({"pool": pool, "fold": mm["fold"] if scope == "insample" else fold, "scope": scope,
                         "mode": m, "arm": arm, "key": k, "is_370": int(is370), "fallback_370": int(fb),
                         "smoke": int(smoke), "case_id": cid, "benchmark": mm["benchmark"],
                         "difficulty": mm["difficulty"], "homolog_in_train": mm["homolog_in_train"],
                         **{c: r.get(c) for c in ("recall", "tc", "precision", "f1", "width_ratio", "cpu_s",
                                                  "wall_s", "error")}})
    return rows


def check_fit(d, allow_smoke):
    run = json.loads((d / "run.json").read_text())
    sel = json.loads((d / "selection.json").read_text())
    if not run.get("registered") and not allow_smoke:
        raise SystemExit(f"{d}: fit used non-registered (smoke) settings; refusing without --allow-smoke")
    if run.get("settings", {}).get("fake_eval") and not allow_smoke:
        raise SystemExit(f"{d}: fake-eval fit; refusing")
    core = FIT.read_ids(d / "cases_core.txt")
    if FIT.ids_hash(core) != run["core_sha256"]:
        raise SystemExit(f"{d}: cases_core.txt does not match run.json")
    return run, sel, core


def limit_ids(ids, n):
    """Smoke only: n evenly spaced ids of the sorted list (all of them when n is 0 or >= len)."""
    ids = sorted(ids)
    if not n or n >= len(ids):
        return ids
    return [ids[int(i * len(ids) / n)] for i in range(n)]


def test_fold(pool, k, out, ws, a, log=FIT.log):
    d = fold_dir(out, pool, k)
    run, sel, core = check_fit(d, a.allow_smoke)
    if run["fold"] != k or sel["fold"] != k:
        raise SystemExit(f"{d}: fit is for fold {run['fold']}, not {k}")
    rows, meta, _, _, test = FIT.load_fold_inputs(pool, k, ws)
    test = limit_ids(test, a.test_limit)
    leak = set(test) & set(core)
    if leak:
        raise FIT.LeakError(f"{pool} fold {k}: {len(leak)} test cases were in the fit's training core")
    arms = arms_of(sel, pool)
    logfh = open(d / "test.log", "a")
    log(f"{pool} fold {k}: {len(test)} test cases x {len({a_[3] for a_ in arms})} distinct presets", logfh)
    ev = _make_evaluator(ws, rows, test, a, logfh)
    out_rows = evaluate_arms(ev, arms, test, meta, pool, k, "oof", not run.get("registered") or bool(a.test_limit))
    write_rows(d / "oof_rows.tsv", out_rows)
    n_err = sum(1 for r in out_rows if r["error"])
    log(f"{pool} fold {k}: {len(out_rows)} rows ({n_err} errors) -> {d / 'oof_rows.tsv'}", logfh)


def final_insample(pool, out, ws, a, log=FIT.log):
    d = fold_dir(out, pool, "all")
    run, sel, core = check_fit(d, a.allow_smoke)
    rows, meta, _, _, _ = FIT.load_fold_inputs(pool, "all", ws)
    ids = limit_ids(rows, a.test_limit)
    arms = arms_of(sel, pool)
    logfh = open(d / "insample.log", "a")
    log(f"{pool} final: in-sample scores on all {len(ids)} cases", logfh)
    ev = _make_evaluator(ws, rows, ids, a, logfh)
    out_rows = evaluate_arms(ev, arms, ids, meta, pool, "all", "insample", not run.get("registered") or bool(a.test_limit))
    write_rows(d / "insample_rows.tsv", out_rows)
    log(f"{pool} final: {len(out_rows)} rows -> {d / 'insample_rows.tsv'}", logfh)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["cv", "run", "test", "final", "insample"])
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--ws", default=str(FIT.default_ws()))
    ap.add_argument("--out")
    ap.add_argument("--build-dir")
    ap.add_argument("--cache")
    ap.add_argument("--slots-dir")
    ap.add_argument("--threads", type=int, default=FIT.MAX_SLOTS)
    ap.add_argument("--allow-smoke", action="store_true")
    ap.add_argument("--test-limit", type=int, default=0,
                    help="SMOKE ONLY: evaluate n evenly spaced test (or in-sample) cases instead of all; rows are "
                         "marked smoke; requires --allow-smoke")
    ap.add_argument("--fit-args", default="", help="extra arguments for recipe.fit (smoke tests only)")
    a = ap.parse_args(argv)
    if a.test_limit and not a.allow_smoke:
        ap.error("--test-limit reduces the test set and is smoke-only: needs --allow-smoke")
    ws = Path(a.ws).expanduser()
    out = Path(a.out) if a.out else default_out(ws)
    folds = [int(f) for f in a.folds.split(",")]
    fit_args = a.fit_args.split() if a.fit_args else []
    for opt in ("build_dir", "cache", "slots_dir"):
        v = getattr(a, opt)
        if v:
            fit_args += [f"--{opt.replace('_', '-')}", v]
    fit_args += ["--threads", str(a.threads)]
    if a.cmd in ("cv", "run"):
        run_fits(a.pool, folds, out, ws, fit_args)
    if a.cmd in ("cv", "test"):
        for k in folds:
            test_fold(a.pool, k, out, ws, a)
    if a.cmd == "final":
        run_fits(a.pool, ["all"], out, ws, fit_args)
        final_insample(a.pool, out, ws, a)
    if a.cmd == "insample":
        final_insample(a.pool, out, ws, a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
