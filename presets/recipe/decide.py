"""Decision per pool x mode, 4-worker re-timing, and the generated preset table.

    python -m recipe.decide retime --pool protein [--scope both|final|oof]
    python -m recipe.decide decide --pool protein
    python -m recipe.decide emit   [--dest presets/out/presets.json]

retime  re-measures CPU at a fixed low concurrency of 4 aligner processes: the slot pool is reduced to 4
        and the other 26 machine-wide slots are held for the duration, so no other engine aligner runs.
        Results go to a SEPARATE cache (WS/cache/retime4.sqlite; the main cache keeps the first, 30-worker
        measurement), same build, same worker path.
          final: the final-fit selection of every mode and the 3.7.0 presets on the final-fit core
                 (-> the budget check of the decision)
          oof:   each fold's selected presets on that fold's test cases and the 3.7.0 presets on all cases
                 (-> reported CPU, column cpu4 of recipe.pool_results --retime)
        -> <out>/<pool>/retime/retime_rows.tsv
decide  per mode, the new preset (= the final-fit selection) ships iff
          (1) out of fold (pooled over the 5 folds, benchmark-weighted means as in the fit objective)
              d_SP >= -0.005 and d_TC >= -0.005 and (d_SP >= +0.005 or d_TC >= +0.005), d = new - 3.7.0;
          (2) within budget: CPU of the final selection on the final-fit core, re-timed at 4 workers,
              <= 1.2 x CPU of the 3.7.0 preset of that mode (protein recall: accurate), same cases, same
              re-timing;
          (3) the final selection differs from 3.7.0 (otherwise there is nothing to ship).
        Otherwise the 3.7.0 preset stays. -> <out>/<pool>/decision.json
emit    presets/out/presets.json: the parameter dict per pool x mode that ships (new or 3.7.0) with provenance
        (decision numbers, keys, build id, sha256 of the recipe files).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

from . import fit as FIT
from . import nested as NEST
from . import presets_370 as P370

MARGIN_LOSS = -0.005
MARGIN_GAIN = 0.005
RETIME_WORKERS = 4
RETIME_COLS = ("pool", "scope", "fold", "mode", "arm", "key", "case_id", "cpu_s", "wall_s", "error")


def recipe_hashes():
    d = Path(__file__).resolve().parent
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(d.glob("*.py"))}


# ------------------------------------------------------------------------------------ retime
def retime(pool, out, ws, a, scope="both", log=FIT.log, folds=None):
    folds = list(range(NEST_FOLDS)) if folds is None else folds
    from .evaluate import Build
    root = Path(out) / pool
    rd = root / "retime"
    rd.mkdir(parents=True, exist_ok=True)
    rows_all, meta, _, _, _ = FIT.load_fold_inputs(pool, "all", ws)
    build = Build.load(a.build_dir or Path(ws) / "kalign" / "960a4f6")
    slots = FIT.SlotPool(a.slots_dir or Path(ws) / "runs" / "slots", FIT.MAX_SLOTS)
    logfh = open(rd / "retime.log", "a")
    jobs = []   # (scope, fold, [(mode, arm, params, key)], case ids)
    if scope in ("both", "final"):
        d = NEST.fold_dir(out, pool, "all")
        _, sel, core = NEST.check_fit(d, a.allow_smoke)
        arms = [(m, arm, p, k) for m, arm, p, k, _, _ in NEST.arms_of(sel, pool)]
        jobs.append(("final_core", "all", arms, sorted(core)))
    if scope in ("both", "oof"):
        refs = [(m, "v370", FIT._py(P370.params_370(pool, m)), None) for m in FIT.MODES[pool]]
        refs = [(m, arm, p, FIT.pkey(p)) for m, arm, p, _ in refs]
        tests = {k: NEST.limit_ids([c for c, mm in meta.items() if mm["fold"] == k], a.test_limit) for k in folds}
        jobs.append(("oof", "all", refs, sorted(c for v in tests.values() for c in v)))
        for k in folds:
            d = NEST.fold_dir(out, pool, k)
            _, sel, _ = NEST.check_fit(d, a.allow_smoke)
            arms = [(m, arm, p, kk) for m, arm, p, kk, _, _ in NEST.arms_of(sel, pool) if arm == "new"]
            jobs.append(("oof", k, arms, tests[k]))
    allowed = set(rows_all)
    ev = FIT.Evaluator(build, rows_all.values(), allowed, cache_path=a.retime_cache or Path(ws) / "cache" /
                       "retime4.sqlite", slots=slots, threads=RETIME_WORKERS, timeout_s=NEST.TEST_TIMEOUT_S, rerun_short_timeouts=True,
                       logfh=logfh)
    log(f"{pool} retime: waiting for exclusive use of the machine (holding {FIT.MAX_SLOTS - RETIME_WORKERS} "
        f"slots, {RETIME_WORKERS} aligners)", logfh)
    slots.hold_others(RETIME_WORKERS)
    out_rows = []
    try:
        for sc, fold, arms, ids in jobs:
            specs = {k: p for _, _, p, k in arms}
            log(f"{pool} retime {sc} fold {fold}: {len(specs)} presets x {len(ids)} cases", logfh)
            res = ev.evaluate(specs, ids)
            for m, arm, p, k in arms:
                for cid in ids:
                    r = res[k][cid]
                    out_rows.append({"pool": pool, "scope": sc, "fold": fold, "mode": m, "arm": arm, "key": k,
                                     "case_id": cid, "cpu_s": r.get("cpu_s"), "wall_s": r.get("wall_s"),
                                     "error": r.get("error")})
    finally:
        slots.release_others()
    path = rd / "retime_rows.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(RETIME_COLS)
        for r in out_rows:
            w.writerow(["" if r.get(c) is None else r.get(c) for c in RETIME_COLS])
    log(f"{pool} retime: {len(out_rows)} rows -> {path}", logfh)
    return path


NEST_FOLDS = 5


def read_retime(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------------------------ decision
def rule(d_sp, d_tc, loss=MARGIN_LOSS, gain=MARGIN_GAIN):
    """The pre-stated decision rule on out-of-fold differences (new - 3.7.0)."""
    return d_sp >= loss and d_tc >= loss and (d_sp >= gain or d_tc >= gain)


def budget_check(pool, retime_rows, sel, budget_factor=FIT.BUDGET_FACTOR):
    """{mode: dict(cpu_new, cpu_ref, ratio, within)} on the final-fit core, 4-worker re-timing."""
    rows = [r for r in retime_rows if r["scope"] == "final_core"]
    tot, err = {}, {}
    for r in rows:
        tot[r["key"]] = tot.get(r["key"], 0.0) + float(r["cpu_s"] or 0.0)
        err[r["key"]] = err.get(r["key"], 0) + (1 if r["error"] else 0)
    out = {}
    for m in FIT.MODES[pool]:
        ref_key = FIT.pkey(FIT._py(P370.params_370(pool, FIT.BUDGET_REF[pool][m])))
        s = sel["selection"][m]
        k = s.get("key") or FIT.pkey(FIT._py(P370.params_370(pool, m)))
        if k not in tot or ref_key not in tot:
            raise SystemExit(f"{pool} {m}: final-core re-timing missing (run `decide retime --scope final`)")
        cpu_new, cpu_ref = tot[k], tot[ref_key]
        out[m] = {"cpu_new_4w": cpu_new, "cpu_ref_4w": cpu_ref, "ref_mode": FIT.BUDGET_REF[pool][m],
                  "ratio": cpu_new / cpu_ref if cpu_ref else float("inf"),
                  "within": cpu_new <= budget_factor * cpu_ref and err[k] == 0,
                  "errors_new": err[k], "budget_factor": budget_factor}
    return out


def decide(pool, out, allow_smoke=False, log=FIT.log):
    root = Path(out) / pool
    res = json.loads((root / "results" / "results.json").read_text())
    if (res.get("smoke") or res.get("partial")) and not allow_smoke:
        raise SystemExit(f"{pool}: pooled results are SMOKE/PARTIAL; refusing without --allow-smoke")
    dfin = NEST.fold_dir(out, pool, "all")
    run, sel, _ = NEST.check_fit(dfin, allow_smoke)
    rt = read_retime(root / "retime" / "retime_rows.tsv")
    bud = budget_check(pool, rt, sel)
    modes = {}
    for m in FIT.MODES[pool]:
        di = res["decision_inputs"].get(m)
        if di is None:
            raise SystemExit(f"{pool} {m}: no out-of-fold differences in results.json")
        bw = di["benchmark"]
        d_sp, d_tc = bw["d_sp"], bw["d_tc"]
        s = sel["selection"][m]
        ref = FIT._py(P370.params_370(pool, m))
        is370 = s.get("key") is None or s["key"] == FIT.pkey(ref)
        ok_rule = rule(d_sp, d_tc)
        ships = ok_rule and bud[m]["within"] and not is370
        reason = ("ships: OOF rule met, within budget" if ships else
                  "3.7.0 stays: final fit selected 3.7.0" if is370 else
                  "3.7.0 stays: OOF rule not met" if not ok_rule else
                  "3.7.0 stays: final selection over budget at 4 workers")
        modes[m] = {"ships_new": ships, "reason": reason,
                    "oof_benchmark_weighted": {k: bw[k] for k in bw},
                    "oof_case_weighted": di["case"], "n_oof": di["n"],
                    "oof_folds_fallback_370": di["folds_fallback"], "oof_folds_selected_370": di["folds_is_370"],
                    "rule": {"d_sp": d_sp, "d_tc": d_tc, "loss_margin": MARGIN_LOSS, "gain_margin": MARGIN_GAIN,
                             "met": ok_rule},
                    "budget": bud[m], "final_selection_key": s.get("key"), "final_is_370": is370,
                    "params": s["params"] if ships else ref, "param_key": s["key"] if ships else FIT.pkey(ref),
                    "final_selection_params": s["params"]}
        log(f"{pool} {m:8s}: dSP={d_sp:+.4f} dTC={d_tc:+.4f} budget ratio {bud[m]['ratio']:.2f} -> {reason}")
    dec = {"pool": pool, "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "rule": "ship iff OOF (benchmark-weighted, pooled) dSP>=-0.005 and dTC>=-0.005 and one >=+0.005, "
                   "final selection within 1.2x 3.7.0 CPU on the final core at 4 workers, and not 3.7.0",
           "smoke": bool(res.get("smoke") or res.get("partial") or not run.get("registered")),
           "build_id": run.get("build_id"), "final_fit_dir": str(dfin), "modes": modes}
    (root / "decision.json").write_text(json.dumps(dec, indent=1) + "\n")
    return dec


def emit(out, dest, log=FIT.log):
    pools = {}
    for pool in ("protein", "nucleotide"):
        p = Path(out) / pool / "decision.json"
        if p.exists():
            pools[pool] = json.loads(p.read_text())
    if not pools:
        raise SystemExit("no decision.json found")
    table = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
             "procedure": "experiments/redesign/docs/PLAN_pooled_refit_v1.md (pooled nested CV, v1)",
             "recipe_sha256": recipe_hashes(), "pools": {}}
    for pool, dec in pools.items():
        table["pools"][pool] = {
            "build_id": dec["build_id"], "smoke": dec["smoke"], "decided": dec["generated"],
            "modes": {m: {"source": "refit-v1" if d["ships_new"] else "3.7.0", "param_key": d["param_key"],
                          "params": d["params"], "reason": d["reason"], "rule": d["rule"],
                          "budget": d["budget"]} for m, d in dec["modes"].items()}}
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(table, indent=1) + "\n")
    log(f"preset table -> {dest}")
    return table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["retime", "decide", "emit"])
    ap.add_argument("--pool", choices=["protein", "nucleotide"])
    ap.add_argument("--ws", default=str(FIT.default_ws()))
    ap.add_argument("--out")
    ap.add_argument("--scope", default="both", choices=["both", "final", "oof"])
    ap.add_argument("--build-dir")
    ap.add_argument("--slots-dir")
    ap.add_argument("--retime-cache")
    ap.add_argument("--dest", default=str(Path(__file__).resolve().parents[1] / "out" / "presets.json"))
    ap.add_argument("--allow-smoke", action="store_true")
    ap.add_argument("--folds", default="0,1,2,3,4", help="folds to re-time (all five except in smoke tests)")
    ap.add_argument("--test-limit", type=int, default=0,
                    help="SMOKE ONLY: same evenly spaced test subset as recipe.nested --test-limit")
    a = ap.parse_args(argv)
    if a.test_limit and not a.allow_smoke:
        ap.error("--test-limit is smoke-only: needs --allow-smoke")
    ws = Path(a.ws).expanduser()
    out = Path(a.out) if a.out else ws / "runs" / "refit"
    if a.cmd in ("retime", "decide") and not a.pool:
        ap.error("--pool is required")
    if a.cmd == "retime":
        retime(a.pool, out, ws, a, a.scope, folds=[int(f) for f in a.folds.split(',')])
    elif a.cmd == "decide":
        decide(a.pool, out, a.allow_smoke)
    else:
        emit(out, a.dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
