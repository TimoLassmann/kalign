"""Re-scoring and selection of one ensemble fit under AMENDMENT_ensemble_4.md (replaces steps 3-5 of recipe.ensfit,
whose own selection.json is kept but superseded).

  1. reference: the fast setting on the screening data and on the core without 10,000 (small core + every training
     family at 1,000);
  2. re-scoring set = the plan's set (screen front, P_small >= fast - 0.01, <= 40) UNION the 40 archive points with
     the highest screening P_1k among those that beat the fast setting on screening precision (small and 1,000) with
     screening SP >= the fast setting's SP (no 0.005 slack), small and 1,000;
  3. every re-scored point on the core small + 1,000 cases, run to completion (no early stop);
  4. eligibility (plan section 6 step 1 without the 10,000 part): no failed case, P_small > fast,
     SP_small >= fast - 0.005, SP_1k >= fast - 0.005; then the knee rule unchanged (recipe.ensfit.knee);
  5. reported, not used: the pick and the best point on the fold's 24 selection families at 10,000 (1,800 s, 1 thread,
     timeouts allowed), the pre-release accurate preset and the pre-release design with the fast setting as member 1
     (seed 1, the declared fallback) on the core.

    python -m recipe.ensfit_select --fit final|fold0..fold4 [--ws WS]
Output: <fit>/selection_a4.json, <fit>/core_a4.json; log lines in <fit>/progress.log tagged [a4].
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from . import cache as C
from . import ensfit as EF
from .fit import _atomic_json, nondominated, pkey

AMENDMENT4 = "AMENDMENT_ensemble_4.md"
A4_MAX = 40


def eligible_a4(c, ref):
    """Plan section 6 step 1 without the 10,000 requirement (AMENDMENT_ensemble_4.md (a))."""
    s = c["core"]
    return bool(s["valid"] and s["p_small"] > ref["p_small"] and s["sp_small"] >= ref["sp_small"] - EF.SP_FLOOR
                and s["sp_1k"] >= ref["sp_1k"] - EF.SP_FLOOR)


def fit_dir(ws, name):
    return Path(ws) / "runs" / "ensfit" / name


def load_archive(d):
    arch = {}
    for line in open(d / "evaluated.jsonl"):
        try:
            e = json.loads(line)
        except json.JSONDecodeError:
            continue
        arch.setdefault(e["key"], e)
    return arch


def rescore_set_a4(archive, ref_screen):
    ents = sorted((e for e in archive.values() if e["screen"]["feasible"]), key=lambda e: e["key"])
    front = [ents[i] for i in nondominated([EF.raw_F(e["screen"]) for e in ents])]
    plan_set, n_ok = EF.rescore_set(front, ref_screen)
    rs = ref_screen
    margin = [e for e in ents if e["screen"]["sp_small"] >= rs["sp_small"] and e["screen"]["sp_1k"] >= rs["sp_1k"]
              and e["screen"]["p_small"] > rs["p_small"] and e["screen"]["p_1k"] > rs["p_1k"]]
    a4_set = sorted(margin, key=lambda e: (-e["screen"]["p_1k"], e["key"]))[:A4_MAX]
    out = {e["key"]: dict(e, why=["plan"]) for e in plan_set}
    for e in a4_set:
        out.setdefault(e["key"], dict(e, why=[]))["why"].append("a4")
    info = {"n_feasible_screen": len(ents), "n_front_screen": len(front), "n_plan_qualified": n_ok,
            "n_plan": len(plan_set), "n_a4_qualified": len(margin), "n_a4": len(a4_set), "n_union": len(out)}
    return sorted(out.values(), key=lambda e: e["key"]), info


def select(ws, name):
    ws = Path(ws)
    d = fit_dir(ws, name)
    run = json.loads((d / "run.json").read_text())
    fold = run["fold"]
    rows, meta, core, screen, test, lockbox, guard = EF.load_ens_fold(ws, fold)
    small_core = sorted(c for c in core if meta[c]["kind"] == "small")
    k1_core = sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)
    k10_sel = sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 10000)
    core_a4 = small_core + k1_core
    ev, b, cache = EF.make_evaluator(ws, rows, core, threads=EF.SLOTS, timeout_s=EF.TIMEOUT_S)
    if b.build_id != run["build"]:
        raise SystemExit(f"{name}: build {b.build_id[:12]} != the fit's {run['build'][:12]}")
    fast, acc = EF.fast_setting(ws), EF.prerelease_accurate(run.get("presets_v2", EF.PRESETS_V2))
    if pkey(fast) != run["key_fast"] or pkey(acc) != run["key_acc"]:
        raise SystemExit(f"{name}: reference settings differ from the fit's run.json")
    f = EF.EnsFit(fold, d, rows, meta, core_a4, screen, test, lockbox, None, ev, EF.Settings(), {}, fast, acc)
    ev.logfh = f.logfh
    tag = "[a4]"
    f.log(f"{tag} {AMENDMENT4}: re-scoring and selection without the 10,000 requirement")

    kf = pkey(fast)
    f.ref = {}
    for scope in ("screen", "core"):
        s = f.score({kf: fast}, scope, raw=True)[kf]
        if not s["valid"]:
            raise RuntimeError(f"the fast setting failed on {s['n_err']} {scope} cases")
        f.ref[scope] = EF.objectives(s, s)
        f.summ[(scope, kf)] = f.ref[scope]
    ref = f.ref["core"]
    f.log(f"{tag} fast on core small+1,000 ({len(core_a4)}): P_small={ref['p_small']:.4f} "
          f"SP_small={ref['sp_small']:.4f} P_1k={ref['p_1k']:.4f} SP_1k={ref['sp_1k']:.4f}")

    archive = load_archive(d)
    seeds = f.space.seed_params()
    design = seeds[0]                                     # pre-release design, fast as member 1 (fallback)
    kd, ka = pkey(design), pkey(acc)
    cand_e, info = rescore_set_a4(archive, f.ref["screen"])
    f.log(f"{tag} re-scoring set: {info}")
    specs = {e["key"]: e["params"] for e in cand_e}
    specs[kd] = design
    specs[ka] = acc
    t0 = time.time()
    core_s = f.score(specs, "core")                       # no early stop: every core case runs to completion
    f.log(f"{tag} core re-scoring done ({len(specs)} settings, {time.time() - t0:.0f}s)")

    cands = [{"key": e["key"], "params": e["params"], "source": f"search gen {e['gen']}", "why": e["why"],
              "screen": e["screen"], "core": core_s[e["key"]]} for e in cand_e]
    EF.eligible = eligible_a4                              # the knee rule reads the module's eligibility test
    elig = [c["key"] for c in cands if eligible_a4(c, ref)]
    k = EF.knee(cands, ref, f.fam_precision(elig))
    by = {c["key"]: c for c in cands}

    def brief(key, why):
        if key is None:
            return None
        c = by[key]
        canon = C.canonical_json(c["params"])
        return {"key": key, "why": why, "source": c["source"], "rescore_set": c["why"], "params": c["params"],
                "params_sha256": hashlib.sha256(canon.encode()).hexdigest(), "short": f.space.short(c["params"]),
                "screen": c["screen"], "core": EF.brief_s(c["core"])}

    # reported: 10,000 on the fold's 24 selection families (timeouts allowed, no early stop)
    rep10 = {}
    keys10 = [x for x in dict.fromkeys([k["pick"], k["b"]]) if x]
    if keys10:
        res = ev.evaluate({**{x: by[x]["params"] for x in keys10}, kf: fast}, k10_sel)
        fr = res[kf]
        for x in keys10:
            r = res[x]
            ok = [c for c in k10_sel if not r[c].get("error") and not fr[c].get("error")]
            rep10[x] = {"n": len(k10_sel), "n_done": sum(1 for c in k10_sel if not r[c].get("error")),
                        "errors": sorted({str(r[c].get("error")) for c in k10_sel if r[c].get("error")}),
                        "n_paired": len(ok),
                        "p_10k": float(np.mean([r[c]["precision"] for c in ok])) if ok else None,
                        "sp_10k": float(np.mean([r[c]["recall"] for c in ok])) if ok else None,
                        "fast_p_10k": float(np.mean([fr[c]["precision"] for c in ok])) if ok else None,
                        "fast_sp_10k": float(np.mean([fr[c]["recall"] for c in ok])) if ok else None,
                        "cpu_10k": float(sum(float(r[c].get("cpu_s") or 0) for c in k10_sel))}
        f.log(f"{tag} 10,000 (reported): " + "; ".join(f"{x[:8]} done {v['n_done']}/{v['n']}" for x, v in rep10.items()))

    out = {"fold": fold, "fit": name, "rule": f"PLAN_ensemble.md section 6 + {AMENDMENT4}",
           "ref_core": EF.brief_s(ref), "ref_screen": EF.brief_s(f.ref["screen"]),
           "rescore_info": info, "n_candidates": len(cands), "n_eligible": k["n_eligible"],
           "pick": brief(k["pick"], f"cheapest of {len(k['within'])} within 1.96 SE of the best"),
           "b": brief(k["b"], "eligible point with the highest P_1k"),
           "knee_rows": k["rows"], "eligible": k["eligible"],
           "reported_10k": rep10,
           "design_fast_m1": {"key": kd, "params": design, "eligible": eligible_a4({"core": core_s[kd]}, ref),
                              "core": EF.brief_s(core_s[kd])},
           "prerelease_accurate": {"key": ka, "core": EF.brief_s(core_s[ka])},
           "candidates": [{"key": c["key"], "source": c["source"], "rescore_set": c["why"],
                           "eligible": c["key"] in elig, "screen": c["screen"], "core": EF.brief_s(c["core"])}
                          for c in sorted(cands, key=lambda c: c["key"])]}
    for dd in EF.SENS_DELTAS:
        out[f"sensitivity_d{dd}"] = brief(k[f"pick_d{dd}"], f"sensitivity: cheapest within {dd} of the best")
    out["written"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    _atomic_json(d / "selection_a4.json", out)
    if out["pick"]:
        pk = out["pick"]["core"]
        f.log(f"{tag} pick {out['pick']['key'][:8]} ({out['pick']['source']}): P_small={pk['p_small']:.4f} "
              f"SP_small={pk['sp_small']:.4f} P_1k={pk['p_1k']:.4f} SP_1k={pk['sp_1k']:.4f} R={pk['R']:.2f} "
              f"| {out['pick']['short']}")
    else:
        f.log(f"{tag} no eligible point among {len(cands)}: no accurate setting from this fit")
    f.log(f"{tag} selection done")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True)
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    EF.use_ws_tmp(a.ws)
    select(Path(a.ws).expanduser(), a.fit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
