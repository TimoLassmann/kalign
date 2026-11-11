"""Step 2 of the nucfit plan's accurate setting (presets/docs/nucfit/PLAN_nucfit.md, section 5): confidence step.

Per fit (folds 0-4 and the final fit), the members of that fit's step-1 pick (runs/nucens/<fit>/selection.json) are
kept; every member takes member 1's guide tree; a grid over the single-linkage tree noise of members 2-5 (SIGMAS) and
the consensus support (SUPPORTS) is scored on the fit's training core (small core + Rfam large + simulations at 1,000;
no input above 5,000 sequences) with build 9137658 and the confidence scorer (recipe.scoring_conf: pair Brier score).
Selection (as protein, recipe.ensfit_conf): eligible = no failed case, SP and TC not below the fit's fast pick (small
and large, point estimates on the core), no large case with SP more than 0.2 below the fast pick; pick = lowest
J = (B_small + B_large) / 2; ties: smaller noise, then smaller support. Large values = mean of the category means.
No eligible point: the fit has no pick (counts as the fast setting out of fold).

  fit    [--fits 0,1,2,3,4,all]   -> runs/nucconf/<fit>/{core.json, selection.json}
  check  identity checks before the fit -> runs/nucconf/identity.json

    python -m recipe.nucconf fit|check [--ws WS]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from . import ensfit as EF
from . import ensfit_conf as EC
from . import evaluate as E
from . import nucens as NE
from . import nucfast as NF
from .fit import LeakError, SlotPool, _atomic_json, _py, log, pkey

SIGMAS = EC.SIGMAS
SUPPORTS = (1, 2, 3)
DROP = EC.DROP
FITS = ("0", "1", "2", "3", "4", "all")
MET = EC.MET


def out_dir(ws, *parts):
    d = Path(ws) / "runs" / "nucconf" / Path(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d


def fit_dir_name(fit):
    return "final" if str(fit) == "all" else f"fold{fit}"


def fold_arg(fit):
    return "all" if str(fit) == "all" else int(fit)


def fast_pick(ws, fit):
    """The fit's nucleotide fast pick exactly as selected (its own gap model): the comparison arm."""
    sel = json.loads((Path(ws) / "runs" / "nucfast" / fit_dir_name(fit) / "selection.json").read_text())
    p = sel["fast"]["params"]
    if pkey(p) != sel["fast"]["key"]:
        raise SystemExit(f"nucfast {fit}: pick params do not hash to their key")
    return sel["fast"]["key"], _py(p)


def ens_pick(ws, fit):
    sel = json.loads((Path(ws) / "runs" / "nucens" / fit_dir_name(fit) / "selection.json").read_text())
    if sel["pick"] is None:
        return None, None
    p = sel["pick"]["params"]
    if pkey(p) != sel["pick"]["key"]:
        raise SystemExit(f"nucens {fit}: pick params do not hash to their key")
    return sel["pick"]["key"], _py(p)


def make_evaluator(ws, rows, allowed, logfh=None):
    E.WORKER_SCRIPT = EC.WORKER
    EF.use_ws_tmp(ws)
    b = E.Build.load(Path(ws) / "kalign" / NE.BUILD)
    slots = SlotPool(Path(ws) / "runs" / "slots", EF.SLOTS)
    ev = EF.EnsEvaluator(b, rows.values(), allowed, cache_path=Path(ws) / "cache" / "nucconf.sqlite", slots=slots,
                         threads=EF.SLOTS, timeout_s=EF.TIMEOUT_S, scorer=EC.SCORER, logfh=logfh)
    return ev, b


def summarise(res, rows, meta, small, large, fast_res=None):
    errs = [c for c in small + large if res[c].get("error")]
    s = {"errors": len(errs), "error_cases": errs[:10],
         "cpu_small": sum(float(res[c].get("cpu_s") or 0) for c in small),
         "cpu_large": sum(float(res[c].get("cpu_s") or 0) for c in large)}
    if errs:
        return s
    bench = sorted({rows[c]["benchmark"] for c in small})
    cats = sorted({meta[c]["k"] for c in large})
    for m in MET:
        s[f"{m}_small"] = float(np.mean([np.mean([EC.val(res[c], m) for c in small if rows[c]["benchmark"] == b])
                                         for b in bench]))
        s[f"{m}_large"] = float(np.mean([np.mean([EC.val(res[c], m) for c in large if meta[c]["k"] == k])
                                         for k in cats]))
    if fast_res is not None:
        d = [EC.val(res[c], "recall") - EC.val(fast_res[c], "recall") for c in large]
        s["drops_large"] = {"below_-0.2": sum(x < -DROP for x in d), "below_-0.1": sum(x < -0.1 for x in d),
                            "above_+0.1": sum(x > 0.1 for x in d)}
    s["J"] = (s["brier_pair_small"] + s["brier_pair_large"]) / 2
    return s


def eligible(s, f):
    return (s["errors"] == 0 and s["recall_small"] >= f["recall_small"] and s["tc_small"] >= f["tc_small"]
            and s["recall_large"] >= f["recall_large"] and s["tc_large"] >= f["tc_large"]
            and s["drops_large"]["below_-0.2"] == 0)


def fit_one(ws, fit):
    od = out_dir(ws, fit_dir_name(fit))
    logfh = open(od / "progress.log", "a")
    rows, meta, core, _, _, lockbox = NE.load_fold(ws, fold_arg(fit))
    small = [c for c in core if meta[c]["kind"] == "small"]
    large = [c for c in core if meta[c]["kind"] == "large"]
    if set(core) & set(lockbox):
        raise LeakError(f"fit {fit}: lockbox case in the core")
    akey, pick = ens_pick(ws, fit)
    fkey, fast = fast_pick(ws, fit)
    if pick is None:
        out = {"fit": fit, "ensemble_pick": None, "pick": None, "n_eligible": 0,
               "note": "step 1 had no eligible point", "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
        _atomic_json(od / "selection.json", out)
        log(f"[nucconf {fit}] step 1 has no pick: no accurate setting from this fit", logfh)
        return out
    g = {(s, m): EC.config(pick, s, m) for s in SIGMAS for m in SUPPORTS}
    specs = {fkey: fast, **{pkey(p): p for p in g.values()}}
    log(f"[nucconf {fit}] {len(small)} small + {len(large)} large; step-1 pick {akey[:8]}, fast {fkey[:8]}; "
        f"{len(g)} grid points", logfh)
    ev, b = make_evaluator(ws, rows, small + large, logfh)
    res = ev.evaluate(specs, small + large)
    fs = summarise(res[fkey], rows, meta, small, large)
    if fs["errors"]:
        raise SystemExit(f"fit {fit}: the fast pick failed on {fs['error_cases']}")
    pts = []
    for (sig, sup), p in g.items():
        k = pkey(p)
        s = summarise(res[k], rows, meta, small, large, res[fkey])
        s.update({"sigma": sig, "min_support": sup, "key": k})
        s["eligible"] = bool(s["errors"] == 0 and eligible(s, fs))
        pts.append(s)
    el = sorted((s for s in pts if s["eligible"]), key=lambda s: (s["J"], s["sigma"], s["min_support"]))
    sel = el[0] if el else None
    _atomic_json(od / "core.json", {"fit": fit, "build": b.build_id, "ensemble_pick": akey, "fast_pick": fkey,
                                    "fast": fs, "points": pts, "n_small": len(small), "n_large": len(large)})
    out = {"fit": fit, "plan": "nucfit/PLAN_nucfit.md section 5 step 2", "build": b.build_id, "ensemble_pick": akey,
           "fast_pick": fkey, "n_eligible": len(el), "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "pick": None if sel is None else {
               "key": sel["key"], "sigma": sel["sigma"], "min_support": sel["min_support"],
               "params": g[(sel["sigma"], sel["min_support"])],
               "sha256": hashlib.sha256(json.dumps(g[(sel["sigma"], sel["min_support"])], sort_keys=True)
                                        .encode()).hexdigest(), "core": sel}}
    _atomic_json(od / "selection.json", out)
    log(f"[nucconf {fit}] eligible {len(el)}/{len(pts)}; pick "
        + ("none" if sel is None else f"sigma {sel['sigma']} support {sel['min_support']} J {sel['J']:.4f}"), logfh)
    for s in sorted(pts, key=lambda s: (s["sigma"], s["min_support"])):
        if s["errors"]:
            log(f"   sigma {s['sigma']:.2f} ms {s['min_support']}: {s['errors']} errors", logfh)
            continue
        log(f"   sigma {s['sigma']:.2f} ms {s['min_support']}: SP {s['recall_small']:.4f}/{s['recall_large']:.4f} "
            f"TC {s['tc_small']:.4f}/{s['tc_large']:.4f} B {s['brier_pair_small']:.4f}/{s['brier_pair_large']:.4f} "
            f"J {s['J']:.4f} drops {s['drops_large']['below_-0.2']} {'*' if s['eligible'] else ''}", logfh)
    log(f"   fast: SP {fs['recall_small']:.4f}/{fs['recall_large']:.4f} TC {fs['tc_small']:.4f}/{fs['tc_large']:.4f} "
        f"B {fs['brier_pair_small']:.4f}/{fs['brier_pair_large']:.4f}", logfh)
    return out


def check(ws, n_small=20, n_large=3):
    """Before any fit (fit data of fold 0 only, never the lockbox), with the pre-release nucleotide presets:
    C1 the confidence worker gives the same alignment sha256 and SP / precision / TC as the plain worker of the fast
       fit (cache.sqlite) for the pre-release fast setting;
    C2 sigma 0 with run_tree_sl_noise == the same dict without it (same alignment sha256), on a 5-member ensemble;
    C3 correct_pairs / ref_pairs_b == the pool scorer's recall (1e-12) on every scored row;
    C4 a single run is scored as hard (conf_source "hard", B = (missed + wrong) / reference pairs) and the ensemble
       carries confidence (conf_source "ensemble")."""
    rows, meta, core, *_ = NE.load_fold(ws, 0)
    small = sorted(c for c in core if meta[c]["kind"] == "small")[:: 37][:n_small]
    large = sorted(c for c in core if meta[c]["kind"] == "large" and int(rows[c]["n_seq"]) <= 600)[:n_large]
    ids = small + large
    from . import fastfit as FF
    fast = _py(NF.with_gap(FF.with_tree(NF.v2_params(NF.NUC_FAST_KEY), FF.TREE_400), 0))
    m1 = NE.fast_member(_py(NF.with_gap(fast, NE.GAP_MODEL)))
    acc = NF.v2_params(NF.NUC_ACC_KEY)
    tree = {k: m1[k] for k in ("tree_method", "tree_leaf", "tree_dist", "tree_alpha", "tree_kmeans_dist")}
    a = NE.single_member(acc, tree)
    space = NE.NucEnsSpace(_py(NF.with_gap(fast, NE.GAP_MODEL)))
    ens = space.wd([m1, a, a, a, a], {k: acc[k] for k in EF.ENS_KEYS} | {"min_support": 2})
    z = EC.config(ens, 0.0, 2)
    z_none = {k: v for k, v in z.items() if k != "run_tree_sl_noise"}
    ev, b = make_evaluator(ws, rows, ids)
    res = ev.evaluate({pkey(p): p for p in (fast, z, z_none)}, ids)
    from .evaluate import Build, evaluate
    plain = {r["case_id"]: r for r in evaluate([rows[c] for c in ids], fast, Build.load(Path(ws) / "kalign" / NE.BUILD),
                                               workers=20, cache_path=Path(ws) / "cache" / "cache.sqlite",
                                               timeout_s=EF.TIMEOUT_S)}
    out = {"C1": [], "C2": [], "C3": [], "C4": []}
    for c in ids:
        n, o = res[pkey(fast)][c], plain[c]
        out["C1"].append({"case": c, "same": o.get("aln_sha256") == n["aln_sha256"]
                          and all(float(o[m]) == float(n[m]) for m in ("recall", "precision", "tc"))})
        out["C2"].append({"case": c, "same": res[pkey(z)][c]["aln_sha256"] == res[pkey(z_none)][c]["aln_sha256"]})
        for p in (fast, z):
            r = res[pkey(p)][c]
            x = r["extra"]
            out["C3"].append({"case": c, "ok": abs(x["correct_pairs"] / x["ref_pairs_b"] - float(r["recall"])) < 1e-12})
        xf, xz = res[pkey(fast)][c]["extra"], res[pkey(z)][c]["extra"]
        hard = (xf["ref_pairs_b"] - xf["correct_pairs"] + xf["wrong_pairs"]) / xf["ref_pairs_b"]
        out["C4"].append({"case": c, "ok": xf["conf_source"] == "hard" and xz["conf_source"] == "ensemble"
                          and abs(xf["brier_pair"] - hard) < 1e-12})
    summ = {k: f"{sum(r.get('same', r.get('ok')) for r in v)}/{len(v)}" for k, v in out.items()}
    _atomic_json(out_dir(ws) / "identity.json", {"summary": summ, **out, "build": b.build_id, "cases": ids})
    print(json.dumps(summ, indent=1))
    return 0 if all(r.get("same", r.get("ok")) for v in out.values() for r in v) else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["fit", "check"])
    ap.add_argument("--fits", default=",".join(FITS))
    ap.add_argument("--ws", default=str(NF.default_ws()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    EF.use_ws_tmp(ws)
    if a.cmd == "check":
        return check(ws)
    for f in a.fits.split(","):
        fit_one(ws, f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
