"""Knee-rule selection of DECISION_fastfit_selection.md (replaces the fast rule of PLAN_fastfit.md section 6).

  delta   per fit: the eligible set E (core, step 1), the best b, and the paired family-bootstrap SE of
          Q_large(b) - Q_large(c) for every c in E; prints the noise level and how many points are within reach.
          Does not print which point the rule picks.
  picks   per fit: the knee pick (the rule as decided), the sensitivity picks at fixed delta 0.005 and 0.01, and the
          registered picks (selection.json), written to <fit>/knee_selection.json; for the final fit also
          runs/fastfit/final_pick.json with the SHA-256 of the pick's canonical parameter JSON (written before any
          lockbox case is scored).

Run on the workstation from ~/tmp/kalign-fastfit/presets with the fit venv:
  python -m recipe.fastfit_knee delta|picks [fit ...]
"""
import argparse
import hashlib
import json
import sqlite3
import time
from pathlib import Path

import numpy as np

from . import cache as C
from .fastfit import load_fold

N_BOOT = 2000
SEED = 1234
Z = 1.96
FITS = ["fold0", "fold1", "fold2", "fold3", "fold4", "final", "final_seed4321"]


def fold_of(fit):
    return "all" if fit.startswith("final") else int(fit[-1])


def eligible(cands, ref):
    """Step 1 on the core data: feasible, small SP and small TC >= 4.0.0 fast."""
    return sorted((c for c in cands if c["core"]["feasible"] and c["core"]["sp_small"] >= ref["sp_small"]
                   and c["core"]["tc_small"] >= ref["tc_small"]), key=lambda c: c["key"])


def large_q(con, build, key, cases):
    """{case_id: (SP+TC)/2} for one candidate on the given extHomFam cases (from the fit's cache)."""
    out = {}
    for i in range(0, len(cases), 500):
        part = cases[i:i + 500]
        q = (f"SELECT case_id, recall, tc, error FROM results WHERE build_id=? AND param_key=? "
             f"AND case_id IN ({','.join('?' * len(part))})")
        for cid, rec, tc, err in con.execute(q, [build, key, *part]):
            out[cid] = 0.0 if err else 0.5 * (float(rec) + float(tc))
    missing = set(cases) - set(out)
    if missing:
        raise SystemExit(f"{key[:12]}: {len(missing)} core cases missing from the cache, e.g. {sorted(missing)[:3]}")
    return out


def boot_q_large(qv, fam_idx, k_is_10k, rng_idx):
    """Q_large per resample: mean over k of the per-k means of the resampled families' cases."""
    out = np.empty(len(rng_idx))
    for r, fams in enumerate(rng_idx):
        sel = np.concatenate([fam_idx[f] for f in fams])
        k1 = sel[~k_is_10k[sel]]
        k10 = sel[k_is_10k[sel]]
        out[r] = 0.5 * (qv[k1].mean() + qv[k10].mean())
    return out


def delta_fit(ws, fit):
    d = ws / "runs" / "fastfit" / fit
    fc = json.load(open(d / "front_core.json"))
    sel = json.load(open(d / "selection.json"))
    ref = sel["ref_core"]
    build = json.load(open(d / "run.json"))["build_new"]
    _, meta, core, _, _, _ = load_fold(ws, fold_of(fit))
    cases = sorted(c for c in core if meta[c]["kind"] == "large")
    fams = sorted({meta[c]["family"] for c in cases})
    fpos = {f: i for i, f in enumerate(fams)}
    fam_idx = [np.array([i for i, c in enumerate(cases) if fpos[meta[c]["family"]] == j]) for j in range(len(fams))]
    k_is_10k = np.array([meta[c]["k"] == 10000 for c in cases])
    rng = np.random.default_rng(SEED)
    rng_idx = rng.integers(0, len(fams), size=(N_BOOT, len(fams)))

    E = eligible(fc["rescored"], ref)
    if not E:
        return {"fit": fit, "n_rescored": len(fc["rescored"]), "n_eligible": 0}
    con = sqlite3.connect(f"file:{ws / 'cache' / 'cache.sqlite'}?mode=ro", uri=True)
    b = max(E, key=lambda c: (c["core"]["q_large"], -c["core"]["R"]))
    qb = large_q(con, build, b["key"], cases)
    vb = np.array([qb[c] for c in cases])
    full_b = 0.5 * (vb[~k_is_10k].mean() + vb[k_is_10k].mean())
    if abs(full_b - b["core"]["q_large"]) > 1e-9:
        raise SystemExit(f"{fit}: recomputed Q_large {full_b:.6f} != core summary {b['core']['q_large']:.6f}")
    bb = boot_q_large(vb, fam_idx, k_is_10k, rng_idx)
    rows = []
    for c in E:
        if c["key"] == b["key"]:
            continue
        qc = large_q(con, build, c["key"], cases)
        vc = np.array([qc[x] for x in cases])
        diff = b["core"]["q_large"] - c["core"]["q_large"]
        se = float(np.std(bb - boot_q_large(vc, fam_idx, k_is_10k, rng_idx), ddof=1))
        rows.append((diff, se, diff <= Z * se, c))
    se = np.array([r[1] for r in rows]) if rows else np.array([np.nan])
    return {"fit": fit, "n_rescored": len(fc["rescored"]), "n_eligible": len(E), "n_families": len(fams),
            "n_cases": len(cases), "se_min": float(np.nanmin(se)), "se_median": float(np.nanmedian(se)),
            "se_max": float(np.nanmax(se)), "delta_median": float(Z * np.nanmedian(se)),
            "n_within": 1 + sum(r[2] for r in rows),
            "_E": E, "_b": b, "_rows": rows, "_sel": sel}


def cheapest(cs):
    """Step 4: the lowest R; ties: higher Q_large, then key."""
    return min(cs, key=lambda c: (c["core"]["R"], -c["core"]["q_large"], c["key"]))


def brief(c, why):
    if c is None:
        return None
    return {"key": c["key"], "why": why, "params": c["params"], "source": c.get("source", f"search gen {c.get('gen')}"),
            "core": c["core"]}


def picks_fit(ws, fit):
    r = delta_fit(ws, fit)
    sel = r.get("_sel") or json.load(open(ws / "runs" / "fastfit" / fit / "selection.json"))
    gm1 = next(c for c in sel["candidates"] if c["source"] == "gm1 fast preset")
    out = {"fit": fit, "rule": "DECISION_fastfit_selection.md", "n_eligible": r["n_eligible"],
           "registered_fast": sel["fast"], "registered_fast_is_fallback": sel["fast_is_fallback"],
           "registered_ultra": sel["ultra"], "gm1": gm1}
    if r["n_eligible"] == 0:
        out.update({"knee": brief(gm1, "fallback: no eligible point"), "knee_d0.005": None, "knee_d0.01": None,
                    "n_within": 0})
        return out
    b, rows = r["_b"], r["_rows"]
    within = [b] + [c for _, _, ok, c in rows if ok]
    out["n_within"] = len(within)
    out["knee"] = brief(cheapest(within), f"cheapest of {len(within)} within 1.96 SE of the best")
    for d in (0.005, 0.01):
        w = [b] + [c for diff, _, _, c in rows if diff <= d]
        out[f"knee_d{d}"] = brief(cheapest(w), f"sensitivity: cheapest of {len(w)} within {d} of the best")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["delta", "picks"])
    ap.add_argument("fits", nargs="*", default=FITS)
    ap.add_argument("--ws", default=str(Path.home() / "tmp" / "kalign-fastfit"))
    a = ap.parse_args(argv)
    ws = Path(a.ws)
    if a.cmd == "picks":
        return picks(ws, a.fits)
    out = []
    for f in a.fits:
        r = delta_fit(ws, f)
        out.append({k: v for k, v in r.items() if not k.startswith("_")})
        if r["n_eligible"] == 0:
            print(f"{f:15s} rescored {r['n_rescored']}, eligible 0")
            continue
        print(f"{f:15s} rescored {r['n_rescored']:2d}, eligible {r['n_eligible']:2d}, families {r['n_families']}, "
              f"SE of the paired difference {r['se_min']:.4f}-{r['se_max']:.4f} (median {r['se_median']:.4f}), "
              f"1.96 x median SE = {r['delta_median']:.4f}, within reach of the best: {r['n_within']}")
    (ws / "runs" / "fastfit" / "knee_delta.json").write_text(json.dumps(out, indent=1))


def picks(ws, fits):
    for f in fits:
        o = picks_fit(ws, f)
        o["written"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        (ws / "runs" / "fastfit" / f / "knee_selection.json").write_text(json.dumps(o, indent=1) + "\n")
        k = o["knee"]
        print(f"{f:15s} eligible {o['n_eligible']:2d}, within reach {o['n_within']:2d}: knee R {k['core']['R']:.3f} "
              f"Q_small {k['core']['q_small']:.4f} Q_large {k['core']['q_large']:.4f} ({k['source']})")
        if f == "final":
            fp = ws / "runs" / "fastfit" / "final_pick.json"
            if fp.exists():
                old = json.loads(fp.read_text())
                if old["key"] != k["key"]:
                    raise SystemExit(f"{fp} exists with a different pick; refusing to overwrite")
                continue
            canon = C.canonical_json(k["params"])
            fp.write_text(json.dumps({"fit": "final", "rule": "DECISION_fastfit_selection.md (knee)", "key": k["key"],
                                      "params": k["params"], "params_sha256": hashlib.sha256(canon.encode()).hexdigest(),
                                      "written": o["written"]}, indent=1) + "\n")
            print(f"  final pick written: {fp} sha256 {hashlib.sha256(canon.encode()).hexdigest()}")
    return 0


if __name__ == "__main__":
    main()
