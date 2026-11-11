"""Pooled out-of-fold results of the nested CV (the numbers of the paper) for one pool.

    python -m recipe.pool_results --pool protein [--out WS/runs/refit] [--retime TSV] [--externals TSV ...]
    -> <out>/<pool>/results/
         oof_summary.tsv     per mode x arm x subset: n, errors, mean and 95% BCa CI of SP, TC, precision, F1,
                             CPU (s per case; also cpu4 = re-timed at 4 workers when --retime is given)
         oof_bw.tsv          benchmark-weighted means (per-case means within benchmark, benchmarks weighted
                             equally = the fit objective) with stratified BCa CIs, overall and per homolog flag
         oof_diff.tsv        paired differences new - 3.7.0 per mode (case-weighted and benchmark-weighted),
                             overall, per benchmark, per homolog_in_train, per difficulty bin, with BCa CIs
         per_fold.tsv        per fold x mode x arm benchmark-weighted means (descriptive)
         selected.tsv        per fold x mode: selected key, 3.7.0 flags, budget, training-core scores
         optimism.tsv        final-fit in-sample minus pooled OOF (final/insample_rows.tsv needed), per mode
         results.json        all of the above + provenance (read by recipe.decide)

Arms: "new" = the preset selected inside the fold (3.7.0 if no candidate was within budget: fallback flag),
"v370" = the 3.7.0 preset of the mode (parameter dict, same worker path), externals = MAFFT L-INS-i, MUSCLE 5,
Clustal Omega on the same cases (comparison only; rows of recipe.externals --out).
A failed alignment (error or timeout) scores 0 on SP, TC, precision and F1 (counted in n_err); its CPU is kept.
CIs: BCa, 10,000 resamples, seed 1234 (the same resample stream for every cell); benchmark-weighted
statistics resample cases within benchmark (stratified bootstrap); jackknife acceleration over cases.
SP = recall of the scorer (sum-of-pairs on the reference's scored columns).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist

import numpy as np

from . import fit as FIT

B_RESAMPLES = 10_000
CI_SEED = 1234
ALPHA = 0.05
METRICS = ("sp", "tc", "precision", "f1", "cpu")
COL = {"sp": "recall", "tc": "tc", "precision": "precision", "f1": "f1", "cpu": "cpu_s"}
N_FOLDS = 5
_ND = NormalDist()


# ------------------------------------------------------------------------------------ bootstrap
def bca(X, strata=None, B=B_RESAMPLES, seed=CI_SEED, alpha=ALPHA):
    """BCa CIs of a (stratum-weighted) mean for every column of X (n x m).

    statistic = mean over strata of the per-stratum column means (strata=None: the plain mean).
    Resampling is within stratum (stratified bootstrap). -> (est, lo, hi), arrays of length m."""
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    n, m = X.shape
    if n == 0:
        nan = np.full(m, np.nan)
        return nan, nan, nan
    labels = np.zeros(n, dtype=int) if strata is None else np.unique(np.asarray(strata), return_inverse=True)[1]
    groups = [X[labels == s] for s in range(labels.max() + 1)]
    K = len(groups)
    means = np.array([g.mean(axis=0) for g in groups])          # K x m
    est = means.mean(axis=0)
    rng = np.random.default_rng(seed)
    theta = np.zeros((B, m))
    for g in groups:
        ng = len(g)
        step = max(1, 4_000_000 // (ng * m))
        for s0 in range(0, B, step):
            s1 = min(B, s0 + step)
            idx = rng.integers(0, ng, size=(s1 - s0, ng))
            theta[s0:s1] += g[idx].mean(axis=1)
    theta /= K
    # jackknife (leave one case out; a stratum of one case is skipped)
    jk = []
    tot = means.sum(axis=0)
    for gi, g in enumerate(groups):
        ng = len(g)
        if ng < 2:
            continue
        loo = (g.sum(axis=0) - g) / (ng - 1)                     # ng x m
        jk.append((tot - means[gi] + loo) / K)
    jk = np.vstack(jk) if jk else np.zeros((0, m))
    lo, hi = np.empty(m), np.empty(m)
    zl, zh = _ND.inv_cdf(alpha / 2), _ND.inv_cdf(1 - alpha / 2)
    for j in range(m):
        t = theta[:, j]
        if np.all(t == t[0]):
            lo[j] = hi[j] = est[j]
            continue
        prop = (np.sum(t < est[j]) + 0.5 * np.sum(t == est[j])) / B
        prop = min(max(prop, 1.0 / (2 * B)), 1 - 1.0 / (2 * B))
        z0 = _ND.inv_cdf(prop)
        a = 0.0
        if len(jk):
            d = jk[:, j].mean() - jk[:, j]
            ss = np.sum(d ** 2)
            a = float(np.sum(d ** 3) / (6 * ss ** 1.5)) if ss > 0 else 0.0
        q = []
        for z in (zl, zh):
            den = 1 - a * (z0 + z)
            q.append(_ND.cdf(z0 + (z0 + z) / den) if den > 0 else (0.0 if z < 0 else 1.0))
        lo[j], hi[j] = np.quantile(t, q[0]), np.quantile(t, q[1])
    return est, lo, hi


# ------------------------------------------------------------------------------------ inputs
def _f(v):
    return float(v) if v not in (None, "") else None


def load_rows(path):
    """oof/insample rows -> list of dicts with floats; failed rows score 0 (quality) and keep CPU."""
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            err = r.get("error") or ""
            d = dict(r)
            d["error"] = err
            for k in METRICS:
                v = _f(r.get(COL[k]))
                d[k] = 0.0 if (v is None or (err and k != "cpu")) else v
            d["homolog_in_train"] = int(r["homolog_in_train"])
            out.append(d)
    return out


def load_externals(paths, pool, meta):
    """recipe.externals --out TSVs -> rows in the same layout (arm = tool, mode '-')."""
    out = []
    for p in paths:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["pool"] != pool or r["case_id"] not in meta:
                    continue
                err = r.get("error") or ""
                mm = meta[r["case_id"]]
                d = {"pool": pool, "fold": mm["fold"], "mode": "-", "arm": r["tool"], "case_id": r["case_id"],
                     "benchmark": mm["benchmark"], "difficulty": mm["difficulty"],
                     "homolog_in_train": mm["homolog_in_train"], "error": err, "key": "", "is_370": "0",
                     "fallback_370": "0"}
                for k in METRICS:
                    v = _f(r.get(COL[k]))
                    d[k] = 0.0 if (v is None or (err and k != "cpu")) else v
                out.append(d)
    return out


def attach_retime(rows, path):
    """Adds 'cpu4' (re-timed at 4 workers, recipe.decide retime) by (key, case_id)."""
    t = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("cpu_s") not in (None, ""):
                t[(r["key"], r["case_id"])] = float(r["cpu_s"])
    n = 0
    for r in rows:
        v = t.get((r.get("key"), r["case_id"]))
        r["cpu4"] = v
        n += v is not None
    return n


# ------------------------------------------------------------------------------------ tables
def subsets(rows):
    """[(subset_type, subset, rows)]"""
    out = [("all", "all", rows)]
    for key in ("benchmark", "homolog_in_train", "difficulty"):
        vals = sorted({str(r[key]) for r in rows})
        for v in vals:
            out.append((key, v, [r for r in rows if str(r[key]) == v]))
    return out


def metric_list(rows):
    return list(METRICS) + (["cpu4"] if rows and all(r.get("cpu4") is not None for r in rows) else [])


def summary_cells(rows, bw=False):
    ms = metric_list(rows)
    X = np.array([[r[m] for m in ms] for r in rows], dtype=float)
    est, lo, hi = bca(X, [r["benchmark"] for r in rows] if bw else None)
    return {m: (float(est[i]), float(lo[i]), float(hi[i])) for i, m in enumerate(ms)}


def paired(rows_new, rows_ref):
    """Join by case id -> list of (case row of new, difference dict)."""
    ref = {r["case_id"]: r for r in rows_ref}
    if set(ref) != {r["case_id"] for r in rows_new}:
        raise SystemExit("paired rows cover different cases")
    out = []
    for r in rows_new:
        q = ref[r["case_id"]]
        d = {k: r[k] - q[k] for k in metric_list([r, q])}
        out.append({**{k: r[k] for k in ("case_id", "benchmark", "homolog_in_train", "difficulty")}, **d})
    return out


def fmt(x):
    return "" if x is None else (f"{x:.6f}" if isinstance(x, float) else str(x))


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([fmt(r.get(c)) for c in cols])


def by_arm(rows):
    g = defaultdict(list)
    for r in rows:
        g[(r["mode"], r["arm"])].append(r)
    return g


def pool_results(pool, out, retime=None, externals=(), allow_partial=False, allow_smoke=False, log=FIT.log):
    out = Path(out)
    root = out / pool
    res_dir = root / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    rows, folds_found, smoke = [], [], False
    for k in range(N_FOLDS):
        p = root / f"fold{k}" / "oof_rows.tsv"
        if p.exists():
            rr = load_rows(p)
            rows += rr
            folds_found.append(k)
            smoke |= any(r.get("smoke") == "1" for r in rr)
    if len(folds_found) < N_FOLDS and not allow_partial:
        raise SystemExit(f"{pool}: OOF rows of folds {folds_found} only; --allow-partial for smoke tests")
    if smoke and not allow_smoke:
        raise SystemExit(f"{pool}: OOF rows come from smoke fits; refusing without --allow-smoke")
    meta = {}
    for r in rows:
        meta.setdefault(r["case_id"], {"fold": r["fold"], "benchmark": r["benchmark"], "difficulty": r["difficulty"],
                                       "homolog_in_train": r["homolog_in_train"]})
    # every case of the pool exactly once per arm
    for (m, arm), rr in by_arm(rows).items():
        ids = [r["case_id"] for r in rr]
        if len(set(ids)) != len(ids):
            raise SystemExit(f"duplicate case rows for {m}/{arm}")
    ext = load_externals(externals, pool, meta) if externals else []
    all_rows = rows + ext
    n_retimed = attach_retime(all_rows, retime) if retime else 0
    if retime:
        log(f"{pool}: {n_retimed} rows with 4-worker CPU")
    flags = {"partial": len(folds_found) < N_FOLDS, "smoke": smoke}
    tag = ("PARTIAL " if flags["partial"] else "") + ("SMOKE " if smoke else "")
    groups = by_arm(all_rows)
    # ---------------- summaries
    summ_rows, bw_rows = [], []
    for (m, arm), rr in sorted(groups.items()):
        for st, sv, sub in subsets(rr):
            cells = summary_cells(sub)
            d = {"pool": pool, "mode": m, "arm": arm, "subset_type": st, "subset": sv, "n": len(sub),
                 "n_err": sum(1 for r in sub if r["error"]), "flags": tag.strip()}
            for met, (e, lo, hi) in cells.items():
                d.update({f"{met}": e, f"{met}_lo": lo, f"{met}_hi": hi})
            d["cpu_total"] = float(sum(r["cpu"] for r in sub))
            summ_rows.append(d)
            if st in ("all", "homolog_in_train"):
                cells = summary_cells(sub, bw=True)
                d2 = {"pool": pool, "mode": m, "arm": arm, "subset_type": st, "subset": sv, "n": len(sub),
                      "n_benchmarks": len({r["benchmark"] for r in sub}), "flags": tag.strip()}
                for met, (e, lo, hi) in cells.items():
                    d2.update({f"{met}": e, f"{met}_lo": lo, f"{met}_hi": hi})
                bw_rows.append(d2)
    mcols = []
    for met in list(METRICS) + ["cpu4"]:
        mcols += [met, f"{met}_lo", f"{met}_hi"]
    write_tsv(res_dir / "oof_summary.tsv", ["pool", "mode", "arm", "subset_type", "subset", "n", "n_err", *mcols,
                                            "cpu_total", "flags"], summ_rows)
    write_tsv(res_dir / "oof_bw.tsv", ["pool", "mode", "arm", "subset_type", "subset", "n", "n_benchmarks", *mcols,
                                       "flags"], bw_rows)
    # ---------------- paired differences new - 3.7.0
    diff_rows, decision_inputs = [], {}
    for m in FIT.MODES[pool]:
        new, ref = groups.get((m, "new")), groups.get((m, "v370"))
        if not new or not ref:
            continue
        pr = paired(new, ref)
        for st, sv, sub in subsets(pr):
            for weighting in ("case", "benchmark"):
                if weighting == "benchmark" and st == "benchmark":
                    continue
                ms = metric_list(sub)
                X = np.array([[r[k] for k in ms] for r in sub], dtype=float)
                est, lo, hi = bca(X, [r["benchmark"] for r in sub] if weighting == "benchmark" else None)
                d = {"pool": pool, "mode": m, "weighting": weighting, "subset_type": st, "subset": sv,
                     "n": len(sub), "flags": tag.strip()}
                for i, k in enumerate(ms):
                    d.update({f"d_{k}": float(est[i]), f"d_{k}_lo": float(lo[i]), f"d_{k}_hi": float(hi[i])})
                diff_rows.append(d)
                if st == "all":
                    decision_inputs.setdefault(m, {})[weighting] = {k: d[k] for k in d if k.startswith("d_")}
        decision_inputs[m]["n"] = len(pr)
        decision_inputs[m]["n_err_new"] = sum(1 for r in new if r["error"])
        decision_inputs[m]["n_err_370"] = sum(1 for r in ref if r["error"])
        decision_inputs[m]["folds_fallback"] = sorted({int(r["fold"]) for r in new if r["fallback_370"] == "1"})
        decision_inputs[m]["folds_is_370"] = sorted({int(r["fold"]) for r in new if r["is_370"] == "1"})
    dcols = []
    for met in list(METRICS) + ["cpu4"]:
        dcols += [f"d_{met}", f"d_{met}_lo", f"d_{met}_hi"]
    write_tsv(res_dir / "oof_diff.tsv", ["pool", "mode", "weighting", "subset_type", "subset", "n", *dcols, "flags"],
              diff_rows)
    # ---------------- per fold (descriptive, benchmark-weighted point estimates)
    pf_rows = []
    for (m, arm), rr in sorted(groups.items()):
        for k in sorted({r["fold"] for r in rr}, key=str):
            sub = [r for r in rr if r["fold"] == k]
            bmeans = defaultdict(list)
            for r in sub:
                bmeans[r["benchmark"]].append(r)
            d = {"pool": pool, "fold": k, "mode": m, "arm": arm, "n": len(sub)}
            for met in ("sp", "tc", "precision", "f1"):
                d[met] = float(np.mean([np.mean([r[met] for r in v]) for v in bmeans.values()]))
            d["cpu_total"] = float(sum(r["cpu"] for r in sub))
            pf_rows.append(d)
    write_tsv(res_dir / "per_fold.tsv", ["pool", "fold", "mode", "arm", "n", "sp", "tc", "precision", "f1",
                                         "cpu_total"], pf_rows)
    # ---------------- selections per fold
    sel_rows = []
    for k in folds_found:
        sel = json.loads((root / f"fold{k}" / "selection.json").read_text())
        for m in FIT.MODES[pool]:
            s = sel["selection"][m]
            c = s.get("core") or {}
            sel_rows.append({"pool": pool, "fold": k, "mode": m, "key": s.get("key"), "is_370": int(s["is_370"]),
                             "fallback_370": int(s["fallback_370"]), "source": s.get("source"),
                             "budget_cpu": s.get("budget_cpu"), "core_sp": c.get("sp"), "core_tc": c.get("tc"),
                             "core_precision": c.get("precision"), "core_cpu": c.get("cpu"),
                             "ref_core_cpu": sel["refs_core"][FIT.BUDGET_REF[pool][m]]["cpu"],
                             "short": s.get("short", "")})
    write_tsv(res_dir / "selected.tsv", ["pool", "fold", "mode", "key", "is_370", "fallback_370", "source",
                                         "budget_cpu", "core_sp", "core_tc", "core_precision", "core_cpu",
                                         "ref_core_cpu", "short"], sel_rows)
    # ---------------- optimism: final-fit selection in-sample minus pooled OOF
    opt_rows = []
    ins_path = root / "final" / "insample_rows.tsv"
    if ins_path.exists():
        ins = by_arm(load_rows(ins_path))
        for m in FIT.MODES[pool]:
            a_in, a_oof = ins.get((m, "new")), groups.get((m, "new"))
            if not a_in or not a_oof:
                continue
            oof_ids = {r["case_id"] for r in a_oof}
            in_ids = {r["case_id"] for r in a_in}
            if not (flags["partial"] or smoke) and in_ids != oof_ids:
                raise SystemExit(f"{pool} {m}: in-sample and OOF rows cover different cases")
            common = oof_ids & in_ids      # = everything except in partial / smoke runs
            pr = paired([r for r in a_in if r["case_id"] in common], [r for r in a_oof if r["case_id"] in common])
            for weighting in ("case", "benchmark"):
                ms = [k for k in metric_list(pr) if k != "cpu4"]
                X = np.array([[r[k] for k in ms] for r in pr], dtype=float)
                est, lo, hi = bca(X, [r["benchmark"] for r in pr] if weighting == "benchmark" else None)
                d = {"pool": pool, "mode": m, "weighting": weighting, "n": len(pr),
                     "final_key": a_in[0]["key"], "flags": tag.strip()}
                for i, k in enumerate(ms):
                    d.update({f"opt_{k}": float(est[i]), f"opt_{k}_lo": float(lo[i]), f"opt_{k}_hi": float(hi[i])})
                opt_rows.append(d)
        ocols = []
        for met in METRICS:
            ocols += [f"opt_{met}", f"opt_{met}_lo", f"opt_{met}_hi"]
        write_tsv(res_dir / "optimism.tsv", ["pool", "mode", "weighting", "n", "final_key", *ocols, "flags"], opt_rows)
    result = {"pool": pool, "generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "folds": folds_found, **flags,
              "n_cases": len(meta), "bootstrap": {"B": B_RESAMPLES, "seed": CI_SEED, "method": "BCa",
                                                  "alpha": ALPHA},
              "retime": str(retime) if retime else None, "externals": [str(p) for p in externals],
              "decision_inputs": decision_inputs, "summary": summ_rows, "benchmark_weighted": bw_rows,
              "diff": diff_rows, "selected": sel_rows, "optimism": opt_rows}
    (res_dir / "results.json").write_text(json.dumps(result, indent=1) + "\n")
    log(f"{pool}: pooled {len(rows)} OOF rows over folds {folds_found} -> {res_dir}")
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--ws", default=str(FIT.default_ws()))
    ap.add_argument("--out")
    ap.add_argument("--retime", help="retime rows TSV of recipe.decide retime (adds cpu4)")
    ap.add_argument("--externals", nargs="*", default=[])
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--allow-smoke", action="store_true")
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    out = Path(a.out) if a.out else ws / "runs" / "refit"
    pool_results(a.pool, out, a.retime, a.externals, a.allow_partial, a.allow_smoke)
    return 0


if __name__ == "__main__":
    sys.exit(main())
