"""Step 0 run: the four 3.7.0 presets once on every case of both pools, then cost/accuracy tables.

    python -m recipe.step0 run     [--workers 30] [--modes fast default accurate recall] [--timeout-s 1800]
    python -m recipe.step0 summary [--out-dir WS/runs/step0]
    python -m recipe.step0 recheck [--mode fast] [--n 20]      # cache check: must be 100% hits, 0 alignments

All results go through recipe.evaluate and the cache. A case whose aligner runs longer than --timeout-s
(handled inside recipe.evaluate: wall seconds of the child, wall >= CPU) is
killed and stored as error "timeout"; it is excluded from the metric means and listed separately.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

from . import cache as C
from . import evaluate as E
from .scoring import read_manifest

WS = Path(os.environ.get("KALIGN_REFIT_WS", Path.home() / "tmp" / "kalign-refit"))
BUILD_DIR = WS / "kalign" / "960a4f6"
OUT_DIR = WS / "runs" / "step0"
MANIFESTS = {"protein": WS / "data" / "cases_protein.tsv", "nucleotide": WS / "data" / "cases_nucleotide.tsv"}
MODES = ("fast", "default", "accurate", "recall")  # recall == accurate in this build (alias), run anyway
TIMEOUT_S = 1800.0


def all_cases():
    rows = []
    for pool, path in MANIFESTS.items():
        rows += read_manifest(path)
    # expensive first: better packing of the worker pool (cost ~ n_seq * ref_len)
    rows.sort(key=lambda r: -int(r["n_seq"]) * int(r["ref_len"]))
    return rows


def cmd_run(a):
    build = E.Build.load(BUILD_DIR)
    cases = all_cases()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"{len(cases)} cases, build {build.build_id[:12]}, workers {a.workers}, timeout {a.timeout_s}s", flush=True)
    for mode in a.modes:
        t0 = time.time()
        rows = E.evaluate(cases, mode, build, a.workers, progress=True, timeout_s=a.timeout_s)
        E.write_rows(rows, OUT_DIR / f"run_{mode}.tsv")
        n_err = sum(1 for r in rows if r.get("error"))
        print(f"[{time.strftime('%H:%M:%S')}] {mode}: {len(rows)} cases, "
              f"{sum(r['cached'] for r in rows)} cached, {n_err} errors, {time.time() - t0:.0f}s wall", flush=True)
    return 0


def cmd_recheck(a):
    build = E.Build.load(BUILD_DIR)
    cases = all_cases()[-a.n:]  # cheapest cases, but the point is only the cache behaviour
    t0 = time.time()
    rows = E.evaluate(cases, a.mode, build, workers=4)
    n_hit = sum(r["cached"] for r in rows)
    msg = f"recheck {a.mode}: {len(rows)} cases, {n_hit} cache hits, {len(rows) - n_hit} aligned, {time.time() - t0:.2f}s"
    print(msg)
    (OUT_DIR / "recheck.txt").write_text(msg + "\n")
    return 0 if n_hit == len(rows) else 1


# ---------------------------------------------------------------- summary
def _q(xs, p):
    xs = sorted(xs)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def _mean(xs):
    return statistics.fmean(xs) if xs else float("nan")


def load_results():
    """-> {mode: [row dict]} read from the cache for the manifest cases (single source of truth)."""
    build = E.Build.load(BUILD_DIR)
    cache = C.Cache()
    cases = all_cases()
    idents = [C.ident(c, E.scorer_version("recipe.scoring:score_case", c)) for c in cases]
    out = {}
    for mode in MODES:
        hits = cache.get_many(build.build_id, C.param_key(mode), idents)
        rows = []
        for c in cases:
            h = hits.get(c["case_id"])
            if h is None:
                continue
            rows.append({"pool": c["pool"], "benchmark": c["benchmark"], "case_id": c["case_id"],
                         "n_seq": int(c["n_seq"]), "ref_len": int(c["ref_len"]), **h})
        out[mode] = rows
    return out


def _fmt(x, nd=4):
    return "" if x is None or x != x else f"{x:.{nd}f}"


def summarise(results):
    tab, share, tails = [], [], []
    cols = ["pool", "benchmark", "mode", "n_cases", "n_scored", "n_errors", "n_timeouts", "SP", "precision", "F1", "TC",
            "cpu_total_s", "cpu_median_s", "cpu_p90_s", "cpu_max_s"]
    for mode, rows in results.items():
        groups = defaultdict(list)
        for r in rows:
            groups[(r["pool"], r["benchmark"])].append(r)
            groups[(r["pool"], "ALL")].append(r)
        for (pool, bm), rs in sorted(groups.items()):
            ok = [r for r in rs if not r["error"]]
            cpu = [r["cpu_s"] for r in rs if r["cpu_s"] is not None]
            tab.append({"pool": pool, "benchmark": bm, "mode": mode, "n_cases": len(rs), "n_scored": len(ok),
                        "n_errors": len(rs) - len(ok),
                        "n_timeouts": sum(1 for r in rs if (r["error"] or "").startswith("timeout")),
                        "SP": _fmt(_mean([r["recall"] for r in ok])), "precision": _fmt(_mean([r["precision"] for r in ok])),
                        "F1": _fmt(_mean([r["f1"] for r in ok])), "TC": _fmt(_mean([r["tc"] for r in ok])),
                        "cpu_total_s": _fmt(sum(cpu), 1), "cpu_median_s": _fmt(_q(cpu, .5), 3),
                        "cpu_p90_s": _fmt(_q(cpu, .9), 3), "cpu_max_s": _fmt(max(cpu) if cpu else float("nan"), 2)})
        for pool in sorted({r["pool"] for r in rows}):
            rs = [r for r in rows if r["pool"] == pool and r["cpu_s"] is not None]
            cpu = sorted((r["cpu_s"] for r in rs), reverse=True)
            tot = sum(cpu)
            k = max(1, round(0.05 * len(cpu)))
            share.append({"pool": pool, "mode": mode, "n_cases": len(cpu), "top5pct_n": k,
                          "cpu_total_s": _fmt(tot, 1), "top5pct_cpu_s": _fmt(sum(cpu[:k]), 1),
                          "top5pct_share": _fmt(sum(cpu[:k]) / tot if tot else float("nan"), 4)})
            for r in sorted(rs, key=lambda r: -r["cpu_s"])[:10]:
                tails.append({"pool": pool, "mode": mode, "case_id": r["case_id"], "n_seq": r["n_seq"],
                              "ref_len": r["ref_len"], "cpu_s": _fmt(r["cpu_s"], 1), "error": r["error"] or ""})
    return cols, tab, share, tails


def _write(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def cmd_summary(a):
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results = load_results()
    cols, tab, share, tails = summarise(results)
    _write(out / "summary_by_benchmark.tsv", cols, tab)
    _write(out / "top5pct_cpu_share.tsv", list(share[0]), share)
    _write(out / "slowest_cases.tsv", list(tails[0]), tails)
    errs = [{"mode": m, "case_id": r["case_id"], "n_seq": r["n_seq"], "cpu_s": _fmt(r["cpu_s"], 1), "error": r["error"]}
            for m, rs in results.items() for r in rs if r["error"]]
    _write(out / "errors.tsv", ["mode", "case_id", "n_seq", "cpu_s", "error"], errs)
    # per-case wide table (one row per case, all modes)
    wide = defaultdict(dict)
    for m, rs in results.items():
        for r in rs:
            w = wide[r["case_id"]]
            w.update(pool=r["pool"], benchmark=r["benchmark"], case_id=r["case_id"], n_seq=r["n_seq"], ref_len=r["ref_len"])
            for k in ("recall", "precision", "f1", "tc", "cpu_s", "wall_s"):
                w[f"{m}_{k}"] = _fmt(r[k], 4)
            w[f"{m}_error"] = r["error"] or ""
    wcols = ["pool", "benchmark", "case_id", "n_seq", "ref_len"] + [f"{m}_{k}" for m in results
                                                                    for k in ("recall", "precision", "f1", "tc", "cpu_s", "wall_s", "error")]
    _write(out / "per_case.tsv", wcols, list(wide.values()))
    print(f"wrote tables to {out}; cases per mode: " + ", ".join(f"{m}={len(r)}" for m, r in results.items()))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=30)
    r.add_argument("--modes", nargs="+", default=list(MODES), choices=MODES)
    r.add_argument("--timeout-s", type=float, default=TIMEOUT_S)
    r.set_defaults(fn=cmd_run)
    s = sub.add_parser("summary")
    s.add_argument("--out-dir", default=str(OUT_DIR))
    s.set_defaults(fn=cmd_summary)
    c = sub.add_parser("recheck")
    c.add_argument("--mode", default="fast", choices=MODES)
    c.add_argument("--n", type=int, default=20)
    c.set_defaults(fn=cmd_recheck)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
