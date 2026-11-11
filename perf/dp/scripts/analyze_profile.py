#!/usr/bin/env python3
"""Summarise kalign_prof timer JSONs (from run_profile.py) into a breakdown
per mode and dataset, and fit how each part scales with N and L.

    uv run --no-project python perf/dp/scripts/analyze_profile.py --results $RES [--md out.md]

Categories (self times of the timed regions, so they add up to the total):
  tree_heur    guide tree without DP: anchors, k-mer/BPM distances, k-means/UPGMA
  tree_dp      DP alignments done to build the guide tree (guide_dist pairwise)
  cons_build   anchor-consistency table build (its pairwise DPs included)
  dp_ss/sp/pp  progressive + refine DP forward+backward per call type
               (with or without consistency bonus)
  dp_meet_tb   Hirschberg meetup + recursion overhead + path conversion
               ("traceback" of the linear-space DP)
  profile      profile construction and update (make_profile_n,
               set_gap_penalties_n, update_n)
  extras       VSM, sparse-bonus construction, SP scoring, realign distances,
               ensemble POAR/scoring/consensus
  other        I/O, setup, make_seq, finalise, aln_mem init
Fits: log t = a + bN log N + bL log L, least squares over the grid points of
one (mode, divergence); bN and bL are the log-log slopes.
"""
import argparse
import glob
import json
import math
import os
import sys
from collections import defaultdict

CATS = ["tree_heur", "tree_dp", "cons_build", "dp_ss", "dp_sp", "dp_pp", "dp_meet_tb",
        "profile", "extras", "other"]


def category(phase, region):
    if phase == "consbuild":
        return "cons_build"
    dp = region.startswith("dp_") or region.startswith("dpc_")
    if phase == "tree":
        if dp or region == "dp_runner" or region == "path_conv":
            return "tree_dp"
        if region in ("prof_make", "prof_gaps", "prof_update", "alnmem_init"):
            return "tree_dp"
        return "tree_heur"
    if dp and region != "dp_runner":
        if region.endswith("_meet"):
            return "dp_meet_tb"
        kind = region.split("_")[1]
        return "dp_" + kind
    if region in ("dp_runner", "path_conv"):
        return "dp_meet_tb"
    if region in ("prof_make", "prof_gaps", "prof_update"):
        return "profile"
    if region in ("vsm", "cons_bonus", "sp_score", "realign_dist", "ens_poar", "ens_score",
                  "ens_consensus", "ensemble_other"):
        return "extras"
    if region.startswith("tree") or region in ("guide_ctx",):
        return "tree_heur"
    return "other"


def load(results):
    rows = []
    for fn in sorted(glob.glob(os.path.join(results, "*.json"))):
        j = json.load(open(fn))
        if "timers" not in j or not j["timers"]:
            continue
        d = j["dataset"]
        tot = j["wall"]["total"]
        cats = defaultdict(float)
        dpcells = defaultdict(float)
        dpns = defaultdict(float)
        fb = defaultdict(float)
        for r in j["timers"]["records"]:
            if r["region"] == "driver":
                cats["other"] += r["self_ns"] / 1e9
                continue
            cats[category(r["phase"], r["region"])] += r["self_ns"] / 1e9
            reg = r["region"]
            if (reg.startswith("dp_") or reg.startswith("dpc_")) and reg != "dp_runner" and not reg.endswith("_meet"):
                kind = reg.split("_")[1]
                key = ("c" if reg.startswith("dpc") else "") + kind
                dpcells[key] += r["cells"]
                dpns[key] += r["self_ns"]
                fb[reg.split("_")[2]] += r["self_ns"] / 1e9
            if reg.endswith("_meet"):
                fb["meet"] += r["self_ns"] / 1e9
        rows.append(dict(name=d["name"], mode=j["mode"], N=d["n_seq"], L=d["mean_len"],
                         div=d["divergence"], bio=d["biotype"], total=tot, cats=dict(cats),
                         dpcells=dict(dpcells), dpns=dict(dpns), fb=dict(fb)))
    return rows


def lstsq(X, y):
    # normal equations, tiny systems
    n = len(X[0])
    A = [[sum(X[k][i] * X[k][j] for k in range(len(X))) for j in range(n)] for i in range(n)]
    b = [sum(X[k][i] * y[k] for k in range(len(X))) for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]
        b[i], b[p] = b[p], b[i]
        for r in range(n):
            if r != i and A[i][i] != 0:
                f = A[r][i] / A[i][i]
                A[r] = [A[r][c] - f * A[i][c] for c in range(n)]
                b[r] -= f * b[i]
    return [b[i] / A[i][i] if A[i][i] else float("nan") for i in range(n)]


def fit(points):
    pts = [(N, L, t) for N, L, t in points if t > 1e-4]
    if len(pts) < 4 or len({p[1] for p in pts}) < 2 or len({p[0] for p in pts}) < 2:
        return None
    X = [[1.0, math.log(N), math.log(L)] for N, L, _ in pts]
    y = [math.log(t) for _, _, t in pts]
    return lstsq(X, y)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--md", default=None, help="also write the tables as Markdown")
    a = ap.parse_args()
    rows = load(a.results)
    out = []
    modes = ["fast", "default", "accurate"]

    out.append("## Breakdown (% of wall time, single thread)\n")
    hdr = "| dataset | mode | total s | " + " | ".join(CATS) + " |"
    out.append(hdr)
    out.append("|" + "---|" * (3 + len(CATS)))
    key = lambda r: (modes.index(r["mode"]) if r["mode"] in modes else 9, r["bio"], r["div"], r["N"], r["L"])
    for r in sorted(rows, key=key):
        cells = [f"{100 * r['cats'].get(c, 0.0) / r['total']:.1f}" for c in CATS]
        out.append(f"| {r['name']} | {r['mode']} | {r['total']:.2f} | " + " | ".join(cells) + " |")

    out.append("\n## Aggregate over the protein grid (sum of times, % of summed wall)\n")
    out.append("| mode | total s | tree (heur+dp) | of which tree_dp | all DP (fwd+bwd+meet, any phase) | fwd | bwd | meet+tb | profile | extras |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for m in modes:
        rs = [r for r in rows if r["mode"] == m and r["bio"] == "protein"]
        if not rs:
            continue
        T = sum(r["total"] for r in rs)
        g = lambda c: sum(r["cats"].get(c, 0.0) for r in rs)
        alldp = sum(r["fb"].get("fwd", 0) + r["fb"].get("bwd", 0) + r["fb"].get("meet", 0) for r in rs)
        alldp += sum(r["cats"].get("dp_meet_tb", 0.0) for r in rs) - sum(r["fb"].get("meet", 0) for r in rs)
        fw = sum(r["fb"].get("fwd", 0) for r in rs)
        bw = sum(r["fb"].get("bwd", 0) for r in rs)
        out.append(f"| {m} | {T:.1f} | {100 * (g('tree_heur') + g('tree_dp')) / T:.1f} | {100 * g('tree_dp') / T:.1f} | "
                   f"{100 * alldp / T:.1f} | {100 * fw / T:.1f} | {100 * bw / T:.1f} | {100 * g('dp_meet_tb') / T:.1f} | "
                   f"{100 * g('profile') / T:.1f} | {100 * g('extras') / T:.1f} |")

    out.append("\n## ns per DP cell inside Kalign runs (fwd+bwd self time / cells; c = with consistency bonus)\n")
    out.append("| mode | ss | sp | pp | css | csp | cpp |")
    out.append("|---|---|---|---|---|---|---|")
    for m in modes:
        rs = [r for r in rows if r["mode"] == m]
        line = []
        for k in ["ss", "sp", "pp", "css", "csp", "cpp"]:
            c = sum(r["dpcells"].get(k, 0) for r in rs)
            ns = sum(r["dpns"].get(k, 0) for r in rs)
            line.append(f"{ns / c:.2f}" if c else "-")
        out.append(f"| {m} | " + " | ".join(line) + " |")

    out.append("\n## Scaling fits: log-log slopes of time vs N and L (protein)\n")
    out.append("| mode | div | part | slope N | slope L | points |")
    out.append("|---|---|---|---|---|---|")
    parts = ["total", "tree_heur", "tree_dp", "cons_build", "dp_all", "profile", "extras"]
    for m in modes:
        for dv in ["lo", "hi"]:
            rs = [r for r in rows if r["mode"] == m and r["div"] == dv and r["bio"] == "protein"]
            for p in parts:
                pts = []
                for r in rs:
                    if p == "total":
                        t = r["total"]
                    elif p == "dp_all":
                        t = sum(r["cats"].get(c, 0.0) for c in ("dp_ss", "dp_sp", "dp_pp", "dp_meet_tb"))
                    else:
                        t = r["cats"].get(p, 0.0)
                    pts.append((r["N"], r["L"], t))
                f = fit(pts)
                if f:
                    out.append(f"| {m} | {dv} | {p} | {f[1]:.2f} | {f[2]:.2f} | {len([q for q in pts if q[2] > 1e-4])} |")
    text = "\n".join(out)
    print(text)
    if a.md:
        open(a.md, "w").write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
