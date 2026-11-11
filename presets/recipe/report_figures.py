"""Report figures (andamentum report look) of a finished re-fit run, plus Addendum A1 (Kalign 3.6.1 / 3.7.0 CLI).

Every chart is drawn from the run's own output files with the andamentum figure style (setup_style("report"),
the four andamentum tones, text kept as text, inline SVG). It runs under the andamentum tool's Python:

    ~/.local/share/uv/tools/andamentum/bin/python -m recipe.report_figures \
        --run RUN --compare V1_RUN --a1 A1_DIR --out OUT

RUN / V1_RUN: run directories with refit/<pool>/<fold>/{stage1,stage2}/hv.json, selection.json and report/tables/.
A1_DIR: Addendum A1 directory with <pool>_kalign{361,370}_<mode>.tsv (recipe.externals rows) and
refit/<pool>/results/oof_bw.tsv (recipe.pool_results with those rows as externals).
Writes <name>.svg (inline-ready), <name>.png (for checking) and <name>.csv (the plotted numbers) per figure.
"""

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402

from andamentum.figures.palettes import ANDAMENTUM_TONES as T  # noqa: E402
from andamentum.figures.style import inline_svg, savefig, setup_style  # noqa: E402

PRESET = setup_style("report")
INK, RULE, ERR = "#1a1a1a", "#d4d0c8", "#777777"
W = 7.5
POOLS = ("protein", "nucleotide")
MODES = ("fast", "default", "recall", "accurate")
MARK = {"fast": "o", "default": "s", "recall": "^", "accurate": "D"}
ARM_C = {"3.6.1": T["grey-light"], "3.7.0": T["slate"], "re-fit": T["rust"], "external": T["ochre"]}
EXT_NAME = {"muscle5": "MUSCLE 5", "mafft_linsi": "MAFFT L-INS-i", "clustalo": "Clustal Omega"}
SHIPS = {"protein": ("default", "recall", "accurate"), "nucleotide": ("recall",)}
_CI = re.compile(r"^\s*([+-]?[\d.]+)\s*\[\s*([+-]?[\d.]+),\s*([+-]?[\d.]+)\s*\]\s*$")


def ci(s):
    return tuple(float(x) for x in _CI.match(str(s)).groups())


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def table(run, name):
    return tsv(Path(run) / "report" / "tables" / f"{name}.tsv")


def fitdir(run, pool, fold):
    return Path(run) / "refit" / pool / ("final" if fold == "final" else f"fold{fold}")


def write_csv(out, name, rows):
    if rows:
        with open(out / f"{name}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


def finish(fig, out, name, draw):
    """draw(fig) builds the figure; called twice (SVG for the page, PNG for checking)."""
    for ext in ("png", "svg"):
        f = plt.figure(figsize=fig)
        draw(f)
        savefig(f, out / f"{name}.{ext}", **({"dpi": 200, "facecolor": "white"} if ext == "png" else {}))
    (out / f"{name}.svg").write_text(inline_svg(out / f"{name}.svg"))


def log_x(ax):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(subs=(1, 3)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())


def log_y(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(subs=(1, 3)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())


# ------------------------------------------------------------------------------------------------ data
def a1_bw(a1, pool):
    """Benchmark-weighted out-of-fold rows of every arm (Addendum A1 pooling): {(arm, mode): row}."""
    out = {}
    for r in tsv(Path(a1) / "refit" / pool / "results" / "oof_bw.tsv"):
        if r["subset_type"] != "all":
            continue
        arm, mode = r["arm"], r["mode"]
        if arm.startswith("kalign3"):
            arm, mode = arm.split("_")
        out[(arm, mode)] = {k: float(r[k]) for k in ("sp", "sp_lo", "sp_hi", "tc", "tc_lo", "tc_hi", "precision",
                                                     "cpu")}
    return out


def cli_cpu(a1, pool):
    """Case-mean single-thread CPU of the Kalign CLI arms (same harness, 8 workers): {(build, mode): (mean, n)}."""
    res = {}
    for b in ("kalign361", "kalign370"):
        for m in MODES:
            rows = tsv(Path(a1) / f"{pool}_{b}_{m}.tsv")
            res[(b, m)] = (sum(float(r["cpu_s"]) for r in rows) / len(rows), len(rows),
                           sum(bool(r["error"]) for r in rows))
    return res


# ------------------------------------------------------------------------------------------------ figures
def fig_tradeoff(a1, out, pool):
    """Out-of-fold SP and TC against CPU: Kalign 3.6.1, 3.7.0, re-fit, external aligners."""
    d = a1_bw(a1, pool)
    pts = []
    for m in MODES:
        pts.append(("3.6.1", m, d[("kalign361", m)]))
        if pool == "protein" and m == "recall":
            continue  # protein recall = accurate in 3.7.0 and in the re-fit; 3.6.1 has a separate recall ensemble
        pts.append(("3.7.0", m, d[("v370", m)]))
        pts.append(("re-fit", m, d[("new", m)]))
    for t, name in EXT_NAME.items():
        pts.append(("external", name, d[(t, "-")]))
    write_csv(out, f"fig_tradeoff_{pool}", [{"arm": a, "mode": m, **v} for a, m, v in pts])

    def draw(f):
        axes = f.subplots(1, 2)
        for ax, met, lab in zip(axes, ("sp", "tc"), ("SP", "TC")):
            for arm, m, v in pts:
                c = ARM_C[arm]
                mk = "*" if arm == "external" else MARK[m]
                filled = arm != "re-fit" or m in SHIPS[pool]
                ax.errorbar(v["cpu"], v[met], yerr=[[v[met] - v[f"{met}_lo"]], [v[f"{met}_hi"] - v[met]]],
                            fmt=mk, ms=11 if mk == "*" else 6.5, color=c, mfc=c if filled else "none", mec=c,
                            ecolor=ERR, elinewidth=0.8, capsize=0, zorder=3)
                if arm == "external":
                    ax.annotate(m, (v["cpu"], v[met]), xytext=(-6, 6), textcoords="offset points", ha="right",
                                fontsize=PRESET.body_font_size, color=INK)
            log_x(ax)
            ax.set_xlabel("CPU per case (s, one thread)")
            ax.set_ylabel(f"{lab}, out of fold")
        hand = [Line2D([], [], marker="o", ls="", color=ARM_C[a], label=lab)
                for a, lab in (("3.6.1", "Kalign 3.6.1"), ("3.7.0", "Kalign 3.7.0"), ("re-fit", "re-fit"))]
        hand += [Line2D([], [], marker=MARK[m], ls="", color=INK, mfc="none", label=m) for m in MODES]
        f.legend(handles=hand, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.04),
                 handletextpad=0.2, columnspacing=1.0)
        f.tight_layout(rect=(0, 0, 1, 0.94))
    finish((W, 3.6), out, f"fig_tradeoff_{pool}", draw)


def fig_speed(a1, out):
    """CPU per case, Kalign 3.6.1 vs 3.7.0 CLI, same harness, by mode and pool."""
    rows = []
    for pool in POOLS:
        c = cli_cpu(a1, pool)
        for m in MODES:
            a, b = c[("kalign361", m)][0], c[("kalign370", m)][0]
            rows.append({"pool": pool, "mode": m, "cpu_361": a, "cpu_370": b, "speedup": a / b,
                         "n": c[("kalign370", m)][1], "errors_361": c[("kalign361", m)][2],
                         "errors_370": c[("kalign370", m)][2]})
    write_csv(out, "fig_speed", rows)

    def draw(f):
        axes = f.subplots(1, 2, sharey=False)
        for ax, pool in zip(axes, POOLS):
            rr = [r for r in rows if r["pool"] == pool]
            x = range(len(rr))
            w = 0.38
            ax.bar([i - w / 2 for i in x], [r["cpu_361"] for r in rr], w, color=ARM_C["3.6.1"], label="Kalign 3.6.1")
            ax.bar([i + w / 2 for i in x], [r["cpu_370"] for r in rr], w, color=ARM_C["3.7.0"], label="Kalign 3.7.0")
            log_y(ax)
            top = max(r["cpu_361"] for r in rr)
            for i, r in zip(x, rr):
                ax.text(i, r["cpu_361"] * 1.15, f"{r['speedup']:.1f}×", ha="center", va="bottom",
                        fontsize=PRESET.body_font_size, color=INK)
            ax.set_ylim(top=top * 3)
            ax.set_xticks(list(x))
            ax.set_xticklabels([r["mode"] for r in rr])
            ax.set_ylabel("CPU per case (s, one thread)")
            ax.set_title(f"{pool}, n = {rr[0]['n']}", loc="left", fontsize=PRESET.label_font_size)
        axes[0].legend(loc="upper left", frameon=False)
        f.tight_layout()
    finish((W, 3.4), out, "fig_speed", draw)


def fig_pareto(run, compare, out, pool):
    sel = json.load(open(fitdir(run, pool, "final") / "selection.json"))
    old = json.load(open(fitdir(compare, pool, "final") / "selection.json"))
    front = sel["front"]
    write_csv(out, f"fig_pareto_{pool}", [{k: p[k] for k in ("source", "sp", "tc", "precision", "cpu")}
                                          for p in front])

    def draw(f):
        axes = f.subplots(1, 2)
        for ax, met, lab in zip(axes, ("tc", "precision"), ("TC", "precision")):
            for src, c in (("stage1", T["grey-light"]), ("stage2", T["slate"])):
                p = [q for q in front if q["source"] == src]
                ax.scatter([q["cpu"] for q in p], [q[met] for q in p], s=7, color=c, lw=0, zorder=2)
            seen = set()
            for m in ("fast", "default", "accurate", "recall"):
                ref = sel["refs_core"][m]
                if ref["key"] not in seen:
                    seen.add(ref["key"])
                    ax.scatter(ref["cpu"], ref[met], marker=MARK[m], s=60, facecolor="none", edgecolor=INK,
                               lw=1.2, zorder=4)
                o = old["selection"][m]["core"]
                ax.scatter(o["cpu"], o[met], marker=MARK[m], s=40, facecolor="none", edgecolor=T["ochre"], lw=1.2,
                           zorder=4)
                s = sel["selection"][m]
                ax.scatter(s["core"]["cpu"], s["core"][met], marker=MARK[m], s=60, color=T["rust"], lw=0, zorder=5)
            log_x(ax)
            ax.set_xlabel("CPU on the training core (s, one thread)")
            ax.set_ylabel(f"{lab}, training core")
        hand = [Line2D([], [], marker="o", ls="", ms=4, color=T["grey-light"], label="front, single run"),
                Line2D([], [], marker="o", ls="", ms=4, color=T["slate"], label="front, ensemble"),
                Line2D([], [], marker="o", ls="", mfc="none", mec=INK, label="3.7.0"),
                Line2D([], [], marker="o", ls="", mfc="none", mec=T["ochre"], label="v1 pick"),
                Line2D([], [], marker="o", ls="", color=T["rust"], label="v1.1 pick")]
        f.legend(handles=hand, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.04),
                 handletextpad=0.2, columnspacing=1.0)
        f.tight_layout(rect=(0, 0, 1, 0.93))
    finish((W, 3.6), out, f"fig_pareto_{pool}", draw)


def fig_convergence(run, out):
    rows = []
    for pool in POOLS:
        for fold in ("0", "1", "2", "3", "4", "final"):
            h = json.load(open(fitdir(run, pool, fold) / "stage2" / "hv.json"))["history"]
            rows += [{"pool": pool, "fit": fold, "generation": g + 1, "hypervolume": v} for g, v in enumerate(h)]
    write_csv(out, "fig_convergence", rows)

    def draw(f):
        axes = f.subplots(1, 2)
        for ax, pool in zip(axes, POOLS):
            for fold in ("0", "1", "2", "3", "4", "final"):
                rr = [r for r in rows if r["pool"] == pool and r["fit"] == fold]
                fin = fold == "final"
                ax.plot([r["generation"] for r in rr], [r["hypervolume"] for r in rr],
                        color=T["rust"] if fin else T["grey-light"], lw=1.8 if fin else 1.0, zorder=3 if fin else 2)
            ax.axvline(40, color=T["ochre"], lw=1.0, ls="--", zorder=1)
            ax.text(42, ax.get_ylim()[0], "v1 limit", color=INK, va="bottom", fontsize=PRESET.body_font_size)
            ax.set_xlim(0, 140)
            ax.set_xlabel("generation (stage 2)")
            ax.set_ylabel("hypervolume")
            ax.set_title(pool, loc="left", fontsize=PRESET.label_font_size)
        hand = [Line2D([], [], color=T["grey-light"], label="fold fits 0 to 4"),
                Line2D([], [], color=T["rust"], lw=1.8, label="final fit")]
        axes[1].legend(handles=hand, loc="lower right", frameon=False)
        f.tight_layout()
    finish((W, 3.2), out, "fig_convergence", draw)


def _dots(f, panels, xlabel, margin=True):
    """panels: [(title, [(label, sp(v,lo,hi), tc(v,lo,hi))])] -> horizontal dot-and-interval plot."""
    axes = f.subplots(1, len(panels), sharex=True, sharey=True)
    axes = axes if len(panels) > 1 else [axes]
    for ax, (title, rows) in zip(axes, panels):
        y = list(range(len(rows)))[::-1]
        for yi, (lab, sp, tc) in zip(y, rows):
            for off, (v, lo, hi), c, mk in ((0.16, sp, T["slate"], "o"), (-0.16, tc, T["rust"], "s")):
                ax.errorbar(v, yi + off, xerr=[[v - lo], [hi - v]], fmt=mk, ms=5, color=c, ecolor=ERR,
                            elinewidth=0.8, capsize=0)
        ax.set_yticks(y)
        ax.set_yticklabels([r[0] for r in rows])
        ax.axvline(0, color=INK, lw=0.8)
        if margin:
            for mm in (-0.005, 0.005):
                ax.axvline(mm, color=RULE, lw=1.0, ls="--")
        ax.set_xlabel(xlabel)
        ax.set_title(title, loc="left", fontsize=PRESET.label_font_size)
    hand = [Line2D([], [], marker="o", ls="", color=T["slate"], label="SP"),
            Line2D([], [], marker="s", ls="", color=T["rust"], label="TC")]
    f.legend(handles=hand, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    f.tight_layout(rect=(0, 0, 1, 0.95))


def fig_breakdown(run, out):
    order = ["<20", "20-30", "30-40", ">=40"]
    panels, csvrows = [], []
    for mode in ("default", "accurate"):
        dec = {r["mode"]: r for r in table(run, "protein_decision")}[mode]
        rows = [("all, benchmark-weighted", ci(dec["dSP (bench-weighted)"]), ci(dec["dTC (bench-weighted)"])),
                ("all, case-weighted", ci(dec["dSP (case)"]), ci(dec["dTC (case)"]))]
        names = {"bali3": "BAliBASE 3", "balifam": "BaliFam100", "ox": "OXBench", "prefab4": "PREFAB 4",
                 "sabre": "SABRE"}
        for r in table(run, "protein_diff_by_benchmark"):
            if r["mode"] == mode:
                rows.append((f"{names[r['benchmark']]} (n = {r['n']})", ci(r["dSP"]), ci(r["dTC"])))
        dd = {r["difficulty"]: r for r in table(run, "protein_diff_by_difficulty") if r["mode"] == mode}
        for k in order:
            r = dd[k]
            rows.append((f"identity {k.replace('>=', '≥ ').replace('<', '< ')}% (n = {r['n']})", ci(r["dSP"]),
                         ci(r["dTC"])))
        panels.append((f"protein {mode}", rows))
        csvrows += [{"mode": mode, "subset": lab, "dSP": s[0], "dSP_lo": s[1], "dSP_hi": s[2], "dTC": t[0],
                     "dTC_lo": t[1], "dTC_hi": t[2]} for lab, s, t in rows]
    write_csv(out, "fig_breakdown", csvrows)

    def draw(f):
        _dots(f, panels, "re-fit minus 3.7.0, out of fold")
    finish((W, 4.6), out, "fig_breakdown", draw)


def fig_optimism(run, compare, out):
    panels, csvrows = [], []
    for pool in POOLS:
        rows = []
        for r in table(run, f"{pool}_optimism"):
            if r["weighting"] != "benchmark" or (pool == "protein" and r["mode"] == "recall"):
                continue
            lab = r["mode"] + (" = recall" if pool == "protein" and r["mode"] == "accurate" else "")
            rows.append((lab, ci(r["dSP"]), ci(r["dTC"])))
            csvrows.append({"pool": pool, "mode": r["mode"], "dSP": r["dSP"], "dTC": r["dTC"]})
        panels.append((pool, rows))
    write_csv(out, "fig_optimism", csvrows)

    def draw(f):
        _dots(f, panels, "in-sample minus out of fold", margin=False)
    finish((W, 2.9), out, "fig_optimism", draw)


def fig_v1_v11(run, compare, out):
    rows, csvrows = [], []
    for pool in POOLS:
        a = {r["mode"]: r for r in table(compare, f"{pool}_decision")}
        b = {r["mode"]: r for r in table(run, f"{pool}_decision")}
        for m in MODES:
            if pool == "protein" and m == "recall":
                continue
            lab = f"{pool} {m}" + (" = recall" if pool == "protein" and m == "accurate" else "")
            rows.append((lab, {k: ci(a[m][f"d{k} (bench-weighted)"]) for k in ("SP", "TC")},
                         {k: ci(b[m][f"d{k} (bench-weighted)"]) for k in ("SP", "TC")}))
            csvrows.append({"row": lab, **{f"v1_d{k}": rows[-1][1][k][0] for k in ("SP", "TC")},
                            **{f"v11_d{k}": rows[-1][2][k][0] for k in ("SP", "TC")}})
    write_csv(out, "fig_v1_v11", csvrows)

    def draw(f):
        axes = f.subplots(1, 2, sharey=True)
        for ax, k in zip(axes, ("SP", "TC")):
            y = list(range(len(rows)))[::-1]
            for yi, (lab, va, vb) in zip(y, rows):
                for off, v, c in ((0.16, va[k], T["ochre"]), (-0.16, vb[k], T["rust"])):
                    ax.errorbar(v[0], yi + off, xerr=[[v[0] - v[1]], [v[2] - v[0]]], fmt="o", ms=5, color=c,
                                ecolor=ERR, elinewidth=0.8, capsize=0)
            ax.set_yticks(y)
            ax.set_yticklabels([r[0] for r in rows])
            ax.axvline(0, color=INK, lw=0.8)
            for mm in (-0.005, 0.005):
                ax.axvline(mm, color=RULE, lw=1.0, ls="--")
            ax.set_xlabel(f"{k}: re-fit minus 3.7.0, out of fold")
        hand = [Line2D([], [], marker="o", ls="", color=T["ochre"], label="v1, stage 2 stopped at 40 generations"),
                Line2D([], [], marker="o", ls="", color=T["rust"], label="v1.1, stage 2 run to convergence")]
        f.legend(handles=hand, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.03))
        f.tight_layout(rect=(0, 0, 1, 0.93))
    finish((W, 3.4), out, "fig_v1_v11", draw)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--compare", required=True, type=Path)
    ap.add_argument("--a1", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    for pool in POOLS:
        fig_tradeoff(a.a1, a.out, pool)
        fig_pareto(a.run, a.compare, a.out, pool)
    fig_speed(a.a1, a.out)
    fig_convergence(a.run, a.out)
    fig_breakdown(a.run, a.out)
    fig_optimism(a.run, a.compare, a.out)
    fig_v1_v11(a.run, a.compare, a.out)
    for pool in POOLS:
        fig_mode_choice(a.run, a.a1, a.out, pool)
    print(f"figures written to {a.out}")


if __name__ == "__main__":
    main()


# ------------------------------------------------------------------------------------------------ mode choice
CHOICE_TIERS = ("fast", "default", "accurate")
PRECISE_SLACK = 0.01   # precise pick: max precision among points within 0.01 (SP+TC)/2 of the SP+TC pick


def core_points(run, a1, pool="protein"):
    """Final-fit front, candidate picks and every other program on the final-fit training core (benchmark-weighted,
    exactly the fit's objective; CPU per case). Returns (front, picks, others, budgets)."""
    fd = fitdir(run, pool, "final")
    core = set((fd / "cases_core.txt").read_text().split())
    bench = {r["case_id"]: r["benchmark"] for r in tsv(fd / "insample_rows.tsv")}
    n = len(core)

    def bw(rows):
        by = {}
        for r in rows:
            by.setdefault(bench[r["case_id"]], []).append(r)
        val = lambda r, k: 0.0 if r.get("error") else float(r[k] or 0)  # noqa: E731
        m = lambda k: sum(sum(val(r, k) for r in v) / len(v) for v in by.values()) / len(by)  # noqa: E731
        assert len(rows) == n, (len(rows), n)
        return {"sp": m("recall"), "tc": m("tc"), "precision": m("precision"),
                "cpu": sum(float(r["cpu_s"]) for r in rows) / n}

    sel = json.load(open(fd / "selection.json"))
    others = {}
    ext = tsv(Path(run) / f"externals_{pool}.tsv")
    for t, name in EXT_NAME.items():
        others[name] = bw([r for r in ext if r["tool"] == t and r["case_id"] in core])
    for m in MODES:
        rows = [r for r in tsv(Path(a1) / f"{pool}_kalign361_{m}.tsv") if r["case_id"] in core]
        others[f"3.6.1 {m}"] = bw(rows)
    for m in CHOICE_TIERS:
        r = sel["refs_core"][m]
        others[f"3.7.0 {m}"] = {"sp": r["sp"], "tc": r["tc"], "precision": r["precision"], "cpu": r["cpu"] / n}
    for v in others.values():
        v["sptc"] = (v["sp"] + v["tc"]) / 2
    front = [dict(p, cpu=p["cpu"] / n, sptc=(p["sp"] + p["tc"]) / 2) for p in sel["front"]]
    fastest_other = min(v["cpu"] for k, v in others.items() if k in EXT_NAME.values())
    picks, limits = knee_modes(front, fastest_other)
    return front, picks, others, limits


KNEE_TOL = 0.002        # do not pay CPU for less than this gain in (SP+TC)/2
DEFAULT_SPLIT = 2.0     # default costs at most half of accurate
FAST_VS_OTHERS = 5.0    # fast is at least 5x faster than the fastest other aligner


def knee_modes(front, fastest_other):
    """Three modes from the (SP+TC)/2 Pareto curve, top down:
    accurate = cheapest point within KNEE_TOL of the best point;
    default  = cheapest point within KNEE_TOL of the best point costing <= CPU(accurate) / DEFAULT_SPLIT;
    fast     = cheapest point within KNEE_TOL of the best point costing <= fastest other aligner / FAST_VS_OTHERS.
    Returns ({("sptc", mode): point}, {mode: CPU limit used})."""
    def pick(limit):
        c = [p for p in front if p["cpu"] <= limit]
        top = max(p["sptc"] for p in c)
        return min((p for p in c if p["sptc"] >= top - KNEE_TOL), key=lambda p: (p["cpu"], -p["sptc"]))
    acc = pick(float("inf"))
    lim = {"accurate": max(p["cpu"] for p in front), "default": acc["cpu"] / DEFAULT_SPLIT,
           "fast": fastest_other / FAST_VS_OTHERS}
    picks = {("sptc", "accurate"): acc, ("sptc", "default"): pick(lim["default"]), ("sptc", "fast"): pick(lim["fast"])}
    return picks, lim


def envelope(points, met):
    """2D Pareto curve: points not beaten on (met up, CPU down), sorted by CPU."""
    pts = sorted(points, key=lambda p: (p["cpu"], -p[met]))
    out, best = [], -1.0
    for p in pts:
        if p[met] > best:
            out.append(p)
            best = p[met]
    return out


def fig_mode_choice(run, a1, out, pool="protein"):
    front, picks, others, budgets = core_points(run, a1, pool)
    rows = [{"point": f"front {p['source']}", **{k: p[k] for k in ("sp", "tc", "sptc", "precision", "cpu")}}
            for p in front]
    rows += [{"point": f"pick {a} {m}", **{k: p[k] for k in ("sp", "tc", "sptc", "precision", "cpu")}}
             for (a, m), p in picks.items()]
    rows += [{"point": k, **{kk: v[kk] for kk in ("sp", "tc", "sptc", "precision", "cpu")}} for k, v in others.items()]
    write_csv(out, f"fig_mode_choice_{pool}", rows)
    fs = PRESET.body_font_size
    ref = list(others.values()) + [p for k, p in picks.items()]
    ylim = {m: (min(v[m] for v in ref) - 0.03, max(max(v[m] for v in ref), max(p[m] for p in front)) + 0.008)
            for m in ("sptc", "precision")}
    lab_off = {"MUSCLE 5": (-8, 6), "MAFFT L-INS-i": (9, -4), "Clustal Omega": (8, -16)}
    lab_off_p = {"MUSCLE 5": (-8, 6), "MAFFT L-INS-i": (8, -16), "Clustal Omega": (8, -16)}

    def draw(f):
        axes = f.subplots(1, 2)
        for ax, met, lab in zip(axes, ("sptc", "precision"), ("(SP + TC) / 2", "precision")):
            lo, hi = ylim[met]
            ax.scatter([p["cpu"] for p in front], [p[met] for p in front], s=6, color=T["grey-light"], lw=0, zorder=2)
            env = envelope(front, met)
            ax.step([p["cpu"] for p in env], [p[met] for p in env], where="post", color=T["grey"], lw=1.0, zorder=3)
            for m in ("fast", "default"):
                ax.axvline(budgets[m], color=RULE, lw=1.0, ls="--", zorder=1)
                if met == "sptc":
                    ax.text(budgets[m], lo + 0.002, f"{m} limit ", rotation=90, ha="right", va="bottom",
                            fontsize=fs * 0.9, color=INK, zorder=6,
                            bbox=dict(boxstyle="square,pad=0.1", fc="#f9f7f4", ec="none"))
            for name, v in others.items():
                ext = name in EXT_NAME.values()
                mk = "*" if ext else ("s" if name.startswith("3.7.0") else "o")
                ax.scatter(v["cpu"], v[met], marker=mk, s=110 if ext else 40, facecolor=T["slate"] if ext else "none",
                           edgecolor=T["slate"], lw=1.2, zorder=4, clip_on=True)
                if ext:
                    off = (lab_off if met == "sptc" else lab_off_p)[name]
                    ax.annotate(name, (v["cpu"], v[met]), xytext=off, textcoords="offset points",
                                ha="right" if off[0] < 0 else "left", fontsize=fs, color=INK)
            for (a, m), p in picks.items():
                if a != "sptc":
                    continue
                ax.scatter(p["cpu"], p[met], marker="o", s=60, color=T["rust"], lw=0, zorder=5)
                if a == "sptc" and met == "sptc":
                    ax.annotate(m, (p["cpu"], p[met]), xytext=(0, 10), textcoords="offset points", ha="center",
                                fontsize=fs, color=INK)
            ax.set_ylim(lo, hi)
            log_x(ax)
            ax.set_xlabel("CPU per case (s, one thread)")
            ax.set_ylabel(f"{lab}, training core")
        hand = [Line2D([], [], color=T["grey"], lw=1.0, marker="o", ms=3, mfc=T["grey-light"], mec=T["grey-light"],
                       label="Pareto front and curve"),
                Line2D([], [], marker="o", ls="", color=T["rust"], label="chosen mode (knee)"),
                Line2D([], [], marker="s", ls="", mfc="none", mec=T["slate"], label="Kalign 3.7.0"),
                Line2D([], [], marker="o", ls="", mfc="none", mec=T["slate"], label="Kalign 3.6.1"),
                Line2D([], [], marker="*", ls="", ms=10, color=T["slate"], label="other programs")]
        f.legend(handles=hand, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.07),
                 handletextpad=0.3, columnspacing=1.2)
        f.tight_layout(rect=(0, 0, 1, 0.88))
    finish((W, 4.2), out, f"fig_mode_choice_{pool}", draw)
    return {"below_axis": sum(p["sptc"] < ylim["sptc"][0] for p in front), "front": len(front)}


# ------------------------------------------------------------------------------------------------ BAliBASE categories
BALI_CATS = ("RV11", "RV12", "RV20", "RV30", "RV40", "RV50")
CAT_ARMS = (("Kalign fast", T["grey-light"]), ("Kalign accurate", T["rust"]), ("Clustal Omega", T["ochre"]),
            ("MAFFT L-INS-i", T["slate"]), ("MUSCLE 5", T["grey"]))
OLD_KALIGN_NAMES = {"kalign31": "Kalign 3.1", "kalign335": "Kalign 3.3.5"}
OLD_ARMS = (("Kalign 3.1", T["grey"]), ("Kalign 3.3.5", T["grey"]))   # drawn as outlined bars


def _boot_ci(vals, n_boot=2000, seed=1234):
    import random
    rng = random.Random(seed)
    n = len(vals)
    means = sorted(sum(vals[rng.randrange(n)] for _ in range(n)) / n for _ in range(n_boot))
    return means[int(0.025 * n_boot)], means[int(0.975 * n_boot) - 1]


def bali_categories(v2_refit, externals_tsv, scope="oof", old_tsvs=()):
    """Per-case SP, TC and CPU on BAliBASE 3 by reference category (BB and BBS together) for the v2 modes and the
    external aligners. scope "oof": each case scored by the presets of the fold that held it out (folds 0-4);
    scope "insample": the final presets on every case (optimistic: the final fit trained on part of them)."""
    ext = [r for r in tsv(externals_tsv) if r["benchmark"] == "bali3"]
    cat = {r["case_id"]: r["subset"] for r in ext}                    # e.g. BB-RV11, BBS-RV50

    def val(r):
        ok = not r["error"]
        return (float(r["recall"]) if ok else 0.0, float(r["tc"]) if ok else 0.0, float(r["cpu_s"]))
    per = {(EXT_NAME[r["tool"]], r["case_id"]): val(r) for r in ext}
    arms = list(CAT_ARMS)
    for f in old_tsvs:
        for r in tsv(f):
            if r["benchmark"] == "bali3":
                per[(OLD_KALIGN_NAMES[r["tool"]], r["case_id"])] = val(r)
                if (OLD_KALIGN_NAMES[r["tool"]], T["grey"]) not in arms:
                    arms.insert(sum(a[0] in ("Kalign 3.1", "Kalign 3.3.5") for a in arms),
                                (OLD_KALIGN_NAMES[r["tool"]], T["grey"]))
    files = ([Path(v2_refit) / "protein" / f"fold{k}" / "oof_rows.tsv" for k in range(5)] if scope == "oof"
             else [Path(v2_refit) / "protein" / "final" / "insample_rows.tsv"])
    for f in files:
        for r in tsv(f):
            if r["benchmark"] == "bali3" and r["arm"] == "new" and r["mode"] in ("fast", "accurate"):
                per[(f"Kalign {r['mode']}", r["case_id"])] = val(r)
    rows = []
    for c in BALI_CATS:
        ids = sorted(i for i, s_ in cat.items() if s_.endswith(c))
        for arm, _ in arms:
            vals = [per[(arm, i)] for i in ids]
            d = {"category": c, "arm": arm, "n": len(ids), "n_bb": sum(cat[i].startswith("BB-") for i in ids),
                 "scope": scope}
            for j, m in enumerate(("sp", "tc", "cpu")):
                v = [x[j] for x in vals]
                lo, hi = _boot_ci(v)
                d.update({m: sum(v) / len(v), f"{m}_lo": lo, f"{m}_hi": hi})
            rows.append(d)
    return rows, per, cat


def fig_bali_categories(v2_refit, externals_tsv, out, scope="oof", old_tsvs=()):
    rows, _, _ = bali_categories(v2_refit, externals_tsv, scope, old_tsvs)
    arms = list(dict.fromkeys(r["arm"] for r in rows))
    colour = dict(CAT_ARMS + OLD_ARMS)
    name = "fig_bali_categories" + ("" if scope == "oof" else "_insample")
    write_csv(out, name, rows)

    def draw(f):
        axes = f.subplots(3, 1, sharex=True, gridspec_kw={"height_ratios": [1, 1, 0.8]})
        w = 0.8 / len(arms)
        for ax, met, lab in zip(axes, ("sp", "tc", "cpu"),
                                ("SP", "TC", "CPU per case (s)")):
            for j, arm in enumerate(arms):
                rr = [r for r in rows if r["arm"] == arm]
                x = [i + (j - (len(arms) - 1) / 2) * w for i in range(len(BALI_CATS))]
                old = arm in ("Kalign 3.1", "Kalign 3.3.5")
                ax.bar(x, [r[met] for r in rr], w, color="none" if old else colour[arm], label=arm,
                       edgecolor=colour[arm] if old else "none", linewidth=1.0 if old else 0,
                       hatch="////" if arm == "Kalign 3.1" else None,
                       yerr=[[r[met] - r[f"{met}_lo"] for r in rr], [r[f"{met}_hi"] - r[met] for r in rr]],
                       ecolor=ERR, error_kw={"elinewidth": 0.7})
            ax.set_ylabel(lab + (", held out" if scope == "oof" and met != "cpu" else
                                 ", in sample" if met != "cpu" else ""))
            if met == "cpu":
                log_y(ax)
            else:
                ax.set_ylim(0, 1)
        n = {r["category"]: r["n"] for r in rows}
        axes[2].set_xticks(range(len(BALI_CATS)))
        axes[2].set_xticklabels([f"{c}\nn = {n[c]}" for c in BALI_CATS])
        axes[0].legend(loc="upper center", ncol=4 if len(arms) > 5 else 5, frameon=False,
                       bbox_to_anchor=(0.5, 1.38 if len(arms) > 5 else 1.25), handlelength=1.0, columnspacing=0.8)
        f.tight_layout()
    finish((W, 7.2 if len(arms) > 5 else 6.8), out, name, draw)
    return rows


def fig_tradeoff_v2(v2_refit, out, pool):
    """Held-out SP and TC against CPU per case (benchmark-weighted): Kalign 3.6.1, 3.7.0, v2 modes, externals."""
    pts = []
    for r in tsv(Path(v2_refit) / pool / "results" / "oof_bw.tsv"):
        if r["subset_type"] != "all":
            continue
        arm, mode = r["arm"], r["mode"]
        if arm.startswith("kalign361"):
            arm, mode = "3.6.1", arm.split("_")[1]
        elif arm == "v370":
            if pool == "protein" and mode == "recall":
                continue
            arm = "3.7.0"
        elif arm == "new":
            if mode not in ("fast", "accurate"):
                continue
            arm = "re-fit"
        elif arm in OLD_KALIGN_NAMES:
            arm, mode = "old", OLD_KALIGN_NAMES[arm]
        elif arm in EXT_NAME:
            arm, mode = "external", EXT_NAME[arm]
        else:
            continue
        pts.append((arm, mode, {k: float(r[k]) for k in ("sp", "sp_lo", "sp_hi", "tc", "tc_lo", "tc_hi",
                                                          "precision", "cpu")}))
    write_csv(out, f"fig_tradeoff_v2_{pool}", [{"arm": a, "mode": m, **v} for a, m, v in pts])
    lab_c = {"3.6.1": "Kalign 3.6.1", "3.7.0": "Kalign 3.7.0", "re-fit": "Kalign new"}

    def draw(f):
        axes = f.subplots(1, 2)
        for ax, met, lab in zip(axes, ("sp", "tc"), ("SP", "TC")):
            for arm, m, v in pts:
                c = T["grey"] if arm == "old" else ARM_C[arm]
                mk = "*" if arm == "external" else (("h" if m.endswith("3.1") else "p") if arm == "old" else MARK[m])
                ax.errorbar(v["cpu"], v[met], yerr=[[v[met] - v[f"{met}_lo"]], [v[f"{met}_hi"] - v[met]]],
                            fmt=mk, ms=12 if mk == "*" else (8 if arm in ("re-fit", "old") else 6), color=c,
                            mfc="none" if arm in ("3.6.1", "old") else c, mec=c, ecolor=ERR, elinewidth=0.8,
                            capsize=0, zorder=5 if arm == "re-fit" else 3)

                if arm == "external":
                    off = {"MUSCLE 5": (-8, 6, "right"), "MAFFT L-INS-i": (8, -16, "left"),
                           "Clustal Omega": (8, -16, "left")}[m] if pool == "nucleotide" else (-6, 6, "right")
                    ax.annotate(m, (v["cpu"], v[met]), xytext=off[:2], textcoords="offset points", ha=off[2],
                                fontsize=PRESET.body_font_size, color=INK)
                if arm == "re-fit":
                    off = (0, 12, "center") if pool == "nucleotide" else (7, -4, "left")
                    ax.annotate(m, (v["cpu"], v[met]), xytext=off[:2], textcoords="offset points", ha=off[2],
                                fontsize=PRESET.body_font_size, color=T["rust"])
            log_x(ax)
            ax.set_xlabel("CPU per case (s, one thread)")
            ax.set_ylabel(f"{lab}, held out")
        hand = [Line2D([], [], marker="o", ls="", color=ARM_C["re-fit"], label="Kalign new"),
                Line2D([], [], marker="o", ls="", color=ARM_C["3.7.0"], label="Kalign 3.7.0"),
                Line2D([], [], marker="o", ls="", mfc="none", mec=ARM_C["3.6.1"], label="Kalign 3.6.1"),
                Line2D([], [], marker="h", ls="", ms=7, mfc="none", mec=T["grey"], label="Kalign 3.1"),
                Line2D([], [], marker="p", ls="", ms=7, mfc="none", mec=T["grey"], label="Kalign 3.3.5"),
                Line2D([], [], marker="*", ls="", ms=10, color=ARM_C["external"], label="other aligners")]
        hand += [Line2D([], [], marker=MARK[m], ls="", color=INK, mfc="none", label=m) for m in MODES]
        f.legend(handles=hand, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.08),
                 handletextpad=0.2, columnspacing=1.0)
        f.tight_layout(rect=(0, 0, 1, 0.88))
    finish((W, 4.0), out, f"fig_tradeoff_v2_{pool}", draw)
