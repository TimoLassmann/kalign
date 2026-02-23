#!/usr/bin/env python3
"""Generate a summary figure from downstream benchmark results.

Produces a 2x2 grid of line plots showing method performance vs difficulty,
with the "true" alignment as a dashed black ceiling line.

- Panel A: Phylo accuracy (nRF) vs tree_depth
- Panel B: Positive selection (F1) vs n_taxa
- Panel C: Alignment accuracy (SP) vs tree_depth
- Panel D: Alignment accuracy (SP) vs indel_rate
"""

import json
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = Path("benchmarks/results")

# Method display order (true last — rendered as dashed black)
METHOD_ORDER = [
    "kalign", "kalign_ens3",
    "mafft", "muscle", "clustalo", "true",
]
METHOD_LABELS = {
    "kalign": "Kalign",
    "kalign_ens3": "Kalign ens3",
    "mafft": "MAFFT",
    "muscle": "MUSCLE",
    "clustalo": "Clustal Omega",
    "true": "True alignment",
}
METHOD_COLORS = {
    "kalign": "#1f77b4",
    "kalign_ens3": "#2ca02c",
    "mafft": "#9467bd",
    "muscle": "#8c564b",
    "clustalo": "#7f7f7f",
    "true": "#000000",
}


def load_cases(pipeline: str) -> list[dict]:
    """Load per-case results from the latest.json for a pipeline."""
    path = RESULTS_DIR / pipeline / "latest.json"
    with open(path) as f:
        data = json.load(f)
    return data.get("cases", [])


def parse_sim_id(sim_id: str) -> dict:
    """Extract parameters from a sim_id string.

    Examples:
        WAG_t16_d2.0_ir0.10_il2.0_r0 -> {model:WAG, n_taxa:16, tree_depth:2.0, ...}
        M8_t16_d0.5_ir0.05_ps0.10_r0 -> {model:M8, n_taxa:16, tree_depth:0.5, ...}
    """
    params = {}
    m = re.match(r"^(\w+)_t(\d+)_d([\d.]+)_ir([\d.]+)", sim_id)
    if m:
        params["model"] = m.group(1)
        params["n_taxa"] = int(m.group(2))
        params["tree_depth"] = float(m.group(3))
        params["indel_rate"] = float(m.group(4))

    # Protein: _il{mean}_r{rep}
    m_il = re.search(r"_il([\d.]+)", sim_id)
    if m_il:
        params["indel_length_mean"] = float(m_il.group(1))

    # Codon: _ps{frac}_r{rep}
    m_ps = re.search(r"_ps([\d.]+)", sim_id)
    if m_ps:
        params["psel_fraction"] = float(m_ps.group(1))

    m_r = re.search(r"_r(\d+)$", sim_id)
    if m_r:
        params["replicate"] = int(m_r.group(1))

    return params


def group_by(cases: list[dict], param_key: str, metric_key: str):
    """Group per-case results by a sim_id parameter.

    Returns {method: {x_value: [metric_values]}}.
    """
    grouped = defaultdict(lambda: defaultdict(list))
    for c in cases:
        if "error" in c:
            continue
        sim_params = parse_sim_id(c.get("sim_id", ""))
        if param_key not in sim_params:
            continue
        x_val = sim_params[param_key]
        method = c["method"]
        val = c.get(metric_key)
        if val is not None and not (isinstance(val, float) and np.isnan(val)):
            grouped[method][x_val].append(val)
    return grouped


def style_ax(ax, hint=""):
    """Apply clean Nature-style aesthetics."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)
    if hint:
        ax.annotate(
            hint, xy=(0.98, 0.95), xycoords="axes fraction",
            ha="right", va="top", fontsize=8, fontstyle="italic", color="gray",
        )


def plot_lines(ax, grouped, methods, xlabel, ylabel):
    """Plot lines with error ribbons for each method.

    Parameters
    ----------
    grouped : dict
        {method: {x_value: [metric_values]}}
    methods : list[str]
        Method names in plot order.
    """
    for method in methods:
        if method not in grouped:
            continue
        data = grouped[method]
        xs = sorted(data.keys())
        means = []
        ses = []
        for x in xs:
            vals = np.asarray(data[x], dtype=float)
            means.append(vals.mean())
            ses.append(vals.std() / max(1, np.sqrt(len(vals))))

        means = np.asarray(means)
        ses = np.asarray(ses)

        color = METHOD_COLORS.get(method, "#333333")
        label = METHOD_LABELS.get(method, method)
        linestyle = "--" if method == "true" else "-"
        linewidth = 1.5 if method != "true" else 2.0
        marker = "o" if method != "true" else ""

        ax.plot(
            xs, means, color=color, linestyle=linestyle,
            linewidth=linewidth, marker=marker, markersize=4,
            label=label, zorder=3 if method != "true" else 2,
        )
        ax.fill_between(
            xs, means - ses, means + ses,
            color=color, alpha=0.12, zorder=1,
        )

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Kalign 3.5 — Downstream Benchmarks by Difficulty",
        fontsize=16, fontweight="bold", y=0.98,
    )

    # Determine which methods are present in results
    available = set()
    for pipeline in ("phylo_accuracy", "positive_selection", "calibration"):
        try:
            cases = load_cases(pipeline)
            for c in cases:
                if "error" not in c:
                    available.add(c["method"])
        except FileNotFoundError:
            pass
    methods = [m for m in METHOD_ORDER if m in available]

    # ── Panel A: Phylo accuracy (nRF) vs tree_depth ──────────────────
    ax = axes[0, 0]
    try:
        cases = load_cases("phylo_accuracy")
        grouped = group_by(cases, "tree_depth", "nrf")
        plot_lines(ax, grouped, methods, "Tree depth", "Normalized RF distance")
        style_ax(ax, "lower = better")
        ax.set_title("A) Phylogenetic tree accuracy", fontweight="bold", loc="left", fontsize=11)
    except FileNotFoundError:
        ax.text(0.5, 0.5, "No phylo_accuracy results", transform=ax.transAxes, ha="center")
        style_ax(ax)

    # ── Panel B: Positive selection (F1) vs n_taxa ───────────────────
    ax = axes[0, 1]
    try:
        cases = load_cases("positive_selection")
        grouped = group_by(cases, "n_taxa", "f1")
        plot_lines(ax, grouped, methods, "Number of taxa", "F1 score")
        style_ax(ax, "higher = better")
        ax.set_title("B) Positive selection detection (FUBAR)", fontweight="bold", loc="left", fontsize=11)
    except FileNotFoundError:
        ax.text(0.5, 0.5, "No positive_selection results", transform=ax.transAxes, ha="center")
        style_ax(ax)

    # ── Panel C: Alignment accuracy (SP) vs tree_depth ───────────────
    ax = axes[1, 0]
    try:
        cases = load_cases("phylo_accuracy")
        grouped = group_by(cases, "tree_depth", "sp_score")
        plot_lines(ax, grouped, methods, "Tree depth", "SP score vs true alignment")
        style_ax(ax, "higher = better")
        ax.set_title("C) Alignment accuracy vs tree depth", fontweight="bold", loc="left", fontsize=11)
    except FileNotFoundError:
        ax.text(0.5, 0.5, "No phylo_accuracy results", transform=ax.transAxes, ha="center")
        style_ax(ax)

    # ── Panel D: Alignment accuracy (SP) vs indel_rate ───────────────
    ax = axes[1, 1]
    try:
        cases = load_cases("phylo_accuracy")
        grouped = group_by(cases, "indel_rate", "sp_score")
        plot_lines(ax, grouped, methods, "Indel rate", "SP score vs true alignment")
        style_ax(ax, "higher = better")
        ax.set_title("D) Alignment accuracy vs indel rate", fontweight="bold", loc="left", fontsize=11)
    except FileNotFoundError:
        ax.text(0.5, 0.5, "No phylo_accuracy results", transform=ax.transAxes, ha="center")
        style_ax(ax)

    # Add legend (shared across panels)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels, loc="lower center",
            ncol=min(len(handles), 7), fontsize=9,
            frameon=False, bbox_to_anchor=(0.5, 0.0),
        )

    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    out = Path("benchmarks/figures/downstream_summary.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved to {out}")
    plt.close()


if __name__ == "__main__":
    main()
