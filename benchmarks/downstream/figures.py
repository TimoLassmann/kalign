"""Publication-quality figure generation for kalign downstream benchmarks.

Each public function produces a PDF figure suitable for inclusion in a
journal manuscript.  All figures read from pipeline result JSON files
via ``_load_latest_json()``.

All figures use:
- matplotlib with a clean, Nature-style aesthetic
- METHOD_COLORS from utils for consistent per-method colouring
- PDF output for vector graphics
- 10 pt base font, 1.5 pt minimum line width
- Single-column width = 3.5 in, double-column = 7 in
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any

from .utils import METHOD_COLORS

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Method ordering (consistent across all figures)
# ---------------------------------------------------------------------------

METHOD_ORDER = ["kalign", "kalign_cons", "kalign_ens3", "mafft", "muscle", "clustalo"]


# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------


def _setup_figure_style() -> None:
    """Configure matplotlib rcParams for publication-quality output."""
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            # Font
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            # Lines
            "lines.linewidth": 1.5,
            "lines.markersize": 5,
            # Axes
            "axes.linewidth": 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            # Grid
            "axes.grid": False,
            # Figure
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            # Legend
            "legend.frameon": False,
            "legend.borderaxespad": 0.5,
        }
    )


def _add_panel_label(ax: Any, label: str) -> None:
    """Add a bold panel label (e.g. 'A', 'B') in the upper-left corner."""
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def _method_color(method: str) -> str:
    """Return the colour for *method*, falling back to dark grey."""
    return METHOD_COLORS.get(method, "#333333")


def _method_label(method: str) -> str:
    """Return a human-friendly label for a method key."""
    labels = {
        "kalign": "Kalign",
        "kalign_cons": "Kalign+cons",
        "kalign_ens3": "Kalign ens3",
        "mafft": "MAFFT",
        "muscle": "MUSCLE",
        "clustalo": "Clustal Omega",
        "true": "True alignment",
    }
    return labels.get(method, method)


def _savefig(fig: Any, output_path: Path) -> None:
    """Save figure as PDF and close it."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), format="pdf")
    import matplotlib.pyplot as plt

    plt.close(fig)
    logger.info("Saved figure: %s", output_path)


def _bootstrap_ci_vec(
    values: list[float], n_bootstrap: int = 5000, alpha: float = 0.05
) -> tuple[float, float, float]:
    """Return (mean, ci_lo, ci_hi) via bootstrap resampling."""
    import numpy as np

    arr = np.array(values, dtype=float)
    n = len(arr)
    if n < 2:
        mu = float(arr.mean())
        return mu, mu, mu
    rng = np.random.default_rng(42)
    boot_means = np.array(
        [float(rng.choice(arr, size=n, replace=True).mean()) for _ in range(n_bootstrap)]
    )
    lo = float(np.percentile(boot_means, 100 * alpha / 2))
    hi = float(np.percentile(boot_means, 100 * (1 - alpha / 2)))
    return float(arr.mean()), lo, hi


def _shared_legend(fig: Any, methods: list[str], ncol: int = 6) -> None:
    """Add a shared legend below all panels."""
    from matplotlib.patches import Patch

    handles = [
        Patch(facecolor=_method_color(m), alpha=0.7, label=_method_label(m))
        for m in methods
    ]
    fig.legend(
        handles=handles, loc="lower center", ncol=ncol,
        fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02),
    )


def _ordered_methods(available: set[str]) -> list[str]:
    """Return methods in canonical order, filtered to those in *available*."""
    return [m for m in METHOD_ORDER if m in available]


def _grouped_boxplot(
    ax: Any,
    data: dict[str, dict[str, list[float]]],
    categories: list[str],
    methods: list[str],
    ylabel: str,
    cat_labels: dict[str, str] | None = None,
) -> None:
    """Draw grouped box plots on *ax*.

    Parameters
    ----------
    data : dict[method -> dict[category -> list[float]]]
    categories : list of category keys
    methods : list of method keys (in order)
    """
    n_methods = len(methods)
    n_cats = len(categories)
    if n_methods == 0 or n_cats == 0:
        return

    width = 0.8 / n_methods
    for j, method in enumerate(methods):
        positions = []
        box_data = []
        color = _method_color(method)
        for i, cat in enumerate(categories):
            vals = data.get(method, {}).get(cat, [])
            box_data.append(vals if vals else [0])
            positions.append(i + j * width - (n_methods - 1) * width / 2)

        bp = ax.boxplot(
            box_data,
            positions=positions,
            widths=width * 0.85,
            patch_artist=True,
            medianprops={"color": "black", "linewidth": 1.2},
            whiskerprops={"linewidth": 0.8},
            capprops={"linewidth": 0.8},
            flierprops={"marker": ".", "markersize": 3, "alpha": 0.5},
            manage_ticks=False,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

    tick_labels = []
    for cat in categories:
        label = (cat_labels or {}).get(cat, cat)
        # Count cases from first method with data
        n = 0
        for m in methods:
            n = len(data.get(m, {}).get(cat, []))
            if n > 0:
                break
        tick_labels.append(f"{label}\n(n={n})" if n > 0 else label)

    ax.set_xticks(range(n_cats))
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(-0.5, n_cats - 0.5)


# ---------------------------------------------------------------------------
# Result loading helpers
# ---------------------------------------------------------------------------


def _load_latest_json(results_dir: Path) -> Any:
    """Load the ``latest.json`` symlink (or most recent run file)."""
    results_dir = Path(results_dir)
    if not results_dir.is_dir():
        return None

    latest = results_dir / "latest.json"
    if latest.exists():
        with open(latest) as fh:
            return json.load(fh)

    json_files = sorted(results_dir.glob("run_*.json"))
    if not json_files:
        return None

    with open(json_files[-1]) as fh:
        return json.load(fh)


def _extract_cases(data: Any) -> list[dict]:
    """Extract the cases list from a pipeline result dict or list."""
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return data.get("cases", [])
    return []


# ---------------------------------------------------------------------------
# Figure 1: BAliBASE SP/TC by category
# ---------------------------------------------------------------------------

_BALIBASE_CATEGORIES = ["RV11", "RV12", "RV20", "RV30", "RV40", "RV50"]


def figure_balibase(cases: list[dict], output_path: Path) -> None:
    """BAliBASE alignment accuracy: SP, Precision, F1, and TC box plots.

    Panel A: SP (recall) per BAliBASE category, grouped by method.
    Panel B: Precision per BAliBASE category, grouped by method.
    Panel C: F1 per BAliBASE category, grouped by method.
    Panel D: TC per BAliBASE category, grouped by method.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _setup_figure_style()

    fig, ((ax_sp, ax_prec), (ax_f1, ax_tc)) = plt.subplots(2, 2, figsize=(14, 10))

    # Filter to balibase cases only
    bb_cases = [c for c in cases if c.get("benchmark") == "balibase" and "error" not in c and not c.get("skipped")]
    if not bb_cases:
        _savefig(fig, output_path)
        return

    available = {c["method"] for c in bb_cases}
    methods = _ordered_methods(available)

    # Build data[method][category] -> list of values
    def _build(metric: str) -> dict[str, dict[str, list[float]]]:
        out: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
        for c in bb_cases:
            m = c["method"]
            ds = c.get("dataset", "")
            val = c.get(metric)
            if val is None or val < 0:
                continue
            for cat in _BALIBASE_CATEGORIES:
                if cat in ds:
                    out[m][cat].append(val)
                    break
        return out

    sp_data = _build("recall")
    prec_data = _build("precision")
    f1_data = _build("f1")
    tc_data = _build("tc")

    _grouped_boxplot(ax_sp, sp_data, _BALIBASE_CATEGORIES, methods, "SP (recall)")
    ax_sp.set_title("Sum-of-pairs score by category", fontsize=10)
    _add_panel_label(ax_sp, "A")

    _grouped_boxplot(ax_prec, prec_data, _BALIBASE_CATEGORIES, methods, "Precision")
    ax_prec.set_title("Precision by category", fontsize=10)
    _add_panel_label(ax_prec, "B")

    _grouped_boxplot(ax_f1, f1_data, _BALIBASE_CATEGORIES, methods, "F1")
    ax_f1.set_title("F1 score by category", fontsize=10)
    _add_panel_label(ax_f1, "C")

    _grouped_boxplot(ax_tc, tc_data, _BALIBASE_CATEGORIES, methods, "TC")
    ax_tc.set_title("Total column score by category", fontsize=10)
    _add_panel_label(ax_tc, "D")

    _shared_legend(fig, methods)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 2: BRAliBASE SP by RNA family
# ---------------------------------------------------------------------------

_RNA_FAMILIES = ["SRP", "tRNA", "rRNA", "g2intron", "U5"]
_RNA_LABELS = {
    "SRP": "SRP",
    "tRNA": "tRNA",
    "rRNA": "rRNA",
    "g2intron": "Group II\nintron",
    "U5": "U5",
}


def figure_bralibase(cases: list[dict], output_path: Path) -> None:
    """BRAliBASE RNA alignment accuracy.

    Panel A: SP (recall) per RNA family, grouped by method.
    Panel B: F1 per RNA family, grouped by method.
    Panel C: Wall time box plots per method (log scale).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_sp, ax_f1, ax_time) = plt.subplots(1, 3, figsize=(18, 5))

    br_cases = [c for c in cases if c.get("benchmark") == "bralibase" and "error" not in c and not c.get("skipped")]
    if not br_cases:
        _savefig(fig, output_path)
        return

    available = {c["method"] for c in br_cases}
    methods = _ordered_methods(available)

    # Build data per family
    def _build_family(metric: str) -> dict[str, dict[str, list[float]]]:
        out: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
        for c in br_cases:
            m = c["method"]
            ds = c.get("dataset", "")
            val = c.get(metric)
            if val is None or val < 0:
                continue
            for fam in _RNA_FAMILIES:
                if fam.lower() in ds.lower() or fam in ds:
                    out[m][fam].append(val)
                    break
        return out

    sp_data = _build_family("recall")
    f1_data = _build_family("f1")

    _grouped_boxplot(ax_sp, sp_data, _RNA_FAMILIES, methods, "SP (recall)", _RNA_LABELS)
    ax_sp.set_title("Sum-of-pairs score by RNA family", fontsize=10)
    _add_panel_label(ax_sp, "A")

    _grouped_boxplot(ax_f1, f1_data, _RNA_FAMILIES, methods, "F1", _RNA_LABELS)
    ax_f1.set_title("F1 score by RNA family", fontsize=10)
    _add_panel_label(ax_f1, "B")

    # Panel C: Wall time
    time_data = []
    time_colors = []
    time_labels = []
    for m in methods:
        times = [c["wall_time"] for c in br_cases if c.get("method") == m and c.get("wall_time", 0) > 0]
        if times:
            time_data.append(times)
            time_colors.append(_method_color(m))
            time_labels.append(_method_label(m))

    if time_data:
        bp = ax_time.boxplot(
            time_data,
            patch_artist=True,
            widths=0.6,
            medianprops={"color": "black", "linewidth": 1.2},
        )
        for patch, c in zip(bp["boxes"], time_colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_time.set_yscale("log")
        ax_time.set_xticks(range(1, len(time_labels) + 1))
        ax_time.set_xticklabels(time_labels, fontsize=8, rotation=45, ha="right")
        ax_time.set_ylabel("Wall time (s)")
        ax_time.set_title("Alignment time", fontsize=10)
    else:
        ax_time.text(0.5, 0.5, "No timing data", ha="center", va="center",
                     transform=ax_time.transAxes, fontsize=10, color="#999999")
    _add_panel_label(ax_time, "C")

    _shared_legend(fig, methods)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 3: Phylogenetic tree accuracy
# ---------------------------------------------------------------------------


def figure_phylo_accuracy(cases: list[dict], output_path: Path) -> None:
    """Phylogenetic tree accuracy.

    Panel A: nRF box plots per method.
    Panel B: nRF vs tree depth with CI bands.
    Panel C: Branch score distance box plots.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_nrf, ax_depth, ax_bsd) = plt.subplots(1, 3, figsize=(14, 5))

    valid = [c for c in cases if "error" not in c]
    available = {c["method"] for c in valid}
    methods = _ordered_methods(available | ({"true"} & available))

    # Panel A: nRF box plots
    if methods:
        box_data = []
        colors = []
        labels = []
        for m in methods:
            vals = [c["nrf"] for c in valid if c.get("method") == m]
            box_data.append(vals if vals else [0])
            colors.append(_method_color(m))
            labels.append(_method_label(m))

        bp = ax_nrf.boxplot(
            box_data, patch_artist=True, widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_nrf.set_xticks(range(1, len(methods) + 1))
        ax_nrf.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax_nrf.set_ylabel("nRF distance")
    _add_panel_label(ax_nrf, "A")

    # Panel B: nRF vs tree depth
    import re as _re
    for m in methods:
        depth_nrf: dict[float, list[float]] = defaultdict(list)
        for c in valid:
            if c.get("method") != m:
                continue
            d = c.get("tree_depth")
            if d is None:
                dm = _re.search(r"_d([\d.]+)_", c.get("sim_id", ""))
                if dm:
                    d = float(dm.group(1))
            nrf = c.get("nrf")
            if d is not None and nrf is not None:
                depth_nrf[d].append(nrf)
        if depth_nrf:
            depths = sorted(depth_nrf.keys())
            means = [float(np.mean(depth_nrf[d])) for d in depths]
            stds = [float(np.std(depth_nrf[d])) for d in depths]
            lo = [max(0, mu - s) for mu, s in zip(means, stds)]
            hi = [mu + s for mu, s in zip(means, stds)]
            color = _method_color(m)
            ax_depth.plot(depths, means, marker="o", color=color, label=_method_label(m), markersize=4)
            ax_depth.fill_between(depths, lo, hi, color=color, alpha=0.07)

    ax_depth.set_xlabel("Tree depth (subs/site)")
    ax_depth.set_ylabel("Mean nRF distance")
    handles, labels = ax_depth.get_legend_handles_labels()
    if handles:
        ax_depth.legend(fontsize=6, ncol=2, loc="upper left")
    _add_panel_label(ax_depth, "B")

    # Panel C: Branch score distance
    if methods:
        bsd_data = []
        colors = []
        labels = []
        for m in methods:
            vals = [c["branch_score_dist"] for c in valid if c.get("method") == m and "branch_score_dist" in c]
            bsd_data.append(vals if vals else [0])
            colors.append(_method_color(m))
            labels.append(_method_label(m))

        bp = ax_bsd.boxplot(
            bsd_data, patch_artist=True, widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_bsd.set_xticks(range(1, len(methods) + 1))
        ax_bsd.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax_bsd.set_ylabel("Branch score distance")
    _add_panel_label(ax_bsd, "C")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 4: Positive selection detection
# ---------------------------------------------------------------------------


def figure_positive_selection(cases: list[dict], output_path: Path) -> None:
    """Positive selection detection.

    Panel A: Precision/Recall/F1 bars per method.
    Panel B: F1 vs tree depth.
    Panel C: F1 vs indel rate.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import re as _re

    _setup_figure_style()

    fig, (ax_bar, ax_depth, ax_indel) = plt.subplots(1, 3, figsize=(14, 5))

    valid = [c for c in cases if "error" not in c]
    by_method: dict[str, list[dict]] = defaultdict(list)
    for c in valid:
        by_method[c.get("method", "")].append(c)

    available = set(by_method.keys())
    methods = _ordered_methods(available | ({"true"} & available))

    # Panel A: bar chart
    if methods:
        n_methods = len(methods)
        x_pos = np.arange(n_methods)
        width = 0.25

        prec_means = [np.mean([c.get("precision", 0) for c in by_method[m]]) for m in methods]
        rec_means = [np.mean([c.get("recall", c.get("tpr", 0)) for c in by_method[m]]) for m in methods]
        f1_means = [np.mean([c.get("f1", 0) for c in by_method[m]]) for m in methods]

        colors_prec = [_method_color(m) for m in methods]
        ax_bar.bar(x_pos - width, prec_means, width, label="Precision",
                   color=colors_prec, alpha=0.8)
        ax_bar.bar(x_pos, rec_means, width, label="Recall",
                   color=colors_prec, alpha=0.55)
        ax_bar.bar(x_pos + width, f1_means, width, label="F1",
                   color=colors_prec, alpha=0.3)
        ax_bar.set_xticks(x_pos)
        ax_bar.set_xticklabels(
            [_method_label(m) for m in methods], rotation=45, ha="right", fontsize=7
        )
        ax_bar.set_ylabel("Score")
        all_vals = list(prec_means) + list(rec_means) + list(f1_means)
        if all_vals:
            ax_bar.set_ylim(0, max(all_vals) * 1.15)
        ax_bar.legend(fontsize=7)
    _add_panel_label(ax_bar, "A")

    # Panel B: F1 vs tree depth
    for m in methods:
        depth_f1: dict[float, list[float]] = defaultdict(list)
        for c in by_method[m]:
            d = c.get("tree_depth")
            if d is None:
                dm = _re.search(r"_d([\d.]+)_", c.get("sim_id", ""))
                if dm:
                    d = float(dm.group(1))
            f1 = c.get("f1")
            if d is not None and f1 is not None:
                depth_f1[d].append(f1)
        if depth_f1:
            depths_sorted = sorted(depth_f1.keys())
            means = [np.mean(depth_f1[d]) for d in depths_sorted]
            ax_depth.plot(depths_sorted, means, marker="o", color=_method_color(m),
                          label=_method_label(m))

    ax_depth.set_xlabel("Tree depth (subs/site)")
    ax_depth.set_ylabel("F1 score")
    handles, _ = ax_depth.get_legend_handles_labels()
    if handles:
        ax_depth.legend(fontsize=6, ncol=2, loc="best")
    _add_panel_label(ax_depth, "B")

    # Panel C: F1 vs indel rate
    for m in methods:
        indel_f1: dict[float, list[float]] = defaultdict(list)
        for c in by_method[m]:
            ir = c.get("indel_rate")
            if ir is None:
                irm = _re.search(r"_ir([\d.]+)_", c.get("sim_id", ""))
                if irm:
                    ir = float(irm.group(1))
            f1 = c.get("f1")
            if ir is not None and f1 is not None:
                indel_f1[ir].append(f1)
        if indel_f1:
            rates_sorted = sorted(indel_f1.keys())
            means = [np.mean(indel_f1[r]) for r in rates_sorted]
            ax_indel.plot(rates_sorted, means, marker="o", color=_method_color(m),
                          label=_method_label(m))

    ax_indel.set_xlabel("Indel rate")
    ax_indel.set_ylabel("F1 score")
    handles, _ = ax_indel.get_legend_handles_labels()
    if handles:
        ax_indel.legend(fontsize=6, ncol=2, loc="best")
    _add_panel_label(ax_indel, "C")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 5: HMMER homology detection
# ---------------------------------------------------------------------------


def figure_hmmer_detection(cases: list[dict], output_path: Path) -> None:
    """HMMER homology detection — single bar chart with bootstrap CIs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    valid = [c for c in cases if "error" not in c]
    by_method: dict[str, list[dict]] = defaultdict(list)
    for c in valid:
        by_method[c.get("method", "")].append(c)

    available = set(by_method.keys())
    methods = _ordered_methods(available)

    def _sens(case: dict) -> float:
        if "sensitivity" in case:
            return case["sensitivity"]
        h = case.get("hits_at_1e3", {})
        return h.get("sens", 0.0) if isinstance(h, dict) else 0.0

    if methods:
        means = []
        ci_lo = []
        ci_hi = []
        colors = []
        for m in methods:
            sens_vals = [_sens(c) for c in by_method[m]]
            mu = float(np.mean(sens_vals)) if sens_vals else 0.0
            means.append(mu)
            if len(sens_vals) >= 2:
                _, lo, hi = _bootstrap_ci_vec(sens_vals, n_bootstrap=2000)
                ci_lo.append(mu - lo)
                ci_hi.append(hi - mu)
            else:
                ci_lo.append(0)
                ci_hi.append(0)
            colors.append(_method_color(m))

        x = np.arange(len(methods))
        ax.bar(
            x, means, color=colors, alpha=0.8,
            yerr=[ci_lo, ci_hi], capsize=3, ecolor="#333333",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(
            [_method_label(m) for m in methods], rotation=45, ha="right", fontsize=8
        )
        ax.set_ylabel("Mean sensitivity (E < 1e-5)")
        if means and max(means) > 0:
            ymin_bar = max(0, min(means) - 0.05)
            ymax = min(1.0, max(mu + h for mu, h in zip(means, ci_hi)) * 1.02)
            ax.set_ylim(ymin_bar, ymax)
        ax.set_title("HMMER homology detection (50 Pfam families)", fontsize=10)

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 6: Speed comparison (aggregated across all pipelines)
# ---------------------------------------------------------------------------


def figure_speed_comparison(all_timed: list[dict], output_path: Path) -> None:
    """Speed and memory comparison aggregated across all pipelines.

    Panel A: Wall time box plots (log scale) per method.
    Panel B: Peak memory box plots per method.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_time, ax_mem) = plt.subplots(1, 2, figsize=(10, 5))

    by_method: dict[str, list[dict]] = defaultdict(list)
    for r in all_timed:
        m = r.get("method")
        if m:
            by_method[m].append(r)

    available = set(by_method.keys())
    methods = _ordered_methods(available)

    # Panel A: Wall time
    if methods:
        time_data = []
        colors = []
        labels = []
        for m in methods:
            # Try multiple time field names
            times = []
            for c in by_method[m]:
                t = c.get("wall_time_align") or c.get("wall_time")
                if t and t > 0:
                    times.append(t)
            time_data.append(times if times else [0.001])
            colors.append(_method_color(m))
            labels.append(_method_label(m))

        bp = ax_time.boxplot(
            time_data, patch_artist=True, widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_time.set_yscale("log")
        ax_time.set_xticks(range(1, len(methods) + 1))
        ax_time.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax_time.set_ylabel("Wall time (s)")
    _add_panel_label(ax_time, "A")

    # Panel B: Peak memory
    if methods:
        mem_data = []
        colors = []
        labels = []
        for m in methods:
            mems = [c.get("peak_memory_mb", 0) for c in by_method[m] if c.get("peak_memory_mb", 0) > 0]
            mem_data.append(mems if mems else [0])
            colors.append(_method_color(m))
            labels.append(_method_label(m))

        bp = ax_mem.boxplot(
            mem_data, patch_artist=True, widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_mem.set_xticks(range(1, len(methods) + 1))
        ax_mem.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax_mem.set_ylabel("Peak memory (MB)")
    _add_panel_label(ax_mem, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 7: Summary heatmap
# ---------------------------------------------------------------------------


def figure_summary_heatmap(all_results: dict, output_path: Path) -> None:
    """Summary heatmap across all pipelines.

    Rows = methods, columns = key metrics per pipeline.
    Colour-coded: green=best, red=worst.  Annotated with raw values.

    Parameters
    ----------
    all_results : dict
        Mapping with optional keys for each pipeline.
        Each value is a dict mapping method names to metric values (float).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    # Columns: (label, data_dict, higher_is_better)
    columns_spec: list[tuple[str, dict[str, float], bool]] = []

    # Alignment accuracy F1
    aln = all_results.get("alignment_accuracy")
    if aln:
        columns_spec.append(("Aln F1", aln, True))

    # Calibration ECE
    cal = all_results.get("calibration")
    if cal:
        columns_spec.append(("ECE", cal, False))

    # Positive selection F1
    psel = all_results.get("positive_selection")
    if psel:
        columns_spec.append(("Sel. F1", psel, True))

    # Phylo nRF
    phylo = all_results.get("phylo_accuracy")
    if phylo:
        columns_spec.append(("nRF", phylo, False))

    # HMMER sensitivity
    hmmer = all_results.get("hmmer_detection")
    if hmmer:
        columns_spec.append(("HMMER sens.", hmmer, True))

    # Speed
    speed = all_results.get("speed")
    if speed:
        columns_spec.append(("Speed (s)", speed, False))

    if not columns_spec:
        logger.warning("figure_summary_heatmap: no data, skipping.")
        return

    # Collect union of methods
    all_methods: set[str] = set()
    for _, mdict, _ in columns_spec:
        all_methods.update(mdict.keys())
    methods = _ordered_methods(all_methods)
    if not methods:
        methods = sorted(all_methods)

    n_methods = len(methods)
    n_cols = len(columns_spec)

    raw = np.full((n_methods, n_cols), np.nan)
    for j, (_, mdict, _) in enumerate(columns_spec):
        for i, m in enumerate(methods):
            if m in mdict:
                raw[i, j] = mdict[m]

    # Normalise each column to [0, 1] (1 = best)
    normed = np.full_like(raw, np.nan)
    for j, (_, _, higher_better) in enumerate(columns_spec):
        col = raw[:, j]
        valid_mask = ~np.isnan(col)
        if valid_mask.sum() < 2:
            normed[valid_mask, j] = 0.5
            continue
        vmin = np.nanmin(col)
        vmax = np.nanmax(col)
        if vmax - vmin < 1e-12:
            normed[valid_mask, j] = 0.5
        else:
            scaled = (col[valid_mask] - vmin) / (vmax - vmin)
            if not higher_better:
                scaled = 1.0 - scaled
            normed[valid_mask, j] = scaled

    fig, ax = plt.subplots(figsize=(max(5, 1.7 * n_cols + 2), 0.65 * n_methods + 1.5))

    cmap = plt.get_cmap("RdYlGn")
    display = np.where(np.isnan(normed), 0.5, normed)
    im = ax.imshow(display, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    for i in range(n_methods):
        for j in range(n_cols):
            val = raw[i, j]
            if np.isnan(val):
                text = "-"
            elif abs(val) < 0.01 or abs(val) >= 1000:
                text = f"{val:.2e}"
            else:
                text = f"{val:.3f}"
            ax.text(
                j, i, text, ha="center", va="center", fontsize=7,
                color="black" if 0.3 < display[i, j] < 0.7 else "white",
            )

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([name for name, _, _ in columns_spec], fontsize=9)
    ax.set_yticks(range(n_methods))
    ax.set_yticklabels([_method_label(m) for m in methods], fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Normalised score (1 = best)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Calibration figure (for calibration pipeline)
# ---------------------------------------------------------------------------


def figure_calibration(results: dict, output_path: Path) -> None:
    """Confidence calibration reliability diagram.

    Panel A: Calibration curve (predicted confidence vs fraction correct).
    Panel B: Histogram of confidence values.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_cal, ax_hist) = plt.subplots(1, 2, figsize=(10, 5))

    ax_cal.plot([0, 1], [0, 1], ls="--", color="#aaaaaa", lw=1.0, label="Perfect")

    for method, data in results.items():
        bins = data.get("bins")
        if not bins:
            continue
        midpoints = np.array([(b["bin_low"] + b["bin_high"]) / 2 for b in bins])
        accuracies = np.array([b.get("accuracy", 0.0) for b in bins])

        ece = data.get("ece")
        brier = data.get("brier_score")
        label_parts = [_method_label(method)]
        if ece is not None:
            label_parts.append(f"ECE={ece:.3f}")
        if brier is not None:
            label_parts.append(f"Brier={brier:.3f}")

        ax_cal.plot(
            midpoints, accuracies, marker="o",
            color=_method_color(method), label=" ".join(label_parts),
        )

    ax_cal.set_xlabel("Predicted confidence")
    ax_cal.set_ylabel("Observed accuracy")
    ax_cal.set_xlim(0, 1)
    ax_cal.set_ylim(0, 1)
    ax_cal.set_aspect("equal")
    ax_cal.legend(loc="lower right", fontsize=7)
    _add_panel_label(ax_cal, "A")

    for method, data in results.items():
        bins = data.get("bins")
        if not bins:
            continue
        midpoints = np.array([(b["bin_low"] + b["bin_high"]) / 2 for b in bins])
        counts = np.array([b.get("n_pairs", 0) for b in bins], dtype=float)
        total = counts.sum()
        if total > 0:
            counts = counts / total
        widths = np.array([b["bin_high"] - b["bin_low"] for b in bins])
        ax_hist.bar(
            midpoints, counts, width=widths * 0.85,
            alpha=0.5, color=_method_color(method), label=_method_label(method),
        )

    ax_hist.set_xlabel("Confidence score")
    ax_hist.set_ylabel("Fraction of residue pairs")
    ax_hist.set_xlim(0, 1)
    ax_hist.legend(loc="upper left", fontsize=7)
    _add_panel_label(ax_hist, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Aggregation helpers for heatmap
# ---------------------------------------------------------------------------


def _aggregate_for_heatmap(
    pipeline_results: list[dict],
    method_key: str,
    metric_key: str,
) -> dict[str, float]:
    """Compute mean of *metric_key* per method from a list of per-case dicts."""
    import numpy as np

    by_method: dict[str, list[float]] = defaultdict(list)
    parts = metric_key.split(".")
    for case in pipeline_results:
        m = case.get(method_key)
        if m is None:
            continue
        val = case
        for p in parts:
            if isinstance(val, dict):
                val = val.get(p)
            else:
                val = None
                break
        if val is not None:
            try:
                by_method[m].append(float(val))
            except (TypeError, ValueError):
                pass

    return {m: float(np.mean(vs)) for m, vs in by_method.items() if vs}


def _adapt_calibration_summary(cal_methods: dict) -> dict:
    """Convert calibration summary format to the bins format expected by figure_calibration."""
    import math

    adapted: dict = {}
    for method, data in cal_methods.items():
        if not isinstance(data, dict):
            continue
        cc = data.get("calibration_curve")
        if not isinstance(cc, dict):
            adapted[method] = data
            continue
        centers = cc.get("bin_centers", [])
        counts = cc.get("bin_counts", [])
        fracs = cc.get("fraction_correct", [])
        if not centers:
            adapted[method] = data
            continue
        bins = []
        for i, center in enumerate(centers):
            if len(centers) > 1:
                width = centers[1] - centers[0] if i == 0 else centers[i] - centers[i - 1]
            else:
                width = 0.1
            acc = fracs[i] if i < len(fracs) else 0.0
            if isinstance(acc, float) and math.isnan(acc):
                acc = 0.0
            bins.append({
                "bin_low": center - width / 2,
                "bin_high": center + width / 2,
                "n_pairs": counts[i] if i < len(counts) else 0,
                "accuracy": acc,
            })
        adapted[method] = {
            "ece": data.get("ece"),
            "brier_score": data.get("brier_score"),
            "bins": bins,
        }
    return adapted


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def generate_all_figures(
    results_dir: str | Path,
    output_dir: str | Path,
) -> None:
    """Generate all publication figures from saved results.

    Reads ``latest.json`` from each pipeline's results sub-directory and
    calls the corresponding figure function.  Skips pipelines whose
    results are not available.
    """
    results_dir = Path(results_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated: list[str] = []

    # --- Fig 1 & 2: Alignment accuracy (BAliBASE + BRAliBASE) ----------------
    aln_data = _load_latest_json(results_dir / "alignment_accuracy")
    aln_cases = _extract_cases(aln_data) if aln_data else []
    aln_valid = [c for c in aln_cases if "error" not in c and not c.get("skipped")]

    if aln_valid:
        path = output_dir / "fig1_balibase_accuracy.pdf"
        try:
            figure_balibase(aln_valid, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate BAliBASE figure: %s", exc)

        path = output_dir / "fig2_bralibase_rna.pdf"
        try:
            figure_bralibase(aln_valid, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate BRAliBASE figure: %s", exc)

    # --- Fig 3: Phylo accuracy -----------------------------------------------
    phylo_data = _load_latest_json(results_dir / "phylo_accuracy")
    phylo_cases = _extract_cases(phylo_data)
    if phylo_cases:
        path = output_dir / "fig3_phylo_accuracy.pdf"
        try:
            figure_phylo_accuracy(phylo_cases, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate phylo accuracy figure: %s", exc)

    # --- Fig 4: Positive selection -------------------------------------------
    psel_data = _load_latest_json(results_dir / "positive_selection")
    psel_cases = _extract_cases(psel_data)
    if psel_cases:
        path = output_dir / "fig4_positive_selection.pdf"
        try:
            figure_positive_selection(psel_cases, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate positive selection figure: %s", exc)

    # --- Fig 5: HMMER detection ----------------------------------------------
    hmmer_data = _load_latest_json(results_dir / "hmmer_detection")
    hmmer_cases = _extract_cases(hmmer_data)
    if hmmer_cases:
        path = output_dir / "fig5_hmmer_detection.pdf"
        try:
            figure_hmmer_detection(hmmer_cases, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate HMMER detection figure: %s", exc)

    # --- Fig 6: Speed comparison (aggregated) --------------------------------
    all_timed: list[dict] = []
    # Collect timing data from alignment_accuracy + downstream pipelines
    for c in aln_valid:
        if "wall_time" in c:
            all_timed.append(c)
    for subdir in ("positive_selection", "phylo_accuracy"):
        sub_data = _load_latest_json(results_dir / subdir)
        sub_cases = _extract_cases(sub_data)
        all_timed.extend(c for c in sub_cases if "wall_time_align" in c or "wall_time" in c)

    if all_timed:
        path = output_dir / "fig6_speed_comparison.pdf"
        try:
            figure_speed_comparison(all_timed, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate speed comparison figure: %s", exc)

    # --- Fig 7: Summary heatmap ----------------------------------------------
    heatmap_data: dict[str, dict[str, float]] = {}

    # Alignment accuracy F1
    if aln_valid:
        agg = _aggregate_for_heatmap(aln_valid, "method", "f1")
        if agg:
            heatmap_data["alignment_accuracy"] = agg

    # Calibration ECE
    cal_data = _load_latest_json(results_dir / "calibration")
    if cal_data is not None:
        cal_src = None
        if isinstance(cal_data, dict):
            cal_src = cal_data.get("methods") or cal_data.get("summary") or cal_data
        if isinstance(cal_src, dict):
            ece_map = {}
            for m, v in cal_src.items():
                if isinstance(v, dict) and "ece" in v:
                    ece_map[m] = v["ece"]
            if ece_map:
                heatmap_data["calibration"] = ece_map

    # Calibration figure (standalone)
    if cal_data is not None:
        cal_methods = None
        if isinstance(cal_data, dict):
            cal_methods = cal_data.get("methods") or cal_data.get("summary")
            if cal_methods is None:
                cal_methods = cal_data
        if isinstance(cal_methods, dict):
            cal_methods = _adapt_calibration_summary(cal_methods)
            if cal_methods:
                path = output_dir / "fig_calibration.pdf"
                try:
                    figure_calibration(cal_methods, path)
                    generated.append(str(path))
                except Exception as exc:
                    logger.error("Failed to generate calibration figure: %s", exc)

    # Positive selection F1
    if psel_cases:
        agg = _aggregate_for_heatmap(psel_cases, "method", "f1")
        if agg:
            heatmap_data["positive_selection"] = agg

    # Phylo nRF
    if phylo_cases:
        agg = _aggregate_for_heatmap(phylo_cases, "method", "nrf")
        if agg:
            heatmap_data["phylo_accuracy"] = agg

    # HMMER sensitivity
    if hmmer_cases:
        agg = _aggregate_for_heatmap(hmmer_cases, "method", "sensitivity")
        if agg:
            heatmap_data["hmmer_detection"] = agg

    # Speed (mean wall time)
    if all_timed:
        # Normalise field name
        for c in all_timed:
            if "wall_time_align" in c and "wall_time" not in c:
                c["wall_time"] = c["wall_time_align"]
        agg = _aggregate_for_heatmap(all_timed, "method", "wall_time")
        if agg:
            heatmap_data["speed"] = agg

    if heatmap_data:
        path = output_dir / "fig7_summary_heatmap.pdf"
        try:
            figure_summary_heatmap(heatmap_data, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate summary heatmap: %s", exc)

    # --- Summary --------------------------------------------------------------
    if generated:
        print(f"Generated {len(generated)} figures:")
        for p in generated:
            print(f"  {p}")
    else:
        print("No figures generated (no result data found).")
