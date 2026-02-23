"""Publication-quality figure generation for kalign downstream benchmarks.

Each public function takes result data (dicts/lists) and an output path,
and produces a PDF figure suitable for inclusion in a journal manuscript.

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
        "kalign_basic": "Kalign basic",
        "kalign_vsm": "Kalign VSM",
        "kalign_ens3": "Kalign ens3",
        "kalign_ens3_m50": "Kalign ens3 m50",
        "kalign_ens3_m70": "Kalign ens3 m70",
        "kalign_ens3_m90": "Kalign ens3 m90",
        "kalign_ens3_wt": "Kalign ens3 wt",
        "kalign_best": "Kalign best",
        "kalign_fast": "Kalign fast",
        "guidance2_mafft": "GUIDANCE2+MAFFT",
        "guidance2_m93": "GUIDANCE2 m93",
        "guidance2_wt": "GUIDANCE2 wt",
        "mafft": "MAFFT",
        "muscle": "MUSCLE",
        "clustalo": "Clustal Omega",
        "true": "True alignment",
        "pfam_original": "Pfam original",
    }
    return labels.get(method, method)


def _legend_if_labeled(ax: Any, **kwargs: Any) -> None:
    """Add a legend only if the axes has artists with labels."""
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(**kwargs)


def _savefig(fig: Any, output_path: Path) -> None:
    """Save figure as PDF and close it."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), format="pdf")
    import matplotlib.pyplot as plt

    plt.close(fig)
    logger.info("Saved figure: %s", output_path)


# ---------------------------------------------------------------------------
# Figure 1: Confidence calibration
# ---------------------------------------------------------------------------


def figure_calibration(results: dict, output_path: str | Path) -> None:
    """Generate the confidence calibration reliability diagram.

    Two panels:
      (A) Calibration curve -- predicted confidence vs fraction correct.
      (B) Histogram of confidence values.

    Parameters
    ----------
    results : dict
        Mapping of method name to calibration data.  Each value is a dict
        with keys ``"bins"`` (list of dicts with ``bin_low``, ``bin_high``,
        ``n_pairs``, ``accuracy``), ``"ece"``, ``"brier_score"``.
    output_path : str or Path
        Destination PDF file.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_cal, ax_hist) = plt.subplots(1, 2, figsize=(10, 5))

    # Panel A -- Calibration curve
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
        label = " ".join(label_parts)

        ax_cal.plot(
            midpoints,
            accuracies,
            marker="o",
            color=_method_color(method),
            label=label,
        )

    ax_cal.set_xlabel("Predicted confidence")
    ax_cal.set_ylabel("Observed accuracy")
    ax_cal.set_xlim(0, 1)
    ax_cal.set_ylim(0, 1)
    ax_cal.set_aspect("equal")
    ax_cal.legend(loc="lower right", fontsize=7)
    _add_panel_label(ax_cal, "A")

    # Panel B -- Confidence histogram
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
            midpoints,
            counts,
            width=widths * 0.85,
            alpha=0.5,
            color=_method_color(method),
            label=_method_label(method),
        )

    ax_hist.set_xlabel("Confidence score")
    ax_hist.set_ylabel("Fraction of residue pairs")
    ax_hist.set_xlim(0, 1)
    ax_hist.legend(loc="upper left", fontsize=7)
    _add_panel_label(ax_hist, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 2: Positive selection detection
# ---------------------------------------------------------------------------


def figure_positive_selection(results: list[dict], output_path: str | Path) -> None:
    """Generate the positive selection detection figure.

    Three panels:
      (A) Bar chart of precision, recall, F1 by method.
      (B) F1 vs tree depth, one line per method.
      (C) F1 vs indel rate, one line per method.

    Parameters
    ----------
    results : list[dict]
        Per-case result dicts.  Each must have at least ``method``,
        ``precision``, ``tpr`` (recall), ``f1``, ``tree_depth``,
        ``indel_rate``.
    output_path : str or Path
        Destination PDF file.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, axes = plt.subplots(1, 3, figsize=(10, 5))
    ax_bar, ax_depth, ax_indel = axes

    # Group results by method
    by_method: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        m = r.get("method")
        if m:
            by_method[m].append(r)

    methods = sorted(by_method.keys())

    # Panel A -- precision / recall / F1 bar chart
    if methods:
        n_methods = len(methods)
        x_pos = np.arange(n_methods)
        width = 0.25

        prec_means = []
        rec_means = []
        f1_means = []
        colors = []
        for m in methods:
            cases = by_method[m]
            prec_means.append(np.mean([c.get("precision", 0) for c in cases]))
            rec_means.append(np.mean([c.get("recall", c.get("tpr", 0)) for c in cases]))
            f1_means.append(np.mean([c.get("f1", 0) for c in cases]))
            colors.append(_method_color(m))

        ax_bar.bar(x_pos - width, prec_means, width, label="Precision", color=colors, alpha=0.8)
        ax_bar.bar(x_pos, rec_means, width, label="Recall", color=colors, alpha=0.55)
        ax_bar.bar(x_pos + width, f1_means, width, label="F1", color=colors, alpha=0.3)
        ax_bar.set_xticks(x_pos)
        ax_bar.set_xticklabels(
            [_method_label(m) for m in methods], rotation=45, ha="right", fontsize=7
        )
        ax_bar.set_ylabel("Score")
        # Auto-scale y-axis with a small margin above the tallest bar
        all_vals = prec_means + rec_means + f1_means
        if all_vals:
            ymax = max(all_vals) * 1.15
            ax_bar.set_ylim(0, ymax)
        ax_bar.legend(fontsize=7)
    _add_panel_label(ax_bar, "A")

    # Panel B -- F1 vs tree depth
    has_depth_data = False
    for m in methods:
        cases = by_method[m]
        depth_f1: dict[float, list[float]] = defaultdict(list)
        for c in cases:
            d = c.get("tree_depth")
            if d is None:
                sim_id = c.get("sim_id", "")
                import re as _re
                dm = _re.search(r"_d([\d.]+)_", sim_id)
                if dm:
                    d = float(dm.group(1))
            f1 = c.get("f1")
            if d is not None and f1 is not None:
                depth_f1[d].append(f1)
        if depth_f1:
            has_depth_data = True
            depths_sorted = sorted(depth_f1.keys())
            means = [np.mean(depth_f1[d]) for d in depths_sorted]
            ax_depth.plot(
                depths_sorted,
                means,
                marker="o",
                color=_method_color(m),
                label=_method_label(m),
            )

    ax_depth.set_xlabel("Tree depth (subs/site)")
    ax_depth.set_ylabel("F1 score")
    if not has_depth_data:
        ax_depth.text(0.5, 0.5, "No data", ha="center", va="center",
                      transform=ax_depth.transAxes, fontsize=10, color="#999999")
    _legend_if_labeled(ax_depth, fontsize=6, ncol=2, loc="best")
    _add_panel_label(ax_depth, "B")

    # Panel C -- F1 vs indel rate
    has_indel_data = False
    for m in methods:
        cases = by_method[m]
        indel_f1: dict[float, list[float]] = defaultdict(list)
        for c in cases:
            ir = c.get("indel_rate")
            if ir is None:
                sim_id = c.get("sim_id", "")
                import re as _re
                irm = _re.search(r"_ir([\d.]+)_", sim_id)
                if irm:
                    ir = float(irm.group(1))
            f1 = c.get("f1")
            if ir is not None and f1 is not None:
                indel_f1[ir].append(f1)
        if indel_f1:
            has_indel_data = True
            rates_sorted = sorted(indel_f1.keys())
            means = [np.mean(indel_f1[r]) for r in rates_sorted]
            ax_indel.plot(
                rates_sorted,
                means,
                marker="o",
                color=_method_color(m),
                label=_method_label(m),
            )

    ax_indel.set_xlabel("Indel rate")
    ax_indel.set_ylabel("F1 score")
    if not has_indel_data:
        ax_indel.text(0.5, 0.5, "No data", ha="center", va="center",
                      transform=ax_indel.transAxes, fontsize=10, color="#999999")
    _legend_if_labeled(ax_indel, fontsize=6, ncol=2, loc="best")
    _add_panel_label(ax_indel, "C")

    # Significance markers -- add if pairwise comparisons are present
    _add_significance_markers(ax_bar, results, methods)

    fig.tight_layout()
    _savefig(fig, output_path)


def _add_significance_markers(
    ax: Any, results: list[dict], methods: list[str]
) -> None:
    """Add significance markers (* / ** / ***) above bars if pairwise
    comparison data is embedded in the results.

    Looks for a top-level ``"pairwise_tests"`` entry (list of dicts with
    ``method_a``, ``method_b``, ``p_value``).  Silently skips if absent.
    """
    # Pairwise tests may be stored at the end of the results list as a
    # special sentinel dict, or not present at all.
    pairwise = None
    for r in results:
        if isinstance(r, dict) and "pairwise_tests" in r:
            pairwise = r["pairwise_tests"]
            break
    if not pairwise:
        return

    method_idx = {m: i for i, m in enumerate(methods)}
    y_max = ax.get_ylim()[1]
    offset = 0.03
    for test in pairwise:
        ma = test.get("method_a", "")
        mb = test.get("method_b", "")
        p = test.get("p_value")
        if ma not in method_idx or mb not in method_idx or p is None:
            continue
        if p < 0.001:
            marker = "***"
        elif p < 0.01:
            marker = "**"
        elif p < 0.05:
            marker = "*"
        else:
            continue
        x1 = method_idx[ma]
        x2 = method_idx[mb]
        y = y_max + offset
        ax.plot([x1, x1, x2, x2], [y, y + 0.01, y + 0.01, y], color="black", lw=0.8)
        ax.text(
            (x1 + x2) / 2,
            y + 0.012,
            marker,
            ha="center",
            va="bottom",
            fontsize=8,
        )
        offset += 0.06


# ---------------------------------------------------------------------------
# Figure 3: Phylogenetic tree accuracy
# ---------------------------------------------------------------------------


def figure_phylo_accuracy(results: list[dict], output_path: str | Path) -> None:
    """Generate the phylogenetic tree accuracy figure.

    Three panels:
      (A) nRF distance box plots, one box per method, sorted by median.
      (B) Mean nRF vs tree depth with CI bands, one line per method.
      (C) Branch score distance box plots.

    Parameters
    ----------
    results : list[dict]
        Per-case result dicts with ``method``, ``nrf``,
        ``branch_score_dist``, ``tree_depth``.
    output_path : str or Path
        Destination PDF file.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_nrf, ax_depth, ax_bsd) = plt.subplots(1, 3, figsize=(10, 5))

    by_method: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        m = r.get("method")
        if m:
            by_method[m].append(r)

    # Sort methods by median nRF for panels A and C
    method_medians = []
    for m, cases in by_method.items():
        nrf_vals = [c.get("nrf", float("nan")) for c in cases]
        nrf_vals = [v for v in nrf_vals if v == v]  # drop NaN
        med = float(np.median(nrf_vals)) if nrf_vals else float("inf")
        method_medians.append((m, med))
    method_medians.sort(key=lambda x: x[1])
    sorted_methods = [m for m, _ in method_medians]

    # Panel A -- nRF box plots
    if sorted_methods:
        box_data = []
        colors = []
        for m in sorted_methods:
            nrf_vals = [c.get("nrf", float("nan")) for c in by_method[m]]
            nrf_vals = [v for v in nrf_vals if v == v]
            box_data.append(nrf_vals if nrf_vals else [0])
            colors.append(_method_color(m))

        bp = ax_nrf.boxplot(
            box_data,
            patch_artist=True,
            widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_nrf.set_xticks(range(1, len(sorted_methods) + 1))
        ax_nrf.set_xticklabels(
            [_method_label(m) for m in sorted_methods],
            rotation=45,
            ha="right",
            fontsize=7,
        )
        ax_nrf.set_ylabel("nRF distance")
    _add_panel_label(ax_nrf, "A")

    # Panel B -- nRF vs tree depth with CI bands
    for m in sorted_methods:
        cases = by_method[m]
        depth_nrf: dict[float, list[float]] = defaultdict(list)
        for c in cases:
            d = c.get("tree_depth")
            if d is None:
                # Extract depth from sim_id like "WAG_t16_d0.5_ir0.02_il2.0_r0"
                sim_id = c.get("sim_id", "")
                import re as _re
                dm = _re.search(r"_d([\d.]+)_", sim_id)
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
    _legend_if_labeled(ax_depth, fontsize=6, ncol=2, loc="upper left")
    _add_panel_label(ax_depth, "B")

    # Panel C -- Branch score distance box plots
    if sorted_methods:
        bsd_data = []
        colors = []
        for m in sorted_methods:
            bsd_vals = [c.get("branch_score_dist", float("nan")) for c in by_method[m]]
            bsd_vals = [v for v in bsd_vals if v == v]
            bsd_data.append(bsd_vals if bsd_vals else [0])
            colors.append(_method_color(m))

        bp = ax_bsd.boxplot(
            bsd_data,
            patch_artist=True,
            widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_bsd.set_xticks(range(1, len(sorted_methods) + 1))
        ax_bsd.set_xticklabels(
            [_method_label(m) for m in sorted_methods],
            rotation=45,
            ha="right",
            fontsize=7,
        )
        ax_bsd.set_ylabel("Branch score distance")
    _add_panel_label(ax_bsd, "C")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 4: HMMER homology detection sensitivity
# ---------------------------------------------------------------------------


def figure_hmmer_detection(results: list[dict], output_path: str | Path) -> None:
    """Generate the HMMER detection sensitivity figure.

    Two panels:
      (A) Mean sensitivity per method (bar chart with bootstrap CI).
      (B) Scatter of family-level sensitivity vs family size.

    Parameters
    ----------
    results : list[dict]
        Per-family result dicts with ``method``, ``n_seed_seqs``, and
        at least ``hits_at_1e3`` (dict with ``sens``).
    output_path : str or Path
        Destination PDF file.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_bar, ax_scatter) = plt.subplots(1, 2, figsize=(10, 5))

    by_method: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        m = r.get("method")
        if m:
            by_method[m].append(r)

    methods = sorted(by_method.keys())

    def _sens(case: dict) -> float:
        """Extract sensitivity from a case dict."""
        # Support both flat 'sensitivity' and nested 'hits_at_1e3.sens'
        if "sensitivity" in case:
            return case["sensitivity"]
        h = case.get("hits_at_1e3", {})
        return h.get("sens", 0.0) if isinstance(h, dict) else 0.0

    # Panel A -- mean sensitivity bar chart
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
                from .utils import bootstrap_ci

                lo, hi = bootstrap_ci(sens_vals, n_bootstrap=2000)
                ci_lo.append(mu - lo)
                ci_hi.append(hi - mu)
            else:
                ci_lo.append(0)
                ci_hi.append(0)
            colors.append(_method_color(m))

        x = np.arange(len(methods))
        ax_bar.bar(
            x,
            means,
            color=colors,
            alpha=0.8,
            yerr=[ci_lo, ci_hi],
            capsize=3,
            ecolor="#333333",
        )
        ax_bar.set_xticks(x)
        ax_bar.set_xticklabels(
            [_method_label(m) for m in methods], rotation=45, ha="right", fontsize=7
        )
        ax_bar.set_ylabel("Mean sensitivity (E < 1e-3)")
        # Auto-scale y-axis
        if means and max(means) > 0:
            ymax = max(m + h for m, h in zip(means, ci_hi)) * 1.15
            ax_bar.set_ylim(0, ymax)
        else:
            ax_bar.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax_bar.transAxes, fontsize=10, color="#999999")
    _add_panel_label(ax_bar, "A")

    # Panel B -- sensitivity vs family size scatter
    has_scatter_data = False
    for m in methods:
        cases = by_method[m]
        sizes = [
            c.get("n_seed_seqs", c.get("true_positives", 0) + c.get("false_negatives", 0))
            for c in cases
        ]
        sens = [_sens(c) for c in cases]
        if any(s > 0 for s in sens):
            has_scatter_data = True
        ax_scatter.scatter(
            sizes,
            sens,
            s=15,
            alpha=0.5,
            color=_method_color(m),
            label=_method_label(m),
        )

    ax_scatter.set_xlabel("Number of seed sequences")
    ax_scatter.set_ylabel("Sensitivity (E < 1e-3)")
    if not has_scatter_data:
        ax_scatter.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax_scatter.transAxes, fontsize=10, color="#999999")
    _legend_if_labeled(ax_scatter, fontsize=6, ncol=2, loc="lower right", markerscale=1.5)
    _add_panel_label(ax_scatter, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 5: Speed and memory comparison
# ---------------------------------------------------------------------------


def figure_speed_comparison(results: list[dict], output_path: str | Path) -> None:
    """Generate the speed and memory comparison figure.

    Two panels:
      (A) Wall time box plots (log scale) per method.
      (B) Peak memory box plots per method.

    Parameters
    ----------
    results : list[dict]
        Per-case result dicts with ``method``, ``wall_time_align``,
        ``peak_memory_mb``.
    output_path : str or Path
        Destination PDF file.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    fig, (ax_time, ax_mem) = plt.subplots(1, 2, figsize=(10, 5))

    by_method: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        m = r.get("method")
        if m:
            by_method[m].append(r)

    # Sort methods by median wall time
    method_medians = []
    for m, cases in by_method.items():
        times = [c.get("wall_time_align", float("nan")) for c in cases]
        times = [t for t in times if t == t and t > 0]
        med = float(np.median(times)) if times else float("inf")
        method_medians.append((m, med))
    method_medians.sort(key=lambda x: x[1])
    sorted_methods = [m for m, _ in method_medians]

    # Panel A -- Wall time box plots (log scale)
    if sorted_methods:
        time_data = []
        colors = []
        for m in sorted_methods:
            times = [c.get("wall_time_align", float("nan")) for c in by_method[m]]
            times = [t for t in times if t == t and t > 0]
            time_data.append(times if times else [0.001])
            colors.append(_method_color(m))

        bp = ax_time.boxplot(
            time_data,
            patch_artist=True,
            widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_time.set_yscale("log")
        ax_time.set_xticks(range(1, len(sorted_methods) + 1))
        ax_time.set_xticklabels(
            [_method_label(m) for m in sorted_methods],
            rotation=45,
            ha="right",
            fontsize=7,
        )
        ax_time.set_ylabel("Wall time (s)")

        # Speed ratio annotation
        guidance_times = []
        kalign_ens_times = []
        for m in sorted_methods:
            vals = [c.get("wall_time_align", float("nan")) for c in by_method[m]]
            vals = [t for t in vals if t == t and t > 0]
            if "guidance2" in m and vals:
                guidance_times.extend(vals)
            if m == "kalign_ens3" and vals:
                kalign_ens_times.extend(vals)
        if guidance_times and kalign_ens_times:
            ratio = float(np.median(guidance_times) / np.median(kalign_ens_times))
            ax_time.text(
                0.98,
                0.98,
                f"GUIDANCE2 / kalign_ens3\n{ratio:.0f}x slower",
                transform=ax_time.transAxes,
                ha="right",
                va="top",
                fontsize=7,
                bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.9, "pad": 3},
            )
    _add_panel_label(ax_time, "A")

    # Panel B -- Peak memory box plots
    if sorted_methods:
        mem_data = []
        colors = []
        for m in sorted_methods:
            mems = [c.get("peak_memory_mb", float("nan")) for c in by_method[m]]
            mems = [v for v in mems if v == v and v > 0]
            mem_data.append(mems if mems else [0])
            colors.append(_method_color(m))

        bp = ax_mem.boxplot(
            mem_data,
            patch_artist=True,
            widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_mem.set_xticks(range(1, len(sorted_methods) + 1))
        ax_mem.set_xticklabels(
            [_method_label(m) for m in sorted_methods],
            rotation=45,
            ha="right",
            fontsize=7,
        )
        ax_mem.set_ylabel("Peak memory (MB)")
    _add_panel_label(ax_mem, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 6: Summary heatmap
# ---------------------------------------------------------------------------


def figure_summary_heatmap(all_results: dict, output_path: str | Path) -> None:
    """Generate a summary heatmap across all pipelines.

    Rows = methods, columns = metrics.  Each cell is normalised to [0, 1]
    where 1 = best, 0 = worst.  Annotated with raw values.

    Parameters
    ----------
    all_results : dict
        Mapping with optional keys ``"calibration"``, ``"positive_selection"``,
        ``"phylo_accuracy"``, ``"hmmer_detection"``, ``"speed"``.
        Each value is a dict mapping method names to metric values
        (float).  For ``calibration`` the metric is ECE (lower is better),
        for ``positive_selection`` it is F1 (higher is better), for
        ``phylo_accuracy`` it is mean nRF (lower is better), for
        ``hmmer_detection`` it is mean sensitivity (higher is better), and
        for ``speed`` it is mean wall time in seconds (lower is better).
    output_path : str or Path
        Destination PDF file.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    _setup_figure_style()

    # Columns: metric name -> (method -> raw value), higher_is_better
    columns_spec: list[tuple[str, dict[str, float], bool]] = []

    cal = all_results.get("calibration")
    if cal:
        columns_spec.append(("ECE", cal, False))

    psel = all_results.get("positive_selection")
    if psel:
        columns_spec.append(("Selection F1", psel, True))

    phylo = all_results.get("phylo_accuracy")
    if phylo:
        columns_spec.append(("nRF", phylo, False))

    hmmer = all_results.get("hmmer_detection")
    if hmmer:
        columns_spec.append(("HMMER sens.", hmmer, True))

    speed = all_results.get("speed")
    if speed:
        columns_spec.append(("Speed (s)", speed, False))

    if not columns_spec:
        logger.warning("figure_summary_heatmap: no data provided, skipping.")
        return

    # Collect the union of methods across all columns
    all_methods: set[str] = set()
    for _, mdict, _ in columns_spec:
        all_methods.update(mdict.keys())
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
        valid = ~np.isnan(col)
        if valid.sum() < 2:
            normed[valid, j] = 0.5
            continue
        vmin = np.nanmin(col)
        vmax = np.nanmax(col)
        if vmax - vmin < 1e-12:
            normed[valid, j] = 0.5
        else:
            scaled = (col[valid] - vmin) / (vmax - vmin)
            if not higher_better:
                scaled = 1.0 - scaled
            normed[valid, j] = scaled

    fig, ax = plt.subplots(figsize=(max(5, 1.7 * n_cols + 2), 0.65 * n_methods + 1.5))

    # Use RdYlGn colourmap (green = good, red = bad)
    cmap = plt.get_cmap("RdYlGn")
    display = np.where(np.isnan(normed), 0.5, normed)
    im = ax.imshow(display, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    # Annotate with raw values
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
                j,
                i,
                text,
                ha="center",
                va="center",
                fontsize=7,
                color="black" if 0.3 < display[i, j] < 0.7 else "white",
            )

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([name for name, _, _ in columns_spec], fontsize=9)
    ax.set_yticks(range(n_methods))
    ax.set_yticklabels([_method_label(m) for m in methods], fontsize=8)

    # Colour bar
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Normalised score (1 = best)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------


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


def _paired_wilcoxon(
    results: list[dict],
    method_a: str,
    method_b: str,
    metric: str,
    family_key: str = "family",
    method_key: str = "method",
) -> float | None:
    """Two-sided Wilcoxon signed-rank test on paired family-level values.

    Returns p-value, or None if insufficient data.
    """
    from collections import defaultdict

    by_fam: dict[str, dict[str, float]] = defaultdict(dict)
    for r in results:
        m = r.get(method_key, "")
        fam = r.get(family_key, "")
        val = r.get(metric)
        if m in (method_a, method_b) and fam and val is not None:
            by_fam[fam][m] = val

    pairs_a = []
    pairs_b = []
    for fam, vals in by_fam.items():
        if method_a in vals and method_b in vals:
            pairs_a.append(vals[method_a])
            pairs_b.append(vals[method_b])

    if len(pairs_a) < 10:
        return None
    try:
        from scipy.stats import wilcoxon

        _, p = wilcoxon(pairs_a, pairs_b, alternative="two-sided")
        return float(p)
    except Exception:
        return None


def _sig_marker(p: float | None) -> str:
    """Return significance marker string for a p-value."""
    if p is None:
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def _draw_sig_bracket(
    ax: Any, x1: float, x2: float, y: float, marker: str, fontsize: int = 7
) -> float:
    """Draw a significance bracket between two bar positions.  Returns the
    y-position above the bracket for stacking."""
    if not marker or marker == "n.s.":
        return y
    h = 0.008
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=0.8)
    ax.text((x1 + x2) / 2, y + h + 0.003, marker, ha="center", va="bottom", fontsize=fontsize)
    return y + h + 0.02


_BALIBASE_CATEGORIES = ["RV11", "RV12", "RV20", "RV30", "RV40", "RV50"]


# ---------------------------------------------------------------------------
# Figure 7: BAliBASE alignment accuracy (with CIs and stats)
# ---------------------------------------------------------------------------

# Methods shown in the BAliBASE bar chart — keys into mumsa_precision data
_BAR_METHODS_7 = [
    ("kalign_baseline", "Kalign\n(default)", "#888888"),
    ("consensus_ms4", "Kalign\ncons. s=4", "#6baed6"),
    ("consensus_ms8", "Kalign\ncons. s=8", "#2166ac"),
    ("mafft", "MAFFT", "#e41a1c"),
    ("muscle", "MUSCLE", "#ff7f00"),
    ("clustalo", "Clustal\u03a9", "#984ea3"),
]


def figure_balibase(
    mumsa_results: list[dict],
    output_path: str | Path,
) -> None:
    """BAliBASE alignment accuracy with bootstrap CIs and significance tests.

    Panel A: F1 bar chart for selected methods with 95 % bootstrap CIs and
             Wilcoxon signed-rank significance brackets.
    Panel B: Precision by BAliBASE category (the key result: consensus
             precision beats external tools).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from collections import defaultdict as _dd

    _setup_figure_style()

    fig, (ax_f1, ax_prec) = plt.subplots(1, 2, figsize=(10, 5.5))

    # --- Group by method ---
    by_method: dict[str, list[dict]] = _dd(list)
    for r in mumsa_results:
        by_method[r.get("method", "")].append(r)

    # --- Panel A: F1 with CIs ---
    bar_keys = [(k, lab, col) for k, lab, col in _BAR_METHODS_7 if k in by_method]
    if bar_keys:
        x = np.arange(len(bar_keys))
        means, ci_lo, ci_hi, colors = [], [], [], []
        for key, _lab, col in bar_keys:
            f1s = [r["f1"] for r in by_method[key] if "f1" in r]
            mu, lo, hi = _bootstrap_ci_vec(f1s)
            means.append(mu)
            ci_lo.append(mu - lo)
            ci_hi.append(hi - mu)
            colors.append(col)

        ax_f1.bar(
            x, means, color=colors, alpha=0.85, edgecolor="white", linewidth=0.5,
            yerr=[ci_lo, ci_hi], capsize=3, ecolor="#333333", error_kw={"linewidth": 1},
        )
        ax_f1.set_xticks(x)
        ax_f1.set_xticklabels([lab for _, lab, _ in bar_keys], fontsize=7)
        ax_f1.set_ylabel("Mean F1 score")
        ax_f1.set_title("BAliBASE F1 (218 families)", fontsize=10)
        if means:
            ax_f1.set_ylim(min(means) * 0.88, max(m + h for m, h in zip(means, ci_hi)) * 1.12)

        # Significance tests: compare consensus_ms4 vs each external tool
        ref_method = "consensus_ms4"
        if ref_method in by_method:
            y_top = max(m + h for m, h in zip(means, ci_hi)) * 1.04
            ref_idx = next(i for i, (k, _, _) in enumerate(bar_keys) if k == ref_method)
            for i, (key, _lab, _col) in enumerate(bar_keys):
                if key in ("mafft", "muscle", "clustalo"):
                    p = _paired_wilcoxon(mumsa_results, ref_method, key, "f1")
                    marker = _sig_marker(p)
                    y_top = _draw_sig_bracket(ax_f1, ref_idx, i, y_top, marker)
    _add_panel_label(ax_f1, "A")

    # --- Panel B: Precision by category ---
    by_method_cat: dict[str, dict[str, list[float]]] = _dd(lambda: _dd(list))
    for r in mumsa_results:
        m = r.get("method", "")
        ds = r.get("dataset", "")
        p = r.get("precision")
        if p is None:
            continue
        for cat in _BALIBASE_CATEGORIES:
            if cat in ds:
                by_method_cat[m][cat].append(p)
                break

    prec_methods = [
        ("kalign_baseline", "Kalign (default)", "#888888"),
        ("consensus_ms4", "Cons. s=4", "#6baed6"),
        ("consensus_ms8", "Cons. s=8", "#2166ac"),
        ("mafft", "MAFFT", "#e41a1c"),
        ("muscle", "MUSCLE", "#ff7f00"),
        ("clustalo", "Clustal\u03a9", "#984ea3"),
    ]
    prec_methods = [(k, l, c) for k, l, c in prec_methods if k in by_method_cat]
    n_cats = len(_BALIBASE_CATEGORIES)
    n_bars = len(prec_methods)
    if n_bars > 0:
        width = 0.8 / n_bars
        for j, (key, label, color) in enumerate(prec_methods):
            cat_means = []
            cat_err_lo = []
            cat_err_hi = []
            for cat in _BALIBASE_CATEGORIES:
                vals = by_method_cat[key].get(cat, [])
                if vals:
                    mu, lo, hi = _bootstrap_ci_vec(vals)
                    cat_means.append(mu)
                    cat_err_lo.append(mu - lo)
                    cat_err_hi.append(hi - mu)
                else:
                    cat_means.append(0)
                    cat_err_lo.append(0)
                    cat_err_hi.append(0)
            positions = np.arange(n_cats) + j * width - (n_bars - 1) * width / 2
            ax_prec.bar(
                positions, cat_means, width, label=label, color=color, alpha=0.85,
                yerr=[cat_err_lo, cat_err_hi], capsize=1.5, ecolor="#555555",
                error_kw={"linewidth": 0.6},
            )
        ax_prec.set_xticks(np.arange(n_cats))
        ax_prec.set_xticklabels(_BALIBASE_CATEGORIES, fontsize=8)
        ax_prec.set_ylabel("Precision")
        ax_prec.set_title("Precision by category", fontsize=10)
        ax_prec.legend(fontsize=5.5, ncol=2, loc="lower right")
        ax_prec.set_ylim(0, 1.05)
    _add_panel_label(ax_prec, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 8: BRAliBASE RNA alignment accuracy (with CIs and stats)
# ---------------------------------------------------------------------------

_RNA_FAMILY_ORDER = ["SRP", "tRNA", "rRNA", "g2intron", "U5"]
_RNA_FAMILY_LABELS = {
    "SRP": "SRP",
    "tRNA": "tRNA",
    "rRNA": "rRNA",
    "g2intron": "Group II\nintron",
    "U5": "U5",
}

_BRALI_METHODS = [
    ("python_api", "Kalign\n(default)", "#2ca02c"),
    ("mafft", "MAFFT", "#e41a1c"),
    ("muscle", "MUSCLE", "#ff7f00"),
    ("clustalo", "Clustal\u03a9", "#984ea3"),
]


def figure_bralibase(results: list[dict], output_path: str | Path) -> None:
    """BRAliBASE RNA alignment accuracy with bootstrap CIs and significance.

    Panel A: F1 bar chart by method with 95 % CIs and Wilcoxon tests.
    Panel B: F1 by RNA family, grouped by method, with error bars.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from collections import defaultdict as _dd

    _setup_figure_style()

    fig, (ax_bar, ax_fam) = plt.subplots(1, 2, figsize=(10, 5.5))

    by_method: dict[str, list[dict]] = _dd(list)
    for r in results:
        m = r.get("method", "")
        if r.get("f1") is not None:
            by_method[m].append(r)

    plot_methods = [(k, l, c) for k, l, c in _BRALI_METHODS if k in by_method]

    # --- Panel A: overall F1 with CIs ---
    if plot_methods:
        x = np.arange(len(plot_methods))
        means, err_lo, err_hi, colors = [], [], [], []
        for key, _lab, col in plot_methods:
            f1s = [r["f1"] for r in by_method[key]]
            mu, lo, hi = _bootstrap_ci_vec(f1s)
            means.append(mu)
            err_lo.append(mu - lo)
            err_hi.append(hi - mu)
            colors.append(col)

        ax_bar.bar(
            x, means, color=colors, alpha=0.85, edgecolor="white", linewidth=0.5,
            yerr=[err_lo, err_hi], capsize=3, ecolor="#333333", error_kw={"linewidth": 1},
        )
        ax_bar.set_xticks(x)
        ax_bar.set_xticklabels([lab for _, lab, _ in plot_methods], fontsize=8)
        ax_bar.set_ylabel("Mean F1 score")
        ax_bar.set_title("BRAliBASE RNA", fontsize=10)
        if means:
            ax_bar.set_ylim(min(means) * 0.88, max(m + h for m, h in zip(means, err_hi)) * 1.12)

        # Significance: compare Kalign vs each external tool
        ref_key = "python_api"
        if ref_key in by_method:
            y_top = max(m + h for m, h in zip(means, err_hi)) * 1.04
            ref_idx = next(i for i, (k, _, _) in enumerate(plot_methods) if k == ref_key)
            for i, (key, _lab, _col) in enumerate(plot_methods):
                if key != ref_key:
                    p = _paired_wilcoxon(results, ref_key, key, "f1")
                    marker = _sig_marker(p)
                    y_top = _draw_sig_bracket(ax_bar, ref_idx, i, y_top, marker)
    _add_panel_label(ax_bar, "A")

    # --- Panel B: F1 by family with error bars ---
    by_method_fam: dict[str, dict[str, list[float]]] = _dd(lambda: _dd(list))
    for r in results:
        m = r.get("method", "")
        f1 = r.get("f1")
        ds = r.get("dataset", "")
        if f1 is None:
            continue
        for fam in _RNA_FAMILY_ORDER:
            if fam.lower() in ds.lower() or fam in ds:
                by_method_fam[m][fam].append(f1)
                break

    n_groups = len(_RNA_FAMILY_ORDER)
    n_bars = len(plot_methods)
    if n_bars > 0 and n_groups > 0:
        width = 0.8 / n_bars
        for j, (key, label, color) in enumerate(plot_methods):
            fam_means, fam_elo, fam_ehi = [], [], []
            for fam in _RNA_FAMILY_ORDER:
                vals = by_method_fam[key].get(fam, [])
                if vals:
                    mu, lo, hi = _bootstrap_ci_vec(vals)
                    fam_means.append(mu)
                    fam_elo.append(mu - lo)
                    fam_ehi.append(hi - mu)
                else:
                    fam_means.append(0)
                    fam_elo.append(0)
                    fam_ehi.append(0)
            positions = np.arange(n_groups) + j * width - (n_bars - 1) * width / 2
            ax_fam.bar(
                positions, fam_means, width, label=label, color=color, alpha=0.85,
                yerr=[fam_elo, fam_ehi], capsize=1.5, ecolor="#555555",
                error_kw={"linewidth": 0.6},
            )
        ax_fam.set_xticks(np.arange(n_groups))
        ax_fam.set_xticklabels(
            [_RNA_FAMILY_LABELS.get(f, f) for f in _RNA_FAMILY_ORDER], fontsize=7
        )
        ax_fam.set_ylabel("Mean F1 score")
        ax_fam.set_title("F1 by RNA family", fontsize=10)
        ax_fam.legend(fontsize=6, loc="lower left")
    _add_panel_label(ax_fam, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 9: Ensemble precision advantage (Pareto curve + category precision)
# ---------------------------------------------------------------------------


def figure_ensemble_precision(
    mumsa_results: list[dict],
    output_path: str | Path,
) -> None:
    """Show that ensemble consensus dominates external tools on precision.

    Panel A: Precision-recall Pareto curve — consensus threshold s=1..8 as a
             curve with bootstrap CI bands; external tools as reference points.
    Panel B: Precision by BAliBASE category at s=4 and s=8 vs external tools.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from collections import defaultdict as _dd

    _setup_figure_style()

    fig, (ax_pareto, ax_cat) = plt.subplots(1, 2, figsize=(10, 5.5))

    by_method: dict[str, list[dict]] = _dd(list)
    for r in mumsa_results:
        by_method[r.get("method", "")].append(r)

    # --- Panel A: Precision-recall Pareto curve ---
    consensus_keys = sorted(
        [k for k in by_method if k.startswith("consensus_ms")],
        key=lambda k: int(k.replace("consensus_ms", "")),
    )

    if consensus_keys:
        rec_means, rec_lo, rec_hi = [], [], []
        prec_means, prec_lo, prec_hi = [], [], []
        support_vals = []
        for key in consensus_keys:
            rs = [r["recall"] for r in by_method[key] if "recall" in r]
            ps = [r["precision"] for r in by_method[key] if "precision" in r]
            mu_r, lo_r, hi_r = _bootstrap_ci_vec(rs)
            mu_p, lo_p, hi_p = _bootstrap_ci_vec(ps)
            rec_means.append(mu_r)
            rec_lo.append(lo_r)
            rec_hi.append(hi_r)
            prec_means.append(mu_p)
            prec_lo.append(lo_p)
            prec_hi.append(hi_p)
            support_vals.append(int(key.replace("consensus_ms", "")))

        # Draw the Pareto curve with CI ellipses (approximated as bands)
        ax_pareto.fill_between(
            rec_means, prec_lo, prec_hi, color="#2166ac", alpha=0.12,
        )
        ax_pareto.plot(
            rec_means, prec_means, "o-", color="#2166ac", linewidth=2,
            markersize=6, zorder=5, label="Kalign consensus",
        )
        # Label selected support values
        for i, s in enumerate(support_vals):
            if s in (1, 3, 5, 8):
                ax_pareto.annotate(
                    f"s={s}", (rec_means[i], prec_means[i]),
                    textcoords="offset points", xytext=(6, 4),
                    fontsize=7, color="#2166ac",
                )

    # Reference points for external tools and baselines
    ref_points = [
        ("kalign_baseline", "Kalign (default)", "s", "#888888"),
        ("ensemble", "Ensemble", "D", "#4daf4a"),
        ("mafft", "MAFFT", "^", "#e41a1c"),
        ("muscle", "MUSCLE", "v", "#ff7f00"),
        ("clustalo", "Clustal\u03a9", "P", "#984ea3"),
    ]
    for key, label, marker, color in ref_points:
        if key not in by_method:
            continue
        rs = [r["recall"] for r in by_method[key] if "recall" in r]
        ps = [r["precision"] for r in by_method[key] if "precision" in r]
        mu_r, lo_r, hi_r = _bootstrap_ci_vec(rs)
        mu_p, lo_p, hi_p = _bootstrap_ci_vec(ps)
        ax_pareto.errorbar(
            mu_r, mu_p,
            xerr=[[mu_r - lo_r], [hi_r - mu_r]],
            yerr=[[mu_p - lo_p], [hi_p - mu_p]],
            fmt=marker, color=color, markersize=7, label=label,
            capsize=2, elinewidth=0.8, markeredgewidth=0.8, zorder=4,
        )

    ax_pareto.set_xlabel("Recall (sensitivity)")
    ax_pareto.set_ylabel("Precision")
    ax_pareto.set_title("Precision-recall tradeoff", fontsize=10)
    ax_pareto.legend(fontsize=6, loc="lower left", markerscale=0.8)
    _add_panel_label(ax_pareto, "A")

    # --- Panel B: Precision by category (consensus beats external tools) ---
    by_method_cat: dict[str, dict[str, list[float]]] = _dd(lambda: _dd(list))
    for r in mumsa_results:
        m = r.get("method", "")
        ds = r.get("dataset", "")
        p = r.get("precision")
        if p is None:
            continue
        for cat in _BALIBASE_CATEGORIES:
            if cat in ds:
                by_method_cat[m][cat].append(p)
                break

    cat_methods = [
        ("kalign_baseline", "Kalign (default)", "#888888"),
        ("consensus_ms4", "Cons. s=4", "#6baed6"),
        ("consensus_ms8", "Cons. s=8", "#2166ac"),
        ("mafft", "MAFFT", "#e41a1c"),
        ("muscle", "MUSCLE", "#ff7f00"),
        ("clustalo", "Clustal\u03a9", "#984ea3"),
    ]
    cat_methods = [(k, l, c) for k, l, c in cat_methods if k in by_method_cat]
    n_cats = len(_BALIBASE_CATEGORIES)
    n_bars = len(cat_methods)
    if n_bars > 0:
        width = 0.8 / n_bars
        for j, (key, label, color) in enumerate(cat_methods):
            cat_means, cat_elo, cat_ehi = [], [], []
            for cat in _BALIBASE_CATEGORIES:
                vals = by_method_cat[key].get(cat, [])
                if vals:
                    mu, lo, hi = _bootstrap_ci_vec(vals)
                    cat_means.append(mu)
                    cat_elo.append(mu - lo)
                    cat_ehi.append(hi - mu)
                else:
                    cat_means.append(0)
                    cat_elo.append(0)
                    cat_ehi.append(0)
            positions = np.arange(n_cats) + j * width - (n_bars - 1) * width / 2
            ax_cat.bar(
                positions, cat_means, width, label=label, color=color, alpha=0.85,
                yerr=[cat_elo, cat_ehi], capsize=1.5, ecolor="#555555",
                error_kw={"linewidth": 0.6},
            )
        ax_cat.set_xticks(np.arange(n_cats))
        ax_cat.set_xticklabels(_BALIBASE_CATEGORIES, fontsize=8)
        ax_cat.set_ylabel("Precision")
        ax_cat.set_title("Precision by category", fontsize=10)
        ax_cat.legend(fontsize=5.5, ncol=2, loc="lower right")
        ax_cat.set_ylim(0, 1.05)
    _add_panel_label(ax_cat, "B")

    fig.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Figure 10: BAliBASE SP/TC comparison (standard baliscore-style)
# ---------------------------------------------------------------------------

_SP_TC_METHODS = [
    ("kalign_default", "Kalign", "#1f77b4"),
    ("kalign_refine", "Kalign\n+refine", "#6baed6"),
    ("kalign_best", "Kalign\nbest", "#2ca02c"),
    ("mafft", "MAFFT", "#e41a1c"),
    ("muscle", "MUSCLE", "#ff7f00"),
    ("clustalo", "Clustal\u03a9", "#984ea3"),
]


def figure_balibase_sp_tc(
    results: list[dict],
    output_path: str | Path,
) -> None:
    """BAliBASE SP and TC by category with violin + strip plots.

    Grid: 2 rows (SP, TC) x 6 columns (RV11..RV50).
    Each cell shows violin + jittered dots per method with median/IQR.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns

    _setup_figure_style()

    plot_methods = [(k, l, c) for k, l, c in _SP_TC_METHODS
                    if any(r.get("method") == k for r in results)]
    key_to_label = {k: l.replace("\n", " ") for k, l, _ in plot_methods}
    palette = {key_to_label[k]: c for k, _, c in plot_methods}
    method_order = [key_to_label[k] for k, _, _ in plot_methods]

    # Build long-form DataFrame with category column
    rows = []
    for r in results:
        m = r.get("method", "")
        if m not in key_to_label or "error" in r:
            continue
        ds = r.get("dataset", "")
        cat = None
        for c in _BALIBASE_CATEGORIES:
            if c in ds:
                cat = c
                break
        if cat is None:
            continue
        rows.append({
            "Method": key_to_label[m],
            "Category": cat,
            "SP": r.get("recall", 0),
            "TC": r.get("tc", 0),
            "family": r.get("family", ""),
        })
    df = pd.DataFrame(rows)

    n_cats = len(_BALIBASE_CATEGORIES)
    fig, axes = plt.subplots(2, n_cats, figsize=(3.2 * n_cats, 8),
                             sharey="row")

    metrics = [("SP", "SP score"), ("TC", "TC score")]
    for row, (metric, ylabel) in enumerate(metrics):
        for col, cat in enumerate(_BALIBASE_CATEGORIES):
            ax = axes[row, col]
            cat_df = df[df["Category"] == cat]
            n_fam = cat_df["family"].nunique()

            if len(cat_df) == 0:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9, color="#999")
                continue

            sns.violinplot(
                data=cat_df, x="Method", y=metric, hue="Method",
                order=method_order, hue_order=method_order,
                palette=palette, inner=None, linewidth=0,
                saturation=0.3, cut=0, legend=False, ax=ax,
            )
            sns.stripplot(
                data=cat_df, x="Method", y=metric, hue="Method",
                order=method_order, hue_order=method_order,
                palette=palette, size=2.5, alpha=0.5,
                jitter=0.25, legend=False, ax=ax,
            )

            # Median + IQR whiskers
            for i, label in enumerate(method_order):
                vals = cat_df.loc[cat_df["Method"] == label, metric].values
                if len(vals) == 0:
                    continue
                q25, med, q75 = np.percentile(vals, [25, 50, 75])
                ax.hlines(med, i - 0.22, i + 0.22,
                          color="black", linewidth=1.5, zorder=5)
                ax.vlines(i, q25, q75, color="black", linewidth=1.0, zorder=4)

            ax.set_ylim(-0.02, 1.05)
            ax.set_xlabel("")
            # Only show x-tick labels on bottom row
            ticks = list(range(len(method_order)))
            ax.set_xticks(ticks)
            if row == 0:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(
                    [l.replace(" ", "\n") for l in method_order],
                    fontsize=5.5, rotation=45, ha="right",
                )
            # y-label only on leftmost column
            if col == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.set_ylabel("")
            # Category title on top row
            if row == 0:
                ax.set_title(f"{cat} (n={n_fam})", fontsize=10, fontweight="bold")

    # Shared legend at the bottom
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=palette[l], label=l) for l in method_order]
    fig.legend(handles=handles, loc="lower center", ncol=len(method_order),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.01))

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    _savefig(fig, output_path)


def _load_balibase_sp_data(results_dir: Path) -> list[dict] | None:
    """Load BAliBASE SP/TC data from baseline, best-mode, and external results.

    Combines:
      - balibase_baseline.json (kalign default + refine)
      - balibase_best.json (kalign best mode)
      - full_comparison.json or mumsa_precision.json (external tools)
    """
    results_dir = Path(results_dir)
    all_results = []

    # Kalign baseline and refine
    baseline_path = results_dir / "balibase_baseline.json"
    if baseline_path.exists():
        with open(baseline_path) as fh:
            data = json.load(fh)
        raw = data.get("results", data) if isinstance(data, dict) else data
        for r in raw:
            if isinstance(r, dict) and r.get("method") == "python_api":
                refine = r.get("refine", "none")
                method = "kalign_refine" if refine == "confident" else "kalign_default"
                all_results.append({
                    "family": r.get("family", ""),
                    "dataset": r.get("dataset", ""),
                    "method": method,
                    "recall": r.get("recall", r.get("sp_score", 0)),
                    "precision": r.get("precision", 0),
                    "f1": r.get("f1", 0),
                    "tc": r.get("tc", 0),
                })

    # Kalign best mode
    best_path = results_dir / "balibase_best.json"
    if best_path.exists():
        with open(best_path) as fh:
            best_data = json.load(fh)
        raw = best_data if isinstance(best_data, list) else best_data.get("results", [])
        for r in raw:
            if isinstance(r, dict):
                all_results.append({
                    "family": r.get("family", ""),
                    "dataset": r.get("dataset", ""),
                    "method": "kalign_best",
                    "recall": r.get("recall", 0),
                    "precision": r.get("precision", 0),
                    "f1": r.get("f1", 0),
                    "tc": r.get("tc", 0),
                })

    # External tools from mumsa_precision.json (most up-to-date)
    mumsa_path = results_dir / "mumsa_precision.json"
    if mumsa_path.exists():
        with open(mumsa_path) as fh:
            mumsa_data = json.load(fh)
        raw = mumsa_data.get("results", []) if isinstance(mumsa_data, dict) else mumsa_data
        for r in raw:
            if isinstance(r, dict) and r.get("method") in ("mafft", "muscle", "clustalo"):
                all_results.append({
                    "family": r.get("family", ""),
                    "dataset": r.get("dataset", ""),
                    "method": r["method"],
                    "recall": r.get("recall", 0),
                    "precision": r.get("precision", 0),
                    "f1": r.get("f1", 0),
                    "tc": r.get("tc", 0),
                })
    elif (results_dir / "full_comparison.json").exists():
        with open(results_dir / "full_comparison.json") as fh:
            fc_data = json.load(fh)
        for r in fc_data.get("results", []):
            if r.get("method") in ("mafft", "muscle", "clustalo"):
                all_results.append({
                    "family": r.get("family", ""),
                    "dataset": r.get("dataset", ""),
                    "method": r["method"],
                    "recall": r.get("recall", r.get("sp_score", 0)),
                    "precision": r.get("precision", 0),
                    "f1": r.get("f1", 0),
                    "tc": r.get("tc", 0),
                })

    return all_results if all_results else None


# ---------------------------------------------------------------------------
# Result loading helpers
# ---------------------------------------------------------------------------


def _load_latest_json(results_dir: Path) -> Any:
    """Load the ``latest.json`` symlink (or most recent run file) from a
    results directory.

    Returns the parsed JSON data, or ``None`` if the directory is empty
    or missing.
    """
    results_dir = Path(results_dir)
    if not results_dir.is_dir():
        logger.warning("Results directory does not exist: %s", results_dir)
        return None

    latest = results_dir / "latest.json"
    if latest.exists():
        with open(latest) as fh:
            return json.load(fh)

    # Fallback: pick the most recent JSON file by name
    json_files = sorted(results_dir.glob("run_*.json"))
    if not json_files:
        logger.warning("No result files found in %s", results_dir)
        return None

    with open(json_files[-1]) as fh:
        return json.load(fh)


def _aggregate_for_heatmap(
    pipeline_results: list[dict],
    method_key: str,
    metric_key: str,
) -> dict[str, float]:
    """Compute mean of *metric_key* per method from a list of per-case dicts.

    Handles nested metric keys like ``"hits_at_1e3.sens"`` by splitting on
    the dot.
    """
    import numpy as np  # noqa: lazy

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
    """Convert current calibration summary format to the bins format
    expected by figure_calibration.

    Current format per method::

        {"ece": float, "brier_score": float,
         "calibration_curve": {"bin_centers": [...], "bin_counts": [...],
                               "fraction_correct": [...]}}

    Expected format per method::

        {"ece": float, "brier_score": float,
         "bins": [{"bin_low": ..., "bin_high": ..., "n_pairs": ...,
                   "accuracy": ...}, ...]}
    """
    import math

    adapted: dict = {}
    for method, data in cal_methods.items():
        if not isinstance(data, dict):
            continue
        cc = data.get("calibration_curve")
        if not isinstance(cc, dict):
            # Already in bins format or no curve data
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
            # Infer bin width from spacing (default 0.1)
            if len(centers) > 1:
                if i == 0:
                    width = centers[1] - centers[0]
                else:
                    width = centers[i] - centers[i - 1]
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

    Loads ``latest.json`` from each pipeline's results sub-directory and
    calls the corresponding figure function.  Skips any pipeline whose
    results are not available.

    Parameters
    ----------
    results_dir : str or Path
        Root results directory (e.g. ``benchmarks/results``).
    output_dir : str or Path
        Directory where PDF figures will be written.
    """
    import numpy as np  # noqa: lazy

    results_dir = Path(results_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated: list[str] = []

    # --- Calibration (Figure 1) -------------------------------------------
    cal_data = _load_latest_json(results_dir / "calibration")
    if cal_data is not None:
        # cal_data may be:
        #  - {"methods": {method -> cal_info}} (legacy)
        #  - {"summary": {method -> cal_info}, "cases": [...]} (current)
        #  - {method -> cal_info} (bare)
        if isinstance(cal_data, dict) and "methods" in cal_data:
            cal_methods = cal_data["methods"]
        elif isinstance(cal_data, dict) and "summary" in cal_data:
            cal_methods = cal_data["summary"]
        elif isinstance(cal_data, dict):
            cal_methods = cal_data
        else:
            cal_methods = None

        # Adapt current summary format to what figure_calibration expects
        if isinstance(cal_methods, dict):
            cal_methods = _adapt_calibration_summary(cal_methods)

        if cal_methods:
            path = output_dir / "fig1_calibration.pdf"
            try:
                figure_calibration(cal_methods, path)
                generated.append(str(path))
            except Exception as exc:
                logger.error("Failed to generate calibration figure: %s", exc)

    # --- Positive selection (Figure 2) ------------------------------------
    psel_data = _load_latest_json(results_dir / "positive_selection")
    if psel_data is not None:
        cases = psel_data if isinstance(psel_data, list) else psel_data.get("cases", [])
        if cases:
            path = output_dir / "fig2_positive_selection.pdf"
            try:
                figure_positive_selection(cases, path)
                generated.append(str(path))
            except Exception as exc:
                logger.error("Failed to generate positive selection figure: %s", exc)

    # --- Phylogenetic accuracy (Figure 3) ---------------------------------
    phylo_data = _load_latest_json(results_dir / "phylo_accuracy")
    if phylo_data is not None:
        cases = phylo_data if isinstance(phylo_data, list) else phylo_data.get("cases", [])
        if cases:
            path = output_dir / "fig3_phylo_accuracy.pdf"
            try:
                figure_phylo_accuracy(cases, path)
                generated.append(str(path))
            except Exception as exc:
                logger.error("Failed to generate phylo accuracy figure: %s", exc)

    # --- HMMER detection (Figure 4) ---------------------------------------
    hmmer_data = _load_latest_json(results_dir / "hmmer_detection")
    if hmmer_data is not None:
        cases = (
            hmmer_data if isinstance(hmmer_data, list) else hmmer_data.get("cases", [])
        )
        if cases:
            path = output_dir / "fig4_hmmer_detection.pdf"
            try:
                figure_hmmer_detection(cases, path)
                generated.append(str(path))
            except Exception as exc:
                logger.error("Failed to generate HMMER detection figure: %s", exc)

    # --- Speed comparison (Figure 5) --------------------------------------
    # Combine all per-case results that have timing data
    all_timed: list[dict] = []
    for subdir in ("positive_selection", "phylo_accuracy", "hmmer_detection"):
        sub_data = _load_latest_json(results_dir / subdir)
        if sub_data is not None:
            sub_cases = (
                sub_data if isinstance(sub_data, list) else sub_data.get("cases", [])
            )
            all_timed.extend(
                c for c in sub_cases if "wall_time_align" in c
            )
    if all_timed:
        path = output_dir / "fig5_speed_comparison.pdf"
        try:
            figure_speed_comparison(all_timed, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate speed comparison figure: %s", exc)

    # --- Summary heatmap (Figure 6) ---------------------------------------
    heatmap_data: dict[str, dict[str, float]] = {}

    # Calibration ECE
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

    # Positive selection F1
    if psel_data is not None:
        psel_cases = (
            psel_data if isinstance(psel_data, list) else psel_data.get("cases", [])
        )
        agg = _aggregate_for_heatmap(psel_cases, "method", "f1")
        if agg:
            heatmap_data["positive_selection"] = agg

    # Phylo nRF
    if phylo_data is not None:
        phylo_cases = (
            phylo_data if isinstance(phylo_data, list) else phylo_data.get("cases", [])
        )
        agg = _aggregate_for_heatmap(phylo_cases, "method", "nrf")
        if agg:
            heatmap_data["phylo_accuracy"] = agg

    # HMMER sensitivity
    if hmmer_data is not None:
        hmmer_cases = (
            hmmer_data if isinstance(hmmer_data, list) else hmmer_data.get("cases", [])
        )
        agg = _aggregate_for_heatmap(hmmer_cases, "method", "sensitivity")
        if agg:
            heatmap_data["hmmer_detection"] = agg

    # Speed (mean wall time)
    if all_timed:
        agg = _aggregate_for_heatmap(all_timed, "method", "wall_time_align")
        if agg:
            heatmap_data["speed"] = agg

    if heatmap_data:
        path = output_dir / "fig6_summary_heatmap.pdf"
        try:
            figure_summary_heatmap(heatmap_data, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate summary heatmap: %s", exc)

    # --- BAliBASE accuracy (Figure 7) ------------------------------------
    mumsa_path = results_dir / "mumsa_precision.json"
    mumsa_cases = None
    if mumsa_path.exists():
        with open(mumsa_path) as fh:
            mumsa_data = json.load(fh)
        mumsa_cases = (
            mumsa_data.get("results", mumsa_data)
            if isinstance(mumsa_data, dict) else mumsa_data
        )
        path = output_dir / "fig7_balibase_accuracy.pdf"
        try:
            figure_balibase(mumsa_cases, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate BAliBASE figure: %s", exc)

    # --- BRAliBASE RNA accuracy (Figure 8) --------------------------------
    brali_path = results_dir / "bralibase_comparison.json"
    if brali_path.exists():
        with open(brali_path) as fh:
            brali_data = json.load(fh)
        brali_cases = (
            brali_data.get("results", brali_data) if isinstance(brali_data, dict) else brali_data
        )
        path = output_dir / "fig8_bralibase_rna.pdf"
        try:
            figure_bralibase(brali_cases, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate BRAliBASE figure: %s", exc)

    # --- Ensemble precision advantage (Figure 9) --------------------------
    if mumsa_cases is not None:
        path = output_dir / "fig9_ensemble_precision.pdf"
        try:
            figure_ensemble_precision(mumsa_cases, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate ensemble precision figure: %s", exc)

    # --- BAliBASE SP/TC comparison (Figure 10) ----------------------------
    balibase_sp_cases = _load_balibase_sp_data(results_dir)
    if balibase_sp_cases:
        path = output_dir / "fig10_balibase_sp_tc.pdf"
        try:
            figure_balibase_sp_tc(balibase_sp_cases, path)
            generated.append(str(path))
        except Exception as exc:
            logger.error("Failed to generate BAliBASE SP/TC figure: %s", exc)

    # --- Summary ----------------------------------------------------------
    if generated:
        print(f"Generated {len(generated)} figures:")
        for p in generated:
            print(f"  {p}")
    else:
        print("No figures generated (no result data found).")
