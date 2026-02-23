"""Generate publication figures for MUMSA consensus precision analysis.

Reads results from benchmarks/results/mumsa_precision.json (produced by
mumsa_precision.py) and generates two figures:

  Figure 1: Precision-recall tradeoff curve showing consensus support
            thresholds vs single-tool baselines.
  Figure 2: Per-category precision bar chart at selected thresholds
            alongside external tools.

Usage:
    python -m benchmarks.mumsa_plots [--output-dir figures/]

Only requires matplotlib + json (no kalign dependency).
"""

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Consistent colours
COLORS = {
    "consensus": "#2166ac",
    "kalign_baseline": "#888888",
    "ensemble": "#4daf4a",
    "mafft": "#e41a1c",
    "muscle": "#ff7f00",
    "clustalo": "#984ea3",
}

TOOL_LABELS = {
    "kalign_baseline": "kalign",
    "ensemble": "kalign ensemble",
    "mafft": "MAFFT",
    "muscle": "MUSCLE",
    "clustalo": "Clustal Omega",
}


def load_results(json_path):
    data = json.load(open(json_path))
    return data["results"], data.get("n_runs", 8)


def aggregate(results, group_key="method"):
    """Group results by method, return {method: {metric: mean_value}}."""
    groups = defaultdict(list)
    for r in results:
        groups[r[group_key]].append(r)
    agg = {}
    for method, entries in groups.items():
        agg[method] = {
            "recall": statistics.mean(e["recall"] for e in entries),
            "precision": statistics.mean(e["precision"] for e in entries),
            "f1": statistics.mean(e["f1"] for e in entries),
            "tc": statistics.mean(e["tc"] for e in entries),
            "n": len(entries),
        }
    return agg


def aggregate_by_category(results):
    """Group results by (category, method)."""
    groups = defaultdict(list)
    for r in results:
        cat = r["dataset"].replace("balibase_", "")
        groups[(cat, r["method"])].append(r)
    agg = {}
    for (cat, method), entries in groups.items():
        agg[(cat, method)] = {
            "recall": statistics.mean(e["recall"] for e in entries),
            "precision": statistics.mean(e["precision"] for e in entries),
            "f1": statistics.mean(e["f1"] for e in entries),
            "tc": statistics.mean(e["tc"] for e in entries),
            "n": len(entries),
        }
    return agg


def figure_precision_vs_support(results, n_runs, output_path):
    """Precision and recall as a function of consensus support threshold,
    with horizontal reference lines for external tools and kalign baseline."""
    agg = aggregate(results)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    # Consensus curves
    consensus_methods = sorted(
        [m for m in agg if m.startswith("consensus_ms")],
        key=lambda m: int(m.replace("consensus_ms", ""))
    )
    support_vals = [int(m.replace("consensus_ms", "")) for m in consensus_methods]
    prec = [agg[m]["precision"] for m in consensus_methods]
    rec = [agg[m]["recall"] for m in consensus_methods]

    ax.plot(support_vals, prec, "o-", color=COLORS["consensus"], linewidth=2.5,
            markersize=7, zorder=5, label="Consensus precision")
    ax.plot(support_vals, rec, "s--", color=COLORS["consensus"], linewidth=1.5,
            markersize=5, alpha=0.5, zorder=4, label="Consensus recall")

    # Reference lines for external tools and baselines (precision only)
    ref_methods = [
        ("kalign_baseline", COLORS["kalign_baseline"], "kalign"),
        ("ensemble", COLORS["ensemble"], "kalign ensemble"),
        ("mafft", COLORS["mafft"], "MAFFT"),
        ("muscle", COLORS["muscle"], "MUSCLE"),
        ("clustalo", COLORS["clustalo"], "Clustal Omega"),
    ]
    for method, color, label in ref_methods:
        if method in agg:
            ax.axhline(agg[method]["precision"], color=color, linestyle=":",
                       linewidth=1.5, alpha=0.8, label=f"{label} precision")

    ax.set_xlabel("Minimum support threshold (s)", fontsize=11)
    ax.set_ylabel("Score", fontsize=11)
    ax.set_title("Consensus precision increases with support threshold", fontsize=11)
    ax.set_xticks(support_vals)
    ax.legend(fontsize=7.5, loc="center left", bbox_to_anchor=(1.02, 0.5),
              framealpha=0.9)
    ax.set_ylim(0.35, 1.0)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")


def figure_precision_by_category(results, n_runs, output_path):
    """Per-category grouped bar chart of precision at selected thresholds."""
    agg = aggregate_by_category(results)

    categories = sorted({r["dataset"].replace("balibase_", "") for r in results})

    # Methods to show: baseline, ms=3, ms=n_runs, ensemble, mafft, muscle, clustalo
    mid = max(1, n_runs // 2)
    bar_methods = [
        ("kalign_baseline", "kalign", COLORS["kalign_baseline"]),
        (f"consensus_ms{mid}", f"consensus s={mid}", "#6baed6"),
        (f"consensus_ms{n_runs}", f"consensus s={n_runs}", COLORS["consensus"]),
        ("ensemble", "ensemble", COLORS["ensemble"]),
        ("mafft", "MAFFT", COLORS["mafft"]),
        ("muscle", "MUSCLE", COLORS["muscle"]),
        ("clustalo", "Clustal\u03a9", COLORS["clustalo"]),
    ]

    n_cats = len(categories)
    n_bars = len(bar_methods)
    bar_width = 0.11
    x = range(n_cats)

    fig, ax = plt.subplots(figsize=(9, 4.5))

    for bi, (method, label, color) in enumerate(bar_methods):
        vals = []
        for cat in categories:
            key = (cat, method)
            vals.append(agg[key]["precision"] if key in agg else 0)
        offsets = [xi + (bi - n_bars / 2 + 0.5) * bar_width for xi in x]
        ax.bar(offsets, vals, bar_width, label=label, color=color,
               edgecolor="white", linewidth=0.3)

    ax.set_xlabel("BAliBASE category", fontsize=11)
    ax.set_ylabel("Precision", fontsize=11)
    ax.set_title("Precision by category and method", fontsize=11)
    ax.set_xticks(list(x))
    ax.set_xticklabels(categories, fontsize=10)
    ax.legend(fontsize=7.5, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.12), framealpha=0.9)
    ax.set_ylim(0, 1.05)
    ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate MUMSA precision figures from saved results."
    )
    parser.add_argument(
        "--input", default="benchmarks/results/mumsa_precision.json",
        help="Path to results JSON (default: benchmarks/results/mumsa_precision.json)"
    )
    parser.add_argument(
        "--output-dir", default="benchmarks/figures",
        help="Directory for output figures (default: benchmarks/figures/)"
    )
    parser.add_argument(
        "--format", default="pdf", choices=["pdf", "png", "svg"],
        help="Output format (default: pdf)"
    )
    args = parser.parse_args()

    results, n_runs = load_results(args.input)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    figure_precision_vs_support(
        results, n_runs,
        out_dir / f"mumsa_precision_vs_support.{args.format}",
    )
    figure_precision_by_category(
        results, n_runs,
        out_dir / f"mumsa_precision_by_category.{args.format}",
    )


if __name__ == "__main__":
    main()
