#!/usr/bin/env python3
"""Interactive Dash app to visualize the Pareto front from a unified optimizer checkpoint.

Usage:
    # View local checkpoint
    uv run python -m benchmarks.view_pareto benchmarks/results/unified_optim/gen_checkpoint.pkl

    # Pull from server and view (auto-refresh every 30s)
    uv run python -m benchmarks.view_pareto --remote tki-workstation:tmp/kalign35/benchmarks/results/unified_optim/gen_checkpoint.pkl

    # Custom port and refresh interval
    uv run python -m benchmarks.view_pareto --port 8051 --refresh 60 checkpoint.pkl
"""

import argparse
import json
import os
import pickle
import subprocess
import time
from pathlib import Path

# Disable all Dash/Plotly telemetry before importing
os.environ["DASH_DISABLE_TELEMETRY"] = "1"
os.environ["PLOTLY_RENDERER"] = "browser"

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, ctx, dash_table, dcc, html

try:
    from kneed import KneeLocator
except ImportError:
    KneeLocator = None

from .optimize_unified import (
    MATRIX_NAMES,
    REFINE_LONG,
    decode_unified_params,
    mode_label,
)


def load_checkpoint(path: str):
    """Load checkpoint, optionally from remote via rsync."""
    with open(path, "rb") as f:
        return pickle.load(f)  # noqa: S301


def sync_remote(remote_path: str, local_path: str):
    """Pull checkpoint from remote server."""
    result = subprocess.run(
        ["rsync", "-az", remote_path, local_path],
        capture_output=True, text=True,
    )
    return result.returncode == 0


def build_pareto_df(ckpt: dict, max_runs: int) -> tuple[pd.DataFrame, dict]:
    """Build a DataFrame from checkpoint population.

    Returns (DataFrame, params_by_idx) where params_by_idx maps idx -> full decoded params.
    Supports mixed_v1 and mixed_v2 (dict-per-individual) checkpoints.
    """
    if ckpt.get("format") not in ("mixed_v1", "mixed_v2"):
        raise ValueError(
            "Old-format checkpoint (float arrays). "
            "Re-run the optimizer to generate a mixed_v1/v2 checkpoint."
        )
    pop_X = ckpt["pop_X"]
    pop_F = ckpt["pop_F"]

    rows = []
    params_by_idx = {}
    for i in range(len(pop_X)):
        params = decode_unified_params(pop_X[i], max_runs)
        params_by_idx[i] = params
        f1 = -pop_F[i, 0]
        tc = -pop_F[i, 1]
        wt = pop_F[i, 2]

        mode = mode_label(params)

        # Per-run details
        run_details = []
        for k in range(params["n_runs"]):
            mat = MATRIX_NAMES.get(params["run_types"][k], "?")
            ref = REFINE_LONG.get(params["run_refine"][k], "?")
            noise = params["run_noise"][k]
            noise_str = f" n={noise:.2f}" if noise > 0 else ""
            run_details.append(
                f"R{k}: gpo={params['run_gpo'][k]:.2f} "
                f"gpe={params['run_gpe'][k]:.2f} "
                f"tgpe={params['run_tgpe'][k]:.2f}{noise_str} {mat} "
                f"vsm={params['run_vsm_amax'][k]:.2f} ref={ref}"
            )

        # Summarize per-run refine/vsm for the table
        refs = "/".join(REFINE_LONG.get(r, "?") for r in params["run_refine"])

        # Per-run summaries for hover
        vsm_summary = "/".join(f"{v:.1f}" for v in params["run_vsm_amax"])
        mat_summary = "/".join(
            MATRIX_NAMES.get(params["run_types"][k], "?")
            for k in range(params["n_runs"])
        )

        rows.append({
            "idx": i,
            "f1": round(f1, 4),
            "tc": round(tc, 4),
            "wall_time": round(wt, 1),
            "mode": mode,
            "n_runs": params["n_runs"],
            "vsm": vsm_summary,
            "refine": refs,
            "matrices": mat_summary,
            "seq_weights": round(params["seq_weights"], 3),
            "consistency": params["consistency"],
            "consistency_weight": round(params["consistency_weight"], 3),
            "realign": params["realign"],
            "min_support": params["min_support"],
            "run_details": "\n".join(run_details),
            # Keep run-0 values for coloring
            "vsm_amax_0": round(params["run_vsm_amax"][0], 3),
            "gpo_0": round(params["run_gpo"][0], 2),
            "gpe_0": round(params["run_gpe"][0], 2),
            "tgpe_0": round(params["run_tgpe"][0], 2),
            "matrix_0": MATRIX_NAMES.get(params["run_types"][0], "?"),
        })

    return pd.DataFrame(rows), params_by_idx


def build_history_df(ckpt: dict, max_runs: int) -> pd.DataFrame:
    """Build DataFrame from evaluation history."""
    history = ckpt.get("history", [])
    if not history:
        return pd.DataFrame()

    rows = []
    for h in history:
        cv = h.get("cv_result", {})
        p = h.get("params", {})
        rows.append({
            "eval": h.get("eval", 0),
            "f1": cv.get("f1", 0),
            "tc": cv.get("tc", 0),
            "wall_time": cv.get("wall_time", 0),
            "mode": mode_label(p) if p else "?",
        })
    return pd.DataFrame(rows)


def _compute_pareto_2d(df: pd.DataFrame, x_col: str, y_col: str) -> pd.DataFrame:
    """Extract the 2D Pareto front from a DataFrame.

    Assumes we want to maximize y_col and minimize x_col (if x is time)
    or maximize both (if both are scores).
    Returns points sorted by x_col for drawing as a line.
    """
    if df.empty:
        return df

    # Determine direction: for wall_time we minimize, for f1/tc we maximize
    # The Pareto front is the "upper-left" boundary (max y, min x)
    # or "upper-right" if both axes are scores to maximize
    x_minimize = (x_col == "wall_time")
    y_minimize = (y_col == "wall_time")

    points = df[[x_col, y_col]].values
    n = len(points)
    is_pareto = [True] * n

    for i in range(n):
        if not is_pareto[i]:
            continue
        for j in range(n):
            if i == j or not is_pareto[j]:
                continue
            # Check if j dominates i
            xi, yi = points[i]
            xj, yj = points[j]

            # "Better" depends on direction
            xj_better = (xj < xi) if x_minimize else (xj > xi)
            xj_equal = abs(xj - xi) < 1e-10
            yj_better = (yj < yi) if y_minimize else (yj > yi)
            yj_equal = abs(yj - yi) < 1e-10

            if (xj_better or xj_equal) and (yj_better or yj_equal) and (xj_better or yj_better):
                is_pareto[i] = False
                break

    pareto = df[is_pareto].copy()
    pareto = pareto.sort_values(x_col)
    return pareto


def _compute_pareto_3d(df: pd.DataFrame) -> pd.DataFrame:
    """Extract the 3D Pareto front: maximize F1, maximize TC, minimize wall_time.

    A point is non-dominated if no other point is better in all three objectives.
    """
    if df.empty:
        return df

    n = len(df)
    f1 = df["f1"].values
    tc = df["tc"].values
    wt = df["wall_time"].values
    is_pareto = [True] * n

    for i in range(n):
        if not is_pareto[i]:
            continue
        for j in range(n):
            if i == j or not is_pareto[j]:
                continue
            # j dominates i if: f1_j >= f1_i AND tc_j >= tc_i AND wt_j <= wt_i
            # with at least one strict inequality
            if (f1[j] >= f1[i] and tc[j] >= tc[i] and wt[j] <= wt[i] and
                    (f1[j] > f1[i] or tc[j] > tc[i] or wt[j] < wt[i])):
                is_pareto[i] = False
                break

    return df[is_pareto].copy()


def auto_select_tiers(pareto_df: pd.DataFrame) -> dict:
    """Auto-select fast/default/accurate tiers from the 2D Pareto front (F1 vs time).

    Uses the kneedle algorithm to find knee points on the F1-vs-time Pareto curve.
    The knee point is where spending more time gives diminishing F1 returns.

    Returns dict with keys "fast", "default", "accurate", each containing
    the DataFrame row index (idx column) or None if not enough data.
    """
    if pareto_df.empty or len(pareto_df) < 3:
        return {"fast": None, "default": None, "accurate": None}

    # Get the 2D Pareto front: maximize F1, minimize time
    front = _compute_pareto_2d(pareto_df, "wall_time", "f1")
    if len(front) < 3:
        # Not enough points for knee detection
        front_sorted = front.sort_values("wall_time")
        return {
            "fast": int(front_sorted.iloc[0]["idx"]),
            "default": int(front_sorted.iloc[len(front_sorted) // 2]["idx"]),
            "accurate": int(front_sorted.iloc[-1]["idx"]),
        }

    # Sort by time (ascending) for kneedle
    front = front.sort_values("wall_time").reset_index(drop=True)
    x = front["wall_time"].values
    y = front["f1"].values

    # Accurate: best F1 on the front
    accurate_iloc = int(np.argmax(y))
    accurate_idx = int(front.iloc[accurate_iloc]["idx"])

    # Fast: fastest Pareto-optimal point (lowest time, already sorted)
    fast_iloc = 0
    fast_idx = int(front.iloc[0]["idx"])

    # Default: use kneedle to find the knee point (diminishing returns)
    default_idx = None
    if KneeLocator is not None and len(x) >= 3:
        try:
            kn = KneeLocator(
                x, y,
                curve="concave",
                direction="increasing",
                S=1.0,
            )
            if kn.knee is not None:
                # Find the closest front point to the knee
                knee_iloc = int(np.argmin(np.abs(x - kn.knee)))
                default_idx = int(front.iloc[knee_iloc]["idx"])
        except Exception:
            pass

    # Fallback: pick the point closest to 2/3 of the time range
    if default_idx is None:
        target_time = x[0] + 0.5 * (x[-1] - x[0])
        default_iloc = int(np.argmin(np.abs(x - target_time)))
        default_idx = int(front.iloc[default_iloc]["idx"])

    # Ensure fast < default < accurate in time and all different
    # Find the iloc positions for ordering checks
    default_iloc = int(front[front["idx"] == default_idx].index[0]) if default_idx is not None else None
    if default_iloc is not None:
        if default_iloc <= fast_iloc and fast_iloc + 1 < len(front):
            default_iloc = fast_iloc + 1
            default_idx = int(front.iloc[default_iloc]["idx"])
        if default_iloc >= accurate_iloc and accurate_iloc > 0:
            default_iloc = accurate_iloc - 1
            default_idx = int(front.iloc[default_iloc]["idx"])
    if default_idx == fast_idx:
        if fast_iloc + 1 < len(front):
            default_idx = int(front.iloc[fast_iloc + 1]["idx"])
    if default_idx == accurate_idx:
        if accurate_iloc > 0:
            default_idx = int(front.iloc[accurate_iloc - 1]["idx"])

    return {"fast": fast_idx, "default": default_idx, "accurate": accurate_idx}


def _format_tier_config(row) -> str:
    """Format a tier selection as copy-pasteable Python code."""
    lines = [
        f"# F1={row['f1']:.4f}  TC={row['tc']:.4f}  Time={row['wall_time']:.1f}s",
        f"# Mode: {row['mode']}",
    ]
    # Per-run params
    if row.get("run_details"):
        for line in row["run_details"].split("\n"):
            lines.append(f"# {line}")
    lines.append(f"config = {{")
    lines.append(f'    "n_runs": {row["n_runs"]},')
    lines.append(f'    "vsm": "{row["vsm"]}",  # per-run')
    lines.append(f'    "seq_weights": {row["seq_weights"]},')
    lines.append(f'    "consistency": {row["consistency"]},')
    lines.append(f'    "consistency_weight": {row["consistency_weight"]},')
    lines.append(f'    "realign": {row["realign"]},')
    lines.append(f'    "refine": "{row["refine"]}",')
    lines.append(f'    "min_support": {row["min_support"]},')
    lines.append(f"}}")
    return "\n".join(lines)


def create_app(ckpt_path: str, remote_path: str = "", refresh_sec: int = 30,
               max_runs: int = 5):
    app = Dash(__name__, title="Kalign Pareto Front")

    app.layout = html.Div([
        dcc.Store(id="ckpt-path", data=ckpt_path),
        dcc.Store(id="remote-path", data=remote_path or ""),
        dcc.Store(id="max-runs", data=max_runs),
        dcc.Interval(id="refresh-interval", interval=refresh_sec * 1000,
                     disabled=(not remote_path)),

        html.H2("Kalign Unified Optimizer - Pareto Front"),
        html.Div(id="status-bar", style={"marginBottom": "10px", "color": "#666"}),

        # Controls row
        html.Div([
            html.Div([
                html.Label("Color by:"),
                dcc.Dropdown(
                    id="color-by",
                    options=[
                        {"label": "Mode", "value": "mode"},
                        {"label": "VSM amax (R0)", "value": "vsm_amax_0"},
                        {"label": "Seq weights", "value": "seq_weights"},
                        {"label": "Consistency", "value": "consistency"},
                        {"label": "Realign", "value": "realign"},
                        {"label": "Refine", "value": "refine"},
                        {"label": "Matrix (R0)", "value": "matrix_0"},
                        {"label": "GPO (R0)", "value": "gpo_0"},
                    ],
                    value="mode",
                    clearable=False,
                ),
            ], style={"width": "200px", "display": "inline-block", "marginRight": "20px"}),
            html.Div([
                html.Label("X axis:"),
                dcc.Dropdown(
                    id="x-axis",
                    options=[
                        {"label": "Wall time (s)", "value": "wall_time"},
                        {"label": "F1", "value": "f1"},
                        {"label": "TC", "value": "tc"},
                    ],
                    value="wall_time",
                    clearable=False,
                ),
            ], style={"width": "200px", "display": "inline-block", "marginRight": "20px"}),
            html.Div([
                html.Label("Y axis:"),
                dcc.Dropdown(
                    id="y-axis",
                    options=[
                        {"label": "F1", "value": "f1"},
                        {"label": "TC", "value": "tc"},
                        {"label": "Wall time (s)", "value": "wall_time"},
                    ],
                    value="f1",
                    clearable=False,
                ),
            ], style={"width": "200px", "display": "inline-block", "marginRight": "20px"}),
            html.Div([
                html.Label("Filter mode:"),
                dcc.Checklist(
                    id="mode-filter",
                    options=[
                        {"label": "single", "value": "single"},
                        {"label": "ens3", "value": "ens3"},
                        {"label": "ens5", "value": "ens5"},
                    ],
                    value=["single", "ens3", "ens5"],
                    inline=True,
                ),
            ], style={"display": "inline-block"}),
        ], style={"marginBottom": "10px"}),

        # Main scatter plot
        dcc.Graph(id="pareto-scatter", style={"height": "600px"}),

        # --- Tier Selection ---
        html.H3("Default Configuration Selection"),
        html.P("Auto-suggested tiers using knee-point detection on the F1-vs-time "
               "Pareto front. Click any point on the 2D scatter plot, then assign "
               "it to a tier with the buttons below.",
               style={"color": "#666", "fontSize": "14px"}),

        html.Div([
            html.Div([
                html.Button("Set as Fast", id="set-fast-btn", n_clicks=0,
                            style={"backgroundColor": "#2196F3", "color": "white",
                                   "border": "none", "padding": "8px 16px",
                                   "marginRight": "10px", "cursor": "pointer"}),
                html.Button("Set as Default", id="set-default-btn", n_clicks=0,
                            style={"backgroundColor": "#FF9800", "color": "white",
                                   "border": "none", "padding": "8px 16px",
                                   "marginRight": "10px", "cursor": "pointer"}),
                html.Button("Set as Accurate", id="set-accurate-btn", n_clicks=0,
                            style={"backgroundColor": "#4CAF50", "color": "white",
                                   "border": "none", "padding": "8px 16px",
                                   "marginRight": "10px", "cursor": "pointer"}),
                html.Button("Auto-select (kneedle)", id="auto-select-btn", n_clicks=0,
                            style={"backgroundColor": "#9E9E9E", "color": "white",
                                   "border": "none", "padding": "8px 16px",
                                   "cursor": "pointer"}),
            ], style={"marginBottom": "10px"}),
        ]),

        # Individual stores for each tier (avoids output=state on same component)
        dcc.Store(id="tier-fast", data=None),
        dcc.Store(id="tier-default", data=None),
        dcc.Store(id="tier-accurate", data=None),
        # Store for last clicked point idx
        dcc.Store(id="last-clicked-idx", data=None),

        html.Div(id="tier-display", style={"marginBottom": "20px"}),

        html.Div([
            html.Button("Save tiers to JSON", id="save-tiers-btn", n_clicks=0,
                        style={"backgroundColor": "#673AB7", "color": "white",
                               "border": "none", "padding": "8px 16px",
                               "cursor": "pointer", "marginRight": "10px"}),
            html.Span(id="save-status", style={"color": "#666", "fontSize": "14px"}),
        ], style={"marginBottom": "20px"}),

        # 3D scatter
        html.H3("3D Pareto Surface (F1 vs TC vs Time)"),
        dcc.Graph(id="pareto-3d", style={"height": "600px"}),

        # Convergence plot
        html.H3("Convergence (best F1 / TC over evaluations)"),
        dcc.Graph(id="convergence-plot", style={"height": "350px"}),

        # Selected point details
        html.H3("Click a point to see full parameters:"),
        html.Pre(id="point-details",
                 style={"backgroundColor": "#f5f5f5", "padding": "15px",
                        "fontSize": "14px", "whiteSpace": "pre-wrap"}),

        # Top solutions table
        html.H3("3D Pareto Front (non-dominated in F1, TC, and time)"),
        dash_table.DataTable(
            id="top-table",
            columns=[
                {"name": "#", "id": "idx"},
                {"name": "Mode", "id": "mode"},
                {"name": "F1", "id": "f1"},
                {"name": "TC", "id": "tc"},
                {"name": "Time", "id": "wall_time"},
                {"name": "VSM", "id": "vsm"},
                {"name": "SW", "id": "seq_weights"},
                {"name": "C", "id": "consistency"},
                {"name": "CW", "id": "consistency_weight"},
                {"name": "Re", "id": "realign"},
                {"name": "Ref", "id": "refine"},
                {"name": "MS", "id": "min_support"},
                {"name": "GPO(R0)", "id": "gpo_0"},
                {"name": "GPE(R0)", "id": "gpe_0"},
                {"name": "TGPE(R0)", "id": "tgpe_0"},
                {"name": "Mat(R0)", "id": "matrix_0"},
            ],
            style_cell={"textAlign": "right", "padding": "4px", "fontSize": "13px"},
            style_header={"fontWeight": "bold"},
            style_data_conditional=[
                {"if": {"column_id": "mode"}, "textAlign": "left"},
                {"if": {"column_id": "refine"}, "textAlign": "center"},
            ],
            sort_action="native",
            row_selectable="single",
            page_size=50,
        ),
    ], style={"maxWidth": "1400px", "margin": "auto", "padding": "20px"})

    # Store for current data; params_by_idx maps idx -> full decoded params dict
    app._df_cache = {"df": None, "hist_df": None, "mtime": 0, "params_by_idx": {}}  # type: ignore[attr-defined]

    def _load_data(ckpt_path_arg, remote_path_arg, max_runs_arg):
        """Load or refresh data from checkpoint."""
        local_path = ckpt_path_arg

        # Sync from remote if configured
        if remote_path_arg:
            sync_remote(remote_path_arg, local_path)

        if not Path(local_path).exists():
            return None, None, 0

        mtime = Path(local_path).stat().st_mtime
        if mtime == app._df_cache["mtime"]:  # type: ignore[attr-defined]
            return app._df_cache["df"], app._df_cache["hist_df"], mtime  # type: ignore[attr-defined]

        try:
            ckpt = load_checkpoint(local_path)
        except Exception:
            return app._df_cache["df"], app._df_cache["hist_df"], app._df_cache["mtime"]  # type: ignore[attr-defined]

        mr = ckpt.get("max_runs", max_runs_arg)
        df, params_by_idx = build_pareto_df(ckpt, mr)
        hist_df = build_history_df(ckpt, mr)
        n_gen = ckpt.get("n_gen_completed", "?")

        f_beta = ckpt.get("f_beta", 1.0)
        app._df_cache = {"df": df, "hist_df": hist_df, "mtime": mtime, "n_gen": n_gen, "params_by_idx": params_by_idx, "f_beta": f_beta}  # type: ignore[attr-defined]
        return df, hist_df, mtime

    @app.callback(
        Output("pareto-scatter", "figure"),
        Output("pareto-3d", "figure"),
        Output("convergence-plot", "figure"),
        Output("top-table", "data"),
        Output("status-bar", "children"),
        Input("refresh-interval", "n_intervals"),
        Input("color-by", "value"),
        Input("x-axis", "value"),
        Input("y-axis", "value"),
        Input("mode-filter", "value"),
        Input("tier-fast", "data"),
        Input("tier-default", "data"),
        Input("tier-accurate", "data"),
        Input("ckpt-path", "data"),
        Input("remote-path", "data"),
        Input("max-runs", "data"),
    )
    def update_all(n_intervals, color_by, x_axis, y_axis, mode_filter,
                   tier_fast, tier_default, tier_accurate,
                   ckpt_path, remote_path, max_runs):
        tier_selections = {"fast": tier_fast, "default": tier_default, "accurate": tier_accurate}
        df, hist_df, mtime = _load_data(ckpt_path, remote_path or None, max_runs)

        empty_fig = go.Figure()
        if df is None or df.empty:
            return empty_fig, empty_fig, empty_fig, [], "No data loaded"

        # Filter by mode
        filtered = df[df["mode"].isin(mode_filter)] if mode_filter else df

        # Status
        n_gen = app._df_cache.get("n_gen", "?")
        f_beta = app._df_cache.get("f_beta", 1.0)
        obj_label = f"F{f_beta}" if f_beta != 1.0 else "F1"
        ago = time.time() - mtime if mtime else 0
        status = (f"Generation {n_gen} | {len(df)} individuals | "
                  f"Best {obj_label}={df['f1'].max():.4f} | Best TC={df['tc'].max():.4f} | "
                  f"Last update: {ago:.0f}s ago")

        # 2D scatter
        hover_data = ["mode", "f1", "tc", "wall_time",
                      "n_runs", "vsm", "refine", "matrices",
                      "seq_weights", "consistency", "realign", "min_support"]
        fig2d = px.scatter(
            filtered, x=x_axis, y=y_axis, color=color_by,
            hover_data=hover_data,
            title=f"Population ({y_axis} vs {x_axis})",
            template="plotly_white",
        )
        fig2d.update_traces(marker=dict(size=8, opacity=0.7))

        # Compute and draw overall 2D Pareto front line
        pareto_line = _compute_pareto_2d(filtered, x_axis, y_axis)
        if len(pareto_line) > 1:
            fig2d.add_trace(go.Scatter(
                x=pareto_line[x_axis].tolist(),
                y=pareto_line[y_axis].tolist(),
                mode="lines+markers",
                name="Pareto front (all)",
                line=dict(color="rgba(0,0,0,0.7)", width=2.5),
                marker=dict(size=10, symbol="diamond", color="rgba(0,0,0,0.7)"),
                hovertext=[
                    f"F1={r['f1']:.4f} TC={r['tc']:.4f} t={r['wall_time']:.0f}s "
                    f"{r['mode']}"
                    for _, r in pareto_line.iterrows()
                ],
                hoverinfo="text",
            ))

        # Per-mode Pareto front lines
        mode_colors = {"single": "rgba(31,119,180,0.6)",
                       "ens3": "rgba(255,127,14,0.6)",
                       "ens5": "rgba(44,160,44,0.6)"}
        for mode_name in sorted(filtered["mode"].unique()):
            mode_df = filtered[filtered["mode"] == mode_name]
            mode_pareto = _compute_pareto_2d(mode_df, x_axis, y_axis)
            if len(mode_pareto) > 1:
                fig2d.add_trace(go.Scatter(
                    x=mode_pareto[x_axis].tolist(),
                    y=mode_pareto[y_axis].tolist(),
                    mode="lines",
                    name=f"Pareto ({mode_name})",
                    line=dict(color=mode_colors.get(mode_name, "rgba(128,128,128,0.5)"),
                              width=1.5, dash="dash"),
                    hoverinfo="skip",
                    showlegend=True,
                ))

        # Tier star markers
        if tier_selections:
            tier_styles = [
                ("fast", "Fast", "#2196F3"),
                ("default", "Default", "#FF9800"),
                ("accurate", "Accurate", "#4CAF50"),
            ]
            for key, label, color in tier_styles:
                idx = tier_selections.get(key)
                if idx is None:
                    continue
                match = filtered[filtered["idx"] == idx]
                if match.empty:
                    match = df[df["idx"] == idx]
                if match.empty:
                    continue
                r = match.iloc[0]
                fig2d.add_trace(go.Scatter(
                    x=[r[x_axis]], y=[r[y_axis]],
                    mode="markers+text",
                    name=f"★ {label}",
                    marker=dict(size=18, symbol="star", color=color,
                                line=dict(width=2, color="black")),
                    text=[label], textposition="top center",
                    textfont=dict(size=12, color=color),
                    hovertext=f"{label}: F1={r['f1']:.4f} TC={r['tc']:.4f} t={r['wall_time']:.0f}s",
                    hoverinfo="text",
                ))

        fig2d.update_layout(height=600)

        # 3D scatter
        fig3d = px.scatter_3d(
            filtered, x="wall_time", y="f1", z="tc", color=color_by,
            hover_data=hover_data,
            title="3D Pareto Surface",
            template="plotly_white",
        )
        fig3d.update_traces(marker=dict(size=5, opacity=0.7))
        fig3d.update_layout(height=600, scene=dict(
            xaxis_title="Wall time (s)",
            yaxis_title="F1",
            zaxis_title="TC",
        ))

        # Convergence
        if hist_df is not None and not hist_df.empty:
            hist_df = hist_df.copy()
            hist_df["best_f1"] = hist_df["f1"].cummax()
            hist_df["best_tc"] = hist_df["tc"].cummax()
            fig_conv = go.Figure()
            fig_conv.add_trace(go.Scatter(
                x=hist_df["eval"], y=hist_df["best_f1"],
                mode="lines", name="Best F1",
            ))
            fig_conv.add_trace(go.Scatter(
                x=hist_df["eval"], y=hist_df["best_tc"],
                mode="lines", name="Best TC",
            ))
            fig_conv.update_layout(
                template="plotly_white", height=350,
                xaxis_title="Evaluation #", yaxis_title="Score",
            )
        else:
            fig_conv = empty_fig

        # Pareto front table (3D: maximize F1, maximize TC, minimize time)
        pareto_3d = _compute_pareto_3d(filtered)
        pareto_3d = pareto_3d.sort_values("f1", ascending=False)
        table_data = pareto_3d.drop(columns=["run_details"]).to_dict("records")

        return fig2d, fig3d, fig_conv, table_data, status

    @app.callback(
        Output("last-clicked-idx", "data"),
        Input("pareto-scatter", "clickData"),
        State("mode-filter", "value"),
        State("x-axis", "value"),
        State("y-axis", "value"),
    )
    def track_click(click_data, mode_filter, x_axis, y_axis):
        """Track the idx of the last clicked point on the 2D scatter."""
        df = app._df_cache.get("df")
        if df is None or click_data is None:
            return None
        pts = click_data.get("points", [])
        if not pts:
            return None
        pt = pts[0]
        x_val = pt.get("x")
        y_val = pt.get("y")
        if x_val is not None and y_val is not None:
            filtered = df[df["mode"].isin(mode_filter)] if mode_filter else df
            if filtered.empty:
                return None
            # Normalize distances by range to handle different scales
            x_range = filtered[x_axis].max() - filtered[x_axis].min()
            y_range = filtered[y_axis].max() - filtered[y_axis].min()
            x_range = max(x_range, 1e-10)
            y_range = max(y_range, 1e-10)
            dists = ((filtered[x_axis] - x_val) / x_range)**2 + \
                    ((filtered[y_axis] - y_val) / y_range)**2
            best = dists.idxmin()
            return int(filtered.loc[best, "idx"])
        return None

    @app.callback(
        Output("tier-fast", "data"),
        Output("tier-default", "data"),
        Output("tier-accurate", "data"),
        Input("set-fast-btn", "n_clicks"),
        Input("set-default-btn", "n_clicks"),
        Input("set-accurate-btn", "n_clicks"),
        Input("auto-select-btn", "n_clicks"),
        State("tier-fast", "data"),
        State("tier-default", "data"),
        State("tier-accurate", "data"),
        State("last-clicked-idx", "data"),
        State("mode-filter", "value"),
    )
    def update_tiers(_fc, _dc, _ac, _auto,
                     cur_fast, cur_default, cur_accurate, last_idx, mode_filter):
        """Update tier selections based on button clicks."""
        triggered_id = ctx.triggered_id
        if triggered_id is None:
            return cur_fast, cur_default, cur_accurate

        if triggered_id == "auto-select-btn":
            df = app._df_cache.get("df")
            if df is not None and not df.empty:
                filtered = df[df["mode"].isin(mode_filter)] if mode_filter else df
                tiers = auto_select_tiers(filtered)
                return tiers.get("fast"), tiers.get("default"), tiers.get("accurate")
            return cur_fast, cur_default, cur_accurate

        if last_idx is None:
            return cur_fast, cur_default, cur_accurate

        if triggered_id == "set-fast-btn":
            return last_idx, cur_default, cur_accurate
        elif triggered_id == "set-default-btn":
            return cur_fast, last_idx, cur_accurate
        elif triggered_id == "set-accurate-btn":
            return cur_fast, cur_default, last_idx

        return cur_fast, cur_default, cur_accurate

    @app.callback(
        Output("tier-display", "children"),
        Input("tier-fast", "data"),
        Input("tier-default", "data"),
        Input("tier-accurate", "data"),
    )
    def render_tiers(tier_fast, tier_default, tier_accurate):
        """Render the selected tier configurations."""
        df = app._df_cache.get("df")
        tiers = {"fast": tier_fast, "default": tier_default, "accurate": tier_accurate}
        if df is None:
            return html.Div("No tiers selected yet.", style={"color": "#999"})

        tier_names = [("fast", "Fast", "#2196F3"),
                      ("default", "Default", "#FF9800"),
                      ("accurate", "Accurate", "#4CAF50")]
        cards = []

        for key, label, color in tier_names:
            idx = tiers.get(key)
            if idx is None:
                cards.append(html.Div([
                    html.H4(f"{label}", style={"color": color, "marginBottom": "5px"}),
                    html.Span("Not set — click a point then press the button",
                              style={"color": "#999", "fontSize": "13px"}),
                ], style={"display": "inline-block", "verticalAlign": "top",
                          "width": "32%", "marginRight": "1%",
                          "padding": "10px", "backgroundColor": "#fafafa",
                          "border": f"2px solid {color}", "borderRadius": "8px"}))
                continue

            match = df[df["idx"] == idx]
            if match.empty:
                continue
            row = match.iloc[0]
            config_text = _format_tier_config(row)
            cards.append(html.Div([
                html.H4(f"{label}", style={"color": color, "marginBottom": "5px"}),
                html.Div(f"F1={row['f1']:.4f}  TC={row['tc']:.4f}  Time={row['wall_time']:.1f}s",
                         style={"fontWeight": "bold", "marginBottom": "5px"}),
                html.Div(f"{row['mode']}", style={"fontSize": "13px", "marginBottom": "8px"}),
                html.Pre(config_text,
                         style={"backgroundColor": "#f0f0f0", "padding": "8px",
                                "fontSize": "12px", "whiteSpace": "pre-wrap",
                                "margin": "0", "borderRadius": "4px"}),
            ], style={"display": "inline-block", "verticalAlign": "top",
                      "width": "32%", "marginRight": "1%",
                      "padding": "10px", "backgroundColor": "#fafafa",
                      "border": f"2px solid {color}", "borderRadius": "8px"}))

        return html.Div(cards, style={"marginBottom": "20px"})

    @app.callback(
        Output("save-status", "children"),
        Input("save-tiers-btn", "n_clicks"),
        State("tier-fast", "data"),
        State("tier-default", "data"),
        State("tier-accurate", "data"),
        State("ckpt-path", "data"),
    )
    def save_tiers(n_clicks, tier_fast, tier_default, tier_accurate, ckpt_path):
        """Save selected tiers to a JSON file next to the checkpoint."""
        tiers = {"fast": tier_fast, "default": tier_default, "accurate": tier_accurate}
        if not n_clicks:
            return ""

        df = app._df_cache.get("df")
        params_by_idx = app._df_cache.get("params_by_idx", {})
        if df is None:
            return "No data loaded."

        any_set = any(tiers.get(k) is not None for k in ("fast", "default", "accurate"))
        if not any_set:
            return "No tiers selected yet."

        output = {}
        for tier_name in ("fast", "default", "accurate"):
            idx = tiers.get(tier_name)
            if idx is None:
                continue
            match = df[df["idx"] == idx]
            if match.empty:
                continue
            row = match.iloc[0]

            # Full decoded params from the optimizer
            full_params = params_by_idx.get(idx, {})

            # Build per-run config list
            runs = []
            for k in range(int(row["n_runs"])):
                run = {
                    "gpo": round(float(full_params["run_gpo"][k]), 4),
                    "gpe": round(float(full_params["run_gpe"][k]), 4),
                    "tgpe": round(float(full_params["run_tgpe"][k]), 4),
                    "matrix": MATRIX_NAMES.get(full_params["run_types"][k], "?"),
                    "vsm_amax": round(float(full_params["run_vsm_amax"][k]), 4),
                    "refine": REFINE_LONG.get(full_params["run_refine"][k], "NONE"),
                }
                if full_params["run_noise"][k] > 0:
                    run["noise"] = round(float(full_params["run_noise"][k]), 4)
                runs.append(run)

            output[tier_name] = {
                "scores": {
                    "f1": float(row["f1"]),
                    "tc": float(row["tc"]),
                    "wall_time": float(row["wall_time"]),
                },
                "params": {
                    "n_runs": int(row["n_runs"]),
                    "seq_weights": float(row["seq_weights"]),
                    "consistency": int(row["consistency"]),
                    "consistency_weight": float(row["consistency_weight"]),
                    "realign": int(row["realign"]),
                    "min_support": int(row["min_support"]),
                },
                "runs": runs,
            }

        # Save next to checkpoint file
        out_path = Path(ckpt_path).parent / "kalign_tiers.json"
        with open(out_path, "w") as f:
            json.dump(output, f, indent=2)

        n_tiers = len(output)
        return f"Saved {n_tiers} tier(s) to {out_path}"

    @app.callback(
        Output("point-details", "children"),
        Input("pareto-scatter", "clickData"),
        Input("pareto-3d", "clickData"),
        Input("top-table", "selected_rows"),
        State("top-table", "data"),
        State("x-axis", "value"),
        State("y-axis", "value"),
    )
    def show_details(click_2d, click_3d, selected_rows, table_data,
                     x_axis, y_axis):
        df = app._df_cache.get("df")
        if df is None or df.empty:
            return "Click a point or select a table row to see details."

        triggered_id = ctx.triggered_id

        # From table selection
        if triggered_id == "top-table" and selected_rows and table_data:
            row = table_data[selected_rows[0]]
            idx = row["idx"]
            match = df[df["idx"] == idx]
            if not match.empty:
                return _format_detail(match.iloc[0])

        # From scatter click — match by coordinates
        click = None
        if triggered_id == "pareto-3d" and click_3d:
            click = click_3d
        elif triggered_id == "pareto-scatter" and click_2d:
            click = click_2d

        if click and click.get("points"):
            pt = click["points"][0]
            x_val = pt.get("x")
            y_val = pt.get("y")
            if x_val is not None and y_val is not None:
                # For 3D clicks, match on wall_time and f1
                if triggered_id == "pareto-3d":
                    x_col, y_col = "wall_time", "f1"
                else:
                    x_col, y_col = x_axis, y_axis
                x_range = max(df[x_col].max() - df[x_col].min(), 1e-10)
                y_range = max(df[y_col].max() - df[y_col].min(), 1e-10)
                dists = ((df[x_col] - x_val) / x_range)**2 + \
                        ((df[y_col] - y_val) / y_range)**2
                best = dists.idxmin()
                return _format_detail(df.loc[best])

        return "Click a point or select a table row to see details."

    return app


def _format_detail(row):
    """Format full details for a selected solution."""
    lines = [
        f"Solution #{row['idx']}",
        f"{'='*50}",
        f"F1 = {row['f1']:.4f}    TC = {row['tc']:.4f}    Time = {row['wall_time']:.1f}s",
        f"Mode: {row['mode']}    n_runs={row['n_runs']}",
        f"",
        f"Per-run parameters:",
        row["run_details"],
        f"",
        f"Shared parameters:",
        f"  seq_weights      = {row['seq_weights']}",
        f"  consistency      = {row['consistency']}",
        f"  consistency_wt   = {row['consistency_weight']}",
        f"  realign          = {row['realign']}",
        f"  min_support      = {row['min_support']}",
        f"  (vsm_amax, refine are per-run — see run details above)",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Interactive Pareto front viewer for kalign optimizer")
    parser.add_argument("checkpoint", help="Path to gen_checkpoint.pkl (local)")
    parser.add_argument("--remote", default=None,
                        help="Remote path (e.g. server:path/gen_checkpoint.pkl) to auto-sync")
    parser.add_argument("--port", type=int, default=8050,
                        help="Dash server port (default: 8050)")
    parser.add_argument("--refresh", type=int, default=30,
                        help="Auto-refresh interval in seconds (default: 30, remote only)")
    parser.add_argument("--max-runs", type=int, default=5,
                        help="Max ensemble runs (must match optimizer, default: 5)")
    args = parser.parse_args()

    # If remote specified, do initial sync
    if args.remote:
        print(f"Syncing from {args.remote}...")
        sync_remote(args.remote, args.checkpoint)

    if not Path(args.checkpoint).exists():
        print(f"Checkpoint not found: {args.checkpoint}")
        return

    app = create_app(
        ckpt_path=args.checkpoint,
        remote_path=args.remote or "",
        refresh_sec=args.refresh,
        max_runs=args.max_runs,
    )

    print(f"Starting Pareto viewer on http://localhost:{args.port}")
    if args.remote:
        print(f"Auto-refreshing from {args.remote} every {args.refresh}s")
    app.run(debug=False, port=args.port)


if __name__ == "__main__":
    main()
