"""Dash app for running kalign benchmarks and viewing results.

Usage:
    python -m benchmarks.app [--port 8050]
"""

import argparse
import json
import sys
import threading
from pathlib import Path

try:
    import dash
    from dash import dcc, html, dash_table, callback_context
    from dash.dependencies import Input, Output, State
    import plotly.express as px
    import pandas as pd
except ImportError:
    print("Dash visualization requires extra dependencies:")
    print("  pip install dash plotly pandas")
    sys.exit(1)

from .datasets import DATASETS, get_cases, download_dataset
from .scoring import run_case

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration presets — the 4 varieties of kalign
# ---------------------------------------------------------------------------

CONFIG_PRESETS = {
    "Kalign (default)": {"refine": "none", "ensemble": 0},
    "Kalign + Refinement": {"refine": "confident", "ensemble": 0},
    "Ensemble (8 runs)": {"refine": "none", "ensemble": 8},
    "Ensemble (12 runs)": {"refine": "none", "ensemble": 12},
    "Clustal Omega": {"method": "clustalo", "refine": "none", "ensemble": 0},
    "MAFFT": {"method": "mafft", "refine": "none", "ensemble": 0},
    "MUSCLE": {"method": "muscle", "refine": "none", "ensemble": 0},
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _available_datasets():
    """Return list of (name, available, n_cases) tuples."""
    info = []
    for name, ds in DATASETS.items():
        avail = ds["is_available"]()
        n = len(ds["cases"]()) if avail else 0
        info.append((name, avail, n))
    return info


def _load_json(path):
    with open(path) as f:
        return json.load(f)


_EXTERNAL_LABELS = {
    "clustalo": "Clustal Omega",
    "mafft": "MAFFT",
    "muscle": "MUSCLE",
}


def _config_label(row):
    """Build a human-readable configuration label from result fields."""
    method = row.get("method", "python_api")
    if method in _EXTERNAL_LABELS:
        return _EXTERNAL_LABELS[method]

    ensemble = row.get("ensemble", 0)
    refine = row.get("refine", "none")

    if ensemble and ensemble > 0:
        label = f"Ensemble ({ensemble})"
    elif refine and refine != "none":
        label = f"Kalign + {refine}"
    else:
        label = "Kalign"
    return label


def _results_to_df(results):
    """Convert list of AlignmentResult dicts to DataFrame."""
    df = pd.DataFrame(results)
    if "error" in df.columns:
        df = df[df["error"].isna() | (df["error"] == "None") | (df["error"].isnull())]
    # Ensure columns exist (backward compat with old JSON files)
    if "ensemble" not in df.columns:
        df["ensemble"] = 0
    df["ensemble"] = df["ensemble"].fillna(0).astype(int)
    for col in ("recall", "precision", "f1", "tc"):
        if col not in df.columns:
            df[col] = 0.0
        df[col] = df[col].fillna(0.0)
    # Add config label for plotting
    df["config"] = df.apply(_config_label, axis=1)
    return df


def _build_figures(df):
    """Build plotly figures from a results DataFrame."""
    figs = []
    if df.empty:
        return figs

    # Filter out errors
    clean = df[df["sp_score"] > 0].copy() if "sp_score" in df.columns else df

    if clean.empty:
        return figs

    # Ensure new columns exist (backward compat with old JSON files)
    for col in ("recall", "precision", "f1", "tc"):
        if col not in clean.columns:
            clean[col] = 0.0

    # Determine color column: use 'config' if multiple configs, else 'method'
    configs = clean["config"].nunique()
    color_col = "config" if configs > 1 else "method"

    # bali_score-compatible metrics
    has_recall = clean["recall"].sum() > 0
    if has_recall:
        # SP score (bali_score compatible) = recall
        fig_sp = px.box(
            clean, x="dataset", y="recall", color=color_col,
            title="SP Score by Dataset (bali_score compatible)",
            labels={"recall": "SP Score", "dataset": "Dataset", "config": "Configuration"},
        )
        fig_sp.update_layout(legend=dict(orientation="h", y=-0.2))
        figs.append(("sp_box", fig_sp))

        # TC score
        fig_tc = px.box(
            clean, x="dataset", y="tc", color=color_col,
            title="TC Score by Dataset",
            labels={"tc": "TC Score", "dataset": "Dataset", "config": "Configuration"},
        )
        fig_tc.update_layout(legend=dict(orientation="h", y=-0.2))
        figs.append(("tc_box", fig_tc))

        # Precision
        fig_prec = px.box(
            clean, x="dataset", y="precision", color=color_col,
            title="Precision by Dataset",
            labels={"precision": "Precision", "dataset": "Dataset", "config": "Configuration"},
        )
        fig_prec.update_layout(legend=dict(orientation="h", y=-0.2))
        figs.append(("precision_box", fig_prec))

        # F1
        fig_f1 = px.box(
            clean, x="dataset", y="f1", color=color_col,
            title="F1 Score by Dataset",
            labels={"f1": "F1", "dataset": "Dataset", "config": "Configuration"},
        )
        fig_f1.update_layout(legend=dict(orientation="h", y=-0.2))
        figs.append(("f1_box", fig_f1))
    else:
        # Fallback: legacy SP score (not bali_score compatible)
        fig_box = px.box(
            clean, x="dataset", y="sp_score", color=color_col,
            title="SP Score Distribution by Dataset (legacy)",
            labels={"sp_score": "SP Score", "dataset": "Dataset", "config": "Configuration"},
        )
        fig_box.update_layout(legend=dict(orientation="h", y=-0.2))
        figs.append(("sp_box", fig_box))

    # Strip plot by family
    hover_cols = ["family", "wall_time", "refine", "ensemble"]
    if has_recall:
        hover_cols.extend(["recall", "precision", "f1", "tc"])
    y_strip = "recall" if has_recall else "sp_score"
    fig_scatter = px.strip(
        clean, x="dataset", y=y_strip, color=color_col,
        hover_data=hover_cols,
        title="SP Score per Family",
        labels={"recall": "SP Score", "sp_score": "SP Score (legacy)", "config": "Configuration"},
    )
    fig_scatter.update_layout(legend=dict(orientation="h", y=-0.2))
    figs.append(("scatter", fig_scatter))

    # Timing
    fig_time = px.box(
        clean, x="dataset", y="wall_time", color=color_col,
        title="Alignment Time by Dataset",
        labels={"wall_time": "Wall Time (s)", "dataset": "Dataset", "config": "Configuration"},
    )
    fig_time.update_layout(legend=dict(orientation="h", y=-0.2))
    figs.append(("timing", fig_time))

    # Summary table by config x dataset
    if configs > 1:
        agg_dict = {
            "mean_sp": ("sp_score", "mean"),
            "median_sp": ("sp_score", "median"),
            "mean_time": ("wall_time", "mean"),
            "n_cases": ("sp_score", "count"),
        }
        if has_recall:
            agg_dict.update({
                "mean_recall": ("recall", "mean"),
                "mean_precision": ("precision", "mean"),
                "mean_f1": ("f1", "mean"),
                "mean_tc": ("tc", "mean"),
            })
        summary = clean.groupby(["config", "dataset"]).agg(**agg_dict).reset_index()
        summary = summary.round(3)
        y_col = "mean_recall" if has_recall else "mean_sp"
        y_label = "Mean SP Score" if has_recall else "Mean SP Score (legacy)"
        fig_summary = px.bar(
            summary, x="dataset", y=y_col, color="config",
            barmode="group",
            title=f"{y_label} by Dataset and Configuration",
            labels={y_col: y_label, "dataset": "Dataset", "config": "Configuration"},
        )
        fig_summary.update_layout(legend=dict(orientation="h", y=-0.2))
        figs.append(("summary_bar", fig_summary))

    return figs


# ---------------------------------------------------------------------------
# State shared between callbacks and the benchmark thread
# ---------------------------------------------------------------------------

_run_state = {
    "running": False,
    "progress": "",
    "results": [],
    "done": False,
}


def _run_in_thread(dataset, method, binary, max_cases, n_threads, configs):
    """Run benchmarks in a background thread so the UI stays responsive."""
    global _run_state
    _run_state["running"] = True
    _run_state["done"] = False
    _run_state["results"] = []
    _run_state["progress"] = f"Loading {dataset} cases..."

    try:
        cases = get_cases(dataset, max_cases=max_cases if max_cases > 0 else None)
        if not cases:
            _run_state["progress"] = f"No cases found for {dataset}. Try downloading first."
            _run_state["running"] = False
            _run_state["done"] = True
            return

        total = len(cases) * len(configs)
        done = 0

        for cfg_name, cfg in configs:
            cfg_method = cfg.get("method", method)
            refine = cfg["refine"]
            ensemble = cfg["ensemble"]

            for case in cases:
                done += 1
                _run_state["progress"] = (
                    f"[{done}/{total}] {cfg_name}: {case.family}..."
                )
                result = run_case(
                    case, method=cfg_method, binary=binary,
                    n_threads=n_threads, refine=refine, ensemble=ensemble,
                )
                _run_state["results"].append(result.to_dict())

        # Auto-save
        import time as _time
        save_path = RESULTS_DIR / f"run_{_time.strftime('%Y%m%d_%H%M%S')}.json"
        data = {
            "timestamp": _time.strftime("%Y-%m-%dT%H:%M:%S"),
            "results": _run_state["results"],
        }
        with open(save_path, "w") as f:
            json.dump(data, f, indent=2)

        _run_state["progress"] = f"Done! {done} alignments scored. Saved to {save_path.name}"
    except Exception as e:
        _run_state["progress"] = f"Error: {e}"
    finally:
        _run_state["running"] = False
        _run_state["done"] = True


# ---------------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------------

def create_app():
    app = dash.Dash(__name__)
    app.title = "Kalign Benchmarks"

    # Discover saved result files
    def _saved_files():
        return sorted(RESULTS_DIR.glob("*.json"), reverse=True)

    ds_info = _available_datasets()

    app.layout = html.Div([
        html.H1("Kalign Benchmark Dashboard"),

        # --- Run controls ---
        html.Div([
            html.H3("Run Benchmark"),
            html.Div([
                # Row 1: Dataset, method, binary
                html.Div([
                    html.Div([
                        html.Label("Dataset"),
                        dcc.Dropdown(
                            id="dataset-dropdown",
                            options=[
                                {"label": f"{name} ({'available' if avail else 'needs download'}, {n} cases)",
                                 "value": name}
                                for name, avail, n in ds_info
                            ],
                            value="balibase",
                        ),
                    ], style={"width": "30%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                    html.Div([
                        html.Label("Method"),
                        dcc.RadioItems(
                            id="method-radio",
                            options=[
                                {"label": " Python API", "value": "python_api"},
                                {"label": " C binary", "value": "cli"},
                            ],
                            value="python_api",
                        ),
                    ], style={"width": "12%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                    html.Div([
                        html.Label("C binary path"),
                        dcc.Input(id="binary-input", type="text", value="build/src/kalign",
                                  style={"width": "160px"}),
                    ], style={"width": "18%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                    html.Div([
                        html.Label("Max cases (0=all)"),
                        dcc.Input(id="max-cases-input", type="number", value=0, min=0, style={"width": "80px"}),
                    ], style={"width": "10%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                    html.Div([
                        html.Label("Threads"),
                        dcc.Input(id="threads-input", type="number", value=1, min=1, style={"width": "60px"}),
                    ], style={"width": "8%", "display": "inline-block", "verticalAlign": "top"}),
                ]),

                html.Hr(style={"margin": "10px 0"}),

                # Row 2: Configuration presets
                html.Div([
                    html.Div([
                        html.Label("Configurations to run"),
                        dcc.Checklist(
                            id="config-checklist",
                            options=[{"label": f" {name}", "value": name}
                                     for name in CONFIG_PRESETS],
                            value=["Kalign (default)"],
                            style={"lineHeight": "2"},
                        ),
                    ], style={"width": "35%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                    html.Div([
                        html.Label("Custom ensemble runs"),
                        dcc.Input(id="custom-ensemble-input", type="number", value=8, min=2, max=32,
                                  style={"width": "80px"}),
                        html.Br(),
                        html.Label("Custom refine mode", style={"marginTop": "8px"}),
                        dcc.Dropdown(
                            id="custom-refine-dropdown",
                            options=[
                                {"label": "none", "value": "none"},
                                {"label": "confident", "value": "confident"},
                                {"label": "all", "value": "all"},
                            ],
                            value="none",
                            style={"width": "140px"},
                        ),
                        html.Button("+ Add custom config", id="add-custom-btn", n_clicks=0,
                                    style={"marginTop": "8px"}),
                    ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                    html.Div([
                        html.Br(),
                        html.Button("Download Dataset", id="download-btn", n_clicks=0,
                                    style={"marginRight": "10px"}),
                        html.Br(), html.Br(),
                        html.Button("Run Benchmark", id="run-btn", n_clicks=0,
                                    style={"backgroundColor": "#4CAF50", "color": "white",
                                           "border": "none", "padding": "12px 24px", "cursor": "pointer",
                                           "fontSize": "14px"}),
                    ], style={"width": "20%", "display": "inline-block", "verticalAlign": "top"}),
                ]),
            ]),
            html.Div(id="progress-text", style={"marginTop": "10px", "fontStyle": "italic"}),
            dcc.Interval(id="progress-interval", interval=1000, n_intervals=0, disabled=True),
            # Hidden store for custom configs
            dcc.Store(id="custom-configs-store", data=[]),
        ], style={"padding": "15px", "backgroundColor": "#f5f5f5", "borderRadius": "8px", "marginBottom": "20px"}),

        # --- Load saved results ---
        html.Div([
            html.H3("View Results"),
            html.Div([
                html.Div([
                    html.Label("Saved result files"),
                    dcc.Dropdown(id="results-dropdown", multi=True),
                ], style={"width": "60%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                html.Div([
                    html.Br(),
                    html.Button("Refresh file list", id="refresh-btn", n_clicks=0,
                                style={"marginRight": "10px"}),
                    html.Button("Load & Plot", id="load-btn", n_clicks=0),
                ], style={"width": "30%", "display": "inline-block", "verticalAlign": "top"}),
            ]),
        ], style={"padding": "15px", "backgroundColor": "#f0f8ff", "borderRadius": "8px", "marginBottom": "20px"}),

        # --- Charts ---
        html.Div(id="charts-container"),

        # --- Table ---
        html.Div(id="table-container"),

    ], style={"maxWidth": "1200px", "margin": "0 auto", "padding": "20px", "fontFamily": "sans-serif"})

    # --- Callbacks ---

    @app.callback(
        Output("results-dropdown", "options"),
        Input("refresh-btn", "n_clicks"),
        Input("progress-interval", "n_intervals"),
    )
    def refresh_file_list(_n1, _n2):
        files = _saved_files()
        return [{"label": f.name, "value": str(f)} for f in files]

    @app.callback(
        Output("custom-configs-store", "data"),
        Output("config-checklist", "options"),
        Output("config-checklist", "value"),
        Input("add-custom-btn", "n_clicks"),
        State("custom-ensemble-input", "value"),
        State("custom-refine-dropdown", "value"),
        State("custom-configs-store", "data"),
        State("config-checklist", "options"),
        State("config-checklist", "value"),
        prevent_initial_call=True,
    )
    def add_custom_config(_n, ensemble_n, refine, custom_cfgs, options, selected):
        """Add a custom configuration to the checklist."""
        ensemble_n = ensemble_n or 0
        refine = refine or "none"

        if ensemble_n > 0:
            name = f"Ensemble ({ensemble_n}) + {refine}"
        else:
            name = f"Custom (refine={refine})"

        cfg = {"refine": refine, "ensemble": ensemble_n}

        # Avoid duplicates
        for existing in custom_cfgs:
            if existing["name"] == name:
                return custom_cfgs, options, selected

        custom_cfgs.append({"name": name, **cfg})
        options.append({"label": f" {name}", "value": name})
        selected.append(name)
        return custom_cfgs, options, selected

    @app.callback(
        Output("progress-text", "children"),
        Output("progress-interval", "disabled"),
        Input("run-btn", "n_clicks"),
        Input("download-btn", "n_clicks"),
        Input("progress-interval", "n_intervals"),
        State("dataset-dropdown", "value"),
        State("method-radio", "value"),
        State("binary-input", "value"),
        State("max-cases-input", "value"),
        State("threads-input", "value"),
        State("config-checklist", "value"),
        State("custom-configs-store", "data"),
        prevent_initial_call=True,
    )
    def handle_run_or_download(run_clicks, dl_clicks, _n_intervals,
                                dataset, method, binary, max_cases, n_threads,
                                selected_configs, custom_cfgs):
        triggered = callback_context.triggered_id

        if triggered == "download-btn" and not _run_state["running"]:
            try:
                download_dataset(dataset)
                return f"Downloaded {dataset}.", True
            except Exception as e:
                return f"Download error: {e}", True

        if triggered == "run-btn" and not _run_state["running"]:
            if not selected_configs:
                return "Select at least one configuration.", True

            # Build config list from presets + custom
            configs = []
            custom_by_name = {c["name"]: c for c in (custom_cfgs or [])}
            for name in selected_configs:
                if name in CONFIG_PRESETS:
                    configs.append((name, CONFIG_PRESETS[name]))
                elif name in custom_by_name:
                    c = custom_by_name[name]
                    configs.append((name, {"refine": c["refine"], "ensemble": c["ensemble"]}))

            if not configs:
                return "No valid configurations selected.", True

            t = threading.Thread(
                target=_run_in_thread,
                args=(dataset, method, binary or "build/src/kalign",
                      max_cases or 0, n_threads or 1, configs),
                daemon=True,
            )
            t.start()
            return f"Starting {len(configs)} configuration(s)...", False

        # Polling progress
        if _run_state["running"]:
            return _run_state["progress"], False
        if _run_state["done"]:
            return _run_state["progress"], True

        return dash.no_update, dash.no_update

    @app.callback(
        Output("charts-container", "children"),
        Output("table-container", "children"),
        Input("load-btn", "n_clicks"),
        Input("progress-interval", "n_intervals"),
        State("results-dropdown", "value"),
        prevent_initial_call=True,
    )
    def update_charts(load_clicks, _n_intervals, selected_files):
        triggered = callback_context.triggered_id

        # If benchmark just finished, show those results immediately
        if triggered == "progress-interval" and _run_state["done"] and _run_state["results"]:
            df = _results_to_df(_run_state["results"])
        elif triggered == "load-btn" and selected_files:
            all_results = []
            for path in selected_files:
                data = _load_json(path)
                for r in data["results"]:
                    r["source"] = Path(path).stem
                all_results.extend(data["results"])
            df = _results_to_df(all_results)
        else:
            return dash.no_update, dash.no_update

        if df.empty:
            return html.P("No results to display."), ""

        charts = []
        for fig_id, fig in _build_figures(df):
            charts.append(dcc.Graph(id=fig_id, figure=fig))

        # Table
        display_cols = [c for c in df.columns if c not in ("error",)]
        table = dash_table.DataTable(
            data=df.to_dict("records"),
            columns=[{"name": c, "id": c} for c in display_cols],
            filter_action="native",
            sort_action="native",
            page_size=50,
            style_table={"overflowX": "auto"},
            style_cell={"textAlign": "left", "padding": "5px"},
            style_header={"fontWeight": "bold"},
        )

        return charts, html.Div([html.H3("Results Table"), table])

    return app


def main():
    parser = argparse.ArgumentParser(
        description="Kalign benchmark dashboard",
        prog="python -m benchmarks.app",
    )
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind to (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8050, help="Port (default: 8050)")
    args = parser.parse_args()

    app = create_app()
    print(f"Starting dashboard at http://{args.host}:{args.port}")
    app.run(debug=False, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
