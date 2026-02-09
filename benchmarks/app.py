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


def _results_to_df(results):
    """Convert list of AlignmentResult dicts to DataFrame."""
    df = pd.DataFrame(results)
    if "error" in df.columns:
        df = df[df["error"].isna() | (df["error"] == "None")]
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

    # Box plot by dataset
    fig_box = px.box(
        clean, x="dataset", y="sp_score", color="method",
        title="SP Score Distribution by Dataset",
        labels={"sp_score": "SP Score", "dataset": "Dataset"},
    )
    fig_box.update_layout(legend=dict(orientation="h", y=-0.2))
    figs.append(("sp_box", fig_box))

    # Scatter by family
    fig_scatter = px.strip(
        clean, x="dataset", y="sp_score", color="method",
        hover_data=["family", "wall_time"],
        title="SP Score per Family",
        labels={"sp_score": "SP Score"},
    )
    fig_scatter.update_layout(legend=dict(orientation="h", y=-0.2))
    figs.append(("sp_scatter", fig_scatter))

    # Timing
    fig_time = px.box(
        clean, x="dataset", y="wall_time", color="method",
        title="Alignment Time by Dataset",
        labels={"wall_time": "Wall Time (s)", "dataset": "Dataset"},
    )
    fig_time.update_layout(legend=dict(orientation="h", y=-0.2))
    figs.append(("timing", fig_time))

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


def _run_in_thread(dataset, methods, binary, max_cases, n_threads):
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

        total = len(cases) * len(methods)
        done = 0

        for case in cases:
            for method in methods:
                done += 1
                _run_state["progress"] = f"[{done}/{total}] {case.family} ({method})..."
                result = run_case(case, method=method, binary=binary, n_threads=n_threads)
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
                ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                html.Div([
                    html.Label("Method"),
                    dcc.Checklist(
                        id="method-checklist",
                        options=[
                            {"label": " Python API", "value": "python_api"},
                            {"label": " C binary", "value": "cli"},
                        ],
                        value=["python_api"],
                    ),
                ], style={"width": "12%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                html.Div([
                    html.Label("C binary path"),
                    dcc.Input(id="binary-input", type="text", value="build/src/kalign",
                              style={"width": "160px"}),
                ], style={"width": "18%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                html.Div([
                    html.Label("Max cases (0 = all)"),
                    dcc.Input(id="max-cases-input", type="number", value=0, min=0, style={"width": "80px"}),
                ], style={"width": "12%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                html.Div([
                    html.Label("Threads"),
                    dcc.Input(id="threads-input", type="number", value=1, min=1, style={"width": "60px"}),
                ], style={"width": "8%", "display": "inline-block", "verticalAlign": "top", "marginRight": "2%"}),
                html.Div([
                    html.Br(),
                    html.Button("Download Dataset", id="download-btn", n_clicks=0,
                                style={"marginRight": "10px"}),
                    html.Button("Run Benchmark", id="run-btn", n_clicks=0,
                                style={"backgroundColor": "#4CAF50", "color": "white",
                                       "border": "none", "padding": "8px 16px", "cursor": "pointer"}),
                ], style={"width": "15%", "display": "inline-block", "verticalAlign": "top"}),
            ]),
            html.Div(id="progress-text", style={"marginTop": "10px", "fontStyle": "italic"}),
            dcc.Interval(id="progress-interval", interval=1000, n_intervals=0, disabled=True),
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
        Output("progress-text", "children"),
        Output("progress-interval", "disabled"),
        Input("run-btn", "n_clicks"),
        Input("download-btn", "n_clicks"),
        Input("progress-interval", "n_intervals"),
        State("dataset-dropdown", "value"),
        State("method-checklist", "value"),
        State("binary-input", "value"),
        State("max-cases-input", "value"),
        State("threads-input", "value"),
        prevent_initial_call=True,
    )
    def handle_run_or_download(run_clicks, dl_clicks, _n_intervals,
                                dataset, methods, binary, max_cases, n_threads):
        triggered = callback_context.triggered_id

        if triggered == "download-btn" and not _run_state["running"]:
            try:
                download_dataset(dataset)
                return f"Downloaded {dataset}.", True
            except Exception as e:
                return f"Download error: {e}", True

        if triggered == "run-btn" and not _run_state["running"]:
            if not methods:
                return "Select at least one method.", True
            t = threading.Thread(
                target=_run_in_thread,
                args=(dataset, methods, binary or "build/src/kalign",
                      max_cases or 0, n_threads or 1),
                daemon=True,
            )
            t.start()
            return "Starting...", False  # enable interval

        # Polling progress
        if _run_state["running"]:
            return _run_state["progress"], False
        if _run_state["done"]:
            return _run_state["progress"], True  # disable interval

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
    parser.add_argument("--port", type=int, default=8050, help="Port (default: 8050)")
    args = parser.parse_args()

    app = create_app()
    print(f"Starting dashboard at http://localhost:{args.port}")
    app.run(debug=False, port=args.port)


if __name__ == "__main__":
    main()
