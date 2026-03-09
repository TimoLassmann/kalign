#!/usr/bin/env python3
"""Multi-objective hyperparameter optimization for kalign using pymoo.

Uses stratified k-fold cross-validation so NSGA-II optimises on held-out
scores, not training scores.  Each evaluation splits BAliBASE families
into k folds (stratified by RV category), aligns k-1 folds, and scores
on the held-out fold.  The objective is the mean held-out F1 / TC.

Objectives (maximized):
    1. Mean held-out F1 across folds (category-averaged within each fold)
    2. Mean held-out TC across folds (category-averaged within each fold)

Decision variables:
    - gpo:                gap open penalty           [2.0, 15.0]
    - gpe:                gap extend penalty          [0.5, 5.0]
    - tgpe:               terminal gap extend         [0.1, 3.0]
    - vsm_amax:           variable scoring matrix     [0.0, 5.0]
    - seq_weights:        profile rebalancing         [0.0, 5.0]
    - consistency:        anchor consistency rounds   {0, 1, 2, 3, 5, 8}
    - consistency_weight: consistency bonus weight    [0.5, 5.0]
    - realign:            tree-rebuild iterations     {0, 1, 2}
    - matrix:             substitution matrix choice  {PFASUM43, PFASUM60, CorBLOSUM66}

Usage:
    # Quick test (small pop, few generations)
    uv run python -m benchmarks.optimize_params --pop-size 20 --n-gen 10

    # Full optimization (5-fold CV, all 218 cases)
    uv run python -m benchmarks.optimize_params --pop-size 60 --n-gen 80 --n-threads 4

    # Faster: 3-fold CV
    uv run python -m benchmarks.optimize_params --n-folds 3 --pop-size 40 --n-gen 50

    # Add wall time as 3rd objective
    uv run python -m benchmarks.optimize_params --optimize-time
"""

import argparse
import os
import pickle
import signal
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from pymoo.algorithms.moo.nsga2 import NSGA2  # type: ignore[import-untyped]
    from pymoo.core.callback import Callback  # type: ignore[import-untyped]
    from pymoo.core.problem import Problem  # type: ignore[import-untyped]
    from pymoo.operators.crossover.sbx import SBX  # type: ignore[import-untyped]
    from pymoo.operators.mutation.pm import PM  # type: ignore[import-untyped]
    from pymoo.operators.sampling.lhs import LHS  # type: ignore[import-untyped]
    from pymoo.optimize import minimize  # type: ignore[import-untyped]
    from pymoo.termination import get_termination  # type: ignore[import-untyped]
except ImportError:
    print("pymoo not installed. Run: uv pip install pymoo")
    sys.exit(1)

from rich.console import Console  # type: ignore[import-untyped]
from rich.layout import Layout  # type: ignore[import-untyped]
from rich.live import Live  # type: ignore[import-untyped]
from rich.panel import Panel  # type: ignore[import-untyped]
from rich.progress import (  # type: ignore[import-untyped]
    BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn,
)
from rich.table import Table  # type: ignore[import-untyped]
from rich.text import Text  # type: ignore[import-untyped]

import kalign  # type: ignore[import-untyped]

from .datasets import BenchmarkCase, balibase_cases, balibase_download, balibase_is_available
from .scoring import score_alignment_detailed

# ---------------------------------------------------------------------------
# Parameter space definition
# ---------------------------------------------------------------------------

# Continuous variables: [gpo, gpe, tgpe, vsm_amax, seq_weights, consistency_weight]
CONT_LOWER = np.array([2.0, 0.5, 0.1, 0.0, 0.0, 0.5])
CONT_UPPER = np.array([15.0, 5.0, 3.0, 5.0, 5.0, 5.0])
CONT_NAMES = ["gpo", "gpe", "tgpe", "vsm_amax", "seq_weights", "consistency_weight"]

# Integer/categorical variables: [consistency, realign, matrix]
INT_LOWER = np.array([0, 0, 0])
INT_UPPER = np.array([6, 2, 2])  # consistency: 0-6 maps via table, realign: 0-2, matrix: 0-2
INT_NAMES = ["consistency", "realign", "matrix"]

CONSISTENCY_MAP = [0, 1, 2, 3, 5, 8, 10]
MATRIX_MAP = ["pfasum43", "pfasum60", "protein"]  # API values for kalign
MATRIX_DISPLAY = {"pfasum43": "PFASUM43", "pfasum60": "PFASUM60", "protein": "CorBLOSUM66"}

N_CONT = len(CONT_LOWER)
N_INT = len(INT_LOWER)
N_VAR = N_CONT + N_INT


def decode_params(x):
    """Decode a parameter vector into a dict of kalign parameters."""
    cont = x[:N_CONT]
    ints = np.round(x[N_CONT:]).astype(int)

    consistency_idx = int(np.clip(ints[0], 0, len(CONSISTENCY_MAP) - 1))
    realign_idx = int(np.clip(ints[1], 0, 2))
    matrix_idx = int(np.clip(ints[2], 0, len(MATRIX_MAP) - 1))

    return {
        "gap_open": float(cont[0]),
        "gap_extend": float(cont[1]),
        "terminal_gap_extend": float(cont[2]),
        "vsm_amax": float(cont[3]),
        "seq_weights": float(cont[4]),
        "consistency_weight": float(cont[5]),
        "consistency": CONSISTENCY_MAP[consistency_idx],
        "realign": realign_idx,
        "seq_type": MATRIX_MAP[matrix_idx],
    }


def encode_params(params):
    """Encode a params dict back into a decision vector (inverse of decode_params)."""
    cont = np.array([
        params["gap_open"], params["gap_extend"], params["terminal_gap_extend"],
        params["vsm_amax"], params["seq_weights"], params["consistency_weight"],
    ])
    consistency_idx = CONSISTENCY_MAP.index(params["consistency"])
    realign_idx = params["realign"]
    matrix_idx = MATRIX_MAP.index(params["seq_type"])
    ints = np.array([consistency_idx, realign_idx, matrix_idx], dtype=float)
    return np.concatenate([cont, ints])


def _matrix_name(params):
    """Human-readable matrix name from seq_type API value."""
    return MATRIX_DISPLAY.get(params["seq_type"], params["seq_type"])


def format_params_short(params):
    """Compact one-line summary of a parameter dict."""
    return (f"gpo={params['gap_open']:.1f} gpe={params['gap_extend']:.2f} "
            f"tgpe={params['terminal_gap_extend']:.2f} vsm={params['vsm_amax']:.1f} "
            f"sw={params['seq_weights']:.1f} c={params['consistency']} "
            f"cw={params['consistency_weight']:.1f} re={params['realign']} "
            f"{_matrix_name(params)}")


def format_params_long(params):
    """Verbose one-line summary."""
    return (f"gpo={params['gap_open']:.2f} gpe={params['gap_extend']:.2f} "
            f"tgpe={params['terminal_gap_extend']:.2f} vsm={params['vsm_amax']:.2f} "
            f"sw={params['seq_weights']:.2f} cons={params['consistency']} "
            f"cw={params['consistency_weight']:.2f} re={params['realign']} "
            f"mat={_matrix_name(params)}")


# ---------------------------------------------------------------------------
# Stratified k-fold CV
# ---------------------------------------------------------------------------

def stratified_kfold(cases: List[BenchmarkCase], k: int, seed: int = 42
                     ) -> List[Tuple[List[BenchmarkCase], List[BenchmarkCase]]]:
    """Split cases into k stratified folds by dataset (RV category).

    Returns list of (train, test) pairs.  Each fold's test set contains
    roughly equal representation from every RV category.
    """
    rng = np.random.RandomState(seed)

    # Group by category
    by_cat: Dict[str, List[BenchmarkCase]] = defaultdict(list)
    for c in cases:
        by_cat[c.dataset].append(c)

    # Shuffle within each category and assign fold indices
    fold_assignments: Dict[str, int] = {}
    for cat_cases in by_cat.values():
        indices = list(range(len(cat_cases)))
        rng.shuffle(indices)
        for rank, idx in enumerate(indices):
            fold_assignments[cat_cases[idx].family] = rank % k

    # Build folds
    folds = []
    for fold_idx in range(k):
        test = [c for c in cases if fold_assignments[c.family] == fold_idx]
        train = [c for c in cases if fold_assignments[c.family] != fold_idx]
        folds.append((train, test))

    return folds


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate_params(params, cases, n_threads=1, quiet=True):
    """Run kalign with given params on all cases, return mean metrics.

    Returns dict with keys: f1, tc, recall, precision, wall_time, per_category.
    """
    results_by_cat: Dict[str, list] = {}
    total_time = 0.0

    for case in cases:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / f"{case.family}_aln.fasta"

            try:
                start = time.perf_counter()
                kalign.align_file_to_file(
                    str(case.unaligned),
                    str(output),
                    format="fasta",
                    seq_type=params["seq_type"],
                    gap_open=params["gap_open"],
                    gap_extend=params["gap_extend"],
                    terminal_gap_extend=params["terminal_gap_extend"],
                    n_threads=n_threads,
                    vsm_amax=params["vsm_amax"],
                    seq_weights=params["seq_weights"],
                    consistency=params["consistency"],
                    consistency_weight=params["consistency_weight"],
                    realign=params["realign"],
                )
                wall_time = time.perf_counter() - start
                total_time += wall_time

                detailed = score_alignment_detailed(case.reference, output)

                cat = case.dataset
                if cat not in results_by_cat:
                    results_by_cat[cat] = []
                results_by_cat[cat].append(detailed)

            except Exception as e:
                if not quiet:
                    print(f"  WARN: {case.family}: {e}", file=sys.stderr)

    if not results_by_cat:
        return {"f1": 0.0, "tc": 0.0, "recall": 0.0, "precision": 0.0,
                "wall_time": total_time, "per_category": {}}

    # Compute per-category means
    per_cat = {}
    for cat, scores in results_by_cat.items():
        per_cat[cat] = {
            "f1": np.mean([s["f1"] for s in scores]),
            "tc": np.mean([s["tc"] for s in scores]),
            "recall": np.mean([s["recall"] for s in scores]),
            "precision": np.mean([s["precision"] for s in scores]),
            "n": len(scores),
        }

    # Overall means (category-averaged so small categories count equally)
    all_f1 = [v["f1"] for v in per_cat.values()]
    all_tc = [v["tc"] for v in per_cat.values()]
    all_recall = [v["recall"] for v in per_cat.values()]
    all_precision = [v["precision"] for v in per_cat.values()]

    return {
        "f1": float(np.mean(all_f1)),
        "tc": float(np.mean(all_tc)),
        "recall": float(np.mean(all_recall)),
        "precision": float(np.mean(all_precision)),
        "wall_time": total_time,
        "per_category": per_cat,
    }


def evaluate_cv(params, folds, n_threads=1, quiet=True):
    """Evaluate params using stratified k-fold CV.

    For each fold, aligns the test cases and scores them.
    Returns the mean held-out F1 and TC across folds,
    plus per-fold details and total wall time.
    """
    fold_f1s = []
    fold_tcs = []
    total_time = 0.0

    for _, test in folds:
        result = evaluate_params(params, test, n_threads, quiet)
        fold_f1s.append(result["f1"])
        fold_tcs.append(result["tc"])
        total_time += result["wall_time"]

    return {
        "f1": float(np.mean(fold_f1s)),
        "tc": float(np.mean(fold_tcs)),
        "f1_std": float(np.std(fold_f1s)),
        "tc_std": float(np.std(fold_tcs)),
        "fold_f1s": fold_f1s,
        "fold_tcs": fold_tcs,
        "wall_time": total_time,
    }


# ---------------------------------------------------------------------------
# Rich live dashboard
# ---------------------------------------------------------------------------

class Dashboard:
    """Rich-based live terminal dashboard for optimization progress."""

    def __init__(self, n_gen: int, pop_size: int, baseline_f1: float, baseline_tc: float,
                 optimize_time: bool = False, baseline_time: float = 0.0):
        self.n_gen = n_gen
        self.pop_size = pop_size
        self.baseline_f1 = baseline_f1
        self.baseline_tc = baseline_tc
        self.optimize_time = optimize_time
        self.baseline_time = baseline_time
        self.console = Console()

        # State
        self.current_gen = 0
        self.eval_count = 0
        self.total_evals = n_gen * pop_size
        self.gen_start_time = time.time()
        self.run_start_time = time.time()
        self.current_eval_params: Optional[dict] = None
        self.current_eval_idx = 0  # within generation

        # Best-ever tracking
        self.best_f1 = 0.0
        self.best_f1_params: Optional[dict] = None
        self.best_tc = 0.0
        self.best_tc_params: Optional[dict] = None
        self.best_time = float("inf")
        self.best_time_params: Optional[dict] = None

        # Pareto front (updated per generation)
        self.pareto_front: List[dict] = []

        # Generation history (best F1 per gen for trend)
        self.gen_history: List[dict] = []

        # Recent evaluations (ring buffer for last 5)
        self.recent_evals: List[dict] = []

        # Progress bar
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(bar_width=40),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
        )
        self.gen_task = self.progress.add_task("Generation", total=n_gen)
        self.eval_task = self.progress.add_task("  Eval", total=pop_size)

        self.live = Live(self._build_layout(), console=self.console, refresh_per_second=2)

    def start(self):
        self.run_start_time = time.time()
        self.live.start()

    def stop(self):
        self.live.stop()

    def _delta_str(self, value: float, baseline: float) -> str:
        """Format a value with delta vs baseline, colored."""
        delta = value - baseline
        if delta > 0.005:
            return f"[bold green]{value:.4f} (+{delta:.4f})[/]"
        elif delta < -0.005:
            return f"[bold red]{value:.4f} ({delta:.4f})[/]"
        else:
            return f"[dim]{value:.4f} ({delta:+.4f})[/]"

    def _build_status_panel(self) -> Panel:
        elapsed = time.time() - self.run_start_time
        gen_elapsed = time.time() - self.gen_start_time

        if self.eval_count > 0:
            avg_per_eval = elapsed / self.eval_count
            remaining = avg_per_eval * (self.total_evals - self.eval_count)
            eta_str = f"{remaining / 60:.0f}m" if remaining > 60 else f"{remaining:.0f}s"
        else:
            eta_str = "..."

        lines = [
            f"Gen [bold]{self.current_gen}[/]/{self.n_gen}  "
            f"Eval [bold]{self.eval_count}[/]/{self.total_evals}  "
            f"Elapsed [bold]{elapsed / 60:.1f}m[/]  "
            f"ETA [bold]{eta_str}[/]  "
            f"Gen time [bold]{gen_elapsed:.1f}s[/]",
        ]

        if self.current_eval_params:
            lines.append(f"[dim]Current: {format_params_short(self.current_eval_params)}[/]")

        return Panel("\n".join(lines), title="Progress", border_style="blue")

    def _build_best_panel(self) -> Panel:
        lines = []
        baseline_str = f"[dim]Baseline:  F1={self.baseline_f1:.4f}  TC={self.baseline_tc:.4f}"
        if self.optimize_time:
            baseline_str += f"  Time={self.baseline_time:.1f}s"
        baseline_str += "[/]"
        lines.append(baseline_str)
        lines.append("")

        if self.best_f1_params:
            lines.append(f"Best F1:   {self._delta_str(self.best_f1, self.baseline_f1)}")
            lines.append(f"  [dim]{format_params_short(self.best_f1_params)}[/]")
        else:
            lines.append("[dim]Best F1:   (no evaluations yet)[/]")

        lines.append("")

        if self.best_tc_params:
            lines.append(f"Best TC:   {self._delta_str(self.best_tc, self.baseline_tc)}")
            lines.append(f"  [dim]{format_params_short(self.best_tc_params)}[/]")
        else:
            lines.append("[dim]Best TC:   (no evaluations yet)[/]")

        if self.optimize_time:
            lines.append("")
            if self.best_time_params:
                delta = self.best_time - self.baseline_time
                if delta < -1.0:
                    time_str = f"[bold green]{self.best_time:.1f}s ({delta:.1f}s)[/]"
                elif delta > 1.0:
                    time_str = f"[bold red]{self.best_time:.1f}s (+{delta:.1f}s)[/]"
                else:
                    time_str = f"[dim]{self.best_time:.1f}s ({delta:+.1f}s)[/]"
                lines.append(f"Fastest:   {time_str}")
                lines.append(f"  [dim]{format_params_short(self.best_time_params)}[/]")

        return Panel("\n".join(lines), title="Best Solutions (CV)", border_style="green")

    def _build_pareto_table(self) -> Panel:
        table = Table(title="Pareto Front", expand=True, show_lines=False, padding=(0, 1))
        table.add_column("#", style="dim", width=3)
        table.add_column("CV F1", justify="right", width=8)
        table.add_column("CV TC", justify="right", width=8)
        if self.optimize_time:
            table.add_column("Time", justify="right", width=7)
        table.add_column("gpo", justify="right", width=5)
        table.add_column("gpe", justify="right", width=5)
        table.add_column("tgpe", justify="right", width=5)
        table.add_column("vsm", justify="right", width=4)
        table.add_column("sw", justify="right", width=4)
        table.add_column("c", justify="right", width=2)
        table.add_column("cw", justify="right", width=4)
        table.add_column("re", justify="right", width=2)
        table.add_column("mat", width=8)

        n_cols = 13 if self.optimize_time else 12

        # Sort by F1 descending
        sorted_front = sorted(self.pareto_front, key=lambda x: -x["f1"])

        for i, sol in enumerate(sorted_front[:10]):  # show top 10
            p = sol["params"]
            f1_style = "bold green" if sol["f1"] > self.baseline_f1 else ""
            tc_style = "bold green" if sol["tc"] > self.baseline_tc else ""
            row = [
                str(i),
                Text(f"{sol['f1']:.4f}", style=f1_style),
                Text(f"{sol['tc']:.4f}", style=tc_style),
            ]
            if self.optimize_time:
                wt = sol.get("wall_time", 0.0)
                row.append(f"{wt:.1f}s" if wt else "?")
            row.extend([
                f"{p['gap_open']:.1f}",
                f"{p['gap_extend']:.2f}",
                f"{p['terminal_gap_extend']:.2f}",
                f"{p['vsm_amax']:.1f}",
                f"{p['seq_weights']:.1f}",
                str(p['consistency']),
                f"{p['consistency_weight']:.1f}",
                str(p['realign']),
                _matrix_name(p),
            ])
            table.add_row(*row)

        if not self.pareto_front:
            table.add_row(*[""] * n_cols)

        return Panel(table, border_style="cyan")

    def _build_trend_panel(self) -> Panel:
        if not self.gen_history:
            return Panel("[dim]No generation data yet[/]", title="Trend", border_style="yellow")

        lines = []
        # Show sparkline-style trend using simple chars
        bar_width = 30
        if len(self.gen_history) > 1:
            f1_vals = [g["best_f1"] for g in self.gen_history]
            lo = min(min(f1_vals), self.baseline_f1) - 0.01
            hi = max(max(f1_vals), self.baseline_f1) + 0.01
            rng = hi - lo if hi > lo else 1.0

            # Baseline marker position
            bl_pos = int((self.baseline_f1 - lo) / rng * bar_width)

            lines.append("[bold]Best CV F1 per generation:[/]")
            for g in self.gen_history[-12:]:  # last 12 gens
                pos = int((g["best_f1"] - lo) / rng * bar_width)
                bar = list("·" * bar_width)
                bar[min(bl_pos, bar_width - 1)] = "│"  # baseline marker
                for j in range(min(pos, bar_width)):
                    bar[j] = "█"
                bar_str = "".join(bar)
                delta = g["best_f1"] - self.baseline_f1
                color = "green" if delta > 0 else "red" if delta < -0.005 else "dim"
                lines.append(f"  Gen {g['gen']:>3d} [{color}]{bar_str} {g['best_f1']:.4f}[/]")

            lines.append(f"  [dim]│ = baseline ({self.baseline_f1:.4f})[/]")
        else:
            g = self.gen_history[0]
            lines.append(f"Gen {g['gen']}: F1={g['best_f1']:.4f} TC={g['best_tc']:.4f}")

        return Panel("\n".join(lines), title="Trend", border_style="yellow")

    def _build_recent_panel(self) -> Panel:
        if not self.recent_evals:
            return Panel("[dim]No evaluations yet[/]", title="Recent", border_style="dim")

        lines = []
        for ev in self.recent_evals[-5:]:
            f1 = ev["f1"]
            tc = ev["tc"]
            delta = f1 - self.baseline_f1
            color = "green" if delta > 0 else "red" if delta < -0.005 else "dim"
            time_str = f" t={ev['wall_time']:.1f}s" if self.optimize_time else ""
            lines.append(
                f"  [{color}]F1={f1:.4f} TC={tc:.4f}{time_str}[/] "
                f"[dim]{format_params_short(ev['params'])}[/]"
            )
        return Panel("\n".join(lines), title="Recent Evaluations", border_style="dim")

    def _build_layout(self):
        layout = Layout()
        layout.split_column(
            Layout(self.progress, name="progress", size=3),
            Layout(name="top", size=7),
            Layout(name="middle"),
            Layout(self._build_recent_panel(), name="bottom", size=8),
        )
        layout["top"].split_row(
            Layout(self._build_status_panel(), name="status"),
            Layout(self._build_best_panel(), name="best"),
        )
        layout["middle"].split_row(
            Layout(self._build_pareto_table(), name="pareto", ratio=3),
            Layout(self._build_trend_panel(), name="trend", ratio=2),
        )
        return layout

    def refresh(self):
        self.live.update(self._build_layout())

    def on_eval_start(self, params: dict, eval_num: int, eval_in_gen: int):
        self.eval_count = eval_num
        self.current_eval_idx = eval_in_gen
        self.current_eval_params = params
        self.progress.update(self.eval_task, completed=eval_in_gen)
        self.refresh()

    def on_eval_end(self, params: dict, cv_result: dict):
        f1 = cv_result["f1"]
        tc = cv_result["tc"]
        wt = cv_result.get("wall_time", 0.0)

        self.recent_evals.append({"params": params, "f1": f1, "tc": tc, "wall_time": wt})
        if len(self.recent_evals) > 5:
            self.recent_evals.pop(0)

        if f1 > self.best_f1:
            self.best_f1 = f1
            self.best_f1_params = params
        if tc > self.best_tc:
            self.best_tc = tc
            self.best_tc_params = params
        if self.optimize_time and wt < self.best_time and f1 > 0.5:
            self.best_time = wt
            self.best_time_params = params

        self.refresh()

    def on_gen_start(self, gen: int):
        self.current_gen = gen
        self.gen_start_time = time.time()
        self.progress.update(self.gen_task, completed=gen)
        self.progress.update(self.eval_task, completed=0)
        self.refresh()

    def on_gen_end(self, gen: int, pareto_front: List[dict]):
        self.pareto_front = pareto_front

        best_f1_in_gen = max((s["f1"] for s in pareto_front), default=0.0)
        best_tc_in_gen = max((s["tc"] for s in pareto_front), default=0.0)
        self.gen_history.append({
            "gen": gen,
            "best_f1": best_f1_in_gen,
            "best_tc": best_tc_in_gen,
            "n_pareto": len(pareto_front),
        })

        self.progress.update(self.gen_task, completed=gen)
        self.refresh()


# ---------------------------------------------------------------------------
# Parallel evaluation helper (must be top-level for pickling)
# ---------------------------------------------------------------------------

def _eval_one(args_tuple):
    """Evaluate one parameter vector.  Top-level function for ProcessPoolExecutor."""
    x, folds, n_threads = args_tuple
    params = decode_params(x)
    cv_result = evaluate_cv(params, folds, n_threads, quiet=True)
    return params, cv_result


def _eval_one_fold(args_tuple):
    """Evaluate one (individual, fold) pair. Finer-grained parallelism."""
    ind_idx, fold_idx, x, test_cases, n_threads = args_tuple
    params = decode_params(x)
    result = evaluate_params(params, test_cases, n_threads, quiet=True)
    return ind_idx, fold_idx, params, result


def _kill_pool(pool: ProcessPoolExecutor) -> None:
    """Forcefully terminate all worker processes in the pool."""
    # Access the internal process map to send SIGTERM to each worker
    for pid in list(pool._processes):  # noqa: SLF001
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    pool.shutdown(wait=False, cancel_futures=True)


# ---------------------------------------------------------------------------
# pymoo Problem + Callback
# ---------------------------------------------------------------------------

class KalignCVProblem(Problem):
    """Multi-objective optimization with stratified CV evaluation."""

    def __init__(self, folds, n_threads=1, n_workers=1, optimize_time=False,
                 dashboard: Optional[Dashboard] = None):
        n_obj = 3 if optimize_time else 2
        xl = np.concatenate([CONT_LOWER, INT_LOWER.astype(float)])
        xu = np.concatenate([CONT_UPPER, INT_UPPER.astype(float)])

        super().__init__(
            n_var=N_VAR,
            n_obj=n_obj,
            xl=xl,
            xu=xu,
        )
        self.folds = folds
        self.n_threads = n_threads
        self.n_workers = n_workers
        self.optimize_time = optimize_time
        self.dashboard = dashboard
        self.eval_count = 0
        self.history: List[dict] = []

    def _evaluate(self, X, out, *args, **kwargs):  # pyright: ignore[reportUnusedVariable]
        F = np.zeros((X.shape[0], self.n_obj))

        if self.n_workers > 1:
            self._evaluate_parallel(X, F)
        else:
            self._evaluate_serial(X, F)

        out["F"] = F

    def _evaluate_serial(self, X, F):
        for i, x in enumerate(X):
            params = decode_params(x)
            self.eval_count += 1

            if self.dashboard:
                self.dashboard.on_eval_start(params, self.eval_count, i)

            cv_result = evaluate_cv(params, self.folds, self.n_threads, quiet=True)
            self._record(i, F, params, cv_result)

    def _evaluate_parallel(self, X, F):
        """Fine-grained parallelism: submit (individual × fold) jobs.

        For pop_size=60, k=5 folds → 300 jobs, much better load balancing
        than 60 coarse-grained jobs where slow individuals block the generation.
        """
        n_folds = len(self.folds)
        n_pop = X.shape[0]

        # Build flat job list: (ind_idx, fold_idx, x, test_cases, n_threads)
        jobs = []
        for i, x in enumerate(X):
            for fold_idx, (_, test) in enumerate(self.folds):
                jobs.append((i, fold_idx, x, test, self.n_threads))

        if self.dashboard:
            self.dashboard.on_eval_start(
                {"gap_open": 0, "gap_extend": 0, "terminal_gap_extend": 0,
                 "vsm_amax": 0, "seq_weights": 0, "consistency_weight": 0,
                 "consistency": 0, "realign": 0, "seq_type": "..."},
                self.eval_count + 1, 0)

        # Collect per-fold results: fold_results[ind_idx] = {fold_idx: result}
        fold_results: Dict[int, Dict[int, dict]] = defaultdict(dict)
        ind_params: Dict[int, dict] = {}

        pool = ProcessPoolExecutor(max_workers=self.n_workers)
        try:
            futures = {pool.submit(_eval_one_fold, j): j[:2] for j in jobs}
            completed_individuals = set()
            completed_jobs = 0

            for future in as_completed(futures):
                ind_idx, fold_idx, params, result = future.result()
                fold_results[ind_idx][fold_idx] = result
                ind_params[ind_idx] = params
                completed_jobs += 1

                # Check if this individual now has all folds complete
                if len(fold_results[ind_idx]) == n_folds:
                    completed_individuals.add(ind_idx)
                    self.eval_count += 1

                    # Aggregate fold results into CV score
                    fold_f1s = [fold_results[ind_idx][fi]["f1"] for fi in range(n_folds)]
                    fold_tcs = [fold_results[ind_idx][fi]["tc"] for fi in range(n_folds)]
                    total_time = sum(fold_results[ind_idx][fi]["wall_time"] for fi in range(n_folds))

                    cv_result = {
                        "f1": float(np.mean(fold_f1s)),
                        "tc": float(np.mean(fold_tcs)),
                        "f1_std": float(np.std(fold_f1s)),
                        "tc_std": float(np.std(fold_tcs)),
                        "fold_f1s": fold_f1s,
                        "fold_tcs": fold_tcs,
                        "wall_time": total_time,
                    }

                    self._record(ind_idx, F, params, cv_result)

                    if self.dashboard:
                        self.dashboard.on_eval_start(
                            params, self.eval_count, len(completed_individuals))

        except KeyboardInterrupt:
            _kill_pool(pool)
            raise
        finally:
            pool.shutdown(wait=False)

    def _record(self, i, F, params, cv_result):
        F[i, 0] = -cv_result["f1"]
        F[i, 1] = -cv_result["tc"]
        if self.optimize_time:
            F[i, 2] = cv_result["wall_time"]

        self.history.append({
            "eval": self.eval_count,
            "params": params,
            "cv_result": cv_result,
        })

        if self.dashboard:
            self.dashboard.on_eval_end(params, cv_result)


class GenerationCallback(Callback):
    """pymoo callback: updates dashboard + saves checkpoint after each generation."""

    def __init__(self, dashboard: Optional[Dashboard] = None,
                 checkpoint_path: Optional[Path] = None,
                 problem: Optional["KalignCVProblem"] = None):
        super().__init__()
        self.dashboard = dashboard
        self.checkpoint_path = checkpoint_path
        self.problem = problem

    def notify(self, algorithm):
        gen = algorithm.n_gen

        if self.dashboard:
            self.dashboard.on_gen_start(gen)

        # Extract Pareto front from algorithm.opt
        pareto = []
        if algorithm.opt is not None and len(algorithm.opt) > 0:
            for ind in algorithm.opt:
                params = decode_params(ind.X)
                entry = {
                    "params": params,
                    "f1": -ind.F[0],
                    "tc": -ind.F[1],
                }
                if len(ind.F) > 2:
                    entry["wall_time"] = ind.F[2]
                pareto.append(entry)

        if self.dashboard:
            self.dashboard.on_gen_end(gen, pareto)

        # Save checkpoint after every generation
        if self.checkpoint_path:
            # Save population state (X, F) — not the algorithm object (has unpicklable locks)
            pop = algorithm.pop
            ckpt = {
                "pop_X": pop.get("X").copy(),
                "pop_F": pop.get("F").copy(),
                "n_gen_completed": gen,
                "history": self.problem.history if self.problem else [],
                "pop_size": len(pop),
            }
            # Write to temp file then rename for atomicity
            tmp = self.checkpoint_path.with_suffix(".tmp")
            with open(tmp, "wb") as f:
                pickle.dump(ckpt, f)
            tmp.rename(self.checkpoint_path)


def load_checkpoint(path: Path):
    """Load a generation checkpoint. Returns (pop_X, pop_F, n_gen_completed, history)."""
    with open(path, "rb") as f:
        ckpt = pickle.load(f)  # noqa: S301
    return ckpt["pop_X"], ckpt["pop_F"], ckpt["n_gen_completed"], ckpt.get("history", [])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Optimize kalign hyperparameters with NSGA-II + stratified k-fold CV")
    parser.add_argument("--pop-size", type=int, default=40,
                        help="Population size (default: 40)")
    parser.add_argument("--n-gen", type=int, default=50,
                        help="Number of generations (default: 50)")
    parser.add_argument("--n-folds", type=int, default=5,
                        help="Number of CV folds (default: 5)")
    parser.add_argument("--n-threads", type=int, default=1,
                        help="OpenMP threads per kalign alignment (default: 1)")
    parser.add_argument("--n-workers", type=int, default=1,
                        help="Parallel worker processes for evaluating individuals (default: 1)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--optimize-time", action="store_true",
                        help="Add wall time as 3rd objective")
    parser.add_argument("--output-dir", type=str, default="benchmarks/results/optim",
                        help="Output directory for results")
    parser.add_argument("--run-name", type=str, default=None,
                        help="Name for this run (creates subdirectory, e.g. 'run2_time')")
    parser.add_argument("--no-dashboard", action="store_true",
                        help="Disable rich dashboard, use plain text output")
    parser.add_argument("--resume", type=str, default=None,
                        help="Resume from a generation checkpoint file (.pkl)")
    args = parser.parse_args()

    console = Console()

    # Setup
    if not balibase_is_available():
        console.print("Downloading BAliBASE...")
        balibase_download()

    cases = balibase_cases()
    console.print(f"Loaded [bold]{len(cases)}[/] BAliBASE cases")

    # Show category distribution
    cats: Dict[str, int] = {}
    for c in cases:
        cats[c.dataset] = cats.get(c.dataset, 0) + 1
    for cat, n in sorted(cats.items()):
        console.print(f"  {cat}: {n} cases")

    # Build stratified folds
    k = args.n_folds
    folds = stratified_kfold(cases, k, seed=args.seed)
    console.print(f"\nStratified [bold]{k}[/]-fold CV:")
    for i, (train, test) in enumerate(folds):
        test_cats: Dict[str, int] = defaultdict(int)
        for c in test:
            test_cats[c.dataset] += 1
        cat_str = ", ".join(f"{cat.replace('balibase_', '')}:{n}"
                            for cat, n in sorted(test_cats.items()))
        console.print(f"  Fold {i}: {len(test)} test / {len(train)} train  ({cat_str})")

    output_dir = Path(args.output_dir)
    if args.run_name:
        output_dir = output_dir / args.run_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Baseline evaluation (CV) --- parallelized across folds
    # Best F1 from first optimization run (pop_size=60, n_gen=100)
    console.print(f"\n[bold]Baseline evaluation[/] (best from run 1, {k}-fold CV)")
    baseline_params = {
        "gap_open": 8.472, "gap_extend": 0.554, "terminal_gap_extend": 0.409,
        "vsm_amax": 1.359, "seq_weights": 3.407, "consistency_weight": 1.167,
        "consistency": 8, "realign": 2, "seq_type": "pfasum60",
    }
    n_baseline_workers = max(1, args.n_workers)
    if n_baseline_workers > 1:
        # Run CV folds + full-dataset eval in parallel
        baseline_jobs = []
        for fi, (_, test) in enumerate(folds):
            baseline_jobs.append((0, fi, encode_params(baseline_params), test, args.n_threads))
        # Add full-dataset as an extra "fold"
        baseline_jobs.append((0, k, encode_params(baseline_params), cases, args.n_threads))
        console.print(f"  Running {len(baseline_jobs)} baseline jobs in parallel ({n_baseline_workers} workers)...")
        with ProcessPoolExecutor(max_workers=min(n_baseline_workers, len(baseline_jobs))) as pool:
            futures = {pool.submit(_eval_one_fold, j): j[1] for j in baseline_jobs}
            fold_results_bl = {}
            for future in as_completed(futures):
                _, fi, _, result = future.result()
                fold_results_bl[fi] = result
        fold_f1s = [fold_results_bl[fi]["f1"] for fi in range(k)]
        fold_tcs = [fold_results_bl[fi]["tc"] for fi in range(k)]
        total_time = sum(fold_results_bl[fi]["wall_time"] for fi in range(k))
        baseline_cv = {
            "f1": float(np.mean(fold_f1s)), "tc": float(np.mean(fold_tcs)),
            "f1_std": float(np.std(fold_f1s)), "tc_std": float(np.std(fold_tcs)),
            "fold_f1s": fold_f1s, "fold_tcs": fold_tcs, "wall_time": total_time,
        }
        baseline_full = fold_results_bl[k]  # the extra "fold" with all cases
    else:
        baseline_cv = evaluate_cv(baseline_params, folds, args.n_threads)
        baseline_full = evaluate_params(baseline_params, cases, args.n_threads)

    console.print(f"  CV F1=[bold]{baseline_cv['f1']:.4f}[/]±{baseline_cv['f1_std']:.4f}  "
                  f"CV TC=[bold]{baseline_cv['tc']:.4f}[/]±{baseline_cv['tc_std']:.4f}  "
                  f"time={baseline_cv['wall_time']:.1f}s")
    for i, (f1, tc) in enumerate(zip(baseline_cv['fold_f1s'], baseline_cv['fold_tcs'])):
        console.print(f"    Fold {i}: F1={f1:.4f} TC={tc:.4f}")

    console.print(f"  Full-data F1=[bold]{baseline_full['f1']:.4f}[/]  "
                  f"TC=[bold]{baseline_full['tc']:.4f}[/]")
    for cat, v in sorted(baseline_full["per_category"].items()):
        console.print(f"    {cat}: F1={v['f1']:.4f} TC={v['tc']:.4f} (n={v['n']})")

    # --- Optimization ---
    n_evals = args.pop_size * args.n_gen
    est_sec_per_eval = baseline_cv["wall_time"]
    parallelism = max(1, args.n_workers)
    est_hours = n_evals * est_sec_per_eval / parallelism / 3600

    console.print(f"\n[bold]Starting NSGA-II[/]: pop_size={args.pop_size}, n_gen={args.n_gen}, "
                  f"{k}-fold CV, {args.n_workers} worker(s) × {args.n_threads} thread(s)")
    console.print(f"Total evaluations: ~{n_evals} ({n_evals * len(cases)} alignments)")
    console.print(f"Estimated time: ~{est_hours:.1f} hours "
                  f"(~{est_sec_per_eval:.0f}s per eval, "
                  f"{parallelism}× parallel)\n")

    # Set up dashboard or plain mode
    use_dashboard = not args.no_dashboard
    dashboard = None

    if use_dashboard:
        dashboard = Dashboard(
            n_gen=args.n_gen,
            pop_size=args.pop_size,
            baseline_f1=baseline_cv["f1"],
            baseline_tc=baseline_cv["tc"],
            optimize_time=args.optimize_time,
            baseline_time=baseline_cv["wall_time"],
        )

    problem = KalignCVProblem(
        folds=folds,
        n_threads=args.n_threads,
        n_workers=args.n_workers,
        optimize_time=args.optimize_time,
        dashboard=dashboard,
    )

    checkpoint_path = output_dir / "gen_checkpoint.pkl"
    callback = GenerationCallback(
        dashboard=dashboard,
        checkpoint_path=checkpoint_path,
        problem=problem,
    )

    # Resume from checkpoint or start fresh
    resumed_gen = 0
    if args.resume:
        resume_path = Path(args.resume)
        if not resume_path.exists():
            console.print(f"[bold red]Checkpoint not found:[/] {resume_path}")
            return
        pop_X, _pop_F, resumed_gen, prev_history = load_checkpoint(resume_path)
        problem.history = prev_history
        console.print(f"[bold green]Resumed[/] from generation {resumed_gen} "
                      f"({len(prev_history)} prior evaluations)")
        remaining = args.n_gen - resumed_gen
        if remaining <= 0:
            console.print(f"[bold yellow]Already completed {resumed_gen} generations "
                          f"(requested {args.n_gen}). Increase --n-gen to continue.[/]")
            return
        termination = get_termination("n_gen", remaining)
        # Reconstruct algorithm with saved population as initial sampling
        algorithm = NSGA2(
            pop_size=len(pop_X),
            sampling=pop_X,
            crossover=SBX(prob=0.9, eta=15),
            mutation=PM(eta=20),
            eliminate_duplicates=True,
        )
    else:
        algorithm = NSGA2(
            pop_size=args.pop_size,
            sampling=LHS(),
            crossover=SBX(prob=0.9, eta=15),
            mutation=PM(eta=20),
            eliminate_duplicates=True,
        )
        termination = get_termination("n_gen", args.n_gen)

    if dashboard:
        dashboard.start()

    try:
        res = minimize(
            problem,
            algorithm,
            termination,
            seed=args.seed,
            verbose=not use_dashboard,  # disable pymoo's own output when dashboard is active
            callback=callback,
        )
    except KeyboardInterrupt:
        if dashboard:
            dashboard.stop()
        console.print("\n[bold yellow]Interrupted![/] Checkpoint was saved after last "
                      f"completed generation to:\n  {checkpoint_path}")
        console.print("Resume with: [bold]--resume " + str(checkpoint_path) + "[/]")
        os._exit(1)  # noqa: SLF001 — skip atexit handlers that hang on worker join
    finally:
        if dashboard:
            dashboard.stop()

    # --- Results ---
    console.print(f"\n[bold]Optimization complete.[/] "
                  f"{len(res.F)} Pareto-optimal solutions found.\n")

    pareto_configs = []
    for i, (x, f) in enumerate(zip(res.X, res.F)):
        params = decode_params(x)
        f1 = -f[0]
        tc = -f[1]
        wt = f[2] if args.optimize_time else None
        pareto_configs.append({"params": params, "f1_cv": f1, "tc_cv": tc, "wall_time": wt})

    # Print Pareto front as a rich table
    table = Table(title="Pareto Front (sorted by CV F1)")
    table.add_column("#", style="dim", width=3)
    table.add_column("CV F1", justify="right")
    table.add_column("CV TC", justify="right")
    table.add_column("Parameters")

    sorted_pareto = sorted(pareto_configs, key=lambda x: -x["f1_cv"])
    for i, cfg in enumerate(sorted_pareto):
        f1_style = "bold green" if cfg["f1_cv"] > baseline_cv["f1"] else ""
        tc_style = "bold green" if cfg["tc_cv"] > baseline_cv["tc"] else ""
        table.add_row(
            str(i),
            Text(f"{cfg['f1_cv']:.4f}", style=f1_style),
            Text(f"{cfg['tc_cv']:.4f}", style=tc_style),
            format_params_short(cfg["params"]),
        )
    console.print(table)

    # --- Re-evaluate best on FULL dataset ---
    best_f1_idx = np.argmin(res.F[:, 0])
    best = pareto_configs[best_f1_idx]

    console.print(f"\n{'='*60}")
    console.print(f"[bold]Best CV F1 solution:[/] CV F1={best['f1_cv']:.4f} CV TC={best['tc_cv']:.4f}")
    console.print(f"  {format_params_long(best['params'])}")

    console.print(f"\n[bold]Full-dataset evaluation[/] (checking for overfit)")
    best_full = evaluate_params(best["params"], cases, args.n_threads)
    console.print(f"  Full F1=[bold]{best_full['f1']:.4f}[/]  Full TC=[bold]{best_full['tc']:.4f}[/]  "
                  f"Recall={best_full['recall']:.4f}  Precision={best_full['precision']:.4f}")
    for cat, v in sorted(best_full["per_category"].items()):
        console.print(f"    {cat}: F1={v['f1']:.4f} TC={v['tc']:.4f} (n={v['n']})")

    gap_f1 = best_full["f1"] - best["f1_cv"]
    gap_tc = best_full["tc"] - best["tc_cv"]
    console.print(f"\n  Overfit check (full - CV):  F1 {gap_f1:+.4f}  TC {gap_tc:+.4f}")
    if gap_f1 > 0.02:
        console.print("  [bold yellow]Warning:[/] Full-data F1 notably higher than CV — possible overfit")
    else:
        console.print("  [bold green]OK:[/] Full-data and CV scores are consistent")

    # Comparison vs baseline
    console.print(f"\n[bold]Improvement vs baseline[/] (full dataset)")
    comp_table = Table()
    comp_table.add_column("Category")
    comp_table.add_column("F1 before", justify="right")
    comp_table.add_column("F1 after", justify="right")
    comp_table.add_column("ΔF1", justify="right")
    comp_table.add_column("TC before", justify="right")
    comp_table.add_column("TC after", justify="right")
    comp_table.add_column("ΔTC", justify="right")

    for cat in sorted(set(list(baseline_full["per_category"].keys()) +
                         list(best_full["per_category"].keys()))):
        b = baseline_full["per_category"].get(cat, {"f1": 0.0, "tc": 0.0})
        o = best_full["per_category"].get(cat, {"f1": 0.0, "tc": 0.0})
        df1 = o["f1"] - b["f1"]
        dtc = o["tc"] - b["tc"]
        f1_style = "green" if df1 > 0.005 else "red" if df1 < -0.005 else ""
        tc_style = "green" if dtc > 0.005 else "red" if dtc < -0.005 else ""
        comp_table.add_row(
            cat.replace("balibase_", ""),
            f"{b['f1']:.4f}", f"{o['f1']:.4f}",
            Text(f"{df1:+.4f}", style=f1_style),
            f"{b['tc']:.4f}", f"{o['tc']:.4f}",
            Text(f"{dtc:+.4f}", style=tc_style),
        )

    df1 = best_full["f1"] - baseline_full["f1"]
    dtc = best_full["tc"] - baseline_full["tc"]
    f1_style = "bold green" if df1 > 0.005 else "bold red" if df1 < -0.005 else "bold"
    tc_style = "bold green" if dtc > 0.005 else "bold red" if dtc < -0.005 else "bold"
    comp_table.add_row(
        "[bold]OVERALL[/]",
        f"{baseline_full['f1']:.4f}", f"{best_full['f1']:.4f}",
        Text(f"{df1:+.4f}", style=f1_style),
        f"{baseline_full['tc']:.4f}", f"{best_full['tc']:.4f}",
        Text(f"{dtc:+.4f}", style=tc_style),
    )
    console.print(comp_table)

    # --- Save ---
    checkpoint = {
        "pareto_configs": pareto_configs,
        "history": problem.history,
        "baseline_cv": baseline_cv,
        "baseline_full": baseline_full,
        "best_full": best_full,
        "folds_info": [(len(tr), len(te)) for tr, te in folds],
        "args": vars(args),
    }
    ckpt_path = output_dir / "optim_checkpoint.pkl"
    with open(ckpt_path, "wb") as f:
        pickle.dump(checkpoint, f)
    console.print(f"\nCheckpoint saved to {ckpt_path}")

    summary_path = output_dir / "pareto_front.txt"
    with open(summary_path, "w") as f:
        f.write("# Pareto-optimal kalign configs (NSGA-II + {}-fold stratified CV)\n".format(k))
        f.write(f"# pop_size={args.pop_size} n_gen={args.n_gen} "
                f"n_cases={len(cases)} seed={args.seed}\n")
        f.write(f"# Baseline CV: F1={baseline_cv['f1']:.4f} TC={baseline_cv['tc']:.4f}\n\n")
        for i, cfg in enumerate(sorted_pareto):
            p = cfg["params"]
            f.write(f"[{i}] CV_F1={cfg['f1_cv']:.4f} CV_TC={cfg['tc_cv']:.4f}\n")
            f.write(f"    gap_open={p['gap_open']:.3f}\n")
            f.write(f"    gap_extend={p['gap_extend']:.3f}\n")
            f.write(f"    terminal_gap_extend={p['terminal_gap_extend']:.3f}\n")
            f.write(f"    vsm_amax={p['vsm_amax']:.3f}\n")
            f.write(f"    seq_weights={p['seq_weights']:.3f}\n")
            f.write(f"    consistency={p['consistency']}\n")
            f.write(f"    consistency_weight={p['consistency_weight']:.3f}\n")
            f.write(f"    realign={p['realign']}\n")
            f.write(f"    matrix={_matrix_name(p)}\n\n")
    console.print(f"Pareto front saved to {summary_path}")


if __name__ == "__main__":
    main()
