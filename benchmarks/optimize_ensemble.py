#!/usr/bin/env python3
"""Multi-objective ensemble hyperparameter optimization for kalign using pymoo.

Optimizes per-run parameters for kalign's ensemble alignment mode, where N
independent alignments (each with different gap penalties, matrices, and tree
noise) are combined via POAR consensus.

Uses stratified k-fold cross-validation so NSGA-II optimises on held-out
scores, not training scores.

Objectives (maximized):
    1. Mean held-out F1 across folds (category-averaged within each fold)
    2. Mean held-out TC across folds (category-averaged within each fold)

Per-run decision variables (× N runs):
    - gpo:    gap open penalty           [2.0, 15.0]
    - gpe:    gap extend penalty          [0.5, 5.0]
    - tgpe:   terminal gap extend         [0.1, 3.0]
    - noise:  tree perturbation sigma     [0.0, 0.5]
    - matrix: substitution matrix         {PFASUM43, PFASUM60, CorBLOSUM66}

Shared decision variables:
    - vsm_amax:           variable scoring matrix     [0.0, 5.0]
    - consistency:        anchor consistency rounds   {0, 1, 2, 3, 5, 8}
    - consistency_weight: consistency bonus weight    [0.5, 5.0]
    - realign:            tree-rebuild iterations     {0, 1, 2}
    - min_support:        POAR consensus threshold    {0, 1, ..., N}

Usage:
    # Quick test
    uv run python -m benchmarks.optimize_ensemble --n-runs 3 --pop-size 20 --n-gen 5

    # Production (3 runs, Threadripper)
    uv run python -m benchmarks.optimize_ensemble \\
        --n-runs 3 --pop-size 100 --n-gen 50 --n-workers 56

    # Production (5 runs)
    uv run python -m benchmarks.optimize_ensemble \\
        --n-runs 5 --pop-size 150 --n-gen 60 --n-workers 56

    # Production (8 runs)
    uv run python -m benchmarks.optimize_ensemble \\
        --n-runs 8 --pop-size 200 --n-gen 80 --n-workers 56

    # Resume after interrupt
    uv run python -m benchmarks.optimize_ensemble \\
        --resume benchmarks/results/ensemble_optim/gen_checkpoint.pkl --n-gen 80
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

from kalign._core import (  # type: ignore[import-untyped]
    PROTEIN, PROTEIN_PFASUM43, PROTEIN_PFASUM60, REFINE_CONFIDENT,
    ensemble_custom_file_to_file,
)
import kalign  # type: ignore[import-untyped]

from .datasets import BenchmarkCase, balibase_cases, balibase_download, balibase_is_available
from .scoring import score_alignment_detailed

# ---------------------------------------------------------------------------
# Parameter space definition
# ---------------------------------------------------------------------------

# Per-run continuous: [gpo, gpe, tgpe, noise] × N_RUNS
PER_RUN_CONT_LOWER = np.array([2.0, 0.5, 0.1, 0.0])
PER_RUN_CONT_UPPER = np.array([15.0, 5.0, 3.0, 0.5])
PER_RUN_CONT_NAMES = ["gpo", "gpe", "tgpe", "noise"]
N_PER_RUN_CONT = len(PER_RUN_CONT_LOWER)

# Per-run integer: [matrix] × N_RUNS
# matrix: 0=PFASUM43, 1=PFASUM60, 2=CorBLOSUM66
N_PER_RUN_INT = 1  # just matrix

# Shared continuous: [vsm_amax, consistency_weight]
SHARED_CONT_LOWER = np.array([0.0, 0.5])
SHARED_CONT_UPPER = np.array([5.0, 5.0])
SHARED_CONT_NAMES = ["vsm_amax", "consistency_weight"]
N_SHARED_CONT = len(SHARED_CONT_LOWER)

# Shared integer: [consistency, realign, min_support]
# consistency: 0-5 maps via CONSISTENCY_MAP
# realign: 0-2
# min_support: 0-N_RUNS (0=auto)
N_SHARED_INT = 3  # consistency, realign, min_support

CONSISTENCY_MAP = [0, 1, 2, 3, 5, 8]
MATRIX_MAP = [PROTEIN_PFASUM43, PROTEIN_PFASUM60, PROTEIN]
MATRIX_NAMES = {PROTEIN_PFASUM43: "P43", PROTEIN_PFASUM60: "P60", PROTEIN: "CB66"}


def get_var_counts(n_runs: int):
    """Return (n_cont, n_int, n_var) for a given number of ensemble runs."""
    n_cont = n_runs * N_PER_RUN_CONT + N_SHARED_CONT
    n_int = n_runs * N_PER_RUN_INT + N_SHARED_INT
    n_var = n_cont + n_int
    return n_cont, n_int, n_var


def get_bounds(n_runs: int):
    """Return (xl, xu) arrays for the full decision vector."""
    # Continuous: [per_run_cont × N_RUNS, shared_cont]
    cont_lower = np.concatenate([np.tile(PER_RUN_CONT_LOWER, n_runs), SHARED_CONT_LOWER])
    cont_upper = np.concatenate([np.tile(PER_RUN_CONT_UPPER, n_runs), SHARED_CONT_UPPER])

    # Integer: [matrix × N_RUNS, consistency, realign, min_support]
    int_lower = np.zeros(n_runs * N_PER_RUN_INT + N_SHARED_INT)
    int_upper = np.concatenate([
        np.full(n_runs, 2.0),  # matrix: 0-2
        [len(CONSISTENCY_MAP) - 1],  # consistency: 0-5
        [2.0],  # realign: 0-2
        [float(n_runs)],  # min_support: 0-N
    ])

    xl = np.concatenate([cont_lower, int_lower])
    xu = np.concatenate([cont_upper, int_upper])
    return xl, xu


def decode_ensemble_params(x, n_runs: int):
    """Decode a decision vector into ensemble parameter dict.

    Returns dict with:
        run_gpo, run_gpe, run_tgpe, run_noise: lists of float (length n_runs)
        run_types: list of int (length n_runs)
        vsm_amax, consistency_weight: float
        consistency, realign, min_support: int
    """
    n_cont, n_int, _ = get_var_counts(n_runs)
    cont = x[:n_cont]
    ints = np.round(x[n_cont:]).astype(int)

    # Per-run continuous
    run_gpo = []
    run_gpe = []
    run_tgpe = []
    run_noise = []
    for k in range(n_runs):
        offset = k * N_PER_RUN_CONT
        run_gpo.append(float(cont[offset + 0]))
        run_gpe.append(float(cont[offset + 1]))
        run_tgpe.append(float(cont[offset + 2]))
        run_noise.append(float(cont[offset + 3]))

    # Shared continuous
    shared_offset = n_runs * N_PER_RUN_CONT
    vsm_amax = float(cont[shared_offset + 0])
    consistency_weight = float(cont[shared_offset + 1])

    # Per-run integer (matrix)
    run_types = []
    for k in range(n_runs):
        matrix_idx = int(np.clip(ints[k], 0, len(MATRIX_MAP) - 1))
        run_types.append(MATRIX_MAP[matrix_idx])

    # Shared integer
    shared_int_offset = n_runs * N_PER_RUN_INT
    consistency_idx = int(np.clip(ints[shared_int_offset + 0], 0, len(CONSISTENCY_MAP) - 1))
    realign = int(np.clip(ints[shared_int_offset + 1], 0, 2))
    min_support = int(np.clip(ints[shared_int_offset + 2], 0, n_runs))

    return {
        "run_gpo": run_gpo,
        "run_gpe": run_gpe,
        "run_tgpe": run_tgpe,
        "run_noise": run_noise,
        "run_types": run_types,
        "vsm_amax": vsm_amax,
        "consistency_weight": consistency_weight,
        "consistency": CONSISTENCY_MAP[consistency_idx],
        "realign": realign,
        "min_support": min_support,
    }


def format_run_short(params, k):
    """Short string for one run's params."""
    mat = MATRIX_NAMES.get(params["run_types"][k], "?")
    return (f"gpo={params['run_gpo'][k]:.1f} gpe={params['run_gpe'][k]:.2f} "
            f"tgpe={params['run_tgpe'][k]:.2f} n={params['run_noise'][k]:.2f} {mat}")


def format_ensemble_short(params):
    """Compact summary of ensemble params."""
    n_runs = len(params["run_gpo"])
    runs = []
    for k in range(n_runs):
        mat = MATRIX_NAMES.get(params["run_types"][k], "?")
        runs.append(f"R{k}:{params['run_gpo'][k]:.1f}/{mat}")
    run_str = " ".join(runs)
    shared_str = (f"vsm={params['vsm_amax']:.1f} c={params['consistency']} "
                  f"cw={params['consistency_weight']:.1f} re={params['realign']} "
                  f"ms={params['min_support']}")
    return f"{run_str} | {shared_str}"


def format_ensemble_long(params):
    """Verbose multi-line summary."""
    lines = []
    n_runs = len(params["run_gpo"])
    for k in range(n_runs):
        lines.append(f"  Run {k}: {format_run_short(params, k)}")
    lines.append(f"  Shared: vsm={params['vsm_amax']:.2f} "
                 f"cons={params['consistency']} cw={params['consistency_weight']:.2f} "
                 f"re={params['realign']} ms={params['min_support']}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Stratified k-fold CV (reused from optimize_params)
# ---------------------------------------------------------------------------

def stratified_kfold(cases: List[BenchmarkCase], k: int, seed: int = 42
                     ) -> List[Tuple[List[BenchmarkCase], List[BenchmarkCase]]]:
    """Split cases into k stratified folds by dataset (RV category)."""
    rng = np.random.RandomState(seed)

    by_cat: Dict[str, List[BenchmarkCase]] = defaultdict(list)
    for c in cases:
        by_cat[c.dataset].append(c)

    fold_assignments: Dict[str, int] = {}
    for cat_cases in by_cat.values():
        indices = list(range(len(cat_cases)))
        rng.shuffle(indices)
        for rank, idx in enumerate(indices):
            fold_assignments[cat_cases[idx].family] = rank % k

    folds = []
    for fold_idx in range(k):
        test = [c for c in cases if fold_assignments[c.family] == fold_idx]
        train = [c for c in cases if fold_assignments[c.family] != fold_idx]
        folds.append((train, test))

    return folds


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate_ensemble(params, cases, n_threads=1, quiet=True):
    """Run ensemble alignment with given params on all cases, return mean metrics."""
    results_by_cat: Dict[str, list] = {}
    total_time = 0.0

    for case in cases:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / f"{case.family}_aln.fasta"

            try:
                start = time.perf_counter()
                ensemble_custom_file_to_file(
                    str(case.unaligned),
                    str(output),
                    run_gpo=params["run_gpo"],
                    run_gpe=params["run_gpe"],
                    run_tgpe=params["run_tgpe"],
                    run_noise=params["run_noise"],
                    run_types=params["run_types"],
                    format="fasta",
                    seq_type=PROTEIN,
                    seed=42,
                    min_support=params["min_support"],
                    refine=REFINE_CONFIDENT,
                    vsm_amax=params["vsm_amax"],
                    realign=params["realign"],
                    seq_weights=-1.0,
                    n_threads=n_threads,
                    consistency_anchors=params["consistency"],
                    consistency_weight=params["consistency_weight"],
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

    per_cat = {}
    for cat, scores in results_by_cat.items():
        per_cat[cat] = {
            "f1": np.mean([s["f1"] for s in scores]),
            "tc": np.mean([s["tc"] for s in scores]),
            "recall": np.mean([s["recall"] for s in scores]),
            "precision": np.mean([s["precision"] for s in scores]),
            "n": len(scores),
        }

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


def evaluate_ensemble_cv(params, folds, n_threads=1, quiet=True):
    """Evaluate ensemble params using stratified k-fold CV."""
    fold_f1s = []
    fold_tcs = []
    total_time = 0.0

    for _, test in folds:
        result = evaluate_ensemble(params, test, n_threads, quiet)
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
    """Rich-based live terminal dashboard for ensemble optimization progress."""

    def __init__(self, n_gen: int, pop_size: int, baseline_f1: float, baseline_tc: float):
        self.n_gen = n_gen
        self.pop_size = pop_size
        self.baseline_f1 = baseline_f1
        self.baseline_tc = baseline_tc
        self.console = Console()

        self.current_gen = 0
        self.eval_count = 0
        self.total_evals = n_gen * pop_size
        self.gen_start_time = time.time()
        self.run_start_time = time.time()
        self.current_eval_params: Optional[dict] = None
        self.current_eval_idx = 0

        self.best_f1 = 0.0
        self.best_f1_params: Optional[dict] = None
        self.best_tc = 0.0
        self.best_tc_params: Optional[dict] = None

        self.pareto_front: List[dict] = []
        self.gen_history: List[dict] = []
        self.recent_evals: List[dict] = []

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
            lines.append(f"[dim]Current: {format_ensemble_short(self.current_eval_params)}[/]")

        return Panel("\n".join(lines), title="Ensemble Optimization", border_style="blue")

    def _build_best_panel(self) -> Panel:
        lines = []
        lines.append(f"[dim]Baseline:  F1={self.baseline_f1:.4f}  TC={self.baseline_tc:.4f}[/]")
        lines.append("")

        if self.best_f1_params:
            lines.append(f"Best F1:   {self._delta_str(self.best_f1, self.baseline_f1)}")
            lines.append(f"  [dim]{format_ensemble_short(self.best_f1_params)}[/]")
        else:
            lines.append("[dim]Best F1:   (no evaluations yet)[/]")

        lines.append("")

        if self.best_tc_params:
            lines.append(f"Best TC:   {self._delta_str(self.best_tc, self.baseline_tc)}")
            lines.append(f"  [dim]{format_ensemble_short(self.best_tc_params)}[/]")
        else:
            lines.append("[dim]Best TC:   (no evaluations yet)[/]")

        return Panel("\n".join(lines), title="Best Solutions", border_style="green")

    def _build_pareto_table(self) -> Panel:
        table = Table(show_header=True, expand=True, padding=(0, 1))
        table.add_column("#", style="dim", width=3)
        table.add_column("CV F1", justify="right", width=8)
        table.add_column("CV TC", justify="right", width=8)
        table.add_column("Parameters")

        for i, cfg in enumerate(self.pareto_front[:8]):
            f1_style = "bold green" if cfg["f1"] > self.baseline_f1 else ""
            tc_style = "bold green" if cfg["tc"] > self.baseline_tc else ""
            table.add_row(
                str(i),
                Text(f"{cfg['f1']:.4f}", style=f1_style),
                Text(f"{cfg['tc']:.4f}", style=tc_style),
                format_ensemble_short(cfg["params"]),
            )

        return Panel(table, title="Pareto Front", border_style="yellow")

    def _build_trend_panel(self) -> Panel:
        if not self.gen_history:
            return Panel("[dim]No data yet[/]", title="Trend", border_style="cyan")

        entries = []
        step = max(1, len(self.gen_history) // 8)
        for i in range(0, len(self.gen_history), step):
            h = self.gen_history[i]
            entries.append(f"Gen {h['gen']:>3}: F1={h['best_f1']:.4f}")
        if self.gen_history[-1] not in [self.gen_history[i] for i in range(0, len(self.gen_history), step)]:
            h = self.gen_history[-1]
            entries.append(f"Gen {h['gen']:>3}: F1={h['best_f1']:.4f}")

        return Panel("  ".join(entries), title="Trend", border_style="cyan")

    def _build_recent_panel(self) -> Panel:
        if not self.recent_evals:
            return Panel("[dim]No evaluations yet[/]", title="Recent", border_style="magenta")

        lines = []
        for ev in self.recent_evals[-5:]:
            lines.append(f"F1={ev['f1']:.4f} TC={ev['tc']:.4f} {format_ensemble_short(ev['params'])}")

        return Panel("\n".join(lines), title="Recent", border_style="magenta")

    def _build_layout(self) -> Layout:
        layout = Layout()
        layout.split_column(
            Layout(self._build_status_panel(), name="status", size=4),
            Layout(self.progress, name="progress", size=3),
            Layout(self._build_best_panel(), name="best", size=9),
            Layout(self._build_pareto_table(), name="pareto", size=12),
            Layout(self._build_trend_panel(), name="trend", size=3),
            Layout(self._build_recent_panel(), name="recent", size=8),
        )
        return layout

    def _refresh(self):
        self.live.update(self._build_layout())

    def on_eval_start(self, params, eval_count, idx_in_gen):
        self.current_eval_params = params
        self.eval_count = eval_count
        self.current_eval_idx = idx_in_gen
        self.progress.update(self.eval_task, completed=idx_in_gen)
        self._refresh()

    def on_eval_end(self, params, cv_result):
        f1 = cv_result["f1"]
        tc = cv_result["tc"]

        if f1 > self.best_f1:
            self.best_f1 = f1
            self.best_f1_params = params
        if tc > self.best_tc:
            self.best_tc = tc
            self.best_tc_params = params

        self.recent_evals.append({"params": params, "f1": f1, "tc": tc})
        if len(self.recent_evals) > 10:
            self.recent_evals = self.recent_evals[-10:]

        self._refresh()

    def on_gen_start(self, gen):
        self.current_gen = gen
        self.gen_start_time = time.time()
        self.progress.update(self.gen_task, completed=gen - 1)
        self.progress.update(self.eval_task, completed=0)
        self._refresh()

    def on_gen_end(self, gen, pareto):
        self.current_gen = gen
        self.progress.update(self.gen_task, completed=gen)
        self.progress.update(self.eval_task, completed=self.pop_size)

        self.pareto_front = sorted(pareto, key=lambda x: -x["f1"])

        best_f1_in_gen = max(p["f1"] for p in pareto) if pareto else 0.0
        self.gen_history.append({"gen": gen, "best_f1": best_f1_in_gen})

        self._refresh()


# ---------------------------------------------------------------------------
# ProcessPoolExecutor helper
# ---------------------------------------------------------------------------

def _kill_pool(pool: ProcessPoolExecutor) -> None:
    """Forcefully terminate all worker processes in the pool."""
    for pid in list(pool._processes):  # noqa: SLF001
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    pool.shutdown(wait=False, cancel_futures=True)


# Top-level function for pickling (ProcessPoolExecutor requires this)
def _eval_one_ensemble(args_tuple):
    """Evaluate one ensemble parameter vector. Top-level for ProcessPoolExecutor."""
    x, n_runs, folds, n_threads = args_tuple
    params = decode_ensemble_params(x, n_runs)
    cv_result = evaluate_ensemble_cv(params, folds, n_threads, quiet=True)
    return params, cv_result


def _eval_one_ensemble_fold(args_tuple):
    """Evaluate one (individual, fold) pair for ensemble. Fine-grained parallelism."""
    ind_idx, fold_idx, x, n_runs, test_cases, n_threads = args_tuple
    params = decode_ensemble_params(x, n_runs)
    result = evaluate_ensemble(params, test_cases, n_threads, quiet=True)
    return ind_idx, fold_idx, params, result


# ---------------------------------------------------------------------------
# pymoo Problem + Callback
# ---------------------------------------------------------------------------

class EnsembleCVProblem(Problem):
    """Multi-objective optimization for ensemble with stratified CV."""

    def __init__(self, n_runs, folds, n_threads=1, n_workers=1,
                 optimize_time=False, dashboard: Optional[Dashboard] = None):
        xl, xu = get_bounds(n_runs)
        n_cont, n_int, n_var = get_var_counts(n_runs)
        n_obj = 3 if optimize_time else 2

        super().__init__(
            n_var=n_var,
            n_obj=n_obj,
            xl=xl,
            xu=xu,
        )
        self.n_runs = n_runs
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
            params = decode_ensemble_params(x, self.n_runs)
            self.eval_count += 1

            if self.dashboard:
                self.dashboard.on_eval_start(params, self.eval_count, i)

            cv_result = evaluate_ensemble_cv(params, self.folds, self.n_threads, quiet=True)
            self._record(i, F, params, cv_result)

    def _evaluate_parallel(self, X, F):
        """Fine-grained parallelism: submit (individual × fold) jobs."""
        n_folds = len(self.folds)

        jobs = []
        for i, x in enumerate(X):
            for fold_idx, (_, test) in enumerate(self.folds):
                jobs.append((i, fold_idx, x, self.n_runs, test, self.n_threads))

        if self.dashboard:
            self.dashboard.on_eval_start(
                decode_ensemble_params(X[0], self.n_runs),
                self.eval_count + 1, 0)

        fold_results: Dict[int, Dict[int, dict]] = defaultdict(dict)
        ind_params: Dict[int, dict] = {}

        pool = ProcessPoolExecutor(max_workers=self.n_workers)
        try:
            futures = {pool.submit(_eval_one_ensemble_fold, j): j[:2] for j in jobs}
            completed_individuals = set()

            for future in as_completed(futures):
                ind_idx, fold_idx, params, result = future.result()
                fold_results[ind_idx][fold_idx] = result
                ind_params[ind_idx] = params

                if len(fold_results[ind_idx]) == n_folds:
                    completed_individuals.add(ind_idx)
                    self.eval_count += 1

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
                 problem: Optional[EnsembleCVProblem] = None):
        super().__init__()
        self.dashboard = dashboard
        self.checkpoint_path = checkpoint_path
        self.problem = problem

    def notify(self, algorithm):
        gen = algorithm.n_gen

        if self.dashboard:
            self.dashboard.on_gen_start(gen)

        pareto = []
        if algorithm.opt is not None and len(algorithm.opt) > 0:
            n_runs = self.problem.n_runs if self.problem else 3
            for ind in algorithm.opt:
                params = decode_ensemble_params(ind.X, n_runs)
                pareto.append({
                    "params": params,
                    "f1": -ind.F[0],
                    "tc": -ind.F[1],
                })

        if self.dashboard:
            self.dashboard.on_gen_end(gen, pareto)

        if self.checkpoint_path:
            pop = algorithm.pop
            ckpt = {
                "pop_X": pop.get("X").copy(),
                "pop_F": pop.get("F").copy(),
                "n_gen_completed": gen,
                "history": self.problem.history if self.problem else [],
                "pop_size": len(pop),
                "n_runs": self.problem.n_runs if self.problem else 3,
            }
            tmp = self.checkpoint_path.with_suffix(".tmp")
            with open(tmp, "wb") as f:
                pickle.dump(ckpt, f)
            tmp.rename(self.checkpoint_path)


def load_checkpoint(path: Path):
    """Load a generation checkpoint."""
    with open(path, "rb") as f:
        ckpt = pickle.load(f)  # noqa: S301
    return (ckpt["pop_X"], ckpt["pop_F"], ckpt["n_gen_completed"],
            ckpt.get("history", []), ckpt.get("n_runs", 3))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Optimize kalign ensemble hyperparameters with NSGA-II + stratified k-fold CV")
    parser.add_argument("--n-runs", type=int, default=3,
                        help="Number of ensemble runs (default: 3)")
    parser.add_argument("--pop-size", type=int, default=100,
                        help="Population size (default: 100)")
    parser.add_argument("--n-gen", type=int, default=50,
                        help="Number of generations (default: 50)")
    parser.add_argument("--n-folds", type=int, default=5,
                        help="Number of CV folds (default: 5)")
    parser.add_argument("--n-threads", type=int, default=1,
                        help="OpenMP threads per alignment (default: 1)")
    parser.add_argument("--n-workers", type=int, default=1,
                        help="Parallel worker processes (default: 1)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--optimize-time", action="store_true",
                        help="Add wall time as 3rd objective")
    parser.add_argument("--output-dir", type=str, default="benchmarks/results/ensemble_optim",
                        help="Output directory")
    parser.add_argument("--run-name", type=str, default=None,
                        help="Name for this run (creates subdirectory, e.g. 'run1_3runs')")
    parser.add_argument("--no-dashboard", action="store_true",
                        help="Disable rich dashboard")
    parser.add_argument("--resume", type=str, default=None,
                        help="Resume from a generation checkpoint file (.pkl)")
    args = parser.parse_args()

    console = Console()
    n_runs = args.n_runs

    # Dimensions
    n_cont, n_int, n_var = get_var_counts(n_runs)
    console.print(f"[bold]Ensemble optimizer[/]: {n_runs} runs, "
                  f"{n_var} decision variables ({n_cont} continuous + {n_int} integer)")

    # Setup
    if not balibase_is_available():
        console.print("Downloading BAliBASE...")
        balibase_download()

    cases = balibase_cases()
    console.print(f"Loaded [bold]{len(cases)}[/] BAliBASE cases")

    cats: Dict[str, int] = {}
    for c in cases:
        cats[c.dataset] = cats.get(c.dataset, 0) + 1
    for cat, n in sorted(cats.items()):
        console.print(f"  {cat}: {n} cases")

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

    # --- Baseline evaluation (current best: ens3+vsm+ref+ra1) ---
    console.print(f"\n[bold]Baseline evaluation[/] (ens{n_runs}+vsm+ref+ra1, {k}-fold CV)")
    console.print("[dim]This uses the original kalign_ensemble with hardcoded scale factors[/]")

    # Evaluate baseline via the standard ensemble API
    baseline_fold_f1s = []
    baseline_fold_tcs = []
    baseline_total_time = 0.0
    for fold_idx, (_, test) in enumerate(folds):
        results_by_cat: Dict[str, list] = {}
        for case in test:
            with tempfile.TemporaryDirectory() as tmpdir:
                output = Path(tmpdir) / f"{case.family}_aln.fasta"
                try:
                    start = time.perf_counter()
                    kalign.align_file_to_file(
                        str(case.unaligned), str(output),
                        format="fasta",
                        ensemble=n_runs,
                        vsm_amax=2.0,
                        refine="confident",
                        realign=1,
                        ensemble_seed=42,
                    )
                    baseline_total_time += time.perf_counter() - start
                    detailed = score_alignment_detailed(case.reference, output)
                    cat = case.dataset
                    if cat not in results_by_cat:
                        results_by_cat[cat] = []
                    results_by_cat[cat].append(detailed)
                except Exception as e:
                    console.print(f"  [yellow]WARN[/]: {case.family}: {e}")

        if results_by_cat:
            cat_f1s = [np.mean([s["f1"] for s in v]) for v in results_by_cat.values()]
            cat_tcs = [np.mean([s["tc"] for s in v]) for v in results_by_cat.values()]
            baseline_fold_f1s.append(float(np.mean(cat_f1s)))
            baseline_fold_tcs.append(float(np.mean(cat_tcs)))

    baseline_cv = {
        "f1": float(np.mean(baseline_fold_f1s)),
        "tc": float(np.mean(baseline_fold_tcs)),
        "f1_std": float(np.std(baseline_fold_f1s)),
        "tc_std": float(np.std(baseline_fold_tcs)),
        "fold_f1s": baseline_fold_f1s,
        "fold_tcs": baseline_fold_tcs,
        "wall_time": baseline_total_time,
    }
    console.print(f"  CV F1=[bold]{baseline_cv['f1']:.4f}[/]±{baseline_cv['f1_std']:.4f}  "
                  f"CV TC=[bold]{baseline_cv['tc']:.4f}[/]±{baseline_cv['tc_std']:.4f}  "
                  f"time={baseline_cv['wall_time']:.1f}s")
    for i, (f1, tc) in enumerate(zip(baseline_cv['fold_f1s'], baseline_cv['fold_tcs'])):
        console.print(f"    Fold {i}: F1={f1:.4f} TC={tc:.4f}")

    # --- Optimization ---
    n_evals = args.pop_size * args.n_gen
    est_sec_per_eval = baseline_cv["wall_time"] / k  # baseline did k folds
    parallelism = max(1, args.n_workers)
    est_hours = n_evals * est_sec_per_eval / parallelism / 3600

    console.print(f"\n[bold]Starting NSGA-II[/]: pop_size={args.pop_size}, n_gen={args.n_gen}, "
                  f"{k}-fold CV, {n_runs} ensemble runs, "
                  f"{args.n_workers} worker(s) × {args.n_threads} thread(s)")
    console.print(f"Total evaluations: ~{n_evals}")
    console.print(f"Estimated time: ~{est_hours:.1f} hours\n")

    use_dashboard = not args.no_dashboard
    dashboard = None

    if use_dashboard:
        dashboard = Dashboard(
            n_gen=args.n_gen,
            pop_size=args.pop_size,
            baseline_f1=baseline_cv["f1"],
            baseline_tc=baseline_cv["tc"],
        )

    problem = EnsembleCVProblem(
        n_runs=n_runs,
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
        pop_X, _pop_F, resumed_gen, prev_history, ckpt_n_runs = load_checkpoint(resume_path)
        if ckpt_n_runs != n_runs:
            console.print(f"[bold red]Checkpoint n_runs={ckpt_n_runs} != --n-runs={n_runs}[/]")
            return
        problem.history = prev_history
        console.print(f"[bold green]Resumed[/] from generation {resumed_gen} "
                      f"({len(prev_history)} prior evaluations)")
        remaining = args.n_gen - resumed_gen
        if remaining <= 0:
            console.print(f"[bold yellow]Already completed {resumed_gen} generations "
                          f"(requested {args.n_gen}). Increase --n-gen to continue.[/]")
            return
        termination = get_termination("n_gen", remaining)
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
            verbose=not use_dashboard,
            callback=callback,
        )
    except KeyboardInterrupt:
        if dashboard:
            dashboard.stop()
        console.print("\n[bold yellow]Interrupted![/] Checkpoint was saved after last "
                      f"completed generation to:\n  {checkpoint_path}")
        console.print("Resume with: [bold]--resume " + str(checkpoint_path) + "[/]")
        os._exit(1)  # skip atexit handlers that hang on worker join
    finally:
        if dashboard:
            dashboard.stop()

    # --- Results ---
    console.print(f"\n[bold]Optimization complete.[/] "
                  f"{len(res.F)} Pareto-optimal solutions found.\n")

    pareto_configs = []
    for i, (x, f) in enumerate(zip(res.X, res.F)):
        params = decode_ensemble_params(x, n_runs)
        f1 = -f[0]
        tc = -f[1]
        wt = f[2] if args.optimize_time else None
        pareto_configs.append({"params": params, "f1_cv": f1, "tc_cv": tc, "wall_time": wt})

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
            format_ensemble_short(cfg["params"]),
        )
    console.print(table)

    # --- Re-evaluate best on FULL dataset ---
    best_f1_idx = np.argmin(res.F[:, 0])
    best = pareto_configs[best_f1_idx]

    console.print(f"\n{'='*60}")
    console.print(f"[bold]Best CV F1 solution:[/] CV F1={best['f1_cv']:.4f} CV TC={best['tc_cv']:.4f}")
    console.print(format_ensemble_long(best["params"]))

    console.print(f"\n[bold]Full-dataset evaluation[/] (checking for overfit)")
    best_full = evaluate_ensemble(best["params"], cases, args.n_threads)
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

    # --- Save ---
    checkpoint = {
        "pareto_configs": pareto_configs,
        "history": problem.history,
        "baseline_cv": baseline_cv,
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
        f.write(f"# Pareto-optimal kalign ensemble configs (NSGA-II + {k}-fold stratified CV)\n")
        f.write(f"# n_runs={n_runs} pop_size={args.pop_size} n_gen={args.n_gen} "
                f"n_cases={len(cases)} seed={args.seed}\n")
        f.write(f"# Baseline CV: F1={baseline_cv['f1']:.4f} TC={baseline_cv['tc']:.4f}\n\n")
        for i, cfg in enumerate(sorted_pareto):
            p = cfg["params"]
            f.write(f"[{i}] CV_F1={cfg['f1_cv']:.4f} CV_TC={cfg['tc_cv']:.4f}\n")
            for rk in range(n_runs):
                mat = MATRIX_NAMES.get(p["run_types"][rk], "?")
                f.write(f"    Run {rk}: gpo={p['run_gpo'][rk]:.3f} gpe={p['run_gpe'][rk]:.3f} "
                        f"tgpe={p['run_tgpe'][rk]:.3f} noise={p['run_noise'][rk]:.3f} "
                        f"matrix={mat}\n")
            f.write(f"    Shared: vsm_amax={p['vsm_amax']:.3f} "
                    f"consistency={p['consistency']} "
                    f"consistency_weight={p['consistency_weight']:.3f} "
                    f"realign={p['realign']} min_support={p['min_support']}\n\n")
    console.print(f"Pareto front saved to {summary_path}")


if __name__ == "__main__":
    main()
