#!/usr/bin/env python3
"""Unified multi-objective hyperparameter optimization for kalign using pymoo.

Searches across kalign's entire operating range: from fast single-run alignment
through consistency-enhanced single-run to multi-run ensemble with POAR consensus.
Always uses three objectives: recall, precision, and wall time. The resulting 3D
Pareto surface reveals the full speed/accuracy trade-off landscape in one run.
Recall and precision are naturally in tension, yielding a richer Pareto front
than F1+TC.

Objectives (always 3):
    1. Maximize recall (category-averaged SP score, held-out CV)
    2. Maximize precision (category-averaged, held-out CV)
    3. Minimize wall time (total CV evaluation time in seconds)

Per-run decision variables (x max_runs slots):
    - gpo:       gap open penalty           [2.0, 15.0]
    - gpe:       gap extend penalty          [0.5, 5.0]
    - tgpe:      terminal gap extend         [0.1, 3.0]
    - noise:     tree perturbation sigma     [0.0, 0.5]
    - matrix:    substitution matrix         {PFASUM43, PFASUM60, CorBLOSUM66}
    - vsm_amax:  variable scoring matrix     [0.0, 5.0]
    - refine:    post-alignment refinement   {NONE, ALL, CONFIDENT, INLINE}

Shared decision variables:
    - n_runs:             {1, 3, 5}       Core mode variable
    - seq_weights:        [0.0, 5.0]     Profile rebalancing
    - consistency:        {0..6}          Anchor consistency rounds
    - consistency_weight: [0.5, 5.0]     Consistency bonus weight
    - realign:            {0, 1, 2}       Tree-rebuild iterations
    - min_support:        {0..max_runs}   POAR consensus threshold

Usage:
    # Quick smoke test (protein/BAliBASE, default)
    uv run python -m benchmarks.optimize_unified --pop-size 20 --n-gen 5

    # Production run on BAliBASE (protein)
    uv run python -m benchmarks.optimize_unified \\
        --pop-size 200 --n-gen 100 --n-workers 56 --n-threads 1

    # Production run on BRAliBASE (RNA)
    uv run python -m benchmarks.optimize_unified --dataset bralibase \\
        --pop-size 200 --n-gen 100 --n-workers 56 --n-threads 1

    # Resume
    uv run python -m benchmarks.optimize_unified \\
        --resume benchmarks/results/unified_optim/balibase/gen_checkpoint.pkl \\
        --n-gen 150 --n-workers 56
"""

import os
# Must set OMP_NUM_THREADS before importing kalign (or any C extension that
# uses OpenMP).  ProcessPoolExecutor uses fork(), and forked children inherit
# the parent's OpenMP thread-pool state — but the pool threads are dead in the
# child.  If OpenMP later tries to use them, the child segfaults.  Setting
# this to "1" prevents the parent from ever creating extra threads, making
# fork safe.  The actual per-alignment thread count is controlled by the
# n_threads parameter passed to kalign at call time.
if "OMP_NUM_THREADS" not in os.environ:
    os.environ["OMP_NUM_THREADS"] = "1"

import argparse
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
    from pymoo.algorithms.moo.nsga3 import NSGA3  # type: ignore[import-untyped]
    from pymoo.core.callback import Callback  # type: ignore[import-untyped]
    from pymoo.core.mixed import (  # type: ignore[import-untyped]
        MixedVariableDuplicateElimination,
        MixedVariableMating,
        MixedVariableSampling,
    )
    from pymoo.core.population import Population  # type: ignore[import-untyped]
    from pymoo.core.problem import Problem  # type: ignore[import-untyped]
    from pymoo.core.variable import Choice, Integer, Real  # type: ignore[import-untyped]
    from pymoo.optimize import minimize  # type: ignore[import-untyped]
    from pymoo.termination import get_termination  # type: ignore[import-untyped]
    from pymoo.util.ref_dirs import get_reference_directions  # type: ignore[import-untyped]
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
    MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66,
    MATRIX_DNA, MATRIX_RNA,
    MATRIX_NUC_1PAM, MATRIX_NUC_20PAM, MATRIX_NUC_200PAM,
    REFINE_NONE, REFINE_ALL,
    REFINE_CONFIDENT, REFINE_INLINE, ensemble_custom_file_to_file,
)

from .datasets import (BenchmarkCase,
                        balibase_cases, balibase_download, balibase_is_available,
                        bralibase_cases, bralibase_download, bralibase_is_available,
                        mdsa_cases, mdsa_download, mdsa_is_available)
from .scoring import score_alignment_detailed

# ---------------------------------------------------------------------------
# Parameter space definition
# ---------------------------------------------------------------------------

# Maps for categorical variables (Choice options)
# These define the actual values; pymoo Choice handles selection directly.

N_RUNS_MAP = [1, 3, 5]  # index -> actual n_runs
CONSISTENCY_MAP = [0, 1, 2, 3, 5, 8, 10]
REFINE_MAP = [REFINE_NONE, REFINE_ALL, REFINE_CONFIDENT, REFINE_INLINE]
REFINE_NAMES = {REFINE_NONE: "N", REFINE_ALL: "A", REFINE_CONFIDENT: "C", REFINE_INLINE: "I"}
REFINE_LONG = {REFINE_NONE: "NONE", REFINE_ALL: "ALL", REFINE_CONFIDENT: "CONFIDENT",
               REFINE_INLINE: "INLINE"}

# --- Dataset-specific parameter space profiles ---

PARAM_PROFILES = {
    "protein": {
        "per_run_cont_lower": np.array([2.0, 0.5, 0.1, 0.0]),
        "per_run_cont_upper": np.array([15.0, 5.0, 3.0, 0.5]),
        "shared_cont_lower": np.array([0.0, 0.0, 0.5]),
        "shared_cont_upper": np.array([5.0, 5.0, 5.0]),
        "matrix_map_int": [MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66],
        "matrix_map_str": ["pfasum43", "pfasum60", "corblosum66"],
        "matrix_names": {MATRIX_PFASUM43: "P43", MATRIX_PFASUM60: "P60", MATRIX_CORBLOSUM66: "CB66"},
        "n_matrices": 3,
        "seq_type_int": MATRIX_PFASUM43,  # default protein matrix for ensemble seq_type
        "seq_type_str": "protein",  # Python API string
        "max_consistency_idx": len(CONSISTENCY_MAP) - 1,  # full range
    },
    "rna": {
        # RNA gap penalties tend to be different — wider search range for gpo
        "per_run_cont_lower": np.array([1.0, 0.2, 0.05, 0.0]),
        "per_run_cont_upper": np.array([20.0, 5.0, 3.0, 0.5]),
        # vsm_amax: 0 by default for RNA, but let optimizer explore [0, 3]
        # seq_weights: 0 by default for RNA, search [0, 3]
        "shared_cont_lower": np.array([0.0, 0.0, 0.5]),
        "shared_cont_upper": np.array([3.0, 3.0, 5.0]),
        # Kimura 1/20/200 PAM matrices (shared with DNA)
        "matrix_map_int": [MATRIX_NUC_1PAM, MATRIX_NUC_20PAM, MATRIX_NUC_200PAM],
        "matrix_map_str": ["nuc_1pam", "nuc_20pam", "nuc_200pam"],
        "matrix_names": {MATRIX_NUC_1PAM: "K1", MATRIX_NUC_20PAM: "K20", MATRIX_NUC_200PAM: "K200"},
        "n_matrices": 3,
        "seq_type_int": MATRIX_NUC_200PAM,
        "seq_type_str": "rna",
        "max_consistency_idx": len(CONSISTENCY_MAP) - 1,  # full range
    },
    "dna": {
        "per_run_cont_lower": np.array([1.0, 0.2, 0.05, 0.0]),
        "per_run_cont_upper": np.array([20.0, 5.0, 3.0, 0.5]),
        "shared_cont_lower": np.array([0.0, 0.0, 0.5]),
        "shared_cont_upper": np.array([3.0, 3.0, 5.0]),
        # Kimura 1/20/200 PAM matrices (shared with RNA)
        "matrix_map_int": [MATRIX_NUC_1PAM, MATRIX_NUC_20PAM, MATRIX_NUC_200PAM],
        "matrix_map_str": ["nuc_1pam", "nuc_20pam", "nuc_200pam"],
        "matrix_names": {MATRIX_NUC_1PAM: "K1", MATRIX_NUC_20PAM: "K20", MATRIX_NUC_200PAM: "K200"},
        "n_matrices": 3,
        "seq_type_int": MATRIX_NUC_200PAM,
        "seq_type_str": "dna",
        "max_consistency_idx": len(CONSISTENCY_MAP) - 1,  # full range
    },
    "nucleotide": {
        "per_run_cont_lower": np.array([1.0, 0.2, 0.05, 0.0]),
        "per_run_cont_upper": np.array([20.0, 5.0, 3.0, 0.5]),
        "shared_cont_lower": np.array([0.0, 0.0, 0.5]),
        "shared_cont_upper": np.array([3.0, 3.0, 5.0]),
        "matrix_map_int": [MATRIX_NUC_1PAM, MATRIX_NUC_20PAM, MATRIX_NUC_200PAM],
        "matrix_map_str": ["nuc_1pam", "nuc_20pam", "nuc_200pam"],
        "matrix_names": {MATRIX_NUC_1PAM: "K1", MATRIX_NUC_20PAM: "K20", MATRIX_NUC_200PAM: "K200"},
        "n_matrices": 3,
        "seq_type_int": MATRIX_NUC_200PAM,
        "seq_type_str": "nucleotide",
        "max_consistency_idx": len(CONSISTENCY_MAP) - 1,
    },
}

# Active profile (set by main() based on --dataset)
_active_profile = PARAM_PROFILES["protein"]

def set_active_profile(profile_name: str):
    global _active_profile
    _active_profile = PARAM_PROFILES[profile_name]

# Convenience accessors for the active profile
def _matrix_map_int(): return _active_profile["matrix_map_int"]
def _matrix_map_str(): return _active_profile["matrix_map_str"]
def _matrix_names(): return _active_profile["matrix_names"]

# Legacy aliases (used in view_pareto.py etc.)
MATRIX_MAP_INT = PARAM_PROFILES["protein"]["matrix_map_int"]
MATRIX_MAP_STR = PARAM_PROFILES["protein"]["matrix_map_str"]
# Combined matrix names across all profiles (used by dashboard)
MATRIX_NAMES = {}
for _prof in PARAM_PROFILES.values():
    MATRIX_NAMES.update(_prof["matrix_names"])

# Old constant names used in view_pareto.py and checkpoint data — keep for loading old checkpoints
_OLD_MATRIX_COMPAT = {3: MATRIX_PFASUM43, 5: MATRIX_PFASUM43, 6: MATRIX_PFASUM60}


def get_vars(max_runs: int) -> dict:
    """Return pymoo mixed-variable space definition.

    Per-run variables (one slot per max_runs):
        gpo, gpe, tgpe, noise, matrix  — gap/tree params
        vsm_amax                       — variable scoring matrix amplitude
        refine                         — post-alignment refinement mode

    Shared variables:
        seq_weights, consistency_weight — scalars applied to all runs
        n_runs, consistency, realign, min_support — integer/categorical
    """
    profile = _active_profile
    lo = profile["per_run_cont_lower"]
    hi = profile["per_run_cont_upper"]
    slo = profile["shared_cont_lower"]
    shi = profile["shared_cont_upper"]

    variables = {}
    for k in range(max_runs):
        variables[f"gpo_{k}"] = Real(bounds=(float(lo[0]), float(hi[0])))
        variables[f"gpe_{k}"] = Real(bounds=(float(lo[1]), float(hi[1])))
        variables[f"tgpe_{k}"] = Real(bounds=(float(lo[2]), float(hi[2])))
        variables[f"noise_{k}"] = Real(bounds=(float(lo[3]), float(hi[3])))
        variables[f"matrix_{k}"] = Choice(options=list(range(profile["n_matrices"])))
        # Per-run VSM amplitude and refinement mode
        variables[f"vsm_amax_{k}"] = Real(bounds=(float(slo[0]), float(shi[0])))
        variables[f"refine_{k}"] = Choice(options=REFINE_MAP)

    variables["seq_weights"] = Real(bounds=(float(slo[1]), float(shi[1])))
    variables["consistency_weight"] = Real(bounds=(float(slo[2]), float(shi[2])))

    consistency_options = CONSISTENCY_MAP[:profile.get("max_consistency_idx", len(CONSISTENCY_MAP) - 1) + 1]
    variables["n_runs"] = Choice(options=N_RUNS_MAP)
    variables["consistency"] = Choice(options=consistency_options)
    variables["realign"] = Integer(bounds=(0, 2))
    variables["min_support"] = Integer(bounds=(0, max_runs))
    variables["consistency_merge"] = Choice(options=[0, 1])
    variables["consistency_merge_weight"] = Real(bounds=(0.5, 10.0))

    return variables


def decode_unified_params(x, max_runs: int):
    """Decode a mixed-variable dict into a unified parameter dict.

    x is a dict with keys like 'gpo_0', 'n_runs', 'consistency', etc.
    Values are native types (float for Real, int for Integer/Choice).

    vsm_amax and refine are per-run (vsm_amax_{k}, refine_{k}).
    If the per-run keys are missing (old checkpoint), falls back to
    shared 'vsm_amax' / 'refine' keys.
    """
    n_runs = int(x["n_runs"])
    consistency = int(x["consistency"])
    realign = int(x["realign"])
    min_support_raw = int(x["min_support"])
    consistency_merge = int(x.get("consistency_merge", 0))
    consistency_merge_weight = float(x.get("consistency_merge_weight", 2.0))

    seq_weights = float(x["seq_weights"])
    consistency_weight = float(x["consistency_weight"])

    run_gpo, run_gpe, run_tgpe, run_noise = [], [], [], []
    run_types, run_matrices = [], []
    run_vsm_amax, run_refine = [], []

    # Detect old checkpoint format (shared vsm_amax / refine)
    has_per_run_vsm = f"vsm_amax_0" in x
    has_per_run_refine = f"refine_0" in x
    shared_vsm = float(x.get("vsm_amax", 0.0)) if not has_per_run_vsm else 0.0
    shared_refine = int(x.get("refine", REFINE_NONE)) if not has_per_run_refine else REFINE_NONE

    for k in range(n_runs):
        run_gpo.append(float(x[f"gpo_{k}"]))
        run_gpe.append(float(x[f"gpe_{k}"]))
        run_tgpe.append(float(x[f"tgpe_{k}"]))
        run_noise.append(float(x[f"noise_{k}"]))
        matrix_idx = int(x[f"matrix_{k}"])
        run_types.append(_matrix_map_int()[matrix_idx])
        run_matrices.append(_matrix_map_str()[matrix_idx])
        run_vsm_amax.append(float(x[f"vsm_amax_{k}"]) if has_per_run_vsm else shared_vsm)
        run_refine.append(int(x[f"refine_{k}"]) if has_per_run_refine else shared_refine)

    # --- Masking rules ---

    # Single-run: no tree noise (deterministic tree), no min_support
    if n_runs == 1:
        run_noise = [0.0]
        min_support = 0
    else:
        min_support = min(min_support_raw, n_runs)

    # When realign > 0, noise is ineffective (alignment-derived tree)
    if realign > 0:
        run_noise = [0.0] * n_runs

    # When consistency == 0, consistency_weight is irrelevant
    if consistency == 0:
        consistency_weight = 1.0

    # Single-run: no consistency merge (needs ensemble)
    if n_runs == 1:
        consistency_merge = 0

    # When consistency_merge == 0, weight is irrelevant
    if consistency_merge == 0:
        consistency_merge_weight = 2.0

    return {
        "n_runs": n_runs,
        "run_gpo": run_gpo,
        "run_gpe": run_gpe,
        "run_tgpe": run_tgpe,
        "run_noise": run_noise,
        "run_types": run_types,
        "run_matrices": run_matrices,
        "run_vsm_amax": run_vsm_amax,
        "run_refine": run_refine,
        "seq_weights": seq_weights,
        "consistency_weight": consistency_weight,
        "consistency": consistency,
        "realign": realign,
        "min_support": min_support,
        "consistency_merge": consistency_merge,
        "consistency_merge_weight": consistency_merge_weight,
    }


def encode_unified_params(params, max_runs: int) -> dict:
    """Encode a unified params dict into a mixed-variable dict."""
    x: dict = {}
    for k in range(max_runs):
        if k < params["n_runs"]:
            x[f"gpo_{k}"] = params["run_gpo"][k]
            x[f"gpe_{k}"] = params["run_gpe"][k]
            x[f"tgpe_{k}"] = params["run_tgpe"][k]
            x[f"noise_{k}"] = params["run_noise"][k]
            x[f"matrix_{k}"] = _matrix_map_int().index(params["run_types"][k])
            x[f"vsm_amax_{k}"] = params["run_vsm_amax"][k]
            x[f"refine_{k}"] = params["run_refine"][k]
        else:
            x[f"gpo_{k}"] = params["run_gpo"][0]
            x[f"gpe_{k}"] = params["run_gpe"][0]
            x[f"tgpe_{k}"] = params["run_tgpe"][0]
            x[f"noise_{k}"] = params["run_noise"][0]
            x[f"matrix_{k}"] = _matrix_map_int().index(params["run_types"][0])
            x[f"vsm_amax_{k}"] = params["run_vsm_amax"][0]
            x[f"refine_{k}"] = params["run_refine"][0]

    x["seq_weights"] = params["seq_weights"]
    x["consistency_weight"] = params["consistency_weight"]
    x["n_runs"] = params["n_runs"]
    x["consistency"] = params["consistency"]
    x["realign"] = params["realign"]
    x["min_support"] = params["min_support"]
    x["consistency_merge"] = params.get("consistency_merge", 0)
    x["consistency_merge_weight"] = params.get("consistency_merge_weight", 2.0)
    return x


def mode_label(params):
    """Short mode label: 'single', 'ens3', 'ens5'."""
    n = params["n_runs"]
    return "single" if n == 1 else f"ens{n}"


def format_unified_short(params):
    """Compact one-line summary."""
    n_runs = params["n_runs"]

    if n_runs == 1:
        mat = _matrix_names().get(params["run_types"][0], "?")
        ref = REFINE_NAMES.get(params["run_refine"][0], "?")
        return (f"{mode_label(params)} {mat} gpo={params['run_gpo'][0]:.1f} "
                f"vsm={params['run_vsm_amax'][0]:.1f} sw={params['seq_weights']:.1f} "
                f"c={params['consistency']} re={params['realign']} ref={ref}")
    else:
        # Show per-run refine modes compactly
        refs = "/".join(REFINE_NAMES.get(r, "?") for r in params["run_refine"])
        vsms = "/".join(f"{v:.1f}" for v in params["run_vsm_amax"])
        cm = "CM" if params.get("consistency_merge", 0) else "POAR"
        return (f"{mode_label(params)} vsm={vsms} "
                f"sw={params['seq_weights']:.1f} c={params['consistency']} "
                f"re={params['realign']} ref={refs} ms={params['min_support']} {cm}")


def format_unified_long(params):
    """Verbose multi-line summary."""
    lines = [f"mode={mode_label(params)} n_runs={params['n_runs']}"]
    for k in range(params["n_runs"]):
        mat = _matrix_names().get(params["run_types"][k], "?")
        ref = REFINE_LONG.get(params["run_refine"][k], "?")
        lines.append(f"  run_{k}: gpo={params['run_gpo'][k]:.3f} "
                     f"gpe={params['run_gpe'][k]:.3f} "
                     f"tgpe={params['run_tgpe'][k]:.3f} "
                     f"noise={params['run_noise'][k]:.3f} {mat} "
                     f"vsm={params['run_vsm_amax'][k]:.3f} ref={ref}")
    lines.append(f"  seq_weights={params['seq_weights']:.3f}")
    lines.append(f"  consistency={params['consistency']} "
                 f"consistency_weight={params['consistency_weight']:.3f}")
    lines.append(f"  realign={params['realign']} "
                 f"min_support={params['min_support']}")
    if params.get("consistency_merge", 0):
        lines.append(f"  consistency_merge=1 "
                     f"consistency_merge_weight={params['consistency_merge_weight']:.3f}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Stratified k-fold CV (shared with optimize_params)
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

def evaluate_unified(params, cases, n_threads=1, quiet=True):
    """Run kalign with unified params on all cases, return mean metrics."""
    results_by_cat: Dict[str, list] = {}
    total_time = 0.0

    for case in cases:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / f"{case.family}_aln.fasta"

            try:
                start = time.perf_counter()

                # Always use ensemble_custom_file_to_file — it handles
                # n_runs=1 just fine and is the only path that accepts
                # all fine-grained optimizer parameters.
                ensemble_custom_file_to_file(
                    str(case.unaligned),
                    str(output),
                    run_gpo=params["run_gpo"],
                    run_gpe=params["run_gpe"],
                    run_tgpe=params["run_tgpe"],
                    run_noise=params["run_noise"],
                    run_types=params["run_types"],
                    format="fasta",
                    seq_type=_active_profile["seq_type_int"],
                    seed=42,
                    min_support=params["min_support"],
                    realign=params["realign"],
                    seq_weights=params["seq_weights"],
                    n_threads=n_threads,
                    consistency_anchors=params["consistency"],
                    consistency_weight=params["consistency_weight"],
                    # Per-run overrides
                    run_vsm_amax=params["run_vsm_amax"],
                    run_refine=params["run_refine"],
                    # MSA consistency merge
                    consistency_merge=params.get("consistency_merge", 0),
                    consistency_merge_weight=params.get("consistency_merge_weight", 2.0),
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


def fbeta_score(precision, recall, beta=1.0):
    """Compute F-beta score: weights recall beta times more than precision."""
    if precision + recall == 0:
        return 0.0
    b2 = beta * beta
    return (1 + b2) * precision * recall / (b2 * precision + recall)


def evaluate_cv(params, folds, n_threads=1, quiet=True):
    """Evaluate unified params using stratified k-fold CV."""
    fold_f1s = []
    fold_tcs = []
    fold_recalls = []
    fold_precisions = []
    total_time = 0.0

    for _, test in folds:
        result = evaluate_unified(params, test, n_threads, quiet)
        fold_f1s.append(result["f1"])
        fold_tcs.append(result["tc"])
        fold_recalls.append(result["recall"])
        fold_precisions.append(result["precision"])
        total_time += result["wall_time"]

    return {
        "f1": float(np.mean(fold_f1s)),
        "tc": float(np.mean(fold_tcs)),
        "recall": float(np.mean(fold_recalls)),
        "precision": float(np.mean(fold_precisions)),
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
    """Rich-based live terminal dashboard for unified optimization."""

    def __init__(self, n_gen: int, pop_size: int,
                 baselines: Dict[str, dict],
                 max_runs: int):
        self.n_gen = n_gen
        self.pop_size = pop_size
        self.baselines = baselines  # {"fast": {...}, "accurate": {...}, "ensemble": {...}}
        self.max_runs = max_runs
        self.console = Console()

        # State
        self.current_gen = 0
        self.eval_count = 0
        self.total_evals = n_gen * pop_size
        self.gen_start_time = time.time()
        self.run_start_time = time.time()
        self.current_eval_idx = 0

        # Best-ever tracking
        self.best_f1 = 0.0
        self.best_f1_entry: Optional[dict] = None
        self.best_tc = 0.0
        self.best_tc_entry: Optional[dict] = None
        self.fastest = float("inf")
        self.fastest_entry: Optional[dict] = None

        # Pareto front (updated per generation)
        self.pareto_front: List[dict] = []

        # Generation history
        self.gen_history: List[dict] = []

        # Recent evaluations (ring buffer)
        self.recent_evals: List[dict] = []

        # Progress bar
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
        )
        self.gen_task = self.progress.add_task("Generations", total=n_gen)
        self.eval_task = self.progress.add_task("Evaluations", total=pop_size)

        self.live = Live(self._build_layout(), console=self.console, refresh_per_second=2)

    def start(self):
        self.run_start_time = time.time()
        self.live.start()

    def stop(self):
        self.live.stop()

    def _build_baselines_panel(self):
        lines = []
        for name, bl in self.baselines.items():
            lines.append(f"  {name:10s} F1={bl['f1']:.4f}  TC={bl['tc']:.4f}  "
                         f"Time={bl['wall_time']:.0f}s")
        return Panel("\n".join(lines) if lines else "(computing...)", title="Baselines")

    def _build_status_panel(self):
        elapsed = time.time() - self.run_start_time
        elapsed_h = elapsed / 3600
        if self.eval_count > 0:
            rate = elapsed / self.eval_count
            remaining = (self.total_evals - self.eval_count) * rate
            eta_h = remaining / 3600
        else:
            rate = 0
            eta_h = 0

        gen_elapsed = time.time() - self.gen_start_time

        text = (f"Gen {self.current_gen}/{self.n_gen}  "
                f"Eval {self.eval_count}/{self.total_evals}\n"
                f"Elapsed {elapsed_h:.1f}h  ETA {eta_h:.1f}h\n"
                f"Gen time {gen_elapsed:.0f}s  Rate {rate:.1f}s/eval")
        return Panel(text, title="Progress")

    def _format_best_entry(self, label, e, delta_str):
        """Format a best-of entry with full parameter details."""
        p = e["params"]
        header = f"{label}  {mode_label(p)}  {e['wall_time']:.0f}s  {delta_str}"
        # Per-run details
        run_parts = []
        for k in range(p["n_runs"]):
            mat = MATRIX_NAMES.get(p["run_types"][k], "?")
            ref = REFINE_NAMES.get(p["run_refine"][k], "?")
            noise_str = f" n={p['run_noise'][k]:.2f}" if p["run_noise"][k] > 0 else ""
            run_parts.append(f"  R{k}: gpo={p['run_gpo'][k]:.2f} "
                             f"gpe={p['run_gpe'][k]:.2f} "
                             f"tgpe={p['run_tgpe'][k]:.2f}{noise_str} {mat} "
                             f"vsm={p['run_vsm_amax'][k]:.2f} ref={ref}")
        shared = (f"  sw={p['seq_weights']:.2f} "
                  f"c={p['consistency']} cw={p['consistency_weight']:.2f} "
                  f"re={p['realign']} ms={p['min_support']}")
        return "\n".join([header] + run_parts + [shared])

    def _build_best_panel(self):
        lines = []
        if self.best_f1_entry:
            e = self.best_f1_entry
            lines.append(self._format_best_entry(
                f"Best Recall: {e.get('recall', 0):.4f}  F1={e['f1']:.4f}  "
                f"P={e.get('precision', 0):.3f}  TC={e.get('tc', 0):.4f}", e, ""))
        if self.best_tc_entry:
            e = self.best_tc_entry
            lines.append(self._format_best_entry(
                f"Best Prec:   {e.get('precision', 0):.4f}  F1={e['f1']:.4f}  "
                f"R={e.get('recall', 0):.3f}  TC={e.get('tc', 0):.4f}", e, ""))
        if self.fastest_entry:
            e = self.fastest_entry
            lines.append(self._format_best_entry(
                f"Fastest: {e['wall_time']:.1f}s  F1={e['f1']:.4f}  "
                f"R={e.get('recall', 0):.3f}  P={e.get('precision', 0):.3f}", e, ""))
        return Panel("\n".join(lines) if lines else "(no data yet)", title="Best by Objective")

    def _build_pareto_table(self):
        table = Table(title="Pareto Front (top 12 by F1)", box=None, padding=(0, 1))
        table.add_column("#", style="dim", width=2)
        table.add_column("Mode", width=6)
        table.add_column("Recall", justify="right", width=6)
        table.add_column("Prec", justify="right", width=6)
        table.add_column("F1", justify="right", width=6)
        table.add_column("Time", justify="right", width=5)
        table.add_column("c", width=2)
        table.add_column("re", width=2)
        table.add_column("ref", width=3)
        table.add_column("Params", no_wrap=True)

        sorted_front = sorted(self.pareto_front, key=lambda x: -x["f1"])[:12]
        bl_recall = max(bl.get("recall", 0) for bl in self.baselines.values()) if self.baselines else 0
        bl_prec = max(bl.get("precision", 0) for bl in self.baselines.values()) if self.baselines else 0

        for i, entry in enumerate(sorted_front):
            p = entry.get("params", {})
            r_style = "bold green" if entry.get("recall", 0) > bl_recall else ""
            p_style = "bold green" if entry.get("precision", 0) > bl_prec else ""
            refs = "/".join(REFINE_NAMES.get(r, "?") for r in p.get("run_refine", [0]))

            # Build detailed params string
            run_strs = []
            for k in range(p.get("n_runs", 1)):
                mat = MATRIX_NAMES.get(p["run_types"][k], "?")
                vsm = p["run_vsm_amax"][k]
                run_strs.append(f"R{k}:{p['run_gpo'][k]:.1f}/{p['run_gpe'][k]:.2f}/"
                                f"{p['run_tgpe'][k]:.2f}/{mat}/v{vsm:.1f}")
            runs = " ".join(run_strs)
            cm = " CM" if p.get("consistency_merge", 0) else ""
            shared = (f"sw={p.get('seq_weights', 0):.1f} "
                      f"ms={p.get('min_support', 0)}{cm}")
            params_str = f"{runs} | {shared}"

            table.add_row(
                str(i),
                mode_label(p),
                Text(f"{entry.get('recall', 0):.4f}", style=r_style),
                Text(f"{entry.get('precision', 0):.4f}", style=p_style),
                f"{entry['f1']:.4f}",
                f"{entry.get('wall_time', 0):.0f}s",
                str(p.get("consistency", 0)),
                str(p.get("realign", 0)),
                refs,
                params_str,
            )
        return table

    def _build_mode_panel(self):
        # Count modes in Pareto front and recent
        pareto_modes: Dict[str, int] = defaultdict(int)
        for entry in self.pareto_front:
            pareto_modes[mode_label(entry.get("params", {}))] += 1
        parts = [f"{m}={n}" for m, n in sorted(pareto_modes.items())]
        return Panel(f"Pareto: {' '.join(parts)}" if parts else "(none)",
                     title="Mode Distribution")

    def _build_trend_panel(self):
        lines = []
        for h in self.gen_history[-6:]:
            lines.append(f"Gen {h['gen']:3d}: F1={h['best_f1']:.4f}  "
                         f"R={h.get('best_recall', 0):.4f}  "
                         f"P={h.get('best_prec', 0):.4f}  "
                         f"n_pareto={h['n_pareto']}")
        return Panel("\n".join(lines) if lines else "(no data)", title="Trend")

    def _build_recent_panel(self):
        lines = []
        for e in self.recent_evals[-5:]:
            p = e.get("params", {})
            ref0 = REFINE_NAMES.get(p["run_refine"][0], "?") if p.get("run_refine") else "?"
            mat = MATRIX_NAMES.get(p["run_types"][0], "?") if p.get("run_types") else "?"
            gpo = p["run_gpo"][0] if p.get("run_gpo") else 0
            vsm0 = p["run_vsm_amax"][0] if p.get("run_vsm_amax") else 0
            lines.append(f"F1={e['f1']:.4f} TC={e['tc']:.4f} "
                         f"t={e['wall_time']:.0f}s {mode_label(p)} "
                         f"gpo={gpo:.1f} {mat} "
                         f"vsm={vsm0:.1f} "
                         f"c={p.get('consistency', 0)} re={p.get('realign', 0)} ref={ref0}")
        return Panel("\n".join(lines) if lines else "(none)", title="Recent")

    def _build_layout(self):
        layout = Layout()
        layout.split_column(
            Layout(name="top", size=6),
            Layout(name="baselines", size=5),
            Layout(name="best", size=24),
            Layout(name="middle", size=16),
            Layout(name="bottom", size=8),
        )
        layout["top"].split_row(
            Layout(self._build_status_panel(), name="status"),
            Layout(self.progress, name="progress"),
        )
        layout["baselines"].update(self._build_baselines_panel())
        layout["best"].update(self._build_best_panel())
        layout["middle"].split_row(
            Layout(self._build_pareto_table(), name="pareto", ratio=3),
            Layout(self._build_trend_panel(), name="trend", ratio=2),
        )
        layout["bottom"].split_row(
            Layout(self._build_mode_panel(), name="modes"),
            Layout(self._build_recent_panel(), name="recent", ratio=2),
        )
        return layout

    def refresh(self):
        self.live.update(self._build_layout())

    def on_eval_start(self, params: dict, eval_num: int, eval_in_gen: int):
        self.eval_count = eval_num
        self.current_eval_idx = eval_in_gen
        self.progress.update(self.eval_task, completed=eval_in_gen)
        self.refresh()

    def on_eval_end(self, params: dict, cv_result: dict):
        recall = cv_result.get("recall", 0.0)
        precision = cv_result.get("precision", 0.0)
        f1 = fbeta_score(precision, recall, 1.0)
        tc = cv_result.get("tc", 0.0)
        wt = cv_result.get("wall_time", 0.0)

        entry = {"params": params, "f1": f1, "tc": tc, "wall_time": wt,
                 "recall": recall, "precision": precision}
        self.recent_evals.append(entry)
        if len(self.recent_evals) > 5:
            self.recent_evals.pop(0)

        if recall > self.best_f1:  # reuse field for best recall
            self.best_f1 = recall
            self.best_f1_entry = entry
        if precision > self.best_tc:  # reuse field for best precision
            self.best_tc = precision
            self.best_tc_entry = entry
        if wt < self.fastest and f1 > 0.5:
            self.fastest = wt
            self.fastest_entry = entry

        self.refresh()

    def on_gen_start(self, gen: int):
        self.current_gen = gen
        self.gen_start_time = time.time()
        self.progress.update(self.gen_task, completed=gen)
        self.progress.update(self.eval_task, completed=0)
        self.refresh()

    def on_gen_end(self, gen: int, pareto_front: List[dict]):
        self.pareto_front = pareto_front

        best_recall = max((s.get("recall", 0) for s in pareto_front), default=0.0)
        best_prec = max((s.get("precision", 0) for s in pareto_front), default=0.0)
        best_f1_in_gen = max((s["f1"] for s in pareto_front), default=0.0)
        self.gen_history.append({
            "gen": gen,
            "best_f1": best_f1_in_gen,
            "best_recall": best_recall,
            "best_prec": best_prec,
            "n_pareto": len(pareto_front),
        })

        self.progress.update(self.gen_task, completed=gen)
        self.refresh()


# ---------------------------------------------------------------------------
# Parallel evaluation helper (must be top-level for pickling)
# ---------------------------------------------------------------------------

def _eval_one_unified_fold(args_tuple):
    """Evaluate one (individual, fold) pair."""
    ind_idx, fold_idx, x, test_cases, n_threads, max_runs = args_tuple
    import faulthandler, sys
    faulthandler.enable(file=sys.stderr)
    params = decode_unified_params(x, max_runs)
    result = evaluate_unified(params, test_cases, n_threads, quiet=True)
    return ind_idx, fold_idx, params, result


def _eval_baseline(args_tuple):
    """Evaluate one baseline configuration on a set of cases."""
    name, fi, bl_params, test_cases, n_threads = args_tuple
    import faulthandler, sys
    faulthandler.enable(file=sys.stderr)
    result = evaluate_unified(bl_params, test_cases, n_threads, quiet=True)
    return name, fi, result


def _kill_pool(pool: ProcessPoolExecutor) -> None:
    """Forcefully terminate all worker processes."""
    for pid in list(pool._processes):  # noqa: SLF001
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    pool.shutdown(wait=False, cancel_futures=True)


# ---------------------------------------------------------------------------
# pymoo Problem + Callback
# ---------------------------------------------------------------------------

class UnifiedCVProblem(Problem):
    """3-objective optimization with stratified CV evaluation."""

    def __init__(self, folds, max_runs: int, n_threads=1, n_workers=1,
                 dashboard: Optional[Dashboard] = None, f_beta: float = 1.0):
        super().__init__(
            vars=get_vars(max_runs),
            n_obj=3,  # always 3: -recall, -precision, time
        )
        self.folds = folds
        self.max_runs = max_runs
        self.n_threads = n_threads
        self.n_workers = n_workers
        self.dashboard = dashboard
        self.f_beta = f_beta
        self.eval_count = 0
        self.history: List[dict] = []

    def _evaluate(self, X, out, *_args, **_kwargs):
        F = np.zeros((len(X), 3))

        if self.n_workers > 1:
            self._evaluate_parallel(X, F)
        else:
            self._evaluate_serial(X, F)

        out["F"] = F

    def _evaluate_serial(self, X, F):
        for i, x in enumerate(X):
            params = decode_unified_params(x, self.max_runs)
            self.eval_count += 1

            if self.dashboard:
                self.dashboard.on_eval_start(params, self.eval_count, i)

            cv_result = evaluate_cv(params, self.folds, self.n_threads, quiet=True)
            self._record(i, F, params, cv_result)

    def _evaluate_parallel(self, X, F):
        """Fine-grained parallelism: submit (individual x fold) jobs."""
        n_folds = len(self.folds)

        jobs = []
        for i, x in enumerate(X):
            for fold_idx, (_, test) in enumerate(self.folds):
                jobs.append((i, fold_idx, x, test, self.n_threads, self.max_runs))

        if self.dashboard:
            self.dashboard.on_eval_start({}, self.eval_count + 1, 0)

        fold_results: Dict[int, Dict[int, dict]] = defaultdict(dict)
        ind_params: Dict[int, dict] = {}

        pool = ProcessPoolExecutor(max_workers=self.n_workers)
        try:
            futures = {pool.submit(_eval_one_unified_fold, j): j[:2] for j in jobs}
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
                    fold_recalls = [fold_results[ind_idx][fi]["recall"] for fi in range(n_folds)]
                    fold_precisions = [fold_results[ind_idx][fi]["precision"] for fi in range(n_folds)]
                    total_time = sum(fold_results[ind_idx][fi]["wall_time"]
                                     for fi in range(n_folds))

                    cv_result = {
                        "f1": float(np.mean(fold_f1s)),
                        "tc": float(np.mean(fold_tcs)),
                        "recall": float(np.mean(fold_recalls)),
                        "precision": float(np.mean(fold_precisions)),
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
        fb = fbeta_score(cv_result["precision"], cv_result["recall"], self.f_beta)
        cv_result["fbeta"] = fb
        F[i, 0] = -cv_result["recall"]
        F[i, 1] = -cv_result["precision"]
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
                 problem: Optional["UnifiedCVProblem"] = None,
                 max_runs: int = 5):
        super().__init__()
        self.dashboard = dashboard
        self.checkpoint_path = checkpoint_path
        self.problem = problem
        self.max_runs = max_runs

    def notify(self, algorithm):
        gen = algorithm.n_gen

        if self.dashboard:
            self.dashboard.on_gen_start(gen)

        # Extract Pareto front
        pareto = []
        if algorithm.opt is not None and len(algorithm.opt) > 0:
            for ind in algorithm.opt:
                params = decode_unified_params(ind.X, self.max_runs)
                recall = -ind.F[0]
                precision = -ind.F[1]
                f1 = fbeta_score(precision, recall, 1.0)
                entry = {
                    "params": params,
                    "recall": recall,
                    "precision": precision,
                    "f1": f1,
                    "tc": 0.0,  # not an objective; filled from cv_result if available
                    "wall_time": ind.F[2],
                }
                pareto.append(entry)

        if self.dashboard:
            self.dashboard.on_gen_end(gen, pareto)

        # Save checkpoint
        if self.checkpoint_path:
            pop = algorithm.pop
            # Deep-copy dicts in X (mixed-variable: array of dicts)
            raw_X = pop.get("X")
            pop_X = np.array([dict(d) for d in raw_X], dtype=object)
            ckpt = {
                "format": "mixed_v2",
                "pop_X": pop_X,
                "pop_F": pop.get("F").copy(),
                "pop_G": pop.get("G"),
                "pop_H": pop.get("H"),
                "n_gen_completed": gen,
                "history": self.problem.history if self.problem else [],
                "pop_size": len(pop),
                "max_runs": self.max_runs,
                "profile": _active_profile.get("seq_type_str", "protein"),
                "f_beta": self.problem.f_beta if self.problem else 1.0,
            }
            tmp = self.checkpoint_path.with_suffix(".tmp")
            with open(tmp, "wb") as f:
                pickle.dump(ckpt, f)
            tmp.rename(self.checkpoint_path)


def load_checkpoint(path: Path):
    """Load a generation checkpoint.

    Backfills missing keys for new decision variables so old checkpoints
    remain compatible with the current variable space.
    """
    with open(path, "rb") as f:
        ckpt = pickle.load(f)  # noqa: S301

    # Backfill new variables added after the checkpoint was saved.
    # Each individual in pop_X is a dict of decision variables.
    # pymoo's crossover directly indexes parent.X[var_name], so ALL
    # variables in get_vars() must exist in every individual dict.
    max_runs = ckpt.get("max_runs", 5)
    pop_X = ckpt.get("pop_X")
    if pop_X is not None:
        for x in pop_X:
            if not isinstance(x, dict):
                continue

            # v1 → v2: shared vsm_amax/refine → per-run arrays
            if "vsm_amax_0" not in x:
                shared_vsm = float(x.get("vsm_amax", 0.0))
                shared_ref = int(x.get("refine", REFINE_NONE))
                for k in range(max_runs):
                    x[f"vsm_amax_{k}"] = shared_vsm
                    x[f"refine_{k}"] = shared_ref

            # MSA consistency merge (added v3.6)
            if "consistency_merge" not in x:
                x["consistency_merge"] = 0
            if "consistency_merge_weight" not in x:
                x["consistency_merge_weight"] = 2.0

    return ckpt


# ---------------------------------------------------------------------------
# Baseline definitions
# ---------------------------------------------------------------------------

BASELINE_CONFIGS_PROTEIN = {
    "fast": {
        "n_runs": 1,
        "run_gpo": [7.0], "run_gpe": [1.25], "run_tgpe": [1.0], "run_noise": [0.0],
        "run_types": [MATRIX_PFASUM60], "run_matrices": ["pfasum60"],
        "run_vsm_amax": [2.0], "run_refine": [REFINE_NONE],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 0, "realign": 0, "min_support": 0,
    },
    "accurate": {
        "n_runs": 1,
        "run_gpo": [8.472], "run_gpe": [0.554], "run_tgpe": [0.409], "run_noise": [0.0],
        "run_types": [MATRIX_PFASUM60], "run_matrices": ["pfasum60"],
        "run_vsm_amax": [1.359], "run_refine": [REFINE_NONE],
        "seq_weights": 3.407, "consistency_weight": 1.167,
        "consistency": 8, "realign": 2, "min_support": 0,
    },
    "ensemble": {
        "n_runs": 3,
        "run_gpo": [7.0, 3.5, 10.5], "run_gpe": [1.25, 2.5, 0.625],
        "run_tgpe": [1.0, 2.0, 0.5], "run_noise": [0.0, 0.15, 0.15],
        "run_types": [MATRIX_PFASUM60, MATRIX_PFASUM60, MATRIX_PFASUM60],
        "run_matrices": ["pfasum60", "pfasum60", "pfasum60"],
        "run_vsm_amax": [2.0, 2.0, 2.0],
        "run_refine": [REFINE_CONFIDENT, REFINE_CONFIDENT, REFINE_CONFIDENT],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 0, "realign": 1, "min_support": 0,
    },
}

BASELINE_CONFIGS_RNA = {
    "fast": {
        "n_runs": 1,
        "run_gpo": [7.0], "run_gpe": [1.25], "run_tgpe": [1.0], "run_noise": [0.0],
        "run_types": [MATRIX_RNA], "run_matrices": ["rna"],
        "run_vsm_amax": [0.0], "run_refine": [REFINE_NONE],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 0, "realign": 0, "min_support": 0,
    },
    "accurate": {
        "n_runs": 1,
        "run_gpo": [7.0], "run_gpe": [1.25], "run_tgpe": [1.0], "run_noise": [0.0],
        "run_types": [MATRIX_RNA], "run_matrices": ["rna"],
        "run_vsm_amax": [0.0], "run_refine": [REFINE_NONE],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 8, "realign": 2, "min_support": 0,
    },
    "ensemble": {
        "n_runs": 3,
        "run_gpo": [7.0, 3.5, 10.5], "run_gpe": [1.25, 2.5, 0.625],
        "run_tgpe": [1.0, 2.0, 0.5], "run_noise": [0.0, 0.15, 0.15],
        "run_types": [MATRIX_RNA, MATRIX_RNA, MATRIX_RNA], "run_matrices": ["rna", "rna", "rna"],
        "run_vsm_amax": [0.0, 0.0, 0.0],
        "run_refine": [REFINE_CONFIDENT, REFINE_CONFIDENT, REFINE_CONFIDENT],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 0, "realign": 1, "min_support": 0,
    },
}

BASELINE_CONFIGS_DNA = {
    "fast": {
        "n_runs": 1,
        "run_gpo": [7.0], "run_gpe": [1.25], "run_tgpe": [1.0], "run_noise": [0.0],
        "run_types": [MATRIX_DNA], "run_matrices": ["dna"],
        "run_vsm_amax": [0.0], "run_refine": [REFINE_NONE],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 0, "realign": 0, "min_support": 0,
    },
    "accurate": {
        "n_runs": 1,
        "run_gpo": [7.0], "run_gpe": [1.25], "run_tgpe": [1.0], "run_noise": [0.0],
        "run_types": [MATRIX_DNA], "run_matrices": ["dna"],
        "run_vsm_amax": [0.0], "run_refine": [REFINE_NONE],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 8, "realign": 2, "min_support": 0,
    },
    "ensemble": {
        "n_runs": 3,
        "run_gpo": [7.0, 3.5, 10.5], "run_gpe": [1.25, 2.5, 0.625],
        "run_tgpe": [1.0, 2.0, 0.5], "run_noise": [0.0, 0.15, 0.15],
        "run_types": [MATRIX_DNA, MATRIX_DNA, MATRIX_DNA], "run_matrices": ["dna", "dna", "dna"],
        "run_vsm_amax": [0.0, 0.0, 0.0],
        "run_refine": [REFINE_CONFIDENT, REFINE_CONFIDENT, REFINE_CONFIDENT],
        "seq_weights": 0.0, "consistency_weight": 2.0,
        "consistency": 8, "realign": 1, "min_support": 0,
    },
}

def get_baseline_configs(dataset: str) -> dict:
    """Return baseline configs appropriate for the dataset."""
    if dataset == "bralibase":
        return BASELINE_CONFIGS_RNA
    if dataset in ("mdsa", "nucleotide"):
        return BASELINE_CONFIGS_DNA
    return BASELINE_CONFIGS_PROTEIN

# Default for backward compat
BASELINE_CONFIGS = BASELINE_CONFIGS_PROTEIN


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Unified kalign hyperparameter optimization with NSGA-II (3 objectives)")
    parser.add_argument("--max-runs", type=int, default=5,
                        help="Max ensemble runs. n_runs choices: {1,3,5} (default: 5)")
    parser.add_argument("--pop-size", type=int, default=200,
                        help="Population size (default: 200)")
    parser.add_argument("--n-gen", type=int, default=100,
                        help="Number of generations (default: 100)")
    parser.add_argument("--n-folds", type=int, default=5,
                        help="Number of CV folds (default: 5)")
    parser.add_argument("--n-threads", type=int, default=1,
                        help="OpenMP threads per kalign alignment (default: 1)")
    parser.add_argument("--n-workers", type=int, default=1,
                        help="Parallel worker processes (default: 1)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--output-dir", type=str, default="benchmarks/results/unified_optim",
                        help="Output directory for results")
    parser.add_argument("--run-name", type=str, default=None,
                        help="Name for this run (creates subdirectory)")
    parser.add_argument("--no-dashboard", action="store_true",
                        help="Disable rich dashboard, use plain text output")
    parser.add_argument("--f-beta", type=float, default=1.0,
                        help="F-beta weight for recall vs precision. "
                             "1.0=standard F1, 1.5=recall-weighted, 2.0=F2 (default: 1.0)")
    parser.add_argument("--dataset", type=str, default="balibase",
                        choices=["balibase", "bralibase", "mdsa", "nucleotide"],
                        help="Benchmark dataset (default: balibase)")
    parser.add_argument("--resume", type=str, default=None,
                        help="Resume from a generation checkpoint file (.pkl)")
    parser.add_argument("--seed-from", type=str, default=None,
                        help="Seed initial population from another checkpoint (e.g. a beta=1.0 run). "
                             "Parameters are kept, fitness is re-evaluated with current settings.")
    args = parser.parse_args()

    console = Console()

    # --- Dataset setup ---
    if args.dataset == "balibase":
        set_active_profile("protein")
        if not balibase_is_available():
            console.print("Downloading BAliBASE...")
            balibase_download()
        cases = balibase_cases()
    elif args.dataset == "bralibase":
        set_active_profile("rna")
        if not bralibase_is_available():
            console.print("Downloading BRAliBASE...")
            bralibase_download()
        cases = bralibase_cases()
    elif args.dataset == "mdsa":
        set_active_profile("dna")
        if not mdsa_is_available():
            console.print("Downloading MDSA...")
            mdsa_download()
        cases = mdsa_cases()
    elif args.dataset == "nucleotide":
        set_active_profile("nucleotide")
        if not bralibase_is_available():
            console.print("Downloading BRAliBASE...")
            bralibase_download()
        if not mdsa_is_available():
            console.print("Downloading MDSA...")
            mdsa_download()
        cases = bralibase_cases() + mdsa_cases()
        console.print(f"  Combined: {len([c for c in cases if c.seq_type == 'rna'])} RNA + "
                      f"{len([c for c in cases if c.seq_type == 'dna'])} DNA cases")
    else:
        console.print(f"[bold red]Unknown dataset: {args.dataset}[/]")
        return

    # --- Run configuration banner ---
    variables = get_vars(args.max_runs)
    var_names = sorted(variables.keys())
    has_cm = "consistency_merge" in variables
    console.print()
    console.print("[bold cyan]=" * 72)
    console.print("[bold cyan]  Kalign Unified Optimizer")
    console.print("[bold cyan]=" * 72)
    console.print(f"  Dataset:     [bold]{args.dataset}[/] ({_active_profile['seq_type_str']})")
    console.print(f"  Pop size:    {args.pop_size}   Generations: {args.n_gen}")
    console.print(f"  Workers:     {args.n_workers}   Threads/eval: {args.n_threads}")
    console.print(f"  Max runs:    {args.max_runs}")
    console.print(f"  Objectives:  recall, precision, time")
    console.print(f"  Variables:   {len(var_names)}")
    console.print(f"  Resume:      {args.resume or 'no'}")
    console.print(f"  Seed from:   {args.seed_from or 'no'}")
    cm_status = "[bold green]ENABLED[/]" if has_cm else "[bold red]DISABLED[/]"
    console.print(f"  Consistency merge: {cm_status}")
    if has_cm:
        console.print(f"    - consistency_merge: Choice([0, 1])")
        console.print(f"    - consistency_merge_weight: Real([0.5, 10.0])")
    console.print("[bold cyan]" + "-" * 72)
    console.print()

    console.print(f"Loaded [bold]{len(cases)}[/] {args.dataset} cases "
                  f"(profile: {_active_profile['seq_type_str']})")

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
    else:
        output_dir = output_dir / args.dataset
    output_dir.mkdir(parents=True, exist_ok=True)

    max_runs = args.max_runs

    # --- Baseline evaluations (parallelized) ---
    bl_configs = get_baseline_configs(args.dataset)
    console.print(f"\n[bold]Baseline evaluations[/] ({k}-fold CV)")
    n_baseline_workers = max(1, args.n_workers)
    baselines: Dict[str, dict] = {}

    # Build all baseline jobs: (name, fold_idx, cases)
    baseline_jobs = []
    for name, bl_params in bl_configs.items():
        for fi, (_, test) in enumerate(folds):
            baseline_jobs.append((name, fi, bl_params, test, args.n_threads))
        # Full-dataset as extra fold
        baseline_jobs.append((name, k, bl_params, cases, args.n_threads))

    if n_baseline_workers > 1:
        console.print(f"  Running {len(baseline_jobs)} baseline jobs in parallel "
                      f"({n_baseline_workers} workers)...")
        bl_fold_results: Dict[str, Dict[int, dict]] = defaultdict(dict)
        with ProcessPoolExecutor(max_workers=min(n_baseline_workers, len(baseline_jobs))) as pool:
            futures = {pool.submit(_eval_baseline, j): j[:2] for j in baseline_jobs}
            for future in as_completed(futures):
                name, fi, result = future.result()
                bl_fold_results[name][fi] = result
    else:
        bl_fold_results = defaultdict(dict)
        for j in baseline_jobs:
            name, fi, result = _eval_baseline(j)
            bl_fold_results[name][fi] = result

    for name in bl_configs:
        fold_f1s = [bl_fold_results[name][fi]["f1"] for fi in range(k)]
        fold_tcs = [bl_fold_results[name][fi]["tc"] for fi in range(k)]
        total_time = sum(bl_fold_results[name][fi]["wall_time"] for fi in range(k))
        baselines[name] = {
            "f1": float(np.mean(fold_f1s)),
            "tc": float(np.mean(fold_tcs)),
            "f1_std": float(np.std(fold_f1s)),
            "tc_std": float(np.std(fold_tcs)),
            "wall_time": total_time,
        }
        bl_full = bl_fold_results[name][k]
        console.print(f"  {name:10s}  CV F1={baselines[name]['f1']:.4f}+-{baselines[name]['f1_std']:.4f}  "
                      f"CV TC={baselines[name]['tc']:.4f}  "
                      f"Time={baselines[name]['wall_time']:.0f}s  "
                      f"Full F1={bl_full['f1']:.4f}")

    # --- Optimization ---
    n_evals = args.pop_size * args.n_gen
    # Use accurate baseline time as estimate
    est_sec_per_eval = baselines.get("accurate", {}).get("wall_time", 100)
    parallelism = max(1, args.n_workers)
    est_hours = n_evals * est_sec_per_eval / parallelism / 3600

    console.print(f"\n[bold]Starting NSGA-III[/]: pop_size={args.pop_size}, n_gen={args.n_gen}, "
                  f"{k}-fold CV, {args.n_workers} worker(s) x {args.n_threads} thread(s)")
    console.print(f"Total evaluations: ~{n_evals}")
    console.print(f"Estimated time: ~{est_hours:.1f} hours "
                  f"(~{est_sec_per_eval:.0f}s per eval, {parallelism}x parallel)\n")

    # Set up dashboard or plain mode
    use_dashboard = not args.no_dashboard
    dashboard = None

    if use_dashboard:
        dashboard = Dashboard(
            n_gen=args.n_gen,
            pop_size=args.pop_size,
            baselines=baselines,
            max_runs=max_runs,
        )

    if args.f_beta != 1.0:
        console.print(f"[bold yellow]Using F-beta objective with β={args.f_beta}[/] "
                      f"(recall weighted {args.f_beta}x more than precision)")

    problem = UnifiedCVProblem(
        folds=folds,
        max_runs=max_runs,
        n_threads=args.n_threads,
        n_workers=args.n_workers,
        dashboard=dashboard,
        f_beta=args.f_beta,
    )

    checkpoint_path = output_dir / "gen_checkpoint.pkl"
    callback = GenerationCallback(
        dashboard=dashboard,
        checkpoint_path=checkpoint_path,
        problem=problem,
        max_runs=max_runs,
    )

    # NSGA-III reference directions for 3 objectives
    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=12)
    n_ref = len(ref_dirs)  # 91 for n_partitions=12
    mixed_mating = MixedVariableMating(
        eliminate_duplicates=MixedVariableDuplicateElimination())
    mixed_dedup = MixedVariableDuplicateElimination()

    # Resume from checkpoint or start fresh
    resumed_gen = 0
    if args.resume:
        resume_path = Path(args.resume)
        if not resume_path.exists():
            console.print(f"[bold red]Checkpoint not found:[/] {resume_path}")
            return
        ckpt = load_checkpoint(resume_path)
        ckpt_fmt = ckpt.get("format")
        if ckpt_fmt not in ("mixed_v1", "mixed_v2"):
            console.print("[bold red]Cannot resume from old-format checkpoint.[/]")
            console.print("Old checkpoints used float arrays; new format uses mixed variables.")
            console.print("Please start a fresh optimization run (remove --resume).")
            return
        if ckpt_fmt == "mixed_v1":
            console.print("[bold yellow]Note:[/] resuming from v1 checkpoint. "
                          "Old shared vsm_amax/refine will be expanded to per-run arrays.")
        if ckpt.get("max_runs") != max_runs:
            console.print(f"[bold red]max_runs mismatch:[/] checkpoint has "
                          f"{ckpt.get('max_runs')}, requested {max_runs}")
            return
        ckpt_profile = ckpt.get("profile", "protein")
        if ckpt_profile != _active_profile["seq_type_str"]:
            console.print(f"[bold red]Profile mismatch:[/] checkpoint has "
                          f"'{ckpt_profile}', current dataset uses "
                          f"'{_active_profile['seq_type_str']}'")
            return
        pop_X = ckpt["pop_X"]
        pop_F = ckpt["pop_F"]
        pop_G = ckpt.get("pop_G")
        pop_H = ckpt.get("pop_H")
        resumed_gen = ckpt["n_gen_completed"]
        problem.history = ckpt.get("history", [])
        problem.eval_count = len(problem.history)
        console.print(f"[bold green]Resumed[/] from generation {resumed_gen} "
                      f"({len(problem.history)} prior evaluations)")
        remaining = args.n_gen - resumed_gen
        if remaining <= 0:
            console.print(f"[bold yellow]Already completed {resumed_gen} generations "
                          f"(requested {args.n_gen}). Increase --n-gen to continue.[/]")
            return
        termination = get_termination("n_gen", remaining)
        # Reconstruct evaluated population so pymoo skips re-evaluation
        pop = Population.new("X", pop_X)
        pop.set("F", pop_F)
        if pop_G is not None:
            pop.set("G", pop_G)
        if pop_H is not None:
            pop.set("H", pop_H)
        for ind in pop:
            ind.evaluated = {"F", "G", "H"}
        algorithm = NSGA3(
            ref_dirs=ref_dirs,
            pop_size=len(pop_X),
            sampling=pop,
            mating=mixed_mating,
            eliminate_duplicates=mixed_dedup,
        )
    else:
        pop_size = args.pop_size
        if pop_size < n_ref:
            console.print(f"[bold yellow]Warning:[/] pop_size ({pop_size}) < reference "
                          f"directions ({n_ref}). Consider --pop-size {n_ref} or larger.")

        # Seed from another checkpoint (parameters only, fitness will be re-evaluated)
        seed_sampling = MixedVariableSampling()
        if args.seed_from:
            seed_path = Path(args.seed_from)
            if not seed_path.exists():
                console.print(f"[bold red]Seed checkpoint not found:[/] {seed_path}")
                return
            seed_ckpt = load_checkpoint(seed_path)
            seed_X = seed_ckpt["pop_X"]
            seed_beta = seed_ckpt.get("f_beta", 1.0)
            console.print(f"[bold green]Seeding[/] initial population with {len(seed_X)} "
                          f"individuals from {seed_path.name} (was β={seed_beta})")
            # Take up to pop_size individuals, sorted by original fitness (best first)
            seed_F = seed_ckpt["pop_F"]
            # Sort by first objective (best = most negative)
            order = np.argsort(seed_F[:, 0])
            seed_X = seed_X[order][:pop_size]

            # Inject diversity for new binary variables: flip ~25% of
            # ensemble individuals to consistency_merge=1 so the optimizer
            # can compare both paths from the start.
            n_flipped = 0
            for idx, x in enumerate(seed_X):
                if not isinstance(x, dict):
                    continue
                if int(x.get("n_runs", 1)) > 1 and idx % 4 == 0:
                    x["consistency_merge"] = 1
                    n_flipped += 1
            if n_flipped:
                console.print(f"  Flipped {n_flipped} individuals to consistency_merge=1")

            # Create unevaluated population — pymoo will evaluate with new f_beta
            seed_pop = Population.new("X", seed_X)
            seed_sampling = seed_pop

        algorithm = NSGA3(
            ref_dirs=ref_dirs,
            pop_size=pop_size,
            sampling=seed_sampling,
            mating=mixed_mating,
            eliminate_duplicates=mixed_dedup,
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
        os._exit(1)  # noqa: SLF001
    finally:
        if dashboard:
            dashboard.stop()

    # --- Results ---
    console.print(f"\n[bold]Optimization complete.[/] "
                  f"{len(res.F)} Pareto-optimal solutions found.\n")

    pareto_configs = []
    for i, (x, f) in enumerate(zip(res.X, res.F)):
        params = decode_unified_params(x, max_runs)
        recall = -f[0]
        precision = -f[1]
        wt = f[2]
        f1 = fbeta_score(precision, recall, 1.0)
        pareto_configs.append({
            "params": params, "recall_cv": recall, "prec_cv": precision,
            "f1_cv": f1, "wall_time": wt,
        })

    # Print Pareto front
    table = Table(title="Pareto Front (sorted by CV F1)")
    table.add_column("#", style="dim", width=3)
    table.add_column("Mode", width=6)
    table.add_column("Recall", justify="right")
    table.add_column("Prec", justify="right")
    table.add_column("F1", justify="right")
    table.add_column("Time", justify="right")
    table.add_column("Parameters")

    sorted_pareto = sorted(pareto_configs, key=lambda x: -x["f1_cv"])
    bl_best_recall = max(bl.get("recall", 0) for bl in baselines.values())
    bl_best_prec = max(bl.get("precision", 0) for bl in baselines.values())

    for i, cfg in enumerate(sorted_pareto[:30]):
        r_style = "bold green" if cfg.get("recall_cv", 0) > bl_best_recall else ""
        p_style = "bold green" if cfg.get("prec_cv", 0) > bl_best_prec else ""
        table.add_row(
            str(i),
            mode_label(cfg["params"]),
            Text(f"{cfg.get('recall_cv', 0):.4f}", style=r_style),
            Text(f"{cfg.get('prec_cv', 0):.4f}", style=p_style),
            f"{cfg['f1_cv']:.4f}",
            f"{cfg['wall_time']:.0f}s",
            format_unified_short(cfg["params"]),
        )
    console.print(table)

    # --- Mode summary ---
    mode_counts: Dict[str, int] = defaultdict(int)
    mode_best_f1: Dict[str, float] = defaultdict(float)
    for cfg in pareto_configs:
        m = mode_label(cfg["params"])
        mode_counts[m] += 1
        mode_best_f1[m] = max(mode_best_f1[m], cfg["f1_cv"])

    console.print("\n[bold]Mode distribution on Pareto front:[/]")
    for m in sorted(mode_counts.keys()):
        console.print(f"  {m}: {mode_counts[m]} solutions, best F1={mode_best_f1[m]:.4f}")

    # --- Re-evaluate top-3 on FULL dataset ---
    console.print(f"\n[bold]Full-dataset evaluation[/] (top 3 by F1, checking for overfit)")
    top3 = sorted_pareto[:3]
    for i, cfg in enumerate(top3):
        full_result = evaluate_unified(cfg["params"], cases, args.n_threads)
        console.print(f"\n  [{i}] {mode_label(cfg['params'])}  "
                      f"CV R={cfg.get('recall_cv', 0):.4f} P={cfg.get('prec_cv', 0):.4f} "
                      f"F1={cfg['f1_cv']:.4f} -> Full F1={full_result['f1']:.4f}")
        gap_f1 = full_result["f1"] - cfg["f1_cv"]
        console.print(f"    Overfit check: F1 {gap_f1:+.4f}")
        for cat, v in sorted(full_result["per_category"].items()):
            console.print(f"      {cat}: F1={v['f1']:.4f} TC={v['tc']:.4f} (n={v['n']})")
        console.print(f"    {format_unified_long(cfg['params'])}")

    # --- Recommended tiers ---
    console.print(f"\n{'='*60}")
    console.print("[bold]Recommended configurations:[/]")

    # Fast: best F1 among solutions under 15s
    fast_candidates = [c for c in sorted_pareto if c["wall_time"] < 15]
    if fast_candidates:
        best_fast = fast_candidates[0]
        console.print(f"\n  [bold]Fast[/] (< 15s): F1={best_fast['f1_cv']:.4f} "
                      f"R={best_fast.get('recall_cv', 0):.4f} P={best_fast.get('prec_cv', 0):.4f} "
                      f"Time={best_fast['wall_time']:.0f}s")
        console.print(f"    {format_unified_short(best_fast['params'])}")

    # Default: best F1 among solutions under 60s
    default_candidates = [c for c in sorted_pareto if c["wall_time"] < 60]
    if default_candidates:
        best_default = default_candidates[0]
        console.print(f"\n  [bold]Default[/] (< 60s): F1={best_default['f1_cv']:.4f} "
                      f"R={best_default.get('recall_cv', 0):.4f} P={best_default.get('prec_cv', 0):.4f} "
                      f"Time={best_default['wall_time']:.0f}s")
        console.print(f"    {format_unified_short(best_default['params'])}")

    # Accurate: best F1 overall
    best_overall = sorted_pareto[0]
    console.print(f"\n  [bold]Accurate[/] (best F1): F1={best_overall['f1_cv']:.4f} "
                  f"R={best_overall.get('recall_cv', 0):.4f} P={best_overall.get('prec_cv', 0):.4f} "
                  f"Time={best_overall['wall_time']:.0f}s")
    console.print(f"    {format_unified_short(best_overall['params'])}")

    # --- Save ---
    results = {
        "pareto_configs": pareto_configs,
        "history": problem.history,
        "baselines": baselines,
        "folds_info": [(len(tr), len(te)) for tr, te in folds],
        "args": vars(args),
        "max_runs": max_runs,
        "profile": _active_profile["seq_type_str"],
        "dataset": args.dataset,
    }
    results_path = output_dir / "optim_results.pkl"
    with open(results_path, "wb") as f:
        pickle.dump(results, f)
    console.print(f"\nResults saved to {results_path}")

    summary_path = output_dir / "pareto_front.txt"
    with open(summary_path, "w") as f:
        f.write(f"# Unified kalign optimization (NSGA-III, 3 objectives: recall, precision, time)\n")
        f.write(f"# pop_size={args.pop_size} n_gen={args.n_gen} max_runs={max_runs} "
                f"n_folds={k} seed={args.seed}\n")
        f.write(f"# Baselines:\n")
        for name, bl in baselines.items():
            f.write(f"#   {name:10s} F1={bl['f1']:.4f} TC={bl['tc']:.4f} "
                    f"Time={bl['wall_time']:.0f}s\n")
        f.write(f"\n")

        for i, cfg in enumerate(sorted_pareto):
            p = cfg["params"]
            f.write(f"[{i}] mode={mode_label(p)} "
                    f"R={cfg.get('recall_cv', cfg.get('f1_cv', 0)):.4f} "
                    f"P={cfg.get('prec_cv', cfg.get('tc_cv', 0)):.4f} "
                    f"F1={cfg.get('f1_cv', 0):.4f} "
                    f"Time={cfg['wall_time']:.0f}s\n")
            f.write(f"    n_runs={p['n_runs']}\n")
            for run_k in range(p["n_runs"]):
                mat = MATRIX_NAMES.get(p["run_types"][run_k], "?")
                ref = REFINE_LONG.get(p["run_refine"][run_k], "?")
                f.write(f"    run_{run_k}: gpo={p['run_gpo'][run_k]:.3f} "
                        f"gpe={p['run_gpe'][run_k]:.3f} "
                        f"tgpe={p['run_tgpe'][run_k]:.3f} "
                        f"noise={p['run_noise'][run_k]:.3f} {mat} "
                        f"vsm={p['run_vsm_amax'][run_k]:.3f} ref={ref}\n")
            f.write(f"    seq_weights={p['seq_weights']:.3f}\n")
            f.write(f"    consistency={p['consistency']} "
                    f"consistency_weight={p['consistency_weight']:.3f}\n")
            f.write(f"    realign={p['realign']} "
                    f"min_support={p['min_support']}\n")
            if p.get("consistency_merge", 0):
                f.write(f"    consistency_merge=1 "
                        f"consistency_merge_weight={p['consistency_merge_weight']:.3f}\n")
            f.write("\n")

    console.print(f"Pareto front saved to {summary_path}")


if __name__ == "__main__":
    main()
