"""One fit of the presets of a pool on one training core: two NSGA-III stages + the selection rule.

    fit(pool, fold, out_dir, core_ids, screen_ids, evaluator, settings) -> selection dict
    python -m recipe.fit --pool protein --fold 0 --out WS/runs/refit/protein/fold0 [--resume]

The procedure (decisions of PLAN_pooled_refit_v1, see docs/build/engine.md):

  0. Leak guard: the fit receives only training case ids (core and screening subset of the fold). It
     asserts that none is a case of the test fold (folds_<pool>.tsv) and the evaluator refuses any case
     outside the core.
  1. The 3.7.0 presets (recipe.presets_370, parameter dicts, same worker path as every candidate) are
     evaluated on the screening subset (-> objective normalisation for the hypervolume) and on the core
     (-> CPU budgets of the selection rule, and candidates of the selection).
  2. Stage 1: single-run presets (protein space.py stage-1 variables; nucleotide nuc_space.py with
     n_runs = 1). NSGA-III, population 56, Das-Dennis p = 5 (56 directions, 4 objectives), seed 1234.
     Initial population = the 3.7.0 single runs (protein fast, default, the 5 runs of accurate;
     nucleotide fast, recall, the runs of default and accurate) + random individuals. Every generation
     is scored on the screening subset. Stop: >= 20 generations and the hypervolume of the stage's
     archive improved < 0.1 % over the last 5 generations, or 100 generations.
  3. The stage-1 front (non-dominated valid individuals of the stage archive on the screening subset,
     screening CPU <= 1.5 x the largest budget) is re-scored on the full core.
  4. Stage 2: ensembles of 3 or 5 runs whose runs are members of the stage-1 front (+ the 3.7.0
     single runs), ensemble-wide single-run parameters from the lead member, ensemble parameters
     searched (protein: min_support, fill_in, fill_anchor_frac, fill_support as space.py; nucleotide:
     min_support, consistency_merge, consistency_merge_weight from nuc_space.py). Seeds: the 3.7.0
     ensembles and the five best (mean z) front members. Same NSGA-III; stop >= 10, < 0.1 % / 5, <= 140 (v1: 40).
  5. The stage-2 front is re-scored on the core (same CPU rule).
  6. Selection per mode on the core: candidates = core-scored stage-1 and stage-2 fronts + the 3.7.0
     presets; non-dominated set on (SP, TC, precision, CPU); within budget = 1.2 x CPU of the 3.7.0
     preset of that mode on the same core, the maximum mean z-score (z over that set) of SP, TC and
     precision; recall additionally needs SP >= the SP of the accurate choice. Ties: lower CPU, key.

Objectives (both pools): SP (= recall of the scorer), TC, precision: per-case means averaged within
benchmark, benchmarks weighted equally (maximised); total single-thread CPU seconds (minimised). F1 is
reported. An individual with any failed case (error or timeout, 600 s per case) is infeasible
(constraint violation). Search candidates (not the 3.7.0 references) are evaluated with two stops that
only ever hit an infeasible or useless individual: at the first failed case, and once the individual's CPU
on the evaluated cases exceeds 3 x the CPU of 3.7.0 accurate on the same cases (above the hypervolume box,
2 x accurate, and above every budget, <= 1.2 x accurate); such an individual is infeasible.
NSGA-III and the hypervolume use the same monotone transform of the objectives (normalised against the
3.7.0 presets on the screening subset: quality axes (1 - q) / (1 - (q_fast - 0.10)), CPU on a log scale
from 0.5 x fast to 2 x default (stage 1) / 2 x accurate (stage 2)); the hypervolume reference point is
(1.1, 1.1, 1.1, 1.1). The transform does not change which points dominate which.

Evaluation: every alignment goes through recipe.evaluate's cache and parameter-set worker
(evaluate.ParamWorker + evaluate._process_case). A machine-wide pool of aligner slots (file locks in
<slots_dir>, default 30 slots) bounds the number of concurrent aligner processes over all fit
processes, so 5 folds (or two pools) running in parallel never exceed 30 aligners.

Resumable: per stage, the pymoo algorithm is pickled after every generation; all evaluations are in
the cache, so a resumed generation costs only cache lookups. Output (out_dir):
  run.json                  settings, fold, case-list hashes, build, space hashes
  refs.json                 3.7.0 presets on screen and core (summaries)
  stage{1,2}/algorithm.pkl  checkpoint;  stage{1,2}/hv.json  hypervolume per generation + stop reason
  stage{1,2}/evaluated.jsonl  every evaluated individual (generation, key, params, screen summary, F)
  stage{1,2}/front_core.json  the re-scored front;  stage1/pool.json  stage-2 member pool
  selection.json            the selection per mode (+ candidate set, budgets, refs)
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import hashlib
import json
import math
import os
import pickle
import random
import sys
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from . import cache as C
from . import presets_370 as P370

# ------------------------------------------------------------------------------------ settings
SEED = 1234
POP_SIZE = 56
REF_PARTITIONS = 5
STAGE_GENS = {1: (20, 100), 2: (10, 140)}   # v1.1 (PLAN Amendment v1.1): stage-2 cap 40 -> 140
HV_WINDOW = 5
HV_TOL = 0.001
HV_REF = 1.1
Q_MARGIN = 0.10
BUDGET_FACTOR = 1.2
RESCORE_MARGIN = 1.5
CPU_CAP_FACTOR = 3.0
FIT_TIMEOUT_S = 600.0
MAX_SLOTS = 30
CHUNK = 20
OBJ = ("sp", "tc", "precision")
METRIC_OF = {"sp": "recall", "tc": "tc", "precision": "precision", "f1": "f1"}
MODES = {"protein": ("fast", "default", "recall", "accurate"),
         "nucleotide": ("fast", "default", "recall", "accurate")}
REF_MODES = {"protein": ("fast", "default", "accurate"), "nucleotide": ("fast", "default", "recall", "accurate")}
BUDGET_REF = {"protein": {"fast": "fast", "default": "default", "recall": "accurate", "accurate": "accurate"},
              "nucleotide": {"fast": "fast", "default": "default", "recall": "recall", "accurate": "accurate"}}
STAGE2_N_RUNS = (3, 5)
PYMOO_VERSION = "0.6.1.6"     # registered runs (PLAN section 7) refuse any other pymoo
STAGE2_MAX = 5


@dataclass
class Settings:
    pop_size: int = POP_SIZE
    ref_partitions: int = REF_PARTITIONS
    seed: int = SEED
    stage1_gens: tuple = STAGE_GENS[1]
    stage2_gens: tuple = STAGE_GENS[2]
    hv_window: int = HV_WINDOW
    hv_tol: float = HV_TOL
    hv_ref: float = HV_REF
    q_margin: float = Q_MARGIN
    budget_factor: float = BUDGET_FACTOR
    rescore_margin: float = RESCORE_MARGIN
    cpu_cap_factor: float = CPU_CAP_FACTOR
    timeout_s: float = FIT_TIMEOUT_S
    fake_eval: bool = False
    note: str = ""

    def gens(self, stage):
        return tuple(self.stage1_gens if stage == 1 else self.stage2_gens)

    def is_registered(self):
        """True iff these are the frozen settings of the procedure (not a smoke test)."""
        d = Settings()
        return all(getattr(self, k) == getattr(d, k) for k in
                   ("pop_size", "ref_partitions", "seed", "hv_window", "hv_tol", "hv_ref", "q_margin",
                    "budget_factor", "rescore_margin", "cpu_cap_factor", "timeout_s", "fake_eval")) and \
            tuple(self.stage1_gens) == STAGE_GENS[1] and tuple(self.stage2_gens) == STAGE_GENS[2] and \
            not self.note   # note is set by every smoke override that is not a field (e.g. --limit-core)

    def to_json(self):
        d = asdict(self)
        d["stage1_gens"], d["stage2_gens"] = list(self.stage1_gens), list(self.stage2_gens)
        return d


class LeakError(RuntimeError):
    """A test-fold case reached a fit."""


def log(msg, fh=None):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
    print(line, flush=True)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


def pkey(p) -> str:
    """Cache param key of a worker dict (identical to what evaluate() stores)."""
    return C.param_key(p)


def ids_hash(ids) -> str:
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()


def _atomic_json(path, obj):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default) + "\n")
    tmp.replace(path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o).__name__)


def _py(x):
    """numpy scalars in a decoded dict -> python scalars (cache keys, JSON)."""
    if isinstance(x, dict):
        return {k: _py(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_py(v) for v in x]
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    return x


# ------------------------------------------------------------------------------------ search spaces
def _mk(kind, b):
    from pymoo.core.variable import Choice, Integer, Real
    return Real(bounds=b) if kind == "real" else Integer(bounds=b) if kind == "int" else Choice(options=b)


class ProteinSpace:
    """Adapter over experiments/redesign/tools/refit/space.py (read-only). Worker dict = decoded params."""
    pool = "protein"

    def __init__(self):
        self.sp = P370.load_space("protein")

    def stage1_vars(self):
        return self.sp.stage1_vars()

    def stage1_decode(self, x):
        return _py(self.sp.stage1_decode(x))

    def stage1_encode(self, p):
        return self.sp.stage1_encode(p)

    def stage1_seeds(self):
        """3.7.0 single runs: fast, default, the five runs of accurate."""
        out = [P370.params_370("protein", "fast"), P370.params_370("protein", "default")]
        out += P370.protein_stage1_members(P370.params_370("protein", "accurate"))
        return [_py(p) for p in out]

    def members_of(self, p):
        """Single-run members of an ensemble worker dict (for the stage-2 pool)."""
        if p["n_runs"] == 1:
            return [p]
        return [_py(q) for q in P370.protein_stage1_members(p)]

    # stage 2: runs = pool members, n_runs in {3, 5}; space.STAGE2_ENSEMBLE searched, STAGE2_FIXED fixed
    def stage2_vars(self, pool_size):
        v = {f"m_{k}": _mk("choice", list(range(pool_size))) for k in range(STAGE2_MAX)}
        v["n_runs"] = _mk("choice", list(STAGE2_N_RUNS))
        v.update({name: _mk(kind, b) for name, kind, b in self.sp.STAGE2_ENSEMBLE})
        return v

    def stage2_decode(self, x, pool):
        sp = self.sp
        n = int(x["n_runs"])
        idx = sp.stage2_members(x, pool)[:n]   # space.py vsm_mode compatibility rule, first n members
        lead = pool[idx[0]]
        per = ["run_gpo", "run_gpe", "run_tgpe", "run_types", "run_vsm_amax", "run_refine", "run_vsm_kappa",
               "run_vsm_dzero", "run_guide_dist"]
        p = {k: [pool[i][k][0] for i in idx] for k in per}
        p.update({
            "n_runs": n, "run_noise": [sp.STAGE2_FIXED["run_noise"]] * n,
            "seq_weights": lead["seq_weights"], "consistency": lead["consistency"],
            "consistency_weight": lead["consistency_weight"], "realign": lead["realign"],
            "vsm_mode": lead["vsm_mode"],
            "consistency_merge": sp.STAGE2_FIXED["consistency_merge"],
            "consistency_merge_weight": sp.STAGE2_FIXED["consistency_merge_weight"],
            "min_support": int(x["min_support"]), "fill_in": int(x["fill_in"]),
            "fill_anchor_frac": float(x["fill_anchor_frac"]), "fill_support": int(x["fill_support"]),
        })
        return _py(sp.mask(p))

    def stage2_ens_params(self, p):
        return {"min_support": p["min_support"], "fill_in": p["fill_in"],
                "fill_anchor_frac": p["fill_anchor_frac"], "fill_support": p["fill_support"]}

    def stage2_seed_ensembles(self):
        return [_py(P370.params_370("protein", "accurate"))]

    def short(self, p):
        return self.sp.short(p)


class NucleotideSpace:
    """Adapter over experiments/redesign/tools/nested_nuc/nuc_space.py (read-only). Worker dict =
    nuc_space.to_worker(decoded) (3.6+ mechanisms at their off values, as the 3.7.0 nucleotide presets)."""
    pool = "nucleotide"
    BASE = ("n_runs", "run_gpo", "run_gpe", "run_tgpe", "run_noise", "run_types", "run_vsm_amax", "run_refine",
            "seq_weights", "consistency", "consistency_weight", "realign", "min_support", "consistency_merge",
            "consistency_merge_weight")

    def __init__(self):
        self.sp = P370.load_space("nucleotide")

    def _base(self, w):
        return {k: w[k] for k in self.BASE}

    def stage1_vars(self):
        sp = self.sp
        return {
            "gpo": _mk("real", (sp.PER_RUN_LO[0], sp.PER_RUN_HI[0])),
            "gpe": _mk("real", (sp.PER_RUN_LO[1], sp.PER_RUN_HI[1])),
            "tgpe": _mk("real", (sp.PER_RUN_LO[2], sp.PER_RUN_HI[2])),
            "matrix": _mk("choice", list(range(len(sp.MATRIX_MAP)))),
            "vsm_amax": _mk("real", (sp.SHARED_LO[0], sp.SHARED_HI[0])),
            "refine": _mk("choice", list(sp.REFINE_MAP)),
            "seq_weights": _mk("real", (sp.SHARED_LO[1], sp.SHARED_HI[1])),
            "consistency": _mk("choice", list(sp.CONSISTENCY_MAP)),
            "consistency_weight": _mk("real", (sp.SHARED_LO[2], sp.SHARED_HI[2])),
            "realign": _mk("int", (0, 2)),
        }

    def stage1_decode(self, x):
        sp = self.sp
        p = {"n_runs": 1, "run_gpo": [float(x["gpo"])], "run_gpe": [float(x["gpe"])],
             "run_tgpe": [float(x["tgpe"])], "run_noise": [0.0], "run_types": [sp.MATRIX_MAP[int(x["matrix"])]],
             "run_vsm_amax": [float(x["vsm_amax"])], "run_refine": [int(x["refine"])],
             "seq_weights": float(x["seq_weights"]), "consistency": int(x["consistency"]),
             "consistency_weight": float(x["consistency_weight"]), "realign": int(x["realign"]),
             "min_support": 0, "consistency_merge": 0, "consistency_merge_weight": 2.0}
        return _py(sp.to_worker(sp.mask(p)))

    def stage1_encode(self, w):
        assert w["n_runs"] == 1
        return {"gpo": w["run_gpo"][0], "gpe": w["run_gpe"][0], "tgpe": w["run_tgpe"][0],
                "matrix": self.sp.MATRIX_MAP.index(w["run_types"][0]), "vsm_amax": w["run_vsm_amax"][0],
                "refine": w["run_refine"][0], "seq_weights": w["seq_weights"], "consistency": w["consistency"],
                "consistency_weight": w["consistency_weight"], "realign": w["realign"]}

    def members_of(self, w):
        """Run k of an ensemble with the ensemble-wide single-run parameters (n_runs = 1)."""
        if w["n_runs"] == 1:
            return [w]
        out = []
        for k in range(w["n_runs"]):
            p = {"n_runs": 1, "run_gpo": [w["run_gpo"][k]], "run_gpe": [w["run_gpe"][k]],
                 "run_tgpe": [w["run_tgpe"][k]], "run_noise": [0.0], "run_types": [w["run_types"][k]],
                 "run_vsm_amax": [w["run_vsm_amax"][k]], "run_refine": [w["run_refine"][k]],
                 "seq_weights": w["seq_weights"], "consistency": w["consistency"],
                 "consistency_weight": w["consistency_weight"], "realign": w["realign"],
                 "min_support": 0, "consistency_merge": 0, "consistency_merge_weight": 2.0}
            out.append(_py(self.sp.to_worker(self.sp.mask(p))))
        return out

    def stage1_seeds(self):
        out = []
        for m in ("fast", "recall", "default", "accurate"):
            out += self.members_of(_py(P370.params_370("nucleotide", m)))
        return out

    def stage2_vars(self, pool_size):
        v = {f"m_{k}": _mk("choice", list(range(pool_size))) for k in range(STAGE2_MAX)}
        v["n_runs"] = _mk("choice", list(STAGE2_N_RUNS))
        v["min_support"] = _mk("int", (0, self.sp.MAX_RUNS))
        v["consistency_merge"] = _mk("choice", [0, 1])
        v["consistency_merge_weight"] = _mk("real", (0.5, 10.0))
        return v

    def stage2_decode(self, x, pool):
        n = int(x["n_runs"])
        idx = [int(x[f"m_{k}"]) for k in range(n)]
        lead = pool[idx[0]]
        p = {k: [pool[i][k][0] for i in idx] for k in
             ("run_gpo", "run_gpe", "run_tgpe", "run_types", "run_vsm_amax", "run_refine")}
        p.update({"n_runs": n, "run_noise": [0.0] * n, "seq_weights": lead["seq_weights"],
                  "consistency": lead["consistency"], "consistency_weight": lead["consistency_weight"],
                  "realign": lead["realign"], "min_support": int(x["min_support"]),
                  "consistency_merge": int(x["consistency_merge"]),
                  "consistency_merge_weight": float(x["consistency_merge_weight"])})
        return _py(self.sp.to_worker(self.sp.mask(p)))

    def stage2_ens_params(self, w):
        return {"min_support": w["min_support"], "consistency_merge": w["consistency_merge"],
                "consistency_merge_weight": w["consistency_merge_weight"]}

    def stage2_seed_ensembles(self):
        return [_py(P370.params_370("nucleotide", m)) for m in ("default", "accurate")]

    def short(self, w):
        return self.sp.short(self._base(w))


def get_space(pool):
    return ProteinSpace() if pool == "protein" else NucleotideSpace()


def stage2_seed_x(space, ens, pool_params):
    """Stage-2 variable dict that decodes to ensemble ``ens`` (members must be in the pool)."""
    keys = [pkey(p) for p in pool_params]
    mem = [keys.index(pkey(m)) for m in space.members_of(ens)]
    x = {f"m_{k}": mem[k] if k < len(mem) else mem[0] for k in range(STAGE2_MAX)}
    x["n_runs"] = ens["n_runs"]
    x.update(space.stage2_ens_params(ens))
    return x


# ------------------------------------------------------------------------------------ objectives
def summarise(rows: dict, meta: dict) -> dict:
    """rows {case_id: cache row}, meta {case_id: {"benchmark": ...}} -> objective summary.

    Means are per-case means within a benchmark, averaged over benchmarks with equal weight. CPU is the
    total over all cases. valid = no case has an error (errors include "timeout")."""
    by_b = {}
    n_err, cpu, case_means = 0, 0.0, {k: [] for k in ("sp", "tc", "precision", "f1")}
    for cid, r in rows.items():
        cpu += float(r.get("cpu_s") or 0.0)
        if r.get("error"):
            n_err += 1
            continue
        b = meta[cid]["benchmark"]
        d = by_b.setdefault(b, {k: [] for k in ("sp", "tc", "precision", "f1")})
        for k in d:
            v = float(r[METRIC_OF[k]])
            d[k].append(v)
            case_means[k].append(v)
    per_b = {b: {k: float(np.mean(v)) for k, v in d.items()} | {"n": len(d["sp"])} for b, d in sorted(by_b.items())}
    s = {k: (float(np.mean([per_b[b][k] for b in per_b])) if per_b else 0.0) for k in ("sp", "tc", "precision", "f1")}
    s.update({"cpu": cpu, "n": len(rows), "n_err": n_err, "valid": n_err == 0 and len(rows) > 0,
              "per_benchmark": per_b,
              "case_mean": {k: (float(np.mean(v)) if v else 0.0) for k, v in case_means.items()}})
    return s


class Normaliser:
    """Monotone transform of the objectives used by NSGA-III and the hypervolume (all minimised)."""

    def __init__(self, refs_screen: dict, stage: int, q_margin=Q_MARGIN):
        fast = refs_screen["fast"]
        self.qlo = {m: fast[m] - q_margin for m in OBJ}
        top = refs_screen["default" if stage == 1 else "accurate"]
        self.c_lo = math.log(0.5 * max(fast["cpu"], 1e-9))
        self.c_hi = math.log(2.0 * max(top["cpu"], 1e-9))
        if self.c_hi <= self.c_lo:
            raise ValueError("CPU normalisation range is empty")

    def F(self, s):
        q = [(1.0 - s[m]) / (1.0 - self.qlo[m]) for m in OBJ]
        c = (math.log(max(s["cpu"], 1e-9)) - self.c_lo) / (self.c_hi - self.c_lo)
        return q + [c]

    def to_json(self):
        return {"qlo": self.qlo, "log_cpu_lo": self.c_lo, "log_cpu_hi": self.c_hi}


def nondominated(F):
    """Indices of the non-dominated rows of F (all objectives minimised)."""
    F = np.asarray(F, dtype=float)
    if len(F) == 0:
        return []
    keep = []
    for i in range(len(F)):
        dom = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
        if not dom.any():
            keep.append(i)
    return keep


def raw_F(s):
    return [-s["sp"], -s["tc"], -s["precision"], s["cpu"]]


def hypervolume(F, ref):
    """Hypervolume of the points of F that dominate the reference point (all minimised)."""
    from pymoo.indicators.hv import HV
    F = np.asarray(F, dtype=float)
    if len(F) == 0:
        return 0.0
    F = F[np.all(F < ref, axis=1)]
    if len(F) == 0:
        return 0.0
    F = F[nondominated(F)]
    return float(HV(ref_point=np.full(F.shape[1], ref))(F))


def hv_stop(history, gmin, gmax, window=HV_WINDOW, tol=HV_TOL):
    """history[i] = hypervolume after generation i+1. -> (stop, reason)."""
    g = len(history)
    if g >= gmax:
        return True, f"max generations ({gmax})"
    if g >= gmin and g > window:
        old, new = history[g - 1 - window], history[g - 1]
        if old > 0 and (new - old) / old < tol:
            return True, f"hypervolume +{(new - old) / old:.5%} over {window} generations (< {tol:.1%})"
    return False, ""


# ------------------------------------------------------------------------------------ selection
def zmean(entries):
    Z = []
    for m in OBJ:
        v = np.array([e[m] for e in entries], dtype=float)
        sd = v.std()
        Z.append((v - v.mean()) / sd if sd > 0 else np.zeros_like(v))
    return (sum(Z) / len(OBJ)) if entries else np.array([])


def select_modes(cands, refs_core, pool, budget_factor=BUDGET_FACTOR):
    """cands: [entry] with key, params, sp, tc, precision, f1, cpu, valid (core scores); refs_core:
    {mode: entry} of the 3.7.0 presets on the same core. -> dict(selection={mode: ...}, front=[...],
    budgets={mode: cpu}). The rule is identical in every fold and in the final fit."""
    valid = [e for e in cands if e["valid"]]
    seen, uniq = set(), []
    for e in sorted(valid, key=lambda e: e["key"]):
        if e["key"] not in seen:
            seen.add(e["key"])
            uniq.append(e)
    front = [uniq[i] for i in nondominated([raw_F(e) for e in uniq])]
    z = zmean(front)
    for e, zz in zip(front, z):
        e["zmean"] = float(zz)
    budgets = {m: budget_factor * refs_core[BUDGET_REF[pool][m]]["cpu"] for m in MODES[pool]}
    sel = {}

    def pick(m, need_sp=None):
        ok = [e for e in front if e["cpu"] <= budgets[m] and (need_sp is None or e["sp"] >= need_sp)]
        if not ok:
            return None
        return sorted(ok, key=lambda e: (-e["zmean"], e["cpu"], e["key"]))[0]

    for m in ("accurate", "fast", "default", "recall"):
        need = None
        if m == "recall":
            acc = sel["accurate"]["entry"]
            need = acc["sp"] if acc is not None else refs_core[BUDGET_REF[pool]["accurate"]]["sp"]
        e = pick(m, need)
        sel[m] = {"entry": e, "budget_cpu": budgets[m], "need_sp": need}
    return {"selection": sel, "front": front, "budgets": budgets}


# ------------------------------------------------------------------------------------ evaluation
class SlotPool:
    """Machine-wide limit on concurrent aligner processes: ``n`` lock files in ``dirpath``; a slot is
    held with an exclusive flock (released automatically if the holder dies). All engine processes use
    the same directory, so together they never run more than ``n`` aligners."""

    def __init__(self, dirpath, n=MAX_SLOTS):
        if n > MAX_SLOTS:
            raise ValueError(f"{n} slots exceed the limit of {MAX_SLOTS} concurrent aligners")
        self.dir = Path(dirpath)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.n = n
        self.paths = [self.dir / f"slot_{i:02d}.lock" for i in range(MAX_SLOTS)]
        self._held = []

    def _try(self, path):
        fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return fd
        except BlockingIOError:
            os.close(fd)
            return None

    @contextlib.contextmanager
    def slot(self):
        rng = random.Random()
        while True:
            order = list(range(self.n))
            rng.shuffle(order)
            for i in order:
                fd = self._try(self.paths[i])
                if fd is not None:
                    try:
                        yield i
                    finally:
                        fcntl.flock(fd, fcntl.LOCK_UN)
                        os.close(fd)
                    return
            time.sleep(0.05 + 0.1 * rng.random())

    def hold_others(self, keep):
        """Block every slot except the first ``keep`` (waits until they are free): used to re-time at a
        fixed low concurrency with no other engine aligner running."""
        for p in self.paths[keep:]:
            while True:
                fd = self._try(p)
                if fd is not None:
                    self._held.append(fd)
                    break
                time.sleep(0.2)
        self.n = keep

    def release_others(self):
        for fd in self._held:
            fcntl.flock(fd, fcntl.LOCK_UN)
            os.close(fd)
        self._held = []


def _spec_cost(p):
    return p["n_runs"] * (1 + p["consistency"]) * (1 + p["realign"]) * (1.5 if any(p["run_refine"]) else 1.0)


def _case_cost(c):
    try:
        return int(c["n_seq"]) * int(c["ref_len"])
    except (KeyError, ValueError, TypeError):
        return 0


class Evaluator:
    """Parameter dicts x cases through the cache and the parameter-set worker, under the slot pool.

    evaluate({key: params}, case_ids) -> {key: {case_id: cache row}}. Refuses any case outside
    ``allowed`` (LeakError)."""

    def __init__(self, build, case_rows, allowed, *, cache_path=None, slots=None, threads=MAX_SLOTS,
                 timeout_s=FIT_TIMEOUT_S, chunk=CHUNK, scorer="recipe.scoring:score_case", logfh=None,
                 rerun_short_timeouts=False):
        from . import evaluate as E
        self.E = E
        self.build = build
        self.rows = {c["case_id"]: c for c in case_rows}
        self.allowed = frozenset(allowed)
        self.cache = C.Cache(cache_path)
        self.slots = slots
        self.threads = min(threads, MAX_SLOTS)
        self.timeout_s = timeout_s
        self.chunk = chunk
        self.scorer = scorer
        self.logfh = logfh
        self._idents = {}
        # out-of-fold tests and re-timing (timeout 1800 s): a cached "timeout" row written under a SHORTER limit
        # (600 s during fitting, possibly by another fold's search) is re-run with this evaluator's limit and
        # replaced, so every reported row had the full limit
        self.rerun_short_timeouts = rerun_short_timeouts
        self._rerun = set()

    def ident(self, cid):
        if cid not in self._idents:
            c = self.rows[cid]
            self._idents[cid] = C.ident(c, self.E.scorer_version(self.scorer, c))
        return self._idents[cid]

    def check(self, case_ids):
        bad = sorted(set(case_ids) - self.allowed)
        if bad:
            raise LeakError(f"evaluation of {len(bad)} case(s) outside the allowed set, e.g. {bad[:3]}")

    def lookup(self, key, case_ids):
        return self.cache.get_many(self.build.build_id, key, [self.ident(c) for c in case_ids])

    def evaluate(self, specs: dict, case_ids, cpu_cap=None, stop_on_error=False):
        """cpu_cap: stop aligning a spec once its CPU over these cases exceeds cpu_cap; stop_on_error: stop
        a spec at its first failed case. The skipped cases get rows with error "skipped: <why>" (not cached);
        both only ever apply to an individual that is infeasible anyway."""
        case_ids = list(case_ids)
        self.check(case_ids)
        for k, p in specs.items():
            if pkey(p) != k:
                raise ValueError(f"spec key mismatch {k[:12]}")
        tasks = []
        self._state = {}
        for k, p in specs.items():
            hits = self.lookup(k, case_ids)
            if self.rerun_short_timeouts:
                short = [c for c, h in hits.items() if h.get("error") == "timeout"
                         and float(h.get("wall_s") or 0.0) < 0.99 * self.timeout_s]
                for c in short:
                    del hits[c]
                    self._rerun.add((k, c))
            st = {"cpu": sum(float(h.get("cpu_s") or 0.0) for h in hits.values()), "abort": None,
                  "cap": cpu_cap, "stop": stop_on_error}
            if stop_on_error and any(h.get("error") for h in hits.values()):
                st["abort"] = "error"
            if cpu_cap is not None and st["cpu"] > cpu_cap:
                st["abort"] = "cpu_cap"
            self._state[k] = st
            if st["abort"]:
                continue
            todo = sorted((c for c in case_ids if c not in hits), key=lambda c: -_case_cost(self.rows[c]))
            # round-robin split keeps chunks of one spec similar in cost
            nch = max(1, math.ceil(len(todo) / self.chunk))
            for j in range(nch):
                ch = todo[j::nch]
                if ch:
                    tasks.append((k, p, ch, _spec_cost(p) * sum(_case_cost(self.rows[c]) for c in ch)))
        if tasks:
            tasks.sort(key=lambda t: -t[3])
            self._run(tasks)
        out = {}
        for k in specs:
            got = self.lookup(k, case_ids)
            missing = [c for c in case_ids if c not in got]
            if missing:
                why = self._state[k]["abort"]
                if not why:
                    raise RuntimeError(f"{len(missing)} results missing after evaluation of {k[:12]}")
                for c in missing:
                    got[c] = {"case_id": c, "error": f"skipped: {why}", "cpu_s": None}
            out[k] = got
        return out

    def _run(self, tasks):
        E, build = self.E, self.build
        n_cases = sum(len(t[2]) for t in tasks)
        t0, done, lock = time.time(), [0], threading.Lock()

        def work(k, p, ch):
            st = self._state[k]
            if st["abort"]:
                return 0
            with (self.slots.slot() if self.slots else contextlib.nullcontext()):
                w = E.ParamWorker(build.python, p)
                try:
                    for cid in ch:
                        if st["abort"]:
                            break
                        r = E._process_case(build, p, self.rows[cid], self.scorer, None, self.cache,
                                            (k, cid) in self._rerun,
                                            lambda inp, out, tmp: w.align(inp, out, self.timeout_s))
                        with lock:
                            done[0] += 1
                            st["cpu"] += float(r.get("cpu_s") or 0.0)
                            if st["stop"] and r.get("error"):
                                st["abort"] = "error"
                            if st["cap"] is not None and st["cpu"] > st["cap"]:
                                st["abort"] = "cpu_cap"
                finally:
                    w.close()
            return len(ch)

        errs = []
        with ThreadPoolExecutor(max_workers=self.threads) as ex:
            futs = [ex.submit(work, k, p, ch) for k, p, ch, _ in tasks]
            last = 0.0
            for f in as_completed(futs):
                try:
                    f.result()
                except Exception:  # noqa: BLE001 - collected, raised after all tasks
                    errs.append(traceback.format_exc()[-800:])
                if time.time() - last > 60:
                    last = time.time()
                    log(f"    aligned {done[0]}/{n_cases} ({time.time() - t0:.0f}s)", self.logfh)
        if errs:
            raise RuntimeError(f"{len(errs)} evaluation tasks failed; first:\n{errs[0]}")


class FakeEvaluator:
    """Deterministic pseudo-results without an aligner (tests and --fake-eval machinery runs).

    Quality rises with ensemble size / consistency / refinement and depends smoothly on the gap
    penalties; CPU grows with the same knobs, so there is a real trade-off. Records every request."""

    def __init__(self, case_rows, allowed, fail_if=None):
        self.rows = {c["case_id"]: c for c in case_rows}
        self.allowed = frozenset(allowed)
        self.fail_if = fail_if
        self.calls = []

    def check(self, case_ids):
        bad = sorted(set(case_ids) - self.allowed)
        if bad:
            raise LeakError(f"evaluation of {len(bad)} case(s) outside the allowed set, e.g. {bad[:3]}")

    def evaluate(self, specs, case_ids, cpu_cap=None, stop_on_error=False):
        case_ids = list(case_ids)
        self.check(case_ids)
        self.calls.append((sorted(specs), case_ids, cpu_cap))
        out = {}
        for k, p in specs.items():
            g = float(np.mean(p["run_gpo"])) / 20.0
            e = float(np.mean(p["run_gpe"])) / 5.0
            base = 0.55 + 0.25 * math.exp(-((g - 0.45) ** 2 + (e - 0.25) ** 2) * 4)
            boost = 0.02 * math.log(p["n_runs"]) + 0.01 * math.log1p(p["consistency"]) + \
                0.01 * (1 if any(p["run_refine"]) else 0) + 0.005 * p["realign"]
            cost = _spec_cost(p) * (1 + 0.3 * g)
            rows = {}
            for cid in case_ids:
                h = int(hashlib.sha256(f"{k}|{cid}".encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
                hb = int(hashlib.sha256(cid.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
                q = min(1.0, max(0.0, base + boost - 0.2 * hb + 0.02 * (h - 0.5)))
                err = "timeout" if self.fail_if and self.fail_if(p, cid) else None
                rows[cid] = {"case_id": cid, "recall": q, "tc": q * 0.8, "precision": q * 0.95,
                             "f1": q * 0.97, "cpu_s": 0.01 * cost * (1 + hb) * (1 + 0.05 * h), "error": err}
            if cpu_cap is not None:
                acc = 0.0
                for cid in case_ids:
                    if acc > cpu_cap:
                        rows[cid] = {"case_id": cid, "error": "skipped: cpu_cap", "cpu_s": None}
                    else:
                        acc += rows[cid]["cpu_s"]
            out[k] = rows
        return out


# ------------------------------------------------------------------------------------ the fit
from pymoo.core.problem import Problem  # noqa: E402


class RefitProblem(Problem):
    """Module level so the pickled algorithm (checkpoint) can be restored."""

    def __init__(self, vars):
        super().__init__(vars=vars, n_obj=4, n_ieq_constr=1)

    def _evaluate(self, X, out, *a, **kw):
        raise RuntimeError("evaluated via ask/tell")


class Fit:
    def __init__(self, pool, fold, out_dir, core_ids, screen_ids, evaluator, case_meta, settings=None,
                 test_ids=frozenset(), extra_meta=None):
        self.pool, self.fold = pool, fold
        self.out = Path(out_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.core = sorted(core_ids)
        self.screen = sorted(screen_ids)
        self.ev = evaluator
        self.meta = case_meta
        self.st = settings or Settings()
        self.space = get_space(pool)
        self.logfh = open(self.out / "progress.log", "a")
        # ---- leak guard: training ids only
        leak = set(self.core) & set(test_ids)
        if leak:
            raise LeakError(f"{len(leak)} test-fold cases in the training core of {pool} fold {fold}, "
                            f"e.g. {sorted(leak)[:3]}")
        if not set(self.screen) <= set(self.core):
            raise ValueError("screening subset not inside the core")
        if len(set(self.core)) != len(self.core):
            raise ValueError("duplicate case in core")
        self.extra_meta = extra_meta or {}
        self.summ = {}           # (scope, key) -> summary

    # ---------------- helpers
    def log(self, msg):
        log(f"[{self.pool} fold {self.fold}] {msg}", self.logfh)

    def score(self, specs: dict, scope, capped=False):
        """Evaluate {key: params} on scope 'screen' | 'core' -> {key: summary}. capped (search candidates):
        an individual is stopped at its first failed case or once its CPU exceeds cpu_cap_factor x the CPU
        of 3.7.0 accurate on the same cases (infeasible either way; beyond the hypervolume box and every
        budget)."""
        ids = self.screen if scope == "screen" else self.core
        need = {k: p for k, p in specs.items() if (scope, k) not in self.summ}
        if need:
            kw = {}
            if capped:
                ref = (self.refs_screen if scope == "screen" else self.refs_core)[BUDGET_REF[self.pool]["accurate"]]
                kw = {"cpu_cap": self.st.cpu_cap_factor * ref["cpu"], "stop_on_error": True}
            res = self.ev.evaluate(need, ids, **kw)
            for k, rows in res.items():
                self.summ[(scope, k)] = summarise(rows, self.meta)
        return {k: self.summ[(scope, k)] for k in specs}

    def write_run_json(self):
        sp = P370.load_space(self.pool)
        if C.sha256_file(sp.__file__) != P370.SPACE_SHA256[self.pool]:
            raise SystemExit(f"{sp.__file__}: search space differs from the frozen one (presets_370.SPACE_SHA256)")
        import pymoo
        if self.st.is_registered() and pymoo.__version__ != PYMOO_VERSION:
            raise SystemExit(f"registered fit needs pymoo {PYMOO_VERSION}, found {pymoo.__version__}")
        meta = {"pool": self.pool, "fold": self.fold, "settings": self.st.to_json(), "pymoo": pymoo.__version__,
                "registered": self.st.is_registered(),
                "n_core": len(self.core), "n_screen": len(self.screen),
                "core_sha256": ids_hash(self.core), "screen_sha256": ids_hash(self.screen),
                "space_file": sp.__file__, "space_sha256": C.sha256_file(sp.__file__),
                **self.extra_meta}
        path = self.out / "run.json"
        if path.exists():
            old = json.loads(path.read_text())
            amended = _v11_amendment(old.get("settings"), meta.get("settings"))
            for k in ("pool", "fold", "settings", "core_sha256", "screen_sha256", "space_sha256"):
                if k == "settings" and amended:
                    continue
                if old.get(k) != meta.get(k):
                    raise SystemExit(f"{path}: {k} differs from the recorded run; refusing to resume")
            meta["amendments"] = old.get("amendments", [])
            if amended:
                meta["amendments"] = meta["amendments"] + [
                    {"amendment": "v1.1", "change": f"stage2_gens {old['settings']['stage2_gens']} -> "
                     f"{meta['settings']['stage2_gens']}", "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
                     "previous_settings": old["settings"]}]
                self.log(f"Amendment v1.1: continuing the recorded run with stage-2 cap "
                         f"{old['settings']['stage2_gens']} -> {meta['settings']['stage2_gens']}")
            meta["started"] = old.get("started")
            meta["resumed"] = old.get("resumed", []) + [time.strftime("%Y-%m-%dT%H:%M:%S")]
        else:
            meta["started"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        _atomic_json(path, meta)
        (self.out / "cases_core.txt").write_text("".join(i + "\n" for i in self.core))
        (self.out / "cases_screen.txt").write_text("".join(i + "\n" for i in self.screen))

    # ---------------- references
    def references(self):
        self.refs = {m: _py(P370.params_370(self.pool, m)) for m in MODES[self.pool]}
        specs = {pkey(p): p for p in self.refs.values()}
        self.log(f"3.7.0 presets ({len(specs)} distinct) on screen ({len(self.screen)}) and core ({len(self.core)})")
        scr = self.score(specs, "screen")
        cor = self.score(specs, "core")
        self.refs_screen = {m: scr[pkey(p)] for m, p in self.refs.items()}
        self.refs_core = {m: cor[pkey(p)] for m, p in self.refs.items()}
        for m in MODES[self.pool]:
            for scope, s in (("screen", self.refs_screen[m]), ("core", self.refs_core[m])):
                if not s["valid"]:
                    raise RuntimeError(f"3.7.0 {m} failed on {s['n_err']} {scope} cases")
            s = self.refs_core[m]
            self.log(f"  3.7.0 {m:8s} core SP={s['sp']:.4f} TC={s['tc']:.4f} P={s['precision']:.4f} "
                     f"CPU={s['cpu']:.1f}s")
        _atomic_json(self.out / "refs.json", {
            "params": self.refs, "keys": {m: pkey(p) for m, p in self.refs.items()},
            "screen": self.refs_screen, "core": self.refs_core})

    # ---------------- one NSGA-III stage
    def run_stage(self, stage, vars_, decode, seeds):
        from pymoo.algorithms.moo.nsga3 import NSGA3
        from pymoo.core.evaluator import Evaluator as PEval
        from pymoo.core.mixed import (MixedVariableDuplicateElimination, MixedVariableMating,
                                      MixedVariableSampling)
        from pymoo.core.population import Population
        from pymoo.problems.static import StaticProblem
        from pymoo.termination import get_termination
        from pymoo.util.ref_dirs import get_reference_directions

        sd = self.out / f"stage{stage}"
        sd.mkdir(exist_ok=True)
        gmin, gmax = self.st.gens(stage)
        norm = Normaliser(self.refs_screen, stage, self.st.q_margin)
        problem = RefitProblem(vars_)
        alg_path, hv_path, ev_path = sd / "algorithm.pkl", sd / "hv.json", sd / "evaluated.jsonl"
        archive = self._load_archive(ev_path)
        hvh = json.loads(hv_path.read_text()) if hv_path.exists() else {"history": []}
        if alg_path.exists():
            algorithm = pickle.load(open(alg_path, "rb"))
            done_gen = algorithm.n_gen - 1
            hvh["history"] = hvh["history"][:done_gen]
            if len(hvh["history"]) < done_gen:     # crash between the hv write and the pickle
                hvh["history"] = [self._hv([e for e in archive.values() if e["gen"] <= g + 1], norm)
                                  for g in range(done_gen)]
            self.log(f"stage {stage}: resumed after generation {done_gen}")
        else:
            ref_dirs = get_reference_directions("das-dennis", 4, n_partitions=self.st.ref_partitions)
            seeds = seeds[:self.st.pop_size]
            n_rand = self.st.pop_size - len(seeds)
            init = list(MixedVariableSampling()(problem, n_rand, random_state=np.random.default_rng(self.st.seed))
                        .get("X")) if n_rand > 0 else []
            algorithm = NSGA3(ref_dirs=ref_dirs, pop_size=self.st.pop_size,
                              sampling=Population.new("X", np.array(seeds + init, dtype=object)),
                              mating=MixedVariableMating(eliminate_duplicates=MixedVariableDuplicateElimination()),
                              eliminate_duplicates=MixedVariableDuplicateElimination())
            algorithm.setup(problem, termination=get_termination("n_gen", gmax + 1000), seed=self.st.seed,
                            verbose=False)
            self.log(f"stage {stage}: NSGA-III {len(ref_dirs)} directions, pop {self.st.pop_size}, "
                     f"{len(seeds)} seeds, generations {gmin}-{gmax}")
        stop, reason = hv_stop(hvh["history"], gmin, gmax, self.st.hv_window, self.st.hv_tol)
        while not stop:
            gen = algorithm.n_gen if algorithm.n_gen is not None else 1
            tg = time.time()
            pop = algorithm.ask()
            X = pop.get("X")
            plist = [decode(x) for x in X]
            specs = {pkey(p): p for p in plist}
            summ = self.score(specs, "screen", capped=True)
            F, G = [], []
            with open(ev_path, "a") as fh:
                for x, p in zip(X, plist):
                    k = pkey(p)
                    s = summ[k]
                    f = norm.F(s) if s["valid"] else [2.0] * 4
                    F.append(f)
                    G.append([-1.0 if s["valid"] else 1.0 + s["n_err"]])
                    if k not in archive:
                        e = {"stage": stage, "gen": gen, "key": k, "params": p, "x": _py(dict(x)),
                             "screen": {kk: s[kk] for kk in ("sp", "tc", "precision", "f1", "cpu", "n_err",
                                                             "valid")}, "F": f}
                        archive[k] = e
                        fh.write(json.dumps(e, default=_json_default) + "\n")
            PEval().eval(StaticProblem(problem, F=np.array(F, float), G=np.array(G, float)), pop)
            algorithm.tell(infills=pop)
            hv = self._hv(archive.values(), norm)
            hvh["history"].append(hv)
            stop, reason = hv_stop(hvh["history"], gmin, gmax, self.st.hv_window, self.st.hv_tol)
            hvh.update({"stop": stop, "reason": reason, "normaliser": norm.to_json()})
            _atomic_json(hv_path, hvh)
            tmp = alg_path.with_suffix(".tmp")
            with open(tmp, "wb") as fh:
                pickle.dump(algorithm, fh)
            tmp.replace(alg_path)
            nv = sum(1 for p in plist if summ[pkey(p)]["valid"])
            self.log(f"stage {stage} gen {gen}: {len(plist)} individuals ({len(specs)} distinct, {nv} valid), "
                     f"archive {len(archive)}, HV {hv:.5f}, {time.time() - tg:.0f}s")
        self.log(f"stage {stage}: stopped after {len(hvh['history'])} generations: {reason}")
        hvh.update({"stop": True, "reason": reason, "n_gen": len(hvh["history"])})
        _atomic_json(hv_path, hvh)
        return archive, norm

    def _load_archive(self, path):
        arch = {}
        if path.exists():
            for line in open(path):
                try:
                    e = json.loads(line)
                except json.JSONDecodeError:
                    continue
                arch.setdefault(e["key"], e)
        return arch

    def _hv(self, entries, norm):
        F = [e["F"] for e in entries if e["screen"]["valid"]]
        return hypervolume(F, self.st.hv_ref)

    def front_of(self, archive):
        """Non-dominated valid archive entries on the screening objectives, CPU-filtered (step 3)."""
        ents = sorted((e for e in archive.values() if e["screen"]["valid"]), key=lambda e: e["key"])
        fr = [ents[i] for i in nondominated([raw_F(e["screen"]) for e in ents])]
        cut = self.st.rescore_margin * self.st.budget_factor * max(s["cpu"] for s in self.refs_screen.values())
        kept = [e for e in fr if e["screen"]["cpu"] <= cut]
        return kept, len(fr), cut

    def rescore(self, stage, archive):
        path = self.out / f"stage{stage}" / "front_core.json"
        fr, n_front, cut = self.front_of(archive)
        self.log(f"stage {stage}: front {n_front} on screen, {len(fr)} with screen CPU <= {cut:.1f}s -> core")
        core = self.score({e["key"]: e["params"] for e in fr}, "core", capped=True)
        out = [{"key": e["key"], "params": e["params"], "stage": stage, "screen": e["screen"],
                "core": core[e["key"]]} for e in fr]
        _atomic_json(path, {"n_front_screen": n_front, "cpu_cut": cut, "front": out})
        return out

    # ---------------- the whole fit
    def run(self):
        self.write_run_json()
        self.log(f"fit start: core {len(self.core)}, screen {len(self.screen)}, settings registered="
                 f"{self.st.is_registered()}")
        self.references()
        sp = self.space
        # ---- stage 1
        seeds1 = [sp.stage1_encode(p) for p in sp.stage1_seeds()]
        arch1, _ = self.run_stage(1, sp.stage1_vars(), sp.stage1_decode, seeds1)
        front1 = self.rescore(1, arch1)
        # ---- stage-2 pool: stage-1 front + 3.7.0 single runs (so the 3.7.0 ensembles are representable)
        members = {e["key"]: e["params"] for e in front1}
        for p in sp.stage1_seeds():
            members.setdefault(pkey(p), p)
        pool_path = self.out / "stage1" / "pool.json"
        pool_keys = sorted(members)
        if pool_path.exists():
            old = json.loads(pool_path.read_text())
            if old["keys"] != pool_keys:
                raise SystemExit(f"{pool_path} differs from the recomputed pool; refusing to resume")
        pool_params = [members[k] for k in pool_keys]
        zs = {}
        if front1:
            z = zmean([e["screen"] for e in front1])
            zs = {e["key"]: float(v) for e, v in zip(front1, z)}
        _atomic_json(pool_path, {"note": "stage-2 member pool: stage-1 front (CPU-filtered) + 3.7.0 single runs, "
                                         "ordered by key; member index i -> params[i]",
                                 "keys": pool_keys, "params": pool_params,
                                 "zmean_screen": [zs.get(k) for k in pool_keys]})
        # ---- stage 2
        seeds2 = [stage2_seed_x(sp, ens, pool_params) for ens in sp.stage2_seed_ensembles()]
        if zs:
            top = sorted(zs, key=lambda k: (-zs[k], k))
            top = (top * STAGE2_MAX)[:STAGE2_MAX]
            x = dict(seeds2[-1])
            for j in range(STAGE2_MAX):
                x[f"m_{j}"] = pool_keys.index(top[j])
            x["n_runs"] = 5
            seeds2.append(x)
        arch2, _ = self.run_stage(2, sp.stage2_vars(len(pool_params)),
                                  lambda x: sp.stage2_decode(x, pool_params), seeds2)
        front2 = self.rescore(2, arch2)
        # ---- selection on the core
        return self.select(front1, front2)

    def select(self, front1, front2):
        cands = []
        for e in front1 + front2:
            cands.append({"key": e["key"], "params": e["params"], "source": f"stage{e['stage']}",
                          **{k: e["core"][k] for k in ("sp", "tc", "precision", "f1", "cpu", "valid", "n_err")}})
        ref_keys = {}
        for m, p in self.refs.items():
            k = pkey(p)
            ref_keys.setdefault(k, []).append(m)
        for k, ms in ref_keys.items():
            s = self.refs_core[ms[0]]
            cands.append({"key": k, "params": self.refs[ms[0]], "source": "3.7.0:" + "/".join(ms),
                          **{kk: s[kk] for kk in ("sp", "tc", "precision", "f1", "cpu", "valid", "n_err")}})
        res = select_modes(cands, self.refs_core, self.pool, self.st.budget_factor)
        out = {"pool": self.pool, "fold": self.fold, "registered": self.st.is_registered(),
               "n_candidates": len(cands), "n_front_core": len(res["front"]), "budgets_cpu": res["budgets"],
               "budget_rule": f"{self.st.budget_factor} x CPU of the 3.7.0 preset of the mode on this core "
                              f"(parameter-set worker path)",
               "refs_core": {m: {k: s[k] for k in ("sp", "tc", "precision", "f1", "cpu")} | {"key": pkey(self.refs[m])}
                             for m, s in self.refs_core.items()},
               "selection": {}, "front": [{k: e[k] for k in ("key", "source", "sp", "tc", "precision", "f1", "cpu",
                                                              "zmean")} for e in res["front"]]}
        for m in MODES[self.pool]:
            r = res["selection"][m]
            e = r["entry"]
            ref_key = pkey(self.refs[m])
            if e is None:
                out["selection"][m] = {"key": None, "fallback_370": True, "params": self.refs[m], "is_370": True,
                                       "budget_cpu": r["budget_cpu"], "need_sp": r["need_sp"],
                                       "note": "no candidate within budget; 3.7.0 stays"}
                continue
            out["selection"][m] = {"key": e["key"], "params": e["params"], "source": e["source"],
                                   "is_370": e["key"] == ref_key, "fallback_370": False,
                                   "zmean": e["zmean"], "budget_cpu": r["budget_cpu"], "need_sp": r["need_sp"],
                                   "core": {k: e[k] for k in ("sp", "tc", "precision", "f1", "cpu")},
                                   "short": self.space.short(e["params"])}
            self.log(f"select {m:8s}: {e['source']} SP={e['sp']:.4f} TC={e['tc']:.4f} P={e['precision']:.4f} "
                     f"CPU={e['cpu']:.1f}s (budget {r['budget_cpu']:.1f}s) z={e['zmean']:.3f}"
                     f"{' = 3.7.0' if e['key'] == ref_key else ''}")
        _atomic_json(self.out / "selection.json", out)
        self.log("fit done")
        return out


def _v11_amendment(old, new):
    """Amendment v1.1 (PLAN): a run recorded with the v1 settings may be continued with the v1.1 settings, which
    differ ONLY in the stage-2 generation cap (10, 40) -> (10, 140). Anything else must match exactly."""
    if not old or not new or old == new:
        return False
    o, n = dict(old), dict(new)
    if tuple(o.pop("stage2_gens", ())) != (10, 40) or tuple(n.pop("stage2_gens", ())) != (10, 140):
        return False
    return o == n


# ------------------------------------------------------------------------------------ data access
def default_ws() -> Path:
    return Path(os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit")).expanduser()


def read_ids(path):
    return [l.strip() for l in open(path) if l.strip()]


def load_fold_inputs(pool, fold, ws):
    """-> (case_rows {id: manifest row}, meta {id: {...}}, core ids, screen ids, test ids)."""
    from .folds import read_folds
    from .scoring import read_manifest
    ws = Path(ws)
    rows = {r["case_id"]: r for r in read_manifest(ws / "data" / f"cases_{pool}.final.tsv")}
    folds = read_folds(pool, ws)
    meta = {cid: {"benchmark": rows[cid]["benchmark"], "fold": int(f["fold"]),
                  "homolog_in_train": int(f["homolog_in_train"]), "difficulty": f["difficulty"]}
            for cid, f in folds.items()}
    if set(meta) != set(rows):
        raise SystemExit(f"folds_{pool}.tsv and cases_{pool}.final.tsv disagree")
    core = read_ids(ws / "data" / "cores" / f"{pool}_fold{fold}_core.txt")
    screen = read_ids(ws / "data" / "cores" / f"{pool}_fold{fold}_screen.txt")
    test = frozenset() if str(fold) == "all" else frozenset(c for c, m in meta.items() if m["fold"] == int(fold))
    return rows, meta, core, screen, test


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--fold", required=True, help="0..4 or all")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ws", default=str(default_ws()))
    ap.add_argument("--build-dir", help="aligner build dir (default WS/kalign/960a4f6)")
    ap.add_argument("--cache", help="cache path (default WS/cache/cache.sqlite)")
    ap.add_argument("--slots-dir", help="machine-wide slot locks (default WS/runs/slots)")
    ap.add_argument("--threads", type=int, default=MAX_SLOTS, help="max concurrent aligners of this process")
    ap.add_argument("--fake-eval", action="store_true", help="pseudo-results, no aligner (machinery test)")
    # smoke-test overrides (recorded; nested/pool/decide refuse non-registered settings without --allow-smoke)
    ap.add_argument("--pop-size", type=int, default=POP_SIZE)
    ap.add_argument("--ref-partitions", type=int, default=REF_PARTITIONS)
    ap.add_argument("--stage1-gens", default=f"{STAGE_GENS[1][0]},{STAGE_GENS[1][1]}")
    ap.add_argument("--stage2-gens", default=f"{STAGE_GENS[2][0]},{STAGE_GENS[2][1]}")
    ap.add_argument("--limit-core", type=int, default=0, help="smoke: first N core ids (screen restricted)")
    ap.add_argument("--timeout-s", type=float, default=FIT_TIMEOUT_S)
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    fold = a.fold if a.fold == "all" else int(a.fold)
    rows, meta, core, screen, test = load_fold_inputs(a.pool, fold, ws)
    note = ""
    if a.limit_core:
        core = sorted(core)[:a.limit_core]
        screen = [c for c in screen if c in set(core)] or core[: max(1, a.limit_core // 2)]
        note = f"limit_core={a.limit_core}"
    st = Settings(pop_size=a.pop_size, ref_partitions=a.ref_partitions,
                  stage1_gens=tuple(int(v) for v in a.stage1_gens.split(",")),
                  stage2_gens=tuple(int(v) for v in a.stage2_gens.split(",")),
                  timeout_s=a.timeout_s, fake_eval=a.fake_eval, note=note)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    extra = {}
    if a.fake_eval:
        ev = FakeEvaluator(rows.values(), core)
    else:
        from .evaluate import Build
        build = Build.load(a.build_dir or ws / "kalign" / "960a4f6")
        slots = SlotPool(a.slots_dir or ws / "runs" / "slots", MAX_SLOTS)
        ev = Evaluator(build, rows.values(), core, cache_path=a.cache or ws / "cache" / "cache.sqlite",
                       slots=slots, threads=a.threads, timeout_s=st.timeout_s)
        extra = {"build_id": build.build_id, "cache": str(ev.cache.path), "threads": a.threads,
                 "slots_dir": str(slots.dir)}
    f = Fit(a.pool, fold, out, core, screen, ev, meta, st, test_ids=test, extra_meta=extra)
    if not a.fake_eval:
        ev.logfh = f.logfh
    f.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
