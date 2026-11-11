"""One fit of the ensemble plan (PLAN_ensemble.md, frozen at tag ens-plan-freeze): Kalign 4 protein accurate mode,
5 members, member 1 = the fast setting (fast-fit pick 0976900c), members 2..5 and the ensemble options searched;
3 objectives (small precision, precision at 1,000, CPU relative to the fast setting), NSGA-III, re-scoring, knee-rule
selection with precision. Reuses the fast fit's machinery (recipe.fastfit / recipe.fit); every choice the plan leaves
open is in docs/IMPL_NOTES_ensfit.md.

    python -m recipe.ensfit --fold 0 --out WS/runs/ensfit/fold0 [--ws WS]           (0..4 or all)

Needs (WS = ~/tmp/kalign-fastfit):
  WS/data/...                       the fast fit's data, folds, cores and screening lists (recipe.fastfit.load_fold)
  WS/kalign/5b93950/build.json      the engine of the plan (section 1)
  WS/runs/fastfit/final_pick.json   the fast setting (key 0976900c) -> member 1 and the reference
  PRESETS_V2                        the pre-release accurate preset (key 2852d1c2) -> seeds and a reported row
  WS/runs/ensfit/select10k.json     per fold: the 24 selection families at 10,000 and the timing-guard family
                                    (recipe.ensfit_prep, written before launch)
Cache: WS/cache/ensfit.sqlite (separate from the fast fit's cache; keys also carry the build id).

Procedure (plan sections 2-6):
  1. the fast setting on screen and core (reference of every objective and constraint);
  2. NSGA-III, population 56, Das-Dennis 9 partitions, seed 1234; 2 seeds + 54 random; screening data
     (small screening subset + 40 training families at 1,000 + the fold's guard family at 10,000);
     stop >= 20 generations and HV gain < 0.1 % over 5, or 60;
  3. re-scoring of <= 40 front points on the core (small core + every training family at 1,000 + 24 at 10,000);
  4. the pre-release accurate preset on screen and core (reported row);
  5. selection on the core: knee rule of DECISION_fastfit_selection.md with precision (plan section 6).

Output (out_dir): run.json, refs.json, algorithm.pkl, hv.json, evaluated.jsonl, front_core.json, selection.json,
progress.log. Resumable (checkpoint per generation; every alignment is cached).
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import json
import math
import os
import pickle
import sys
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from . import cache as C
from . import fit as FIT
from .fit import (Evaluator, LeakError, SlotPool, _atomic_json, _json_default, _py, hv_stop, hypervolume, log,
                  nondominated, pkey)

SLOTS = 60                      # machine-wide aligner limit (as the fast fit)
FIT.MAX_SLOTS = SLOTS

PLAN = "PLAN_ensemble.md"
PLAN_SHA256 = "77d2f1577f8a092f75a4494ec1aa2c2e911e5e838ea44db4e30cf69030619a23"
SEED = 1234
POP_SIZE = 56
REF_PARTITIONS = 9
GENS = (20, 30)                 # AMENDMENT_ensemble_1.md (c): cap 30 (plan: 60)
AMENDMENT = "AMENDMENT_ensemble_1.md + AMENDMENT_ensemble_2.md"
HV_WINDOW, HV_TOL, HV_REF = 5, 0.001, 1.1
P_MARGIN = 0.10                 # precision axes normalised from (P_ref - 0.10)
R_LO, R_HI = 30.0, 50.0         # AMENDMENT_ensemble_2.md (c): observed R range of generations 1-2 (plan: 0.25-3)
CPU_SMALL_MAX, CPU_1K_MAX = 7.1, 53.6
SP_FLOOR = 0.005                # recall floor: SP >= fast setting - 0.005 (small and 1,000)
GUARD_CPU_MAX = 600.0           # 10,000 timing guard
TIMEOUT_S = 1800.0
RESCORE_MAX, RESCORE_SMALL_MARGIN = 40, 0.01
N_SELECT_10K = 24
N_BOOT, BOOT_SEED, Z = 2000, 1234, 1.96
SENS_DELTAS = (0.005, 0.01)
PYMOO_VERSION = "0.6.1.6"
BUILD = "5b93950"
K = 5

FAST_KEY = "0976900cd9cc212935e8c848cbb3a377e77ec8ee1578efb35b64f49193e5dae0"
ACC_KEY = "2852d1c29232d035cf5f065e560b338321afd3d9de459ab9fa4e62c0d8e9ae32"
PRESETS_V2 = Path("~/tmp/kalign-rel/5b93950/src/presets/results/v2/presets_v2.json").expanduser()

# aligner constants (kalign.h)
MATRIX_MAP = [1, 2, 3]                  # PFASUM43, PFASUM60, CorBLOSUM66
MATRIX_NAMES = {1: "P43", 2: "P60", 3: "CB66"}
REFINE_MAP = [0, 1, 2, 3]               # none, all, confident, inline
CONSISTENCY_MAP = [0, 1, 2, 3, 5, 8, 10]
FILL_MAP = [0, 1, 2, 3]                 # off, realign, project, attach
TREE_METHODS = (0, 3)                   # k-means (+ WPGMA leaves), hybrid (k-means + Prim leaves)
TREE_DISTS = (0, 1, 2, 3)               # Myers raw, Myers normalised, LCS indel, 1 - LCS / shorter
LEAF_LOG = (math.log10(50), math.log10(5000))
MAXLEN = 1024
REALIGN_TREE = 1                        # not searched (plan section 2); engine default = the fast setting's value
GAP_MODEL = 1                           # member gap model in every member (plan section 1)
GUIDE_DIST = 0
TREE_400 = {"method": 0, "dist": 0, "alpha": 0, "leaf": 50, "kd": 0}   # the shipped tree as explicit values
CANON = {"consistency_weight": 1.0, "vsm_kappa": 1.0, "vsm_dzero": 1.0, "vsm_amax": 0.0,
         "consistency_merge_weight": 2.0, "fill_anchor_frac": 0.0, "fill_support": 1, "min_support": 0}

# per-member worker lists (one entry per member); the scalars of the old layout are member 1's values (overridden
# by the lists in _core, present because build_kwargs reads them)
RUN_LISTS = ("run_gpo", "run_gpe", "run_tgpe", "run_noise", "run_types", "run_vsm_amax", "run_refine",
             "run_vsm_kappa", "run_vsm_dzero", "run_guide_dist", "run_tree_method", "run_tree_dist",
             "run_tree_alpha", "run_tree_maxlen", "run_tree_leaf", "run_tree_kmeans_dist", "run_realign_tree",
             "run_gap_model", "run_seq_weights", "run_realign", "run_consistency_anchors", "run_consistency_weight",
             "run_vsm_mode")
SCALAR_OF = {"run_seq_weights": "seq_weights", "run_realign": "realign", "run_consistency_anchors": "consistency",
             "run_consistency_weight": "consistency_weight", "run_vsm_mode": "vsm_mode"}
ENS_KEYS = ("min_support", "consistency_merge", "consistency_merge_weight", "fill_in", "fill_anchor_frac",
            "fill_support")


# ------------------------------------------------------------------------------------ the two fixed settings
def sha256_file(p):
    return C.sha256_file(p)


def fast_setting(ws):
    """The fast setting (fast-fit final pick 0976900c) as a 1-member worker dict + the member gap model."""
    d = json.loads((Path(ws) / "runs" / "fastfit" / "final_pick.json").read_text())
    if d["key"] != FAST_KEY or pkey(d["params"]) != FAST_KEY:
        raise SystemExit(f"final_pick.json: key {d['key'][:12]} / {pkey(d['params'])[:12]} != {FAST_KEY[:12]}")
    p = copy.deepcopy(d["params"])
    p["run_gap_model"] = [GAP_MODEL]
    return _py(p)


def prerelease_accurate(path=PRESETS_V2):
    """The pre-release accurate preset 2852d1c2 exactly as the release header carries it (gap model 0, shipped tree,
    scalar member settings)."""
    d = json.loads(Path(path).read_text())["presets"]["protein"]["accurate"]
    if d["key"] != ACC_KEY:
        raise SystemExit(f"{path}: protein accurate key {d['key'][:12]} != {ACC_KEY[:12]}")
    return _py(copy.deepcopy(d["params"]))


def fast_member(fast):
    """Member dict of the fast setting (the fixed member 1)."""
    return {"gpo": fast["run_gpo"][0], "gpe": fast["run_gpe"][0], "tgpe": fast["run_tgpe"][0],
            "matrix": fast["run_types"][0], "vsm_amax": fast["run_vsm_amax"][0],
            "vsm_kappa": fast["run_vsm_kappa"][0], "vsm_dzero": fast["run_vsm_dzero"][0],
            "vsm_mode": fast["vsm_mode"], "seq_weights": fast["seq_weights"], "realign": fast["realign"],
            "consistency": fast["consistency"], "consistency_weight": fast["consistency_weight"],
            "refine": fast["run_refine"][0], "noise": fast["run_noise"][0],
            "tree_method": fast["run_tree_method"][0], "tree_leaf": fast["run_tree_leaf"][0],
            "tree_dist": fast["run_tree_dist"][0], "tree_alpha": fast["run_tree_alpha"][0],
            "tree_kmeans_dist": fast["run_tree_kmeans_dist"][0], "tree_maxlen": fast["run_tree_maxlen"][0],
            "realign_tree": fast["run_realign_tree"][0], "gap_model": fast["run_gap_model"][0],
            "guide_dist": fast["run_guide_dist"][0]}


def acc_members(acc, tree):
    """Members 2..5 of the pre-release accurate preset with the member gap model and the given (explicit) tree."""
    out = []
    for k in range(1, K):
        out.append({"gpo": acc["run_gpo"][k], "gpe": acc["run_gpe"][k], "tgpe": acc["run_tgpe"][k],
                    "matrix": acc["run_types"][k], "vsm_amax": acc["run_vsm_amax"][k],
                    "vsm_kappa": acc["run_vsm_kappa"][k], "vsm_dzero": acc["run_vsm_dzero"][k],
                    "vsm_mode": acc["vsm_mode"], "seq_weights": acc["seq_weights"], "realign": acc["realign"],
                    "consistency": acc["consistency"], "consistency_weight": acc["consistency_weight"],
                    "refine": acc["run_refine"][k], "noise": acc["run_noise"][k],
                    "tree_method": tree["method"], "tree_leaf": tree["leaf"], "tree_dist": tree["dist"],
                    "tree_alpha": tree["alpha"], "tree_kmeans_dist": tree["kd"]})
    return out


# ------------------------------------------------------------------------------------ search space
def mask_member(m):
    """Canonical values for member parameters that cannot change the output (one cache key per distinct output)."""
    m = dict(m)
    if int(m["consistency"]) == 0:
        m["consistency_weight"] = CANON["consistency_weight"]
    if int(m["vsm_mode"]) == 0:
        m["vsm_kappa"], m["vsm_dzero"] = CANON["vsm_kappa"], CANON["vsm_dzero"]
    else:
        m["vsm_amax"] = CANON["vsm_amax"]
    if int(m["tree_method"]) == 0 and int(m["tree_kmeans_dist"]) == 0:
        m["tree_dist"], m["tree_alpha"] = 0, 0
    m.update({"tree_maxlen": MAXLEN, "realign_tree": REALIGN_TREE, "gap_model": GAP_MODEL, "guide_dist": GUIDE_DIST})
    return m


def mask_ensemble(e):
    e = dict(e)
    if int(e["consistency_merge"]) == 1:
        # the consistency-merge path never reads min_support or the fill-in options (ensemble.c)
        e.update({"min_support": CANON["min_support"], "fill_in": 0,
                  "fill_anchor_frac": CANON["fill_anchor_frac"], "fill_support": CANON["fill_support"]})
    else:
        e["consistency_merge_weight"] = CANON["consistency_merge_weight"]
    if int(e["fill_in"]) not in (1, 2):
        e["fill_anchor_frac"] = CANON["fill_anchor_frac"]
    if int(e["fill_in"]) != 3:
        e["fill_support"] = CANON["fill_support"]
    else:
        e["fill_support"] = max(1, min(int(e["fill_support"]), K))
    e["min_support"] = min(int(e["min_support"]), K)
    return e


def worker_dict(members, ens):
    """5 member dicts (member 1 first) + ensemble options -> worker dict (recipe/_align_worker.py layout)."""
    ms = [mask_member(m) for m in members]
    e = mask_ensemble(ens)
    col = {"run_gpo": ("gpo", float), "run_gpe": ("gpe", float), "run_tgpe": ("tgpe", float),
           "run_noise": ("noise", float), "run_types": ("matrix", int), "run_vsm_amax": ("vsm_amax", float),
           "run_refine": ("refine", int), "run_vsm_kappa": ("vsm_kappa", float),
           "run_vsm_dzero": ("vsm_dzero", float), "run_guide_dist": ("guide_dist", int),
           "run_tree_method": ("tree_method", int), "run_tree_dist": ("tree_dist", int),
           "run_tree_alpha": ("tree_alpha", int), "run_tree_maxlen": ("tree_maxlen", int),
           "run_tree_leaf": ("tree_leaf", int), "run_tree_kmeans_dist": ("tree_kmeans_dist", int),
           "run_realign_tree": ("realign_tree", int), "run_gap_model": ("gap_model", int),
           "run_seq_weights": ("seq_weights", float), "run_realign": ("realign", int),
           "run_consistency_anchors": ("consistency", int), "run_consistency_weight": ("consistency_weight", float),
           "run_vsm_mode": ("vsm_mode", int)}
    p = {"n_runs": len(ms)}
    for k, (name, f) in col.items():
        p[k] = [f(m[name]) for m in ms]
    for rk, sk in SCALAR_OF.items():
        p[sk] = p[rk][0]
    p.update({"min_support": int(e["min_support"]), "consistency_merge": int(e["consistency_merge"]),
              "consistency_merge_weight": float(e["consistency_merge_weight"]), "fill_in": int(e["fill_in"]),
              "fill_anchor_frac": float(e["fill_anchor_frac"]), "fill_support": int(e["fill_support"])})
    return _py(p)


def members_of(p):
    """Inverse of worker_dict: -> (member dicts, ensemble options)."""
    names = {"run_gpo": "gpo", "run_gpe": "gpe", "run_tgpe": "tgpe", "run_noise": "noise", "run_types": "matrix",
             "run_vsm_amax": "vsm_amax", "run_refine": "refine", "run_vsm_kappa": "vsm_kappa",
             "run_vsm_dzero": "vsm_dzero", "run_guide_dist": "guide_dist", "run_tree_method": "tree_method",
             "run_tree_dist": "tree_dist", "run_tree_alpha": "tree_alpha", "run_tree_maxlen": "tree_maxlen",
             "run_tree_leaf": "tree_leaf", "run_tree_kmeans_dist": "tree_kmeans_dist",
             "run_realign_tree": "realign_tree", "run_gap_model": "gap_model", "run_seq_weights": "seq_weights",
             "run_realign": "realign", "run_consistency_anchors": "consistency",
             "run_consistency_weight": "consistency_weight", "run_vsm_mode": "vsm_mode"}
    ms = [{v: p[k][j] for k, v in names.items()} for j in range(p["n_runs"])]
    return ms, {k: p[k] for k in ENS_KEYS}


class Space:
    """Members 2..5 (each searched on its own) and the ensemble options; member 1 = the fast setting."""

    MEMBER_VARS = ("gpo", "gpe", "tgpe", "matrix", "vsm_amax", "vsm_kappa", "vsm_dzero", "vsm_mode", "seq_weights",
                   "realign", "consistency", "consistency_weight", "refine", "noise", "tree_method", "tree_leaf_log",
                   "tree_dist", "tree_alpha", "tree_kmeans_dist")

    def __init__(self, fast, acc):
        self.fast, self.acc = fast, acc
        self.m1 = fast_member(fast)

    def vars(self):
        from pymoo.core.variable import Choice, Integer, Real
        v = {}
        for j in range(2, K + 1):
            v.update({f"gpo_{j}": Real(bounds=(2.0, 15.0)), f"gpe_{j}": Real(bounds=(0.5, 5.0)),
                      f"tgpe_{j}": Real(bounds=(0.1, 3.0)), f"matrix_{j}": Choice(options=list(range(3))),
                      f"vsm_amax_{j}": Real(bounds=(0.0, 5.0)), f"vsm_kappa_{j}": Real(bounds=(0.0, 1.0)),
                      f"vsm_dzero_{j}": Real(bounds=(0.5, 10.0)), f"vsm_mode_{j}": Choice(options=[0, 1]),
                      f"seq_weights_{j}": Real(bounds=(0.0, 5.0)), f"realign_{j}": Integer(bounds=(0, 2)),
                      f"consistency_{j}": Choice(options=list(CONSISTENCY_MAP)),
                      f"consistency_weight_{j}": Real(bounds=(0.5, 5.0)),
                      f"refine_{j}": Choice(options=list(REFINE_MAP)), f"noise_{j}": Real(bounds=(0.0, 0.5)),
                      f"tree_method_{j}": Choice(options=list(TREE_METHODS)), f"tree_leaf_log_{j}": Real(bounds=LEAF_LOG),
                      f"tree_dist_{j}": Choice(options=list(TREE_DISTS)), f"tree_alpha_{j}": Choice(options=[0, 1]),
                      f"tree_kmeans_dist_{j}": Choice(options=[0, 1])})
        v.update({"min_support": Integer(bounds=(0, K)), "consistency_merge": Choice(options=[0, 1]),
                  "consistency_merge_weight": Real(bounds=(0.5, 10.0)), "fill_in": Choice(options=list(FILL_MAP)),
                  "fill_anchor_frac": Real(bounds=(0.0, 1.0)), "fill_support": Integer(bounds=(1, K))})
        return v

    def decode(self, x):
        ms = [self.m1]
        for j in range(2, K + 1):
            ms.append({"gpo": float(x[f"gpo_{j}"]), "gpe": float(x[f"gpe_{j}"]), "tgpe": float(x[f"tgpe_{j}"]),
                       "matrix": MATRIX_MAP[int(x[f"matrix_{j}"])], "vsm_amax": float(x[f"vsm_amax_{j}"]),
                       "vsm_kappa": float(x[f"vsm_kappa_{j}"]), "vsm_dzero": float(x[f"vsm_dzero_{j}"]),
                       "vsm_mode": int(x[f"vsm_mode_{j}"]), "seq_weights": float(x[f"seq_weights_{j}"]),
                       "realign": int(x[f"realign_{j}"]), "consistency": int(x[f"consistency_{j}"]),
                       "consistency_weight": float(x[f"consistency_weight_{j}"]), "refine": int(x[f"refine_{j}"]),
                       "noise": float(x[f"noise_{j}"]), "tree_method": int(x[f"tree_method_{j}"]),
                       "tree_leaf": int(round(10 ** float(x[f"tree_leaf_log_{j}"]))),
                       "tree_dist": int(x[f"tree_dist_{j}"]), "tree_alpha": int(x[f"tree_alpha_{j}"]),
                       "tree_kmeans_dist": int(x[f"tree_kmeans_dist_{j}"])})
        ens = {k: x[k] for k in ENS_KEYS}
        p = worker_dict(ms, ens)
        if p["run_gpo"][0] != self.fast["run_gpo"][0]:
            raise RuntimeError("member 1 is not the fast setting")
        return p

    def encode(self, p):
        ms, ens = members_of(p)
        x = {}
        for j in range(2, K + 1):
            m = ms[j - 1]
            x.update({f"gpo_{j}": m["gpo"], f"gpe_{j}": m["gpe"], f"tgpe_{j}": m["tgpe"],
                      f"matrix_{j}": MATRIX_MAP.index(m["matrix"]), f"vsm_amax_{j}": m["vsm_amax"],
                      f"vsm_kappa_{j}": m["vsm_kappa"], f"vsm_dzero_{j}": m["vsm_dzero"], f"vsm_mode_{j}": m["vsm_mode"],
                      f"seq_weights_{j}": m["seq_weights"], f"realign_{j}": m["realign"],
                      f"consistency_{j}": m["consistency"], f"consistency_weight_{j}": m["consistency_weight"],
                      f"refine_{j}": m["refine"], f"noise_{j}": m["noise"], f"tree_method_{j}": m["tree_method"],
                      f"tree_leaf_log_{j}": math.log10(m["tree_leaf"]), f"tree_dist_{j}": m["tree_dist"],
                      f"tree_alpha_{j}": m["tree_alpha"], f"tree_kmeans_dist_{j}": m["tree_kmeans_dist"]})
        x.update(ens)
        return x

    def seed_params(self):
        """Plan section 5: (1) the pre-release accurate members 2..5 and ensemble options with the fast setting as
        member 1 (shipped tree as explicit values, member gap model); (2) the same with the fast setting's guide tree
        in every member."""
        ens = {k: self.acc[k] for k in ENS_KEYS}
        ft = {"method": self.m1["tree_method"], "dist": self.m1["tree_dist"], "alpha": self.m1["tree_alpha"],
              "leaf": self.m1["tree_leaf"], "kd": self.m1["tree_kmeans_dist"]}
        return [worker_dict([self.m1] + acc_members(self.acc, TREE_400), ens),
                worker_dict([self.m1] + acc_members(self.acc, ft), ens)]

    @staticmethod
    def short(p):
        ms, e = members_of(p)
        tn = {0: "km", 3: "km+Prim"}
        parts = []
        for j, m in enumerate(ms[1:], 2):
            vsm = f"vsm0 a={m['vsm_amax']:.2f}" if m["vsm_mode"] == 0 else \
                f"vsm1 k={m['vsm_kappa']:.2f} d0={m['vsm_dzero']:.1f}"
            parts.append(f"m{j}:{MATRIX_NAMES[m['matrix']]} {m['gpo']:.2f}/{m['gpe']:.2f}/{m['tgpe']:.2f} {vsm} "
                         f"sw={m['seq_weights']:.2f} re={m['realign']} c={m['consistency']} ref={m['refine']} "
                         f"nz={m['noise']:.2f} {tn[m['tree_method']]} leaf={m['tree_leaf']} d={m['tree_dist']} "
                         f"a={m['tree_alpha']} kd={m['tree_kmeans_dist']}")
        ens = (f"ms={e['min_support']} cm={e['consistency_merge']}"
               f"{'/' + format(e['consistency_merge_weight'], '.2f') if e['consistency_merge'] else ''} "
               f"fill={e['fill_in']} frac={e['fill_anchor_frac']:.3f} fs={e['fill_support']}")
        return ens + " | " + " | ".join(parts)


# ------------------------------------------------------------------------------------ objectives
def summarise(rows, meta, guard):
    """rows {case_id: cache row} -> summary. Small: per-case means within benchmark, benchmarks weighted equally;
    1,000 and 10,000: per-family means. CPU per part; the guard row separately."""
    small, large = {}, {1000: [], 10000: []}
    cpu = {"small": 0.0, 1000: 0.0, 10000: 0.0}
    n_err, n_err_part = 0, {"small": 0, 1000: 0, 10000: 0}
    g = rows.get(guard) if guard else None
    for cid, r in rows.items():
        m = meta[cid]
        part = "small" if m["kind"] == "small" else m["k"]
        cpu[part] += float(r.get("cpu_s") or 0.0)
        if r.get("error"):
            n_err += 1
            n_err_part[part] += 1
            continue
        v = (float(r["precision"]), float(r["recall"]), float(r["tc"]), float(r["width_ratio"]))
        if part == "small":
            small.setdefault(m["benchmark"], []).append(v)
        else:
            large[part].append(v)

    def mean(vs, i):
        return float(np.mean([v[i] for v in vs])) if vs else float("nan")
    per_b = {b: {"p": mean(vs, 0), "sp": mean(vs, 1), "tc": mean(vs, 2), "wr": mean(vs, 3), "n": len(vs)}
             for b, vs in sorted(small.items())}

    def bw(name):
        return float(np.mean([v[name] for v in per_b.values()])) if per_b else float("nan")
    s = {"p_small": bw("p"), "sp_small": bw("sp"), "tc_small": bw("tc"), "wr_small": bw("wr"),
         "p_1k": mean(large[1000], 0), "sp_1k": mean(large[1000], 1), "tc_1k": mean(large[1000], 2),
         "wr_1k": mean(large[1000], 3), "n_1k": len(large[1000]),
         "p_10k": mean(large[10000], 0), "sp_10k": mean(large[10000], 1), "tc_10k": mean(large[10000], 2),
         "wr_10k": mean(large[10000], 3), "n_10k": len(large[10000]),
         "cpu_small": cpu["small"], "cpu_1k": cpu[1000], "cpu_10k": cpu[10000],
         "guard": guard, "guard_cpu": float(g.get("cpu_s") or 0.0) if g else None,
         "guard_err": (g.get("error") if g else "missing"),
         "n": len(rows), "n_err": n_err, "n_err_small": n_err_part["small"], "n_err_1k": n_err_part[1000],
         "n_err_10k": n_err_part[10000], "valid": n_err == 0 and len(rows) > 0, "per_benchmark": per_b}
    return s


def objectives(s, ref):
    """CPU ratios, R and the constraints of plan section 4 against the fast setting on the same cases."""
    rs = s["cpu_small"] / max(ref["cpu_small"], 1e-9)
    r1 = s["cpu_1k"] / max(ref["cpu_1k"], 1e-9)
    out = dict(s)
    floor_small = ref["sp_small"] - SP_FLOOR
    floor_1k = ref["sp_1k"] - SP_FLOOR
    guard_ok = s["guard_err"] is None and s["guard_cpu"] is not None and s["guard_cpu"] <= GUARD_CPU_MAX
    out.update({"cpu_ratio_small": rs, "cpu_ratio_1k": r1, "R": math.sqrt(max(rs, 1e-12) * max(r1, 1e-12)),
                "sp_floor_small": floor_small, "sp_floor_1k": floor_1k, "guard_ok": guard_ok})
    # AMENDMENT_ensemble_1.md (a): the guard constraint is removed (guard_ok is only reported)
    # AMENDMENT_ensemble_2.md (a): no CPU constraint (R stays objective 3)
    out["feasible"] = bool(s["valid"]
                           and s["sp_small"] >= floor_small and s["sp_1k"] >= floor_1k)
    return out


class Normaliser:
    """PLAN_fastfit section 5 transform with precision in place of Q and the fast setting as the reference."""

    def __init__(self, ref):
        self.lo_s = ref["p_small"] - P_MARGIN
        self.lo_l = ref["p_1k"] - P_MARGIN

    def F(self, s):
        return [(1.0 - s["p_small"]) / (1.0 - self.lo_s), (1.0 - s["p_1k"]) / (1.0 - self.lo_l),
                (math.log(max(s["R"], 1e-9)) - math.log(R_LO)) / (math.log(R_HI) - math.log(R_LO))]

    def to_json(self):
        return {"p_small_lo": self.lo_s, "p_1k_lo": self.lo_l, "R_lo": R_LO, "R_hi": R_HI}


def raw_F(s):
    return [-s["p_small"], -s["p_1k"], s["R"]]


def G(s):
    """Inequality constraints (<= 0 feasible): recall floors, CPU limits, failed cases (no guard: amendment 1)."""
    if not s["valid"]:
        return [1.0] * 2 + [1.0 + s["n_err"]]
    return [s["sp_floor_small"] - s["sp_small"], s["sp_floor_1k"] - s["sp_1k"], -1.0]   # amendment 2: no CPU limits


def rescore_set(front, ref_screen):
    """PLAN_fastfit section 5 with precision: front points with P_small >= P_small,fast - 0.01; at most 40: the 20 with
    the highest P_1k, then the 20 with the lowest R among the rest, then by P_1k."""
    ok = [e for e in front if e["screen"]["p_small"] >= ref_screen["p_small"] - RESCORE_SMALL_MARGIN]
    if len(ok) <= RESCORE_MAX:
        return sorted(ok, key=lambda e: e["key"]), len(ok)
    by_l = sorted(ok, key=lambda e: (-e["screen"]["p_1k"], e["key"]))
    pick = by_l[:20]
    rest = [e for e in ok if e not in pick]
    pick += sorted(rest, key=lambda e: (e["screen"]["R"], e["key"]))[:20]
    rest = [e for e in by_l if e not in pick]
    pick += rest[:RESCORE_MAX - len(pick)]
    return sorted(pick, key=lambda e: e["key"]), len(ok)


def eligible(c, ref):
    """Plan section 6 step 1 on the core summaries (c and ref: core summaries with objectives)."""
    s = c["core"]
    return bool(s["feasible"] and s["p_small"] > ref["p_small"] and s["sp_small"] >= ref["sp_small"] - SP_FLOOR
                and s["sp_1k"] >= ref["sp_1k"] - SP_FLOOR and s["n_err_10k"] == 0
                and s["n_10k"] == ref["n_10k"] and s["p_10k"] >= ref["p_10k"]
                and s["sp_10k"] >= ref["sp_10k"] - SP_FLOOR)


def knee(cands, ref, fam_prec):
    """Plan section 6 / DECISION_fastfit_selection.md with precision. cands: [{key, params, core, ...}] (core
    summaries), ref: the fast setting's core summary, fam_prec: {key: np.array of per-family precision at 1,000 in one
    fixed family order}. -> dict."""
    E = sorted((c for c in cands if eligible(c, ref)), key=lambda c: c["key"])
    out = {"n_candidates": len(cands), "n_eligible": len(E), "eligible": [c["key"] for c in E]}
    if not E:
        out.update({"pick": None, "b": None, "within": [], "rows": []})
        for d in SENS_DELTAS:
            out[f"pick_d{d}"] = None
        return out
    b = sorted(E, key=lambda c: (-c["core"]["p_1k"], c["core"]["R"], c["key"]))[0]
    vb = fam_prec[b["key"]]
    if abs(float(vb.mean()) - b["core"]["p_1k"]) > 1e-9:
        raise RuntimeError(f"recomputed P_1k {vb.mean():.6f} != core summary {b['core']['p_1k']:.6f}")
    rng = np.random.default_rng(BOOT_SEED)
    idx = rng.integers(0, len(vb), size=(N_BOOT, len(vb)))
    bb = vb[idx].mean(axis=1)
    rows = []
    for c in E:
        if c["key"] == b["key"]:
            rows.append({"key": c["key"], "diff": 0.0, "se": 0.0, "within": True})
            continue
        vc = fam_prec[c["key"]]
        se = float(np.std(bb - vc[idx].mean(axis=1), ddof=1))
        diff = b["core"]["p_1k"] - c["core"]["p_1k"]
        rows.append({"key": c["key"], "diff": diff, "se": se, "within": bool(diff <= Z * se)})
    by_key = {c["key"]: c for c in E}

    def cheapest(keys):
        cs = [by_key[k] for k in keys]
        return min(cs, key=lambda c: (c["core"]["R"], -c["core"]["p_1k"], c["key"]))
    within = [r["key"] for r in rows if r["within"]]
    out.update({"b": b["key"], "within": within, "rows": rows, "pick": cheapest(within)["key"]})
    for d in SENS_DELTAS:
        out[f"pick_d{d}"] = cheapest([r["key"] for r in rows if r["diff"] <= d])["key"]
    return out


# ------------------------------------------------------------------------------------ evaluation with a CPU limit
def _proc_cpu(pid):
    """user+sys CPU seconds of a running process (Linux /proc)."""
    with open(f"/proc/{pid}/stat") as fh:
        f = fh.read().rsplit(")", 1)[1].split()
    return (int(f[11]) + int(f[12])) / os.sysconf("SC_CLK_TCK")


CPU_POLL_S = 2.0


def make_limited_worker(E):
    class LimitedWorker(E.ParamWorker):
        """ParamWorker whose align() also stops the aligner once the case has used more than cpu_limit CPU seconds
        (error "cpu_limit"). Used only for the 10,000 timing guard of search candidates: a candidate stopped there
        is infeasible exactly as it would be after finishing above the limit."""

        def align(self, inp, out, timeout_s=E.DEFAULT_TIMEOUT_S, cpu_limit=None):
            if cpu_limit is None:
                return super().align(inp, out, timeout_s)
            t0 = time.perf_counter()
            try:
                if self.proc is None:
                    self._start(timeout_s)
                base = _proc_cpu(self.proc.pid)
                os.write(self.proc.stdin.fileno(), (json.dumps({"in": str(inp), "out": str(out)}) + "\n").encode())
                deadline = time.monotonic() + timeout_s
                while True:
                    try:
                        line = self._readline(min(deadline, time.monotonic() + CPU_POLL_S))
                        break
                    except TimeoutError:
                        if time.monotonic() >= deadline:
                            raise
                        if _proc_cpu(self.proc.pid) - base > cpu_limit:
                            used = self._reap(kill=True) - self.cpu_reported
                            return time.perf_counter() - t0, max(used, 0.0), f"cpu_limit: > {cpu_limit:.0f} s CPU"
                r = json.loads(line)
            except TimeoutError:
                used = self._reap(kill=True) - self.cpu_reported
                return time.perf_counter() - t0, max(used, 0.0), "timeout"
            except (EOFError, BrokenPipeError, OSError, ValueError, RuntimeError) as e:
                tail = self._stderr_tail()
                used = self._reap(kill=True) - self.cpu_reported
                return time.perf_counter() - t0, max(used, 0.0), f"worker exit ({type(e).__name__}): {tail}"
            self.n += 1
            self.cpu_reported += r["cpu"]
            if self.recycle and self.n >= self.recycle:
                self.close()
            return r["wall"], r["cpu"], r["error"]
    return LimitedWorker


class EnsEvaluator(Evaluator):
    """recipe.fit.Evaluator (cache, slot pool, leak guard, early stops) + a per-case CPU limit (the 10,000 timing
    guard of search candidates). A cached "cpu_limit" row is re-run when the same case is requested without a
    limit, so every row used for anything else had the full 1,800 s."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.Worker = make_limited_worker(self.E)
        self._limits = {}

    def evaluate(self, specs, case_ids, cpu_cap=None, stop_on_error=False, cpu_limits=None):
        case_ids = list(case_ids)
        self.check(case_ids)
        self._limits = dict(cpu_limits or {})
        for k, p in specs.items():
            if pkey(p) != k:
                raise ValueError(f"spec key mismatch {k[:12]}")
        tasks = []
        self._state = {}
        self._rerun = set()
        for k, p in specs.items():
            hits = self.lookup(k, case_ids)
            for c in [c for c, h in hits.items() if str(h.get("error") or "").startswith("cpu_limit")
                      and c not in self._limits]:
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
            todo = sorted((c for c in case_ids if c not in hits), key=lambda c: -FIT._case_cost(self.rows[c]))
            nch = max(1, math.ceil(len(todo) / self.chunk))
            for j in range(nch):
                ch = todo[j::nch]
                if ch:
                    tasks.append((k, p, ch, FIT._spec_cost(p) * sum(FIT._case_cost(self.rows[c]) for c in ch)))
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
                w = self.Worker(build.python, p)
                try:
                    for cid in ch:
                        if st["abort"]:
                            break
                        lim = self._limits.get(cid)
                        r = E._process_case(build, p, self.rows[cid], self.scorer, None, self.cache,
                                            (k, cid) in self._rerun,
                                            lambda inp, out, tmp, lim=lim: w.align(inp, out, self.timeout_s, lim))
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


# ------------------------------------------------------------------------------------ data
def select10k_path(ws):
    return Path(ws) / "runs" / "ensfit" / "select10k.json"


def load_ens_fold(ws, fold, select=None):
    """-> rows, meta, core, screen, test, lockbox, guard (ids). core = small core + every training family at 1,000 +
    the fold's 24 selection families at 10,000; screen = small screening subset + the fast fit's 40 screening families
    at 1,000 + the fold's guard family at 10,000."""
    from .fastfit import load_fold
    rows, meta, core0, screen0, test, lockbox = load_fold(ws, fold)
    sel = select or json.loads(select10k_path(ws).read_text())
    f = sel["folds"][str(fold)]
    sel10 = [f"ehf.{fam}.k10000" for fam in f["families"]]
    guard = f"ehf.{f['guard']}.k10000"
    core0 = set(core0)
    bad = [c for c in sel10 if c not in core0]
    if bad:
        raise LeakError(f"selection families at 10,000 outside the training data of fold {fold}: {bad[:3]}")
    if guard not in sel10:
        raise ValueError("guard family not among the selection families")
    small_core = sorted(c for c in core0 if meta[c]["kind"] == "small")
    k1_core = sorted(c for c in core0 if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)
    small_screen = sorted(c for c in screen0 if meta[c]["kind"] == "small")
    k1_screen = sorted(c for c in screen0 if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)
    core = small_core + k1_core + sorted(sel10)
    screen = small_screen + k1_screen          # AMENDMENT_ensemble_1.md (a): no 10,000 input in the search
    return rows, meta, core, screen, test, lockbox, guard


# ------------------------------------------------------------------------------------ settings
@dataclass
class Settings:
    pop_size: int = POP_SIZE
    ref_partitions: int = REF_PARTITIONS
    seed: int = SEED
    gens: tuple = GENS
    timeout_s: float = TIMEOUT_S
    note: str = ""

    def is_registered(self):
        d = Settings()
        return (self.pop_size, self.ref_partitions, self.seed, tuple(self.gens), self.timeout_s, self.note) == \
            (d.pop_size, d.ref_partitions, d.seed, tuple(d.gens), d.timeout_s, "")

    def to_json(self):
        d = asdict(self)
        d["gens"] = list(self.gens)
        return d


from pymoo.core.problem import Problem  # noqa: E402


class EnsProblem(Problem):
    """Module level so the pickled algorithm can be restored."""

    def __init__(self, vars):
        super().__init__(vars=vars, n_obj=3, n_ieq_constr=3)

    def _evaluate(self, X, out, *a, **kw):
        raise RuntimeError("evaluated via ask/tell")


SUMMARY_KEYS = ("p_small", "sp_small", "tc_small", "wr_small", "p_1k", "sp_1k", "tc_1k", "wr_1k", "n_1k", "p_10k",
                "sp_10k", "tc_10k", "wr_10k", "n_10k", "cpu_small", "cpu_1k", "cpu_10k", "guard_cpu", "guard_err",
                "cpu_ratio_small", "cpu_ratio_1k", "R", "n", "n_err", "n_err_small", "n_err_1k", "n_err_10k",
                "valid", "guard_ok", "feasible")


def brief_s(s):
    return {k: s.get(k) for k in SUMMARY_KEYS}


# ------------------------------------------------------------------------------------ the fit
class EnsFit:
    def __init__(self, fold, out_dir, rows, meta, core, screen, test, lockbox, guard, ev, settings, extra,
                 fast, acc):
        self.fold = fold
        self.out = Path(out_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.meta = meta
        self.core, self.screen = sorted(core), sorted(screen)
        self.guard = guard
        self.ev = ev
        self.st = settings
        self.extra = extra
        self.fast, self.acc = fast, acc
        self.space = Space(fast, acc)
        self.logfh = open(self.out / "progress.log", "a")
        for name, bad in (("test-fold", set(self.core) & set(test)), ("lockbox", set(self.core) & set(lockbox))):
            if bad:
                raise LeakError(f"{len(bad)} {name} cases in the training core of fold {fold}, e.g. {sorted(bad)[:3]}")
        if not set(self.screen) <= set(self.core):
            raise ValueError("screening data not inside the core")
        if len(set(self.core)) != len(self.core):
            raise ValueError("duplicate case in core")
        self.summ = {}

    def log(self, msg):
        log(f"[ensfit fold {self.fold}] {msg}", self.logfh)

    def _ids(self, scope):
        return self.screen if scope == "screen" else self.core

    def cap(self, scope):
        """Early stop (only ever an infeasible candidate): CPU over the evaluated cases above the sum of every CPU
        limit (small and 1,000 limits; 600 s for the guard on screen, 1,800 s per 10,000 case on the core)."""
        r = self.ref[scope]
        lim10 = TIMEOUT_S * r["n_10k"]           # no 10,000 input on the screening data (amendment 1)
        return CPU_SMALL_MAX * r["cpu_small"] + CPU_1K_MAX * r["cpu_1k"] + lim10

    def score(self, specs, scope, capped=False, raw=False):
        ids = self._ids(scope)
        need = {k: p for k, p in specs.items() if (scope, k) not in self.summ}
        if need:
            kw = {}
            if capped:
                kw = {"stop_on_error": True}      # amendment 2: no CPU-based early stop; real failures only
            res = self.ev.evaluate(need, ids, **kw)
            for k, rows in res.items():
                s = summarise(rows, self.meta, self.guard)
                self.summ[(scope, k)] = s if raw else objectives(s, self.ref[scope])
        return {k: self.summ[(scope, k)] for k in specs}

    def write_run_json(self):
        import pymoo
        if self.st.is_registered() and pymoo.__version__ != PYMOO_VERSION:
            raise SystemExit(f"registered fit needs pymoo {PYMOO_VERSION}, found {pymoo.__version__}")
        meta = {"plan": PLAN, "plan_sha256": PLAN_SHA256, "amendment": AMENDMENT, "fold": self.fold, "settings": self.st.to_json(),
                "registered": self.st.is_registered(), "pymoo": pymoo.__version__,
                "n_core": len(self.core), "n_screen": len(self.screen), "guard": self.guard,
                "core_sha256": FIT.ids_hash(self.core), "screen_sha256": FIT.ids_hash(self.screen),
                "key_fast": pkey(self.fast), "key_acc": pkey(self.acc),
                "ensfit_sha256": sha256_file(__file__), **self.extra}
        path = self.out / "run.json"
        if path.exists():
            old = json.loads(path.read_text())
            for k in ("fold", "settings", "core_sha256", "screen_sha256", "guard", "build", "key_fast",
                      "select10k_sha256"):
                if old.get(k) != meta.get(k):
                    raise SystemExit(f"{path}: {k} differs from the recorded run; refusing to resume")
            meta["started"] = old.get("started")
            meta["resumed"] = old.get("resumed", []) + [time.strftime("%Y-%m-%dT%H:%M:%S")]
            meta["ensfit_sha256_history"] = old.get("ensfit_sha256_history", [old.get("ensfit_sha256")]) + \
                [meta["ensfit_sha256"]]
        else:
            meta["started"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        _atomic_json(path, meta)
        (self.out / "cases_core.txt").write_text("".join(i + "\n" for i in self.core))
        (self.out / "cases_screen.txt").write_text("".join(i + "\n" for i in self.screen))

    def _write_refs(self):
        d = {"fast": self.fast, "key_fast": pkey(self.fast), "acc": self.acc, "key_acc": pkey(self.acc),
             "ref": {sc: brief_s(s) | {"per_benchmark": s["per_benchmark"]} for sc, s in self.ref.items()}}
        if getattr(self, "acc_s", None):
            d["acc_summary"] = {sc: brief_s(s) | {"per_benchmark": s["per_benchmark"]} for sc, s in self.acc_s.items()}
        _atomic_json(self.out / "refs.json", d)

    def references(self):
        kf = pkey(self.fast)
        self.log(f"fast setting on screen ({len(self.screen)}) and core ({len(self.core)})")
        self.ref = {}
        for scope in ("screen", "core"):
            s = self.score({kf: self.fast}, scope, raw=True)[kf]
            if not s["valid"]:
                raise RuntimeError(f"the fast setting failed on {s['n_err']} {scope} cases")
            self.ref[scope] = objectives(s, s)
            self.summ[(scope, kf)] = self.ref[scope]
            r = self.ref[scope]
            self.log(f"  {scope}: P_small={r['p_small']:.4f} SP_small={r['sp_small']:.4f} P_1k={r['p_1k']:.4f} "
                     f"SP_1k={r['sp_1k']:.4f} CPU small={r['cpu_small']:.1f}s 1k={r['cpu_1k']:.1f}s "
                     f"10k={r['cpu_10k']:.1f}s")
        self._write_refs()

    def acc_reference(self):
        """The pre-release accurate preset (reported row; not a candidate: its member 1 is not the fast setting)."""
        ka = pkey(self.acc)
        self.acc_s = {}
        for scope in ("screen", "core"):
            self.acc_s[scope] = self.score({ka: self.acc}, scope)[ka]
            a = self.acc_s[scope]
            self.log(f"  pre-release accurate {scope}: P_small={a['p_small']:.4f} P_1k={a['p_1k']:.4f} "
                     f"R={a['R']:.2f} errors={a['n_err']} feasible={a['feasible']}")
        self._write_refs()

    def search(self):
        from pymoo.algorithms.moo.nsga3 import NSGA3
        from pymoo.core.evaluator import Evaluator as PEval
        from pymoo.core.mixed import MixedVariableDuplicateElimination, MixedVariableMating, MixedVariableSampling
        from pymoo.core.population import Population
        from pymoo.problems.static import StaticProblem
        from pymoo.termination import get_termination
        from pymoo.util.ref_dirs import get_reference_directions

        gmin, gmax = self.st.gens
        norm = Normaliser(self.ref["screen"])
        problem = EnsProblem(self.space.vars())
        alg_path, hv_path, ev_path = self.out / "algorithm.pkl", self.out / "hv.json", self.out / "evaluated.jsonl"
        archive = {}
        if ev_path.exists():
            for line in open(ev_path):
                try:
                    e = json.loads(line)
                except json.JSONDecodeError:
                    continue
                archive.setdefault(e["key"], e)
        hvh = json.loads(hv_path.read_text()) if hv_path.exists() else {"history": []}
        if alg_path.exists():
            algorithm = pickle.load(open(alg_path, "rb"))
            done = algorithm.n_gen - 1
            hvh["history"] = hvh["history"][:done]
            if len(hvh["history"]) < done:
                hvh["history"] = [self._hv([e for e in archive.values() if e["gen"] <= g + 1]) for g in range(done)]
            self.log(f"resumed after generation {done}")
        else:
            ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=self.st.ref_partitions)
            sp = self.space.seed_params()
            seeds = [self.space.encode(p) for p in sp][:self.st.pop_size]
            for p, x in zip(sp, seeds):
                if pkey(self.space.decode(x)) != pkey(p):
                    raise RuntimeError("seed does not round-trip through encode/decode")
            n_rand = self.st.pop_size - len(seeds)
            init = list(MixedVariableSampling()(problem, n_rand, random_state=np.random.default_rng(self.st.seed))
                        .get("X")) if n_rand > 0 else []
            algorithm = NSGA3(ref_dirs=ref_dirs, pop_size=self.st.pop_size,
                              sampling=Population.new("X", np.array(seeds + init, dtype=object)),
                              mating=MixedVariableMating(eliminate_duplicates=MixedVariableDuplicateElimination()),
                              eliminate_duplicates=MixedVariableDuplicateElimination())
            algorithm.setup(problem, termination=get_termination("n_gen", gmax + 1000), seed=self.st.seed,
                            verbose=False)
            self.log(f"NSGA-III {len(ref_dirs)} directions, pop {self.st.pop_size}, {len(seeds)} seeds, "
                     f"generations {gmin}-{gmax}")
        stop, reason = hv_stop(hvh["history"], gmin, gmax, HV_WINDOW, HV_TOL)
        while not stop:
            gen = algorithm.n_gen if algorithm.n_gen is not None else 1
            tg = time.time()
            pop = algorithm.ask()
            X = pop.get("X")
            plist = [self.space.decode(x) for x in X]
            specs = {pkey(p): p for p in plist}
            summ = self.score(specs, "screen", capped=True)
            F, Gs = [], []
            with open(ev_path, "a") as fh:
                for x, p in zip(X, plist):
                    k = pkey(p)
                    s = summ[k]
                    f = norm.F(s) if s["valid"] else [2.0] * 3
                    F.append(f)
                    Gs.append(G(s))
                    if k not in archive:
                        e = {"gen": gen, "key": k, "params": p, "x": _py(dict(x)), "F": f, "G": G(s),
                             "screen": brief_s(s)}
                        archive[k] = e
                        fh.write(json.dumps(e, default=_json_default) + "\n")
            PEval().eval(StaticProblem(problem, F=np.array(F, float), G=np.array(Gs, float)), pop)
            algorithm.tell(infills=pop)
            hv = self._hv(archive.values())
            hvh["history"].append(hv)
            stop, reason = hv_stop(hvh["history"], gmin, gmax, HV_WINDOW, HV_TOL)
            hvh.update({"stop": stop, "reason": reason, "normaliser": norm.to_json()})
            _atomic_json(hv_path, hvh)
            tmp = alg_path.with_suffix(".tmp")
            with open(tmp, "wb") as fh:
                pickle.dump(algorithm, fh)
            tmp.replace(alg_path)
            nf = sum(1 for p in plist if summ[pkey(p)]["feasible"])
            cpu = sum((summ[k]["cpu_small"] + summ[k]["cpu_1k"] + summ[k]["cpu_10k"]) for k in specs)
            self.log(f"gen {gen}: {len(plist)} individuals ({len(specs)} distinct, {nf} feasible), archive "
                     f"{len(archive)}, HV {hv:.5f}, CPU {cpu:.0f}s, {time.time() - tg:.0f}s")
        self.log(f"search stopped after {len(hvh['history'])} generations: {reason}")
        hvh.update({"stop": True, "reason": reason, "n_gen": len(hvh["history"])})
        _atomic_json(hv_path, hvh)
        return archive

    def _hv(self, entries):
        return hypervolume([e["F"] for e in entries if e["screen"]["feasible"]], HV_REF)

    def rescore(self, archive):
        ents = sorted((e for e in archive.values() if e["screen"]["feasible"]), key=lambda e: e["key"])
        front = [ents[i] for i in nondominated([raw_F(e["screen"]) for e in ents])]
        pick, n_ok = rescore_set(front, self.ref["screen"])
        self.log(f"front {len(front)} on screen, {n_ok} with P_small >= fast - {RESCORE_SMALL_MARGIN}, "
                 f"{len(pick)} re-scored on the core")
        core = self.score({e["key"]: e["params"] for e in pick}, "core", capped=True)
        out = [{"key": e["key"], "params": e["params"], "gen": e["gen"], "screen": e["screen"],
                "core": core[e["key"]]} for e in pick]
        _atomic_json(self.out / "front_core.json",
                     {"n_front_screen": len(front), "n_qualified": n_ok, "front_screen_keys": [e["key"] for e in front],
                      "rescored": [dict(o, core=brief_s(o["core"])) for o in out]})
        return out

    def fam_precision(self, keys):
        """{key: per-family precision at 1,000 on the core, families in id order} (from the cache)."""
        ids = sorted(c for c in self.core if self.meta[c]["kind"] == "large" and self.meta[c]["k"] == 1000)
        out = {}
        for k in keys:
            got = self.ev.lookup(k, ids)
            if set(got) != set(ids) or any(got[c].get("error") for c in ids):
                raise RuntimeError(f"{k[:12]}: incomplete 1,000 core rows")
            out[k] = np.array([float(got[c]["precision"]) for c in ids])
        return out

    def run(self):
        self.write_run_json()
        self.log(f"fit start: core {len(self.core)}, screen {len(self.screen)}, guard {self.guard}, "
                 f"registered={self.st.is_registered()}")
        self.references()
        archive = self.search()
        resc = self.rescore(archive)
        self.log("pre-release accurate preset (reported row)")
        self.acc_reference()
        cands = [{"key": e["key"], "params": e["params"], "source": f"search gen {e['gen']}", "core": e["core"]}
                 for e in resc]
        elig_keys = [c["key"] for c in cands if eligible(c, self.ref["core"])]
        k = knee(cands, self.ref["core"], self.fam_precision(elig_keys))
        by = {c["key"]: c for c in cands}

        def brief(key, why):
            if key is None:
                return None
            c = by[key]
            canon = C.canonical_json(c["params"])
            return {"key": key, "why": why, "source": c["source"], "params": c["params"],
                    "params_sha256": hashlib.sha256(canon.encode()).hexdigest(), "short": self.space.short(c["params"]),
                    "core": brief_s(c["core"])}
        out = {"fold": self.fold, "registered": self.st.is_registered(), "rule": "PLAN_ensemble.md section 6",
               "ref_core": brief_s(self.ref["core"]), "acc_core": brief_s(self.acc_s["core"]),
               "n_candidates": len(cands), "n_eligible": k["n_eligible"],
               "pick": brief(k["pick"], f"cheapest of {len(k['within'])} within 1.96 SE of the best"),
               "b": brief(k["b"], "eligible point with the highest P_1k"),
               "knee_rows": k["rows"], "eligible": k["eligible"],
               "candidates": [{"key": c["key"], "source": c["source"], "eligible": c["key"] in elig_keys,
                               "core": brief_s(c["core"])} for c in sorted(cands, key=lambda c: c["key"])]}
        for d in SENS_DELTAS:
            out[f"sensitivity_d{d}"] = brief(k[f"pick_d{d}"], f"sensitivity: cheapest within {d} of the best")
        out["written"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        _atomic_json(self.out / "selection.json", out)
        if out["pick"]:
            pk = out["pick"]["core"]
            self.log(f"pick: {out['pick']['source']} P_small={pk['p_small']:.4f} P_1k={pk['p_1k']:.4f} "
                     f"P_10k={pk['p_10k']:.4f} R={pk['R']:.3f} | {out['pick']['short']}")
        else:
            self.log("no eligible point: no accurate setting from this fit")
        self.log("fit done")
        return out


def default_ws():
    return Path("~/tmp/kalign-fastfit").expanduser()


def use_ws_tmp(ws):
    """Every temporary file under WS/tmp (never /tmp): worker dirs, scorer files, stderr captures."""
    import tempfile
    tmp = Path(ws) / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    os.environ["KALIGN_REFIT_TMP"] = os.environ.get("KALIGN_REFIT_TMP") or str(tmp)
    os.environ["TMPDIR"] = os.environ.get("TMPDIR") if str(os.environ.get("TMPDIR", "")).startswith(
        str(Path.home())) else str(tmp)
    tempfile.tempdir = os.environ["TMPDIR"]


def make_evaluator(ws, rows, allowed, threads=SLOTS, timeout_s=TIMEOUT_S, cache=None):
    from .evaluate import Build
    use_ws_tmp(ws)
    b = Build.load(Path(ws) / "kalign" / BUILD)
    slots = SlotPool(Path(ws) / "runs" / "slots", SLOTS)
    cache = cache or Path(ws) / "cache" / "ensfit.sqlite"
    ev = EnsEvaluator(b, rows.values(), allowed, cache_path=cache, slots=slots, threads=threads, timeout_s=timeout_s)
    return ev, b, cache


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fold", required=True, help="0..4 or all")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ws", default=str(default_ws()))
    ap.add_argument("--threads", type=int, default=SLOTS)
    ap.add_argument("--presets-v2", default=str(PRESETS_V2))
    # smoke / pilot overrides (recorded; the run is then not registered)
    ap.add_argument("--pop-size", type=int, default=POP_SIZE)
    ap.add_argument("--gens", default=f"{GENS[0]},{GENS[1]}")
    ap.add_argument("--note", default="", help="recorded; makes the run unregistered (pilot, smoke)")
    ap.add_argument("--limit", type=int, default=0, help="smoke: at most N small and N/10 1,000 core cases")
    ap.add_argument("--search-only", action="store_true", help="pilot: stop after the search")
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    fold = a.fold if a.fold == "all" else int(a.fold)
    sel_path = select10k_path(ws)
    rows, meta, core, screen, test, lockbox, guard = load_ens_fold(ws, fold)
    note = [a.note] if a.note else []
    if a.limit:
        sm = sorted(c for c in core if meta[c]["kind"] == "small")[:a.limit]
        k1 = sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)[:max(2, a.limit // 10)]
        k10 = sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 10000)
        core = sm + k1 + k10
        screen = [c for c in screen if c in set(sm)] + k1[:2]
        note.append(f"limit={a.limit}")
    st = Settings(pop_size=a.pop_size, gens=tuple(int(v) for v in a.gens.split(",")), note=",".join(note))
    ev, b, cache = make_evaluator(ws, rows, core, threads=a.threads, timeout_s=st.timeout_s)
    fast, acc = fast_setting(ws), prerelease_accurate(a.presets_v2)
    extra = {"build": b.build_id, "build_cli": b.cli, "build_python": b.python, "cache": str(cache),
             "threads": a.threads, "slots": SLOTS, "select10k": str(sel_path),
             "select10k_sha256": sha256_file(sel_path),
             "final_pick_sha256": sha256_file(ws / "runs" / "fastfit" / "final_pick.json"),
             "presets_v2": a.presets_v2, "presets_v2_sha256": sha256_file(a.presets_v2)}
    f = EnsFit(fold, a.out, rows, meta, core, screen, test, lockbox, guard, ev, st, extra, fast, acc)
    ev.logfh = f.logfh
    if a.search_only:
        f.write_run_json()
        f.references()
        f.search()
        f.log("search-only run: stopped after the search")
        return 0
    f.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
