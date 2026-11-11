"""One fit of the fast-fit plan (PLAN_fastfit.md): single-run protein settings with the guide tree searched,
3 objectives (small accuracy, large accuracy, CPU relative to 4.0.0 fast), NSGA-III, re-scoring, selection.

    python -m recipe.fastfit --fold 0 --out WS/runs/fastfit/fold0 [--ws WS]           (0..4 or all)

Data (WS/data): small = the gm1 refit's protein pool, folds and cores (cases_protein.final.tsv, folds_protein.tsv,
cores/protein_fold<f>_{core,screen}.txt); large = recipe.largefam output (cases_exthomfam.tsv, folds_exthomfam.tsv,
cores/exthomfam_fold<f>_{core,screen}.txt). Builds: WS/kalign/99aecb2 (candidates; member gap model, tree
options), WS/kalign/daff955 (Kalign 4.0.0: the reference "4.0.0 fast").

Procedure (the plan's sections 4-6):
  1. references on screen and core: 4.0.0 fast (4.0.0 build, its fast parameters) and the gm1 fast preset (new build,
     4.0.0 tree values explicit).
  2. NSGA-III, 3 objectives, population 56, Das-Dennis 9 partitions, seed 1234; 6 seeds + random; screening data;
     stop >= 20 generations and HV gain < 0.1 % over 5, or 100.
  3. re-score <= 40 front points on the core (plan section 5 rule).
  4. selection on the core: fast pick, ultra-fast pick (plan section 6).

Output (out_dir): run.json, refs.json, algorithm.pkl, hv.json, evaluated.jsonl (every evaluated individual),
front_core.json, selection.json, progress.log. Resumable (checkpoint per generation; every alignment is cached).
"""
from __future__ import annotations

import argparse
import json
import math
import pickle
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from . import cache as C
from . import fit as FIT
from . import presets_370 as P370
from .fit import (Evaluator, LeakError, SlotPool, _atomic_json, _json_default, _py, hv_stop, hypervolume, log,
                  nondominated, pkey)

# up to 60 concurrent aligners (plan section 4; 64 hardware threads). SlotPool and Evaluator read FIT.MAX_SLOTS when
# they are called, so this raises their limit for this module only.
SLOTS = 60
FIT.MAX_SLOTS = SLOTS

SEED = 1234
POP_SIZE = 56
REF_PARTITIONS = 9            # Das-Dennis, 3 objectives -> 55 directions
GENS = (20, 100)
HV_WINDOW, HV_TOL, HV_REF = 5, 0.001, 1.1
Q_MARGIN = 0.10               # accuracy axes normalised from (Q_ref - 0.10)
R_LO, R_HI = 0.25, 3.0        # CPU axis: log R mapped from 0.25 to 3
CPU_SMALL_MAX, CPU_LARGE_MAX = 1.5, 2.0
CPU_STOP = 3.0                # a candidate is stopped once its CPU exceeds 3 x the reference on the same cases
TIMEOUT_S = 600.0
RESCORE_MAX, RESCORE_SMALL_MARGIN = 40, 0.01
FAST_SMALL_MARGIN, ULTRA_SMALL_MARGIN = 0.002, 0.01
PYMOO_VERSION = "0.6.1.6"
NEW_BUILD, REF_BUILD = "99aecb2", "daff955"

# Kalign 4.0.0 fast preset (single run; lib/src/aln_wrap.c of 4.0.0, as used by the explorations' BASE dict).
# Identity I1 checks that the 4.0.0 build's worker with this dict reproduces its CLI preset.
P400_FAST = {"n_runs": 1, "run_gpo": [6.190432441439749], "run_gpe": [0.8092411597874494],
             "run_tgpe": [0.5996379665635143], "run_noise": [0.0], "run_types": [1], "run_vsm_amax": [0.0],
             "run_refine": [0], "run_vsm_kappa": [1.0], "run_vsm_dzero": [1.0], "run_guide_dist": [0],
             "seq_weights": 2.8012471630813027, "consistency": 3, "consistency_weight": 2.891914131308884,
             "realign": 1, "min_support": 0, "consistency_merge": 0, "consistency_merge_weight": 2.0, "vsm_mode": 0,
             "fill_in": 0, "fill_anchor_frac": 0.0, "fill_support": 1}
# gm1 fast preset (runs/gm1_v2/refit/protein/final/selection.json, selection.fast.params)
GM1_FAST = {"n_runs": 1, "run_gpo": [7.833692432618236], "run_gpe": [0.9463291267641571],
            "run_tgpe": [2.697449842806514], "run_noise": [0.0], "run_types": [1], "run_vsm_amax": [0.7154252391041576],
            "run_refine": [0], "run_vsm_kappa": [1.0], "run_vsm_dzero": [1.0], "run_guide_dist": [0],
            "seq_weights": 3.6221996245704107, "consistency": 2, "consistency_weight": 0.6357651696858546,
            "realign": 1, "min_support": 0, "consistency_merge": 0, "consistency_merge_weight": 2.0, "vsm_mode": 0,
            "fill_in": 0, "fill_anchor_frac": 0.0, "fill_support": 1}

# ------------------------------------------------------------------------------------ guide-tree variables
TREE_METHODS = (0, 3)         # 0 bisecting k-means + WPGMA leaves (4.0.0), 3 hybrid: k-means + Prim leaves
TREE_DISTS = (0, 1, 2, 3)     # Myers raw, Myers normalised, LCS indel, 1 - LCS / shorter
LEAF_LOG = (math.log10(50), math.log10(5000))
MAXLEN = 1024
TREE_400 = {"method": 0, "dist": 0, "alpha": 0, "leaf": 50, "kd": 0, "rt": 1}
TREE_E21 = {"method": 3, "dist": 2, "alpha": 1, "leaf": 2000, "kd": 1, "rt": 1}


def with_tree(p, t, realign=None):
    """Worker dict p with the guide-tree options t (masked as the decoder does)."""
    q = dict(p)
    if realign is not None:
        q["realign"] = realign
    dist, alpha = (0, 0) if (t["method"] == 0 and t["kd"] == 0) else (t["dist"], t["alpha"])
    rt = 1 if q["realign"] == 0 else t["rt"]
    q.update({"run_tree_method": [t["method"]], "run_tree_dist": [dist], "run_tree_alpha": [alpha],
              "run_tree_maxlen": [MAXLEN], "run_tree_leaf": [int(t["leaf"])], "run_tree_kmeans_dist": [t["kd"]],
              "run_realign_tree": [rt]})
    return q


class Space:
    """Stage-1 single-run variables of the gm1 refit's protein space (read-only, sha256-pinned), without guide_dist,
    plus the guide-tree variables."""

    def __init__(self):
        self.sp = P370.load_space("protein")

    def vars(self):
        from pymoo.core.variable import Choice, Real
        v = {k: x for k, x in self.sp.stage1_vars().items() if k != "guide_dist"}
        v.update({"tree_method": Choice(options=list(TREE_METHODS)), "tree_leaf_log": Real(bounds=LEAF_LOG),
                  "tree_dist": Choice(options=list(TREE_DISTS)), "tree_alpha": Choice(options=[0, 1]),
                  "tree_kmeans_dist": Choice(options=[0, 1]), "realign_tree": Choice(options=[0, 1])})
        return v

    def decode(self, x):
        x = dict(x)
        x["guide_dist"] = 0
        p = _py(self.sp.stage1_decode(x))
        t = {"method": int(x["tree_method"]), "dist": int(x["tree_dist"]), "alpha": int(x["tree_alpha"]),
             "leaf": int(round(10 ** float(x["tree_leaf_log"]))), "kd": int(x["tree_kmeans_dist"]),
             "rt": int(x["realign_tree"])}
        return _py(with_tree(p, t))

    def encode(self, p):
        base = {k: v for k, v in p.items() if not k.startswith("run_tree_") and k != "run_realign_tree"}
        x = self.sp.stage1_encode(base)
        x.pop("guide_dist", None)
        x.update({"tree_method": p["run_tree_method"][0], "tree_leaf_log": math.log10(p["run_tree_leaf"][0]),
                  "tree_dist": p["run_tree_dist"][0], "tree_alpha": p["run_tree_alpha"][0],
                  "tree_kmeans_dist": p["run_tree_kmeans_dist"][0], "realign_tree": p["run_realign_tree"][0]})
        return x

    def seeds(self):
        out = []
        for p in (P400_FAST, GM1_FAST):
            out += [with_tree(p, TREE_400), with_tree(p, TREE_E21), with_tree(p, TREE_E21, realign=0)]
        return [_py(p) for p in out]

    def short(self, p):
        names = {0: "kmeans+WPGMA", 3: "kmeans+Prim"}
        dn = {0: "myers", 1: "myers_norm", 2: "lcs", 3: "lcs_frac"}
        base = {k: v for k, v in p.items() if not k.startswith("run_tree_") and k != "run_realign_tree"}
        return (f"{self.sp.short(base)} | tree={names[p['run_tree_method'][0]]} leaf={p['run_tree_leaf'][0]} "
                f"dist={dn[p['run_tree_dist'][0]]} alpha={'full' if p['run_tree_alpha'][0] else 'red'} "
                f"kmdist={'tree' if p['run_tree_kmeans_dist'][0] else 'bpm'} rt={p['run_realign_tree'][0]}")


# ------------------------------------------------------------------------------------ objectives
def summarise(rows, meta):
    """rows {case_id: cache row} -> summary: Q_small (per-case (SP+TC)/2, benchmark-weighted), SP_small, Q_large (mean
    of the per-k family means), per-k means, CPU per part, errors."""
    small, large = {}, {}
    cpu = {"small": 0.0, "large": 0.0}
    n_err = 0
    for cid, r in rows.items():
        m = meta[cid]
        cpu[m["kind"]] += float(r.get("cpu_s") or 0.0)
        if r.get("error"):
            n_err += 1
            continue
        q = 0.5 * (float(r["recall"]) + float(r["tc"]))
        if m["kind"] == "small":
            small.setdefault(m["benchmark"], []).append((q, float(r["recall"]), float(r["tc"])))
        else:
            large.setdefault(m["k"], []).append((q, float(r["recall"]), float(r["tc"])))
    per_b = {b: {"q": float(np.mean([v[0] for v in vs])), "sp": float(np.mean([v[1] for v in vs])),
                 "tc": float(np.mean([v[2] for v in vs])), "n": len(vs)} for b, vs in sorted(small.items())}
    per_k = {str(k): {"q": float(np.mean([v[0] for v in vs])), "sp": float(np.mean([v[1] for v in vs])),
                      "tc": float(np.mean([v[2] for v in vs])), "n": len(vs)} for k, vs in sorted(large.items())}
    s = {"q_small": float(np.mean([v["q"] for v in per_b.values()])) if per_b else 0.0,
         "sp_small": float(np.mean([v["sp"] for v in per_b.values()])) if per_b else 0.0,
         "tc_small": float(np.mean([v["tc"] for v in per_b.values()])) if per_b else 0.0,
         "q_large": float(np.mean([v["q"] for v in per_k.values()])) if per_k else 0.0,
         "cpu_small": cpu["small"], "cpu_large": cpu["large"], "n": len(rows), "n_err": n_err,
         "valid": n_err == 0 and len(rows) > 0, "per_benchmark": per_b, "per_k": per_k}
    return s


def objectives(s, ref):
    """Add the CPU ratios, R and feasibility (relative to the reference summary on the same cases)."""
    rs = s["cpu_small"] / max(ref["cpu_small"], 1e-9)
    rl = s["cpu_large"] / max(ref["cpu_large"], 1e-9)
    out = dict(s)
    out.update({"cpu_ratio_small": rs, "cpu_ratio_large": rl, "R": math.sqrt(max(rs, 1e-12) * max(rl, 1e-12)),
                "feasible": bool(s["valid"] and rs <= CPU_SMALL_MAX and rl <= CPU_LARGE_MAX)})
    return out


class Normaliser:
    """Monotone transform used by NSGA-III and the hypervolume (all minimised)."""

    def __init__(self, ref):
        self.lo_s = ref["q_small"] - Q_MARGIN
        self.lo_l = ref["q_large"] - Q_MARGIN

    def F(self, s):
        return [(1.0 - s["q_small"]) / (1.0 - self.lo_s), (1.0 - s["q_large"]) / (1.0 - self.lo_l),
                (math.log(max(s["R"], 1e-9)) - math.log(R_LO)) / (math.log(R_HI) - math.log(R_LO))]

    def to_json(self):
        return {"q_small_lo": self.lo_s, "q_large_lo": self.lo_l, "R_lo": R_LO, "R_hi": R_HI}


def raw_F(s):
    return [-s["q_small"], -s["q_large"], s["R"]]


def G(s):
    """Inequality constraints (<= 0 feasible): CPU caps and failed cases."""
    if not s["valid"]:
        return [1.0, 1.0, 1.0 + s["n_err"]]
    return [s["cpu_ratio_small"] - CPU_SMALL_MAX, s["cpu_ratio_large"] - CPU_LARGE_MAX, -1.0]


def rescore_set(front, ref_screen):
    """Plan section 5: front points with Q_small >= Q_small,ref - 0.01; at most 40: the 20 with the highest Q_large,
    then the 20 with the lowest R among the rest, then by Q_large."""
    ok = [e for e in front if e["screen"]["q_small"] >= ref_screen["q_small"] - RESCORE_SMALL_MARGIN]
    if len(ok) <= RESCORE_MAX:
        return sorted(ok, key=lambda e: e["key"]), len(ok)
    by_l = sorted(ok, key=lambda e: (-e["screen"]["q_large"], e["key"]))
    pick = by_l[:20]
    rest = [e for e in ok if e not in pick]
    pick += sorted(rest, key=lambda e: (e["screen"]["R"], e["key"]))[:20]
    rest = [e for e in by_l if e not in pick]
    pick += rest[:RESCORE_MAX - len(pick)]
    return sorted(pick, key=lambda e: e["key"]), len(ok)


def select(cands, ref, fallback):
    """Plan section 6 on core summaries. cands: [{key, params, source, core}], ref: 4.0.0 fast core summary,
    fallback: the gm1 fast candidate."""
    feas = sorted((c for c in cands if c["core"]["feasible"]), key=lambda c: c["key"])
    fast = [c for c in feas if c["core"]["q_small"] >= ref["q_small"] - FAST_SMALL_MARGIN]
    fast = sorted(fast, key=lambda c: (-c["core"]["q_large"], c["core"]["R"], c["key"]))
    ultra = [c for c in feas if c["core"]["q_small"] >= ref["q_small"] - ULTRA_SMALL_MARGIN
             and c["core"]["q_large"] >= ref["q_large"]]
    ultra = sorted(ultra, key=lambda c: (c["core"]["R"], -c["core"]["q_large"], c["key"]))
    return {"fast": fast[0] if fast else None, "fast_fallback": fallback if not fast else None,
            "ultra": ultra[0] if ultra else None, "n_fast_eligible": len(fast), "n_ultra_eligible": len(ultra)}


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


class FastProblem(Problem):
    """Module level so the pickled algorithm can be restored."""

    def __init__(self, vars):
        super().__init__(vars=vars, n_obj=3, n_ieq_constr=3)

    def _evaluate(self, X, out, *a, **kw):
        raise RuntimeError("evaluated via ask/tell")


# ------------------------------------------------------------------------------------ data
def read_ids(path):
    return [l.strip() for l in open(path) if l.strip()]


def load_fold(ws, fold):
    """-> rows {id: manifest row}, meta {id: {...}}, core, screen, test, lockbox (ids)."""
    import csv
    from .scoring import read_manifest
    ws = Path(ws)
    rows, meta = {}, {}
    for r in read_manifest(ws / "data" / "cases_protein.final.tsv"):
        rows[r["case_id"]] = r
    for r in csv.DictReader(open(ws / "data" / "folds_protein.tsv"), delimiter="\t"):
        meta[r["case_id"]] = {"kind": "small", "benchmark": rows[r["case_id"]]["benchmark"], "fold": r["fold"],
                              "group": r["group_id"]}
    for r in read_manifest(ws / "data" / "cases_exthomfam.tsv"):
        rows[r["case_id"]] = r
    for r in csv.DictReader(open(ws / "data" / "folds_exthomfam.tsv"), delimiter="\t"):
        meta[r["case_id"]] = {"kind": "large", "benchmark": rows[r["case_id"]]["benchmark"], "fold": r["fold"],
                              "group": "ehf." + r["family"], "family": r["family"], "k": int(r["k"])}
    if set(meta) != set(rows):
        raise SystemExit("folds and manifests disagree")
    cores = ws / "data" / "cores"
    core = read_ids(cores / f"protein_fold{fold}_core.txt") + read_ids(cores / f"exthomfam_fold{fold}_core.txt")
    screen = read_ids(cores / f"protein_fold{fold}_screen.txt") + read_ids(cores / f"exthomfam_fold{fold}_screen.txt")
    lockbox = frozenset(c for c, m in meta.items() if m["fold"] == "lockbox")
    test = frozenset() if str(fold) == "all" else frozenset(c for c, m in meta.items() if m["fold"] == str(fold))
    return rows, meta, core, screen, test, lockbox


# ------------------------------------------------------------------------------------ the fit
class FastFit:
    def __init__(self, fold, out_dir, rows, meta, core, screen, test, lockbox, ev_new, ev_ref, settings, extra):
        self.fold = fold
        self.out = Path(out_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.meta = meta
        self.core, self.screen = sorted(core), sorted(screen)
        self.ev, self.ev_ref = ev_new, ev_ref
        self.st = settings
        self.extra = extra
        self.space = Space()
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
        log(f"[fastfit fold {self.fold}] {msg}", self.logfh)

    def _ids(self, scope):
        return self.screen if scope == "screen" else self.core

    def score(self, specs, scope, capped=False, ref_build=False, raw=False):
        """raw: summaries without the CPU ratios (only for the 4.0.0 reference itself, before self.ref exists)."""
        ids = self._ids(scope)
        tag = "ref" if ref_build else "new"
        need = {k: p for k, p in specs.items() if (tag, scope, k) not in self.summ}
        if need:
            kw = {}
            if capped:
                r = self.ref[scope]
                kw = {"cpu_cap": CPU_STOP * (r["cpu_small"] + r["cpu_large"]), "stop_on_error": True}
            res = (self.ev_ref if ref_build else self.ev).evaluate(need, ids, **kw)
            for k, rows in res.items():
                s = summarise(rows, self.meta)
                self.summ[(tag, scope, k)] = s if raw else objectives(s, self.ref[scope])
        return {k: self.summ[(tag, scope, k)] for k in specs}

    def write_run_json(self):
        import pymoo
        sp = self.space.sp
        if C.sha256_file(sp.__file__) != P370.SPACE_SHA256["protein"]:
            raise SystemExit(f"{sp.__file__}: search space differs from the pinned one")
        if self.st.is_registered() and pymoo.__version__ != PYMOO_VERSION:
            raise SystemExit(f"registered fit needs pymoo {PYMOO_VERSION}, found {pymoo.__version__}")
        meta = {"plan": "PLAN_fastfit.md", "fold": self.fold, "settings": self.st.to_json(),
                "registered": self.st.is_registered(), "pymoo": pymoo.__version__,
                "n_core": len(self.core), "n_screen": len(self.screen),
                "core_sha256": FIT.ids_hash(self.core), "screen_sha256": FIT.ids_hash(self.screen),
                "space_sha256": C.sha256_file(sp.__file__), **self.extra}
        path = self.out / "run.json"
        if path.exists():
            old = json.loads(path.read_text())
            for k in ("fold", "settings", "core_sha256", "screen_sha256", "space_sha256", "build_new", "build_ref"):
                if old.get(k) != meta.get(k):
                    raise SystemExit(f"{path}: {k} differs from the recorded run; refusing to resume")
            meta["started"] = old.get("started")
            meta["resumed"] = old.get("resumed", []) + [time.strftime("%Y-%m-%dT%H:%M:%S")]
        else:
            meta["started"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        _atomic_json(path, meta)
        (self.out / "cases_core.txt").write_text("".join(i + "\n" for i in self.core))
        (self.out / "cases_screen.txt").write_text("".join(i + "\n" for i in self.screen))

    def references(self):
        """4.0.0 fast (4.0.0 build) and the gm1 fast preset (new build, 4.0.0 tree values explicit)."""
        k400 = pkey(P400_FAST)
        self.log(f"4.0.0 fast (build {REF_BUILD}) on screen ({len(self.screen)}) and core ({len(self.core)})")
        self.ref = {}
        for scope in ("screen", "core"):
            s = self.score({k400: P400_FAST}, scope, ref_build=True, raw=True)[k400]
            if not s["valid"]:
                raise RuntimeError(f"4.0.0 fast failed on {s['n_err']} {scope} cases")
            self.ref[scope] = s
        for scope in ("screen", "core"):
            self.ref[scope] = objectives(self.ref[scope], self.ref[scope])
            self.summ[("ref", scope, k400)] = self.ref[scope]
        self.gm1 = _py(with_tree(GM1_FAST, TREE_400))
        kg = pkey(self.gm1)
        self.gm1_s = {scope: self.score({kg: self.gm1}, scope)[kg] for scope in ("screen", "core")}
        for scope in ("screen", "core"):
            r, g = self.ref[scope], self.gm1_s[scope]
            self.log(f"  {scope}: 4.0.0 fast Q_small={r['q_small']:.4f} Q_large={r['q_large']:.4f} "
                     f"CPU small={r['cpu_small']:.0f}s large={r['cpu_large']:.0f}s | gm1 fast Q_small="
                     f"{g['q_small']:.4f} Q_large={g['q_large']:.4f} R={g['R']:.3f}")
        _atomic_json(self.out / "refs.json", {"p400_fast": P400_FAST, "key_400": k400, "build_400": REF_BUILD,
                                              "gm1_fast": self.gm1, "key_gm1": kg, "ref": self.ref,
                                              "gm1": self.gm1_s})

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
        problem = FastProblem(self.space.vars())
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
            seeds = [self.space.encode(p) for p in self.space.seeds()][:self.st.pop_size]
            for p, x in zip(self.space.seeds(), seeds):     # the encoding must round-trip
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
                             "screen": {kk: s[kk] for kk in ("q_small", "sp_small", "tc_small", "q_large", "per_k",
                                                             "cpu_small", "cpu_large", "cpu_ratio_small",
                                                             "cpu_ratio_large", "R", "n_err", "valid", "feasible")}}
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
            self.log(f"gen {gen}: {len(plist)} individuals ({len(specs)} distinct, {nf} feasible), archive "
                     f"{len(archive)}, HV {hv:.5f}, {time.time() - tg:.0f}s")
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
        self.log(f"front {len(front)} on screen, {n_ok} with Q_small >= ref - {RESCORE_SMALL_MARGIN}, "
                 f"{len(pick)} re-scored on the core")
        core = self.score({e["key"]: e["params"] for e in pick}, "core", capped=True)
        out = [{"key": e["key"], "params": e["params"], "gen": e["gen"], "screen": e["screen"], "core": core[e["key"]]}
               for e in pick]
        _atomic_json(self.out / "front_core.json", {"n_front_screen": len(front), "n_qualified": n_ok,
                                                    "front_screen_keys": [e["key"] for e in front], "rescored": out})
        return out

    def run(self):
        self.write_run_json()
        self.log(f"fit start: core {len(self.core)}, screen {len(self.screen)}, registered={self.st.is_registered()}")
        self.references()
        archive = self.search()
        resc = self.rescore(archive)
        kg = pkey(self.gm1)
        cands = [{"key": e["key"], "params": e["params"], "source": f"search gen {e['gen']}", "core": e["core"]}
                 for e in resc]
        gm1c = {"key": kg, "params": self.gm1, "source": "gm1 fast preset", "core": self.gm1_s["core"]}
        if kg not in {c["key"] for c in cands}:
            cands.append(gm1c)
        sel = select(cands, self.ref["core"], gm1c)

        def brief(c):
            if c is None:
                return None
            s = c["core"]
            return {"key": c["key"], "source": c["source"], "params": c["params"], "short": self.space.short(c["params"]),
                    "core": {k: s[k] for k in ("q_small", "sp_small", "tc_small", "q_large", "per_k", "cpu_small",
                                               "cpu_large", "cpu_ratio_small", "cpu_ratio_large", "R", "feasible")}}
        out = {"fold": self.fold, "registered": self.st.is_registered(), "n_candidates": len(cands),
               "ref_core": {k: self.ref["core"][k] for k in ("q_small", "sp_small", "tc_small", "q_large", "per_k",
                                                             "cpu_small", "cpu_large")},
               "fast": brief(sel["fast"] or sel["fast_fallback"]), "fast_is_fallback": sel["fast"] is None,
               "ultra": brief(sel["ultra"]), "n_fast_eligible": sel["n_fast_eligible"],
               "n_ultra_eligible": sel["n_ultra_eligible"],
               "candidates": [brief(c) for c in sorted(cands, key=lambda c: c["key"])]}
        _atomic_json(self.out / "selection.json", out)
        f = out["fast"]
        self.log(f"fast pick: {f['source']}{' (fallback)' if out['fast_is_fallback'] else ''} Q_small="
                 f"{f['core']['q_small']:.4f} Q_large={f['core']['q_large']:.4f} R={f['core']['R']:.3f} | "
                 f"{f['short']}")
        if out["ultra"]:
            u = out["ultra"]
            self.log(f"ultra-fast pick: {u['source']} Q_small={u['core']['q_small']:.4f} Q_large="
                     f"{u['core']['q_large']:.4f} R={u['core']['R']:.3f} | {u['short']}")
        else:
            self.log("ultra-fast pick: none eligible")
        self.log("fit done")
        return out


def default_ws():
    return Path("~/tmp/kalign-fastfit").expanduser()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fold", required=True, help="0..4 or all")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ws", default=str(default_ws()))
    ap.add_argument("--threads", type=int, default=SLOTS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--random-init", action="store_true", help="no seeds (seed-robustness repeat, plan section 7)")
    # smoke-test overrides (recorded; the run is then not registered)
    ap.add_argument("--pop-size", type=int, default=POP_SIZE)
    ap.add_argument("--ref-partitions", type=int, default=REF_PARTITIONS)
    ap.add_argument("--gens", default=f"{GENS[0]},{GENS[1]}")
    ap.add_argument("--limit", type=int, default=0, help="smoke: at most N small and N/10 large core cases")
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    fold = a.fold if a.fold == "all" else int(a.fold)
    rows, meta, core, screen, test, lockbox = load_fold(ws, fold)
    note = []
    if a.limit:
        sm = sorted(c for c in core if meta[c]["kind"] == "small")[:a.limit]
        lg = sorted(c for c in core if meta[c]["kind"] == "large")[:max(2, a.limit // 10)]
        core = sm + lg
        screen = [c for c in screen if c in set(core) and meta[c]["kind"] == "small"] + lg[:2]
        note.append(f"limit={a.limit}")
    if a.seed != SEED:
        note.append(f"seed={a.seed}")
    if a.random_init:
        note.append("random_init")
    st = Settings(pop_size=a.pop_size, ref_partitions=a.ref_partitions, seed=a.seed,
                  gens=tuple(int(v) for v in a.gens.split(",")), note=",".join(note))
    from .evaluate import Build
    b_new, b_ref = Build.load(ws / "kalign" / NEW_BUILD), Build.load(ws / "kalign" / REF_BUILD)
    slots = SlotPool(ws / "runs" / "slots", SLOTS)
    allowed = set(core)
    cache = ws / "cache" / "cache.sqlite"
    ev_new = Evaluator(b_new, rows.values(), allowed, cache_path=cache, slots=slots, threads=a.threads,
                       timeout_s=st.timeout_s)
    ev_ref = Evaluator(b_ref, rows.values(), allowed, cache_path=cache, slots=slots, threads=a.threads,
                       timeout_s=st.timeout_s)
    extra = {"build_new": b_new.build_id, "build_ref": b_ref.build_id, "cache": str(cache), "threads": a.threads,
             "slots": SLOTS}
    f = FastFit(fold, a.out, rows, meta, core, screen, test, lockbox, ev_new, ev_ref, st, extra)
    if a.random_init:
        f.space.seeds = lambda: []
    ev_new.logfh = ev_ref.logfh = f.logfh
    f.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
