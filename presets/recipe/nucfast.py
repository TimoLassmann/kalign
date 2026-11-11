"""One fit of the nucfit plan's fast setting (presets/docs/nucfit/PLAN_nucfit.md, section 4): single-run nucleotide
settings with the guide tree and the gap model searched, 3 objectives (small accuracy, large accuracy, CPU relative
to the pre-release nucleotide fast setting 4ff7e4d1), NSGA-III, re-scoring, knee selection.

    python -m recipe.nucfast --fold 0 --out WS/runs/nucfast/fold0 [--ws WS]       (0..4 or all)

Everything not listed here is recipe.fastfit unchanged (the protein fast fit: NSGA-III settings, stop rule,
objectives, normalisation, constraints, re-scoring rule, cache, resumable).

Differences from recipe.fastfit:
  * data: recipe.nucfit_data (cases_nucfit.tsv, cores/nucsmall_*, cores/nuclarge_*); small = BRAliBASE, MDSA, Rfam
    (benchmarks weighted equally); large = Rfam seed families > 200 sequences, simulated DNA at 1,000 and 10,000
    (Q_large = mean of the three category means, as the protein per-k means).
  * search space: the pinned nucleotide stage-1 space (recipe.fit.NucleotideSpace over nested_nuc/nuc_space.py)
    + the protein fast fit's guide-tree variables + the gap model (0 Kalign rule, 1 member model).
  * one build for candidates and reference: WS/kalign/9137658 (the release build; its nucleotide presets are the
    reference, so the reference is scored with the same binary as the candidates).
  * reference: the pre-release nucleotide fast setting (presets/results/v2/presets_v2.json, key 4ff7e4d1) with its
    4.0.0 guide tree written explicitly; seeds: pre-release fast and accurate x guide trees {4.0.0, E21, protein
    fast} x gap model {0, 1}.
  * selection: the knee rule of DECISION_fastfit_selection.md, written in from the start: eligible = feasible,
    small SP and small TC >= the reference (core); b = highest Q_large; keep the points whose Q_large is within
    1.96 paired bootstrap SE of b (2,000 resamples of large cases, stratified by category, seed 1234); pick the
    lowest R (ties: higher Q_large, then key). No eligible point: the reference ships (fallback).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import fastfit as FF
from . import fit as FIT
from . import presets_370 as P370
from .fit import Evaluator, SlotPool, _atomic_json, _py, pkey

BUILD = "9137658"
N_BOOT, Z = 2000, 1.96

# PLAN_nucfit.md amendment A3 (2026-11-03, user: "relax the CPU limits"; after generation 2 of the A1 launch, before
# any selection; only feasibility counts were looked at): at generation 2, 0-1 of 56 candidates per fit met the
# limits (1.5x small, 2x large), so the search was constrained by speed, not accuracy. New limits, set here for the
# nucleotide fast fit only (recipe.fastfit reads these module globals at call time; the protein fit is untouched):
NUC_CPU_SMALL_MAX = NUC_CPU_LARGE_MAX = 10.0   # feasible: CPU <= 10x the reference (small and large)
NUC_CPU_STOP = 15.0                            # a candidate is stopped beyond 15x the reference CPU
NUC_R_HI = 10.0                                # CPU axis (log R) normalised over [0.25, 10] so the wider range counts
NUC_TIMEOUT_S = 1800.0                         # per alignment, as the accurate fit
FF.CPU_SMALL_MAX, FF.CPU_LARGE_MAX = NUC_CPU_SMALL_MAX, NUC_CPU_LARGE_MAX
FF.CPU_STOP, FF.R_HI = NUC_CPU_STOP, NUC_R_HI


@dataclass
class NucSettings(FF.Settings):
    timeout_s: float = NUC_TIMEOUT_S

    def is_registered(self):
        d = NucSettings()
        return (self.pop_size, self.ref_partitions, self.seed, tuple(self.gens), self.timeout_s, self.note) == \
            (d.pop_size, d.ref_partitions, d.seed, tuple(d.gens), d.timeout_s, "")
GAP_MODELS = (0, 1)

# presets/results/v2/presets_v2.json, presets.nucleotide.{fast,accurate}.params (keys 4ff7e4d1, bd03cdb7)
NUC_FAST = {"n_runs": 1, "run_gpo": [15.584104422293285], "run_gpe": [0.6147504027876776],
            "run_tgpe": [2.9040783529268577], "run_noise": [0.0], "run_types": [12],
            "run_vsm_amax": [0.2943446472264315], "run_refine": [0], "seq_weights": 0.427509434153184,
            "consistency": 1, "consistency_weight": 1.1977573259542555, "realign": 0, "min_support": 0,
            "consistency_merge": 0, "consistency_merge_weight": 2.0, "vsm_mode": 0, "run_vsm_kappa": [1.0],
            "run_vsm_dzero": [1.0], "run_guide_dist": [0], "fill_in": 0, "fill_anchor_frac": 0.0,
            "fill_support": 1}
NUC_FAST_KEY = "4ff7e4d1ccaafdc4fe9b6602af04f23f91479e43d8c718a5ce3f328cb20aadd2"
NUC_ACC_KEY = "bd03cdb7dcff5dde16e9be6462d0ccddb7143959523e997039704c8eb05dd9d7"
TREE_PFAST = {"method": 3, "dist": 2, "alpha": 0, "leaf": 1163, "kd": 1, "rt": 1}   # protein fast 0976900c


def v2_params(key):
    """Worker dict of a pre-release preset from presets_v2.json, checked against its key."""
    d = json.loads((Path(__file__).resolve().parents[1] / "results" / "v2" / "presets_v2.json").read_text())
    for mode in d["presets"]["nucleotide"].values():
        if mode["key"] == key:
            if pkey(mode["params"]) != key:
                raise SystemExit(f"presets_v2.json: params do not hash to {key[:8]}")
            return _py(mode["params"])
    raise SystemExit(f"{key[:8]} not in presets_v2.json")


def with_gap(p, g):
    q = dict(p)
    q["run_gap_model"] = [int(g)]
    return q


class NucSpace:
    """Nucleotide stage-1 variables (pinned space) + guide tree (as recipe.fastfit) + gap model."""

    def __init__(self):
        self.ns = FIT.NucleotideSpace()
        self.sp = self.ns.sp          # the pinned module (sha256 checked in write_run_json)

    def vars(self):
        from pymoo.core.variable import Choice, Real
        v = dict(self.ns.stage1_vars())
        v.update({"tree_method": Choice(options=list(FF.TREE_METHODS)), "tree_leaf_log": Real(bounds=FF.LEAF_LOG),
                  "tree_dist": Choice(options=list(FF.TREE_DISTS)), "tree_alpha": Choice(options=[0, 1]),
                  "tree_kmeans_dist": Choice(options=[0, 1]), "realign_tree": Choice(options=[0, 1]),
                  "gap_model": Choice(options=list(GAP_MODELS))})
        return v

    def decode(self, x):
        p = self.ns.stage1_decode(x)
        t = {"method": int(x["tree_method"]), "dist": int(x["tree_dist"]), "alpha": int(x["tree_alpha"]),
             "leaf": int(round(10 ** float(x["tree_leaf_log"]))), "kd": int(x["tree_kmeans_dist"]),
             "rt": int(x["realign_tree"])}
        return _py(with_gap(FF.with_tree(p, t), x["gap_model"]))

    def encode(self, p):
        base = {k: v for k, v in p.items()
                if not k.startswith("run_tree_") and k not in ("run_realign_tree", "run_gap_model")}
        x = self.ns.stage1_encode(base)
        x.update({"tree_method": p["run_tree_method"][0], "tree_leaf_log": math.log10(p["run_tree_leaf"][0]),
                  "tree_dist": p["run_tree_dist"][0], "tree_alpha": p["run_tree_alpha"][0],
                  "tree_kmeans_dist": p["run_tree_kmeans_dist"][0], "realign_tree": p["run_realign_tree"][0],
                  "gap_model": p["run_gap_model"][0]})
        return x

    def seeds(self):
        out = []
        for p in (v2_params(NUC_FAST_KEY), v2_params(NUC_ACC_KEY)):
            for t in (FF.TREE_400, FF.TREE_E21, TREE_PFAST):
                for g in GAP_MODELS:
                    out.append(_py(with_gap(FF.with_tree(p, t), g)))
        return out

    def short(self, p):
        names = {0: "kmeans+WPGMA", 3: "kmeans+Prim"}
        dn = {0: "myers", 1: "myers_norm", 2: "lcs", 3: "lcs_frac"}
        base = {k: v for k, v in p.items()
                if not k.startswith("run_tree_") and k not in ("run_realign_tree", "run_gap_model")}
        return (f"{self.ns.short(base)} | gap_model={p.get('run_gap_model', [0])[0]} "
                f"tree={names[p['run_tree_method'][0]]} leaf={p['run_tree_leaf'][0]} "
                f"dist={dn[p['run_tree_dist'][0]]} alpha={'full' if p['run_tree_alpha'][0] else 'red'} "
                f"kmdist={'tree' if p['run_tree_kmeans_dist'][0] else 'bpm'} rt={p['run_realign_tree'][0]}")


def load_fold(ws, fold):
    """-> rows, meta, core, screen, test, lockbox (the recipe.fastfit.load_fold contract)."""
    from .scoring import read_manifest
    ws = Path(ws)
    rows, meta = {}, {}
    extra = {r["case_id"]: r for r in csv.DictReader(open(ws / "data" / "cases_nucfit.tsv"), delimiter="\t")}
    for r in read_manifest(ws / "data" / "cases_nucfit.tsv"):
        e = extra[r["case_id"]]
        rows[r["case_id"]] = r
        # role "excluded" (amendment A4): cases removed to split an oversized group; in no core and no test fold
        fold_s = e["role"] if e["role"] in ("lockbox", "excluded") else str(e["fold"])
        if e["tier"] == "small":
            meta[r["case_id"]] = {"kind": "small", "benchmark": e["benchmark"], "fold": fold_s,
                                  "group": e["group_id"], "seq_kind": e["kind"]}
        else:
            k = "rfam" if e["benchmark"] == "rfam_large" else e["size_cat"]
            meta[r["case_id"]] = {"kind": "large", "benchmark": e["benchmark"], "fold": fold_s,
                                  "group": e["group_id"], "family": e["family"], "k": k, "seq_kind": e["kind"]}
    cores = ws / "data" / "cores"
    core = FF.read_ids(cores / f"nucsmall_fold{fold}_core.txt") + FF.read_ids(cores / f"nuclarge_fold{fold}_core.txt")
    screen = (FF.read_ids(cores / f"nucsmall_fold{fold}_screen.txt")
              + FF.read_ids(cores / f"nuclarge_fold{fold}_screen.txt"))
    # PLAN_nucfit.md amendment A1 (2026-11-03): no simulation at 10,000 in the search (as the protein ensemble fit,
    # AMENDMENT_ensemble_1 (a)). They took 98 % of the search CPU (median 209 s per alignment, every timeout); they
    # stay in the core, so re-scoring, Q_large and the selection still include them.
    screen = [c for c in screen if meta[c].get("k") != "10000"]
    lockbox = frozenset(c for c, m in meta.items() if m["fold"] == "lockbox")
    excluded = {c for c, m in meta.items() if m["fold"] == "excluded"}
    if excluded & (set(core) | set(screen)):
        raise FF.LeakError(f"fold {fold}: excluded case in the core or screen")
    test = frozenset() if str(fold) == "all" else frozenset(c for c, m in meta.items() if m["fold"] == str(fold))
    return rows, meta, core, screen, test, lockbox


# ------------------------------------------------------------------------------------ knee selection
def large_values(ev, key, params, cases):
    """{case_id: (SP+TC)/2} of one candidate on the core's large cases (cached rows; 0 for a failed case)."""
    res = ev.evaluate({key: params}, cases)[key]
    return {c: 0.0 if r.get("error") else 0.5 * (float(r["recall"]) + float(r["tc"])) for c, r in res.items()}


def knee(cands, ref, ev, meta, core):
    """-> (pick or None, info). cands: [{key, params, core, ...}] with core summaries."""
    E = sorted((c for c in cands if c["core"]["feasible"] and c["core"]["sp_small"] >= ref["sp_small"]
                and c["core"]["tc_small"] >= ref["tc_small"]), key=lambda c: c["key"])
    if not E:
        return None, {"n_eligible": 0}
    cases = sorted(c for c in core if meta[c]["kind"] == "large")
    cats = sorted({meta[c]["k"] for c in cases})
    idx = {k: np.array([i for i, c in enumerate(cases) if meta[c]["k"] == k]) for k in cats}
    rng = np.random.default_rng(FF.SEED)
    draws = [{k: v[rng.integers(0, len(v), len(v))] for k, v in idx.items()} for _ in range(N_BOOT)]

    def boot(vals):
        return np.array([np.mean([vals[d[k]].mean() for k in cats]) for d in draws])

    b = max(E, key=lambda c: (c["core"]["q_large"], -c["core"]["R"]))
    vb = np.array([large_values(ev, b["key"], b["params"], cases)[c] for c in cases])
    full = float(np.mean([vb[idx[k]].mean() for k in cats]))
    if abs(full - b["core"]["q_large"]) > 1e-9:
        raise SystemExit(f"recomputed Q_large {full:.6f} != core summary {b['core']['q_large']:.6f}")
    bb = boot(vb)
    within, rows = [b], []
    for c in E:
        if c["key"] == b["key"]:
            continue
        vc = np.array([large_values(ev, c["key"], c["params"], cases)[x] for x in cases])
        diff = b["core"]["q_large"] - c["core"]["q_large"]
        se = float(np.std(bb - boot(vc), ddof=1))
        rows.append({"key": c["key"], "diff": diff, "se": se, "within": diff <= Z * se})
        if diff <= Z * se:
            within.append(c)
    pick = min(within, key=lambda c: (c["core"]["R"], -c["core"]["q_large"], c["key"]))
    return pick, {"n_eligible": len(E), "best": b["key"], "n_within": len(within), "comparisons": rows}


class NucFastFit(FF.FastFit):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.space = NucSpace()

    def log(self, msg):
        FIT.log(f"[nucfast fold {self.fold}] {msg}", self.logfh)

    def write_run_json(self):
        import pymoo
        from . import cache as C
        if C.sha256_file(self.space.sp.__file__) != P370.SPACE_SHA256["nucleotide"]:
            raise SystemExit(f"{self.space.sp.__file__}: search space differs from the pinned one")
        if self.st.is_registered() and pymoo.__version__ != FF.PYMOO_VERSION:
            raise SystemExit(f"registered fit needs pymoo {FF.PYMOO_VERSION}, found {pymoo.__version__}")
        meta = {"plan": "nucfit/PLAN_nucfit.md", "fold": self.fold, "settings": self.st.to_json(),
                "registered": self.st.is_registered(), "pymoo": pymoo.__version__,
                "n_core": len(self.core), "n_screen": len(self.screen),
                "core_sha256": FIT.ids_hash(self.core), "screen_sha256": FIT.ids_hash(self.screen),
                "space_sha256": C.sha256_file(self.space.sp.__file__), **self.extra}
        path = self.out / "run.json"
        if path.exists():
            old = json.loads(path.read_text())
            for k in ("fold", "settings", "core_sha256", "screen_sha256", "space_sha256", "build_new", "build_ref"):
                if old.get(k) != meta.get(k):
                    raise SystemExit(f"{path}: {k} differs from the recorded run; refusing to resume")
            meta["started"] = old.get("started")
            meta["resumed"] = old.get("resumed", []) + [FIT.time.strftime("%Y-%m-%dT%H:%M:%S")]
        else:
            meta["started"] = FIT.time.strftime("%Y-%m-%dT%H:%M:%S")
        _atomic_json(path, meta)
        (self.out / "cases_core.txt").write_text("".join(i + "\n" for i in self.core))
        (self.out / "cases_screen.txt").write_text("".join(i + "\n" for i in self.screen))

    def references(self):
        """The pre-release nucleotide fast setting, with its 4.0.0 guide tree and gap model 0 written explicitly
        (the same alignments as the plain dict; checked by the identity step)."""
        self.refp = _py(with_gap(FF.with_tree(v2_params(NUC_FAST_KEY), FF.TREE_400), 0))
        k = pkey(self.refp)
        self.log(f"reference 4ff7e4d1 (build {BUILD}) on screen ({len(self.screen)}) and core ({len(self.core)})")
        self.ref = {}
        for scope in ("screen", "core"):
            s = self.score({k: self.refp}, scope, raw=True)[k]
            if not s["valid"]:
                raise RuntimeError(f"reference failed on {s['n_err']} {scope} cases")
            self.ref[scope] = FF.objectives(s, s)
            self.summ[("new", scope, k)] = self.ref[scope]
            r = self.ref[scope]
            self.log(f"  {scope}: Q_small={r['q_small']:.4f} SP={r['sp_small']:.4f} TC={r['tc_small']:.4f} "
                     f"Q_large={r['q_large']:.4f} CPU small={r['cpu_small']:.0f}s large={r['cpu_large']:.0f}s")
        _atomic_json(self.out / "refs.json", {"ref_params": self.refp, "ref_key": k, "v2_key": NUC_FAST_KEY,
                                              "build": BUILD, "ref": self.ref})

    def run(self):
        self.write_run_json()
        self.log(f"fit start: core {len(self.core)}, screen {len(self.screen)}, registered={self.st.is_registered()}")
        self.references()
        archive = self.search()
        resc = self.rescore(archive)
        kref = pkey(self.refp)
        cands = [{"key": e["key"], "params": e["params"], "source": f"search gen {e['gen']}", "core": e["core"]}
                 for e in resc]
        refc = {"key": kref, "params": self.refp, "source": "reference 4ff7e4d1", "core": self.ref["core"]}
        if kref not in {c["key"] for c in cands}:
            cands.append(refc)
        pick, info = knee(cands, self.ref["core"], self.ev, self.meta, self.core)

        def brief(c):
            s = c["core"]
            return {"key": c["key"], "source": c["source"], "params": c["params"], "short": self.space.short(c["params"]),
                    "core": {k: s[k] for k in ("q_small", "sp_small", "tc_small", "q_large", "per_k", "per_benchmark",
                                               "cpu_small", "cpu_large", "cpu_ratio_small", "cpu_ratio_large", "R",
                                               "feasible")}}
        out = {"fold": self.fold, "registered": self.st.is_registered(), "rule": "knee (PLAN_nucfit.md section 4)",
               "n_candidates": len(cands), "knee": info,
               "ref_core": {k: self.ref["core"][k] for k in ("q_small", "sp_small", "tc_small", "q_large", "per_k",
                                                             "per_benchmark", "cpu_small", "cpu_large")},
               "fast": brief(pick or refc), "fast_is_fallback": pick is None,
               "candidates": [brief(c) for c in sorted(cands, key=lambda c: c["key"])]}
        _atomic_json(self.out / "selection.json", out)
        f = out["fast"]
        self.log(f"fast pick: {f['source']}{' (fallback)' if pick is None else ''} Q_small={f['core']['q_small']:.4f} "
                 f"Q_large={f['core']['q_large']:.4f} R={f['core']['R']:.3f} | {f['short']}")
        self.log("fit done")
        return out


def default_ws():
    return Path("~/tmp/kalign-nucfit").expanduser()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fold", required=True, help="0..4 or all")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ws", default=str(default_ws()))
    ap.add_argument("--threads", type=int, default=FF.SLOTS)
    ap.add_argument("--seed", type=int, default=FF.SEED)
    ap.add_argument("--random-init", action="store_true", help="no seeds (seed-robustness repeat)")
    ap.add_argument("--pop-size", type=int, default=FF.POP_SIZE)
    ap.add_argument("--ref-partitions", type=int, default=FF.REF_PARTITIONS)
    ap.add_argument("--gens", default=f"{FF.GENS[0]},{FF.GENS[1]}")
    ap.add_argument("--limit", type=int, default=0, help="smoke: at most N small and N/10 large core cases")
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    fold = a.fold if a.fold == "all" else int(a.fold)
    rows, meta, core, screen, test, lockbox = load_fold(ws, fold)
    note = []
    if a.limit:
        sm = sorted(c for c in core if meta[c]["kind"] == "small")[:a.limit]
        lg = sorted(c for c in core if meta[c]["kind"] == "large")[:max(3, a.limit // 10)]
        core = sm + lg
        screen = [c for c in screen if c in set(core)] or core[:10]
        note.append(f"limit={a.limit}")
    if a.seed != FF.SEED:
        note.append(f"seed={a.seed}")
    if a.random_init:
        note.append("random_init")
    st = NucSettings(pop_size=a.pop_size, ref_partitions=a.ref_partitions, seed=a.seed,
                     gens=tuple(int(v) for v in a.gens.split(",")), note=",".join(note))
    from .evaluate import Build
    b = Build.load(ws / "kalign" / BUILD)
    slots = SlotPool(ws / "runs" / "slots", FF.SLOTS)
    cache = ws / "cache" / "cache.sqlite"
    cache.parent.mkdir(parents=True, exist_ok=True)
    ev = Evaluator(b, rows.values(), set(core), cache_path=cache, slots=slots, threads=a.threads,
                   timeout_s=st.timeout_s, rerun_short_timeouts=True)   # A3: 600 s timeouts re-run at 1,800 s
    extra = {"build_new": b.build_id, "build_ref": b.build_id, "cache": str(cache), "threads": a.threads,
             "slots": FF.SLOTS}
    f = NucFastFit(fold, a.out, rows, meta, core, screen, test, lockbox, ev, ev, st, extra)
    if a.random_init:
        f.space.seeds = lambda: []
    ev.logfh = f.logfh
    f.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
