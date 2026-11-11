"""Step 1 of the nucfit plan's accurate setting (presets/docs/nucfit/PLAN_nucfit.md, section 5): 5-member nucleotide
ensemble, member 1 = the fold's own nucleotide fast pick (runs/nucfast/<fit>/selection.json), members 2-5 and the
ensemble options searched; NSGA-III, re-scoring, knee selection with precision.

    python -m recipe.nucens --fold 0 --out WS/runs/nucens/fold0 [--ws WS]          (0..4 or all)

Everything not listed here is recipe.ensfit unchanged (the protein ensemble fit: NSGA-III settings, stop rule,
precision objectives, SP floors, re-scoring rule, evaluator with CPU limits, resumable).

Differences from recipe.ensfit:
  * data: recipe.nucfast.load_fold; inputs above 5,000 sequences are left out (accurate refuses them, as protein),
    so large = Rfam large families + simulations at 1,000; the large precision / SP / TC are the means of the two
    category means (fields *_1k of recipe.ensfit carry the large values; *_10k are unused).
  * member 1 = the fast pick of the same fit (fold k uses fold k's pick; protein used the final pick in every fold);
    every member uses the member gap model (1) and member 1's realign-tree flag.
  * search space per member 2-5: the pinned nucleotide ranges (gpo [1, 20], gpe [0.2, 5], tgpe [0.05, 3], matrices
    Kimura 1 / 20 / 200 PAM, vsm_amax [0, 3], seq_weights [0, 3], consistency_weight [0.5, 5], refine, realign 0-2,
    consistency anchors) + member noise [0, 0.5] + the protein ensemble's guide-tree variables; vsm_mode 0 (the
    nucleotide space has no other); ensemble options as protein.
  * seeds (the pre-release nucleotide accurate preset is a single run): (A) members 2-5 = the pre-release nucleotide
    accurate setting; (B) members 2-5 = accurate, fast, accurate, fast (pre-release nucleotide); both with member 1's
    guide tree and the pre-release ensemble options (min_support 0, no consistency merge, no fill-in).
  * CPU axis normalised over R in [1, 100] (no observed range exists before the fit).
  * knee: bootstrap stratified by large category (as recipe.nucfast), on precision.
  * build 9137658 (the release build), cache WS/cache/nucens.sqlite.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from . import cache as C
from . import ensfit as EF
from . import fit as FIT
from . import nucfast as NF
from .fit import SlotPool, _atomic_json, _py, pkey

BUILD = NF.BUILD
MAX_SEQS = 5000
LARGE_K = ("rfam", "1000")
MATRIX_MAP = [10, 11, 12]                    # Kimura 1, 20, 200 PAM
MATRIX_NAMES = {10: "K1", 11: "K20", 12: "K200"}
R_LO, R_HI = 1.0, 100.0
N_BOOT, Z = 2000, 1.96
GAP_MODEL = 1


def nucfast_dir(fit):
    return "final" if str(fit) == "all" else f"fold{fit}"


def fast_setting(ws, fit):
    """The nucleotide fast pick of the same fit, with the member gap model (member 1 of every ensemble)."""
    sel = json.loads((Path(ws) / "runs" / "nucfast" / nucfast_dir(fit) / "selection.json").read_text())
    p = copy.deepcopy(sel["fast"]["params"])
    if pkey(p) != sel["fast"]["key"]:
        raise SystemExit(f"nucfast {fit}: pick params do not hash to their key")
    p["run_gap_model"] = [GAP_MODEL]
    return _py(p), sel["fast"]["key"]


def fast_member(fast):
    m = EF.fast_member({**fast, "run_vsm_kappa": fast.get("run_vsm_kappa", [1.0]),
                        "run_vsm_dzero": fast.get("run_vsm_dzero", [1.0])})
    return m


def single_member(p, tree):
    """A single-run pre-release nucleotide preset as an ensemble member with the given tree."""
    return {"gpo": p["run_gpo"][0], "gpe": p["run_gpe"][0], "tgpe": p["run_tgpe"][0], "matrix": p["run_types"][0],
            "vsm_amax": p["run_vsm_amax"][0], "vsm_kappa": 1.0, "vsm_dzero": 1.0, "vsm_mode": 0,
            "seq_weights": p["seq_weights"], "realign": p["realign"], "consistency": p["consistency"],
            "consistency_weight": p["consistency_weight"], "refine": p["run_refine"][0], "noise": 0.0,
            "tree_method": tree["tree_method"], "tree_leaf": tree["tree_leaf"], "tree_dist": tree["tree_dist"],
            "tree_alpha": tree["tree_alpha"], "tree_kmeans_dist": tree["tree_kmeans_dist"]}


class NucEnsSpace:
    MEMBER_VARS = ("gpo", "gpe", "tgpe", "matrix", "vsm_amax", "seq_weights", "realign", "consistency",
                   "consistency_weight", "refine", "noise", "tree_method", "tree_leaf_log", "tree_dist", "tree_alpha",
                   "tree_kmeans_dist")

    def __init__(self, fast):
        self.fast = fast
        self.m1 = fast_member(fast)
        self.rt = int(fast["run_realign_tree"][0])

    def wd(self, members, ens):
        p = EF.worker_dict(members, ens)
        p["run_realign_tree"] = [self.rt] * p["n_runs"]      # member 1's flag in every member
        p["run_gap_model"] = [GAP_MODEL] * p["n_runs"]
        p["run_tree_maxlen"] = [int(self.fast["run_tree_maxlen"][0])] * p["n_runs"]
        return _py(p)

    def vars(self):
        from pymoo.core.variable import Choice, Integer, Real
        v = {}
        for j in range(2, EF.K + 1):
            v.update({f"gpo_{j}": Real(bounds=(1.0, 20.0)), f"gpe_{j}": Real(bounds=(0.2, 5.0)),
                      f"tgpe_{j}": Real(bounds=(0.05, 3.0)), f"matrix_{j}": Choice(options=list(range(3))),
                      f"vsm_amax_{j}": Real(bounds=(0.0, 3.0)), f"seq_weights_{j}": Real(bounds=(0.0, 3.0)),
                      f"realign_{j}": Integer(bounds=(0, 2)),
                      f"consistency_{j}": Choice(options=list(EF.CONSISTENCY_MAP)),
                      f"consistency_weight_{j}": Real(bounds=(0.5, 5.0)),
                      f"refine_{j}": Choice(options=list(EF.REFINE_MAP)), f"noise_{j}": Real(bounds=(0.0, 0.5)),
                      f"tree_method_{j}": Choice(options=list(EF.TREE_METHODS)),
                      f"tree_leaf_log_{j}": Real(bounds=EF.LEAF_LOG),
                      f"tree_dist_{j}": Choice(options=list(EF.TREE_DISTS)), f"tree_alpha_{j}": Choice(options=[0, 1]),
                      f"tree_kmeans_dist_{j}": Choice(options=[0, 1])})
        v.update({"min_support": Integer(bounds=(0, EF.K)), "consistency_merge": Choice(options=[0, 1]),
                  "consistency_merge_weight": Real(bounds=(0.5, 10.0)), "fill_in": Choice(options=list(EF.FILL_MAP)),
                  "fill_anchor_frac": Real(bounds=(0.0, 1.0)), "fill_support": Integer(bounds=(1, EF.K))})
        return v

    def decode(self, x):
        ms = [self.m1]
        for j in range(2, EF.K + 1):
            ms.append({"gpo": float(x[f"gpo_{j}"]), "gpe": float(x[f"gpe_{j}"]), "tgpe": float(x[f"tgpe_{j}"]),
                       "matrix": MATRIX_MAP[int(x[f"matrix_{j}"])], "vsm_amax": float(x[f"vsm_amax_{j}"]),
                       "vsm_kappa": 1.0, "vsm_dzero": 1.0, "vsm_mode": 0,
                       "seq_weights": float(x[f"seq_weights_{j}"]), "realign": int(x[f"realign_{j}"]),
                       "consistency": int(x[f"consistency_{j}"]),
                       "consistency_weight": float(x[f"consistency_weight_{j}"]), "refine": int(x[f"refine_{j}"]),
                       "noise": float(x[f"noise_{j}"]), "tree_method": int(x[f"tree_method_{j}"]),
                       "tree_leaf": int(round(10 ** float(x[f"tree_leaf_log_{j}"]))),
                       "tree_dist": int(x[f"tree_dist_{j}"]), "tree_alpha": int(x[f"tree_alpha_{j}"]),
                       "tree_kmeans_dist": int(x[f"tree_kmeans_dist_{j}"])})
        p = self.wd(ms, {k: x[k] for k in EF.ENS_KEYS})
        if p["run_gpo"][0] != self.fast["run_gpo"][0] or p["run_types"][0] != self.fast["run_types"][0]:
            raise RuntimeError("member 1 is not the fast setting")
        return p

    def encode(self, p):
        ms, ens = EF.members_of(p)
        x = {}
        for j in range(2, EF.K + 1):
            m = ms[j - 1]
            x.update({f"gpo_{j}": m["gpo"], f"gpe_{j}": m["gpe"], f"tgpe_{j}": m["tgpe"],
                      f"matrix_{j}": MATRIX_MAP.index(m["matrix"]), f"vsm_amax_{j}": m["vsm_amax"],
                      f"seq_weights_{j}": m["seq_weights"], f"realign_{j}": m["realign"],
                      f"consistency_{j}": m["consistency"], f"consistency_weight_{j}": m["consistency_weight"],
                      f"refine_{j}": m["refine"], f"noise_{j}": m["noise"], f"tree_method_{j}": m["tree_method"],
                      f"tree_leaf_log_{j}": math.log10(m["tree_leaf"]), f"tree_dist_{j}": m["tree_dist"],
                      f"tree_alpha_{j}": m["tree_alpha"], f"tree_kmeans_dist_{j}": m["tree_kmeans_dist"]})
        x.update(ens)
        return x

    def seed_params(self):
        acc, fast0 = NF.v2_params(NF.NUC_ACC_KEY), NF.v2_params(NF.NUC_FAST_KEY)
        ens = {k: acc[k] for k in EF.ENS_KEYS}
        a, f = single_member(acc, self.m1), single_member(fast0, self.m1)
        return [self.wd([self.m1, a, a, a, a], ens), self.wd([self.m1, a, f, a, f], ens)]

    @staticmethod
    def short(p):
        ms, e = EF.members_of(p)
        tn = {0: "km", 3: "km+Prim"}
        parts = [f"m{j}:{MATRIX_NAMES[m['matrix']]} {m['gpo']:.2f}/{m['gpe']:.2f}/{m['tgpe']:.2f} "
                 f"vsm={m['vsm_amax']:.2f} sw={m['seq_weights']:.2f} re={m['realign']} c={m['consistency']} "
                 f"ref={m['refine']} nz={m['noise']:.2f} {tn[m['tree_method']]} leaf={m['tree_leaf']} "
                 f"d={m['tree_dist']} a={m['tree_alpha']} kd={m['tree_kmeans_dist']}"
                 for j, m in enumerate(ms[1:], 2)]
        return (f"ms={e['min_support']} cm={e['consistency_merge']} fill={e['fill_in']} "
                f"frac={e['fill_anchor_frac']:.3f} fs={e['fill_support']} | " + " | ".join(parts))


# ------------------------------------------------------------------------------------ objectives
def summarise(rows, meta, guard=None):
    """Small: benchmark-weighted means; large: mean of the category means (Rfam large, simulations at 1,000), carried
    in the *_1k fields of recipe.ensfit."""
    small, large = {}, {}
    cpu = {"small": 0.0, "large": 0.0}
    n_err = {"small": 0, "large": 0}
    for cid, r in rows.items():
        m = meta[cid]
        part = m["kind"]
        cpu[part] += float(r.get("cpu_s") or 0.0)
        if r.get("error"):
            n_err[part] += 1
            continue
        v = (float(r["precision"]), float(r["recall"]), float(r["tc"]), float(r["width_ratio"]))
        (small.setdefault(m["benchmark"], []) if part == "small" else large.setdefault(m["k"], [])).append(v)

    def means(d):
        return {b: {"p": float(np.mean([v[0] for v in vs])), "sp": float(np.mean([v[1] for v in vs])),
                    "tc": float(np.mean([v[2] for v in vs])), "wr": float(np.mean([v[3] for v in vs])), "n": len(vs)}
                for b, vs in sorted(d.items())}
    pb, pk = means(small), means(large)

    def avg(d, f):
        return float(np.mean([v[f] for v in d.values()])) if d else float("nan")
    ne = n_err["small"] + n_err["large"]
    return {"p_small": avg(pb, "p"), "sp_small": avg(pb, "sp"), "tc_small": avg(pb, "tc"), "wr_small": avg(pb, "wr"),
            "p_1k": avg(pk, "p"), "sp_1k": avg(pk, "sp"), "tc_1k": avg(pk, "tc"), "wr_1k": avg(pk, "wr"),
            "n_1k": sum(v["n"] for v in pk.values()),
            "p_10k": float("nan"), "sp_10k": float("nan"), "tc_10k": float("nan"), "wr_10k": float("nan"),
            "n_10k": 0, "cpu_small": cpu["small"], "cpu_1k": cpu["large"], "cpu_10k": 0.0,
            "guard": None, "guard_cpu": None, "guard_err": None,
            "n": len(rows), "n_err": ne, "n_err_small": n_err["small"], "n_err_1k": n_err["large"], "n_err_10k": 0,
            "valid": ne == 0 and len(rows) > 0, "per_benchmark": pb, "per_k": pk}


class Normaliser(EF.Normaliser):
    def F(self, s):
        return [(1.0 - s["p_small"]) / (1.0 - self.lo_s), (1.0 - s["p_1k"]) / (1.0 - self.lo_l),
                (math.log(max(s["R"], 1e-9)) - math.log(R_LO)) / (math.log(R_HI) - math.log(R_LO))]

    def to_json(self):
        return {"p_small_lo": self.lo_s, "p_large_lo": self.lo_l, "R_lo": R_LO, "R_hi": R_HI}


def eligible(c, ref):
    s = c["core"]
    return bool(s["feasible"] and s["p_small"] > ref["p_small"] and s["sp_small"] >= ref["sp_small"] - EF.SP_FLOOR
                and s["sp_1k"] >= ref["sp_1k"] - EF.SP_FLOOR)


def knee(cands, ref, ev, meta, core):
    """Knee rule on large precision: b = highest P_large among eligible; keep points within 1.96 paired bootstrap SE
    (2,000 resamples stratified by category, seed 1234); pick the lowest R."""
    E = sorted((c for c in cands if eligible(c, ref)), key=lambda c: c["key"])
    out = {"n_candidates": len(cands), "n_eligible": len(E), "eligible": [c["key"] for c in E]}
    if not E:
        out.update({"pick": None, "b": None, "within": [], "rows": []})
        return out
    cases = sorted(c for c in core if meta[c]["kind"] == "large")
    cats = sorted({meta[c]["k"] for c in cases})
    idx = {k: np.array([i for i, c in enumerate(cases) if meta[c]["k"] == k]) for k in cats}
    rng = np.random.default_rng(EF.BOOT_SEED)
    draws = [{k: v[rng.integers(0, len(v), len(v))] for k, v in idx.items()} for _ in range(N_BOOT)]

    def vals(c):
        got = ev.lookup(c["key"], cases)
        if set(got) != set(cases) or any(got[x].get("error") for x in cases):
            raise RuntimeError(f"{c['key'][:12]}: incomplete large core rows")
        return np.array([float(got[x]["precision"]) for x in cases])

    def boot(v):
        return np.array([np.mean([v[d[k]].mean() for k in cats]) for d in draws])
    b = sorted(E, key=lambda c: (-c["core"]["p_1k"], c["core"]["R"], c["key"]))[0]
    vb = vals(b)
    full = float(np.mean([vb[idx[k]].mean() for k in cats]))
    if abs(full - b["core"]["p_1k"]) > 1e-9:
        raise RuntimeError(f"recomputed P_large {full:.6f} != core summary {b['core']['p_1k']:.6f}")
    bb = boot(vb)
    rows = []
    for c in E:
        if c["key"] == b["key"]:
            rows.append({"key": c["key"], "diff": 0.0, "se": 0.0, "within": True})
            continue
        se = float(np.std(bb - boot(vals(c)), ddof=1))
        diff = b["core"]["p_1k"] - c["core"]["p_1k"]
        rows.append({"key": c["key"], "diff": diff, "se": se, "within": bool(diff <= Z * se)})
    by = {c["key"]: c for c in E}
    within = [r["key"] for r in rows if r["within"]]
    pick = min((by[k] for k in within), key=lambda c: (c["core"]["R"], -c["core"]["p_1k"], c["key"]))
    out.update({"b": b["key"], "within": within, "rows": rows, "pick": pick["key"]})
    return out


# ------------------------------------------------------------------------------------ data
def load_fold(ws, fold):
    rows, meta, core, screen, test, lockbox = NF.load_fold(ws, fold)

    def keep(c):
        m = meta[c]
        return m["kind"] == "small" or (m["k"] in LARGE_K and int(rows[c]["n_seq"]) <= MAX_SEQS)
    return rows, meta, [c for c in core if keep(c)], [c for c in screen if keep(c)], test, lockbox


# ------------------------------------------------------------------------------------ the fit
class NucEnsFit(EF.EnsFit):
    def __init__(self, fold, out_dir, rows, meta, core, screen, test, lockbox, ev, settings, extra, fast, fast_key):
        acc = NF.v2_params(NF.NUC_ACC_KEY)
        super().__init__(fold, out_dir, rows, meta, core, screen, test, lockbox, None, ev, settings, extra, fast, acc)
        self.space = NucEnsSpace(fast)
        self.fast_key = fast_key
        self.rows = rows

    def log(self, msg):
        FIT.log(f"[nucens fold {self.fold}] {msg}", self.logfh)

    def score(self, specs, scope, capped=False, raw=False):
        ids = self._ids(scope)
        need = {k: p for k, p in specs.items() if (scope, k) not in self.summ}
        if need:
            res = self.ev.evaluate(need, ids, **({"stop_on_error": True} if capped else {}))
            for k, rows in res.items():
                s = summarise(rows, self.meta)
                self.summ[(scope, k)] = s if raw else EF.objectives(s, self.ref[scope])
        return {k: self.summ[(scope, k)] for k in specs}

    def write_run_json(self):
        import pymoo
        if self.st.is_registered() and pymoo.__version__ != EF.PYMOO_VERSION:
            raise SystemExit(f"registered fit needs pymoo {EF.PYMOO_VERSION}, found {pymoo.__version__}")
        meta = {"plan": "nucfit/PLAN_nucfit.md section 5 step 1", "fold": self.fold, "settings": self.st.to_json(),
                "registered": self.st.is_registered(), "pymoo": pymoo.__version__,
                "n_core": len(self.core), "n_screen": len(self.screen),
                "core_sha256": FIT.ids_hash(self.core), "screen_sha256": FIT.ids_hash(self.screen),
                "fast_pick": self.fast_key, "key_member1": pkey(self.fast), "nucens_sha256": C.sha256_file(__file__),
                **self.extra}
        path = self.out / "run.json"
        if path.exists():
            old = json.loads(path.read_text())
            for k in ("fold", "settings", "core_sha256", "screen_sha256", "build", "fast_pick"):
                if old.get(k) != meta.get(k):
                    raise SystemExit(f"{path}: {k} differs from the recorded run; refusing to resume")
            meta["started"] = old.get("started")
            meta["resumed"] = old.get("resumed", []) + [time.strftime("%Y-%m-%dT%H:%M:%S")]
        else:
            meta["started"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        _atomic_json(path, meta)
        (self.out / "cases_core.txt").write_text("".join(i + "\n" for i in self.core))
        (self.out / "cases_screen.txt").write_text("".join(i + "\n" for i in self.screen))

    def search(self):
        orig = EF.Normaliser
        EF.Normaliser = Normaliser           # recipe.ensfit.search builds its normaliser by module name
        try:
            return super().search()
        finally:
            EF.Normaliser = orig

    def run(self):
        self.write_run_json()
        self.log(f"fit start: core {len(self.core)}, screen {len(self.screen)}, member 1 = nucfast pick "
                 f"{self.fast_key[:8]}, registered={self.st.is_registered()}")
        self.references()
        archive = self.search()
        resc = self.rescore(archive)
        self.log("pre-release nucleotide accurate preset (reported row, single run)")
        self.acc_reference()
        cands = [{"key": e["key"], "params": e["params"], "source": f"search gen {e['gen']}", "core": e["core"]}
                 for e in resc]
        k = knee(cands, self.ref["core"], self.ev, self.meta, self.core)
        by = {c["key"]: c for c in cands}

        def brief(key, why):
            if key is None:
                return None
            c = by[key]
            return {"key": key, "why": why, "source": c["source"], "params": c["params"],
                    "params_sha256": hashlib.sha256(C.canonical_json(c["params"]).encode()).hexdigest(),
                    "short": self.space.short(c["params"]), "core": EF.brief_s(c["core"]) | {
                        "per_k": c["core"]["per_k"], "per_benchmark": c["core"]["per_benchmark"]}}
        out = {"fold": self.fold, "registered": self.st.is_registered(), "rule": "PLAN_nucfit.md section 5 step 1",
               "fast_pick": self.fast_key, "ref_core": EF.brief_s(self.ref["core"]),
               "acc_core": EF.brief_s(self.acc_s["core"]), "n_candidates": len(cands), "n_eligible": k["n_eligible"],
               "pick": brief(k["pick"], f"cheapest of {len(k['within'])} within 1.96 SE of the best"),
               "b": brief(k["b"], "eligible point with the highest P_large"), "knee_rows": k["rows"],
               "eligible": k["eligible"],
               "candidates": [{"key": c["key"], "source": c["source"], "eligible": c["key"] in k["eligible"],
                               "core": EF.brief_s(c["core"])} for c in sorted(cands, key=lambda c: c["key"])],
               "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
        _atomic_json(self.out / "selection.json", out)
        if out["pick"]:
            pk = out["pick"]["core"]
            self.log(f"pick: {out['pick']['source']} P_small={pk['p_small']:.4f} P_large={pk['p_1k']:.4f} "
                     f"SP_small={pk['sp_small']:.4f} SP_large={pk['sp_1k']:.4f} R={pk['R']:.2f} | {out['pick']['short']}")
        else:
            self.log("no eligible point: no accurate setting from this fit")
        self.log("fit done")
        return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fold", required=True, help="0..4 or all")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ws", default=str(NF.default_ws()))
    ap.add_argument("--threads", type=int, default=EF.SLOTS)
    ap.add_argument("--pop-size", type=int, default=EF.POP_SIZE)
    ap.add_argument("--gens", default=f"{EF.GENS[0]},{EF.GENS[1]}")
    ap.add_argument("--note", default="")
    ap.add_argument("--limit", type=int, default=0, help="smoke: at most N small and N/10 large core cases")
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    fold = a.fold if a.fold == "all" else int(a.fold)
    rows, meta, core, screen, test, lockbox = load_fold(ws, fold)
    note = [a.note] if a.note else []
    if a.limit:
        sm = sorted(c for c in core if meta[c]["kind"] == "small")[:a.limit]
        lg = sorted(c for c in core if meta[c]["kind"] == "large")[:max(2, a.limit // 10)]
        core = sm + lg
        screen = [c for c in screen if c in set(core)] or core[:10]
        note.append(f"limit={a.limit}")
    st = EF.Settings(pop_size=a.pop_size, gens=tuple(int(v) for v in a.gens.split(",")), note=",".join(note))
    from .evaluate import Build
    EF.use_ws_tmp(ws)
    b = Build.load(ws / "kalign" / BUILD)
    slots = SlotPool(ws / "runs" / "slots", EF.SLOTS)
    cache = ws / "cache" / "nucens.sqlite"
    ev = EF.EnsEvaluator(b, rows.values(), set(core), cache_path=cache, slots=slots, threads=a.threads,
                         timeout_s=st.timeout_s)
    fast, fkey = fast_setting(ws, fold)
    extra = {"build": b.build_id, "build_cli": b.cli, "build_python": b.python, "cache": str(cache),
             "threads": a.threads, "slots": EF.SLOTS}
    f = NucEnsFit(fold, a.out, rows, meta, core, screen, test, lockbox, ev, st, extra, fast, fkey)
    ev.logfh = f.logfh
    f.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
