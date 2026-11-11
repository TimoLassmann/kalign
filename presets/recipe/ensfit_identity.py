"""Identity checks before the ensemble fit is launched (PLAN_ensemble.md section 11, task checks a-c, plus d).

    python -m recipe.ensfit_identity run [--ws WS]     -> WS/runs/ensfit/identity/{rows.jsonl,summary.json}

Cases (fit data only, training cases of fold "all", never a lockbox input): S20 = the first 4 cases by id of each of
the 5 benchmarks in the fold-"all" small screening subset; F2 = the first 2 families by id of the fold-"all" screening
families at 1,000.

  a  a 1-member run of the fast setting through the fit's evaluator (build 5b93950, member gap model) is byte-identical
     to `kalign --mode fast -n 1` of the f739560 release CLI (S20 + F2);
  b  the evaluator with the pre-release accurate parameters (2852d1c2) is byte-identical to `kalign --mode accurate -n 1`
     of the 5b93950 CLI (S20);
  c  per-member values: 5 members that are all the fast setting ("equal") against the same 5 members where only member
     k (2..5) has a different value of one setting ("member k differs"): sequence weights 0 (fast: 2.86), realign 2
     (fast: 0), consistency 8 rounds / weight 3.558 (fast: 0), VSM mode per-merge (fast: legacy). Two ensemble option
     sets: the pre-release accurate options (min_support 2, fill-in realign 0.0796) and min_support 5 with fill-in off.
     Outputs through the evaluator (S20 + F2); member scores and the winning member from ens_probe (the same call in C
     with logging on; its output must equal the evaluator's);
  d  seed 1 of the search (shipped tree as explicit values) is byte-identical to the same design with the shipped tree
     left at the engine default (-1) (S20 + F2).
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from . import ensfit as EF
from . import evaluate as E
from .fit import SlotPool, pkey

CLI_FAST = Path("~/tmp/kalign-rel/f739560/src/build/src/kalign").expanduser()
CLI_ACC = Path("~/tmp/kalign-rel/5b93950/src/build/src/kalign").expanduser()
PROBE_REL = Path("runs/ens_premeasure/ens_probe")
SETTINGS = {
    "seq_weights": {"seq_weights": 0.0},
    "realign": {"realign": 2},
    "consistency": {"consistency": 8, "consistency_weight": 3.5579641521153507},
    "vsm_mode": {"vsm_mode": 1},
}
ENS_SETS = {"acc_options": None, "ms5_nofill": {"min_support": 5, "consistency_merge": 0,
                                                "consistency_merge_weight": 2.0, "fill_in": 0,
                                                "fill_anchor_frac": 0.0, "fill_support": 1}}
EXPECT = {"MATRIX_PFASUM43": 1, "MATRIX_PFASUM60": 2, "MATRIX_CORBLOSUM66": 3, "REFINE_NONE": 0, "REFINE_ALL": 1,
          "REFINE_CONFIDENT": 2, "REFINE_INLINE": 3, "FILL_OFF": 0, "FILL_REALIGN": 1, "FILL_PROJECT": 2,
          "FILL_ATTACH": 3}


def cases(ws):
    from .fastfit import load_fold
    rows, meta, core, screen, test, lockbox = load_fold(ws, "all")
    small = sorted(c for c in screen if meta[c]["kind"] == "small")
    by_b = defaultdict(list)
    for c in small:
        by_b[meta[c]["benchmark"]].append(c)
    s20 = [c for b in sorted(by_b) for c in by_b[b][:4]]
    f2 = sorted(c for c in screen if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)[:2]
    ids = s20 + f2
    if set(ids) & set(lockbox) or not set(ids) <= set(core):
        raise EF.LeakError("identity cases must be fit (training) cases")
    return rows, meta, s20, f2


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else None


def probe_config(p, seed=42):
    """ens_probe CONFIG from a worker dict (same values the Python path receives)."""
    kw = p
    L = [f"ens {kw['min_support']} {kw['consistency_merge']} {float(kw['consistency_merge_weight']):.17g} "
         f"{kw['fill_in']} {float(kw['fill_anchor_frac']):.17g} {kw['fill_support']} {seed}"]
    n = kw["n_runs"]
    g = lambda k, i, d: kw[k][i] if k in kw else d  # noqa: E731
    for i in range(n):
        vals = [kw["run_types"][i], kw["run_gpo"][i], kw["run_gpe"][i], kw["run_tgpe"][i], kw["run_noise"][i],
                kw["run_vsm_amax"][i], g("run_seq_weights", i, kw["seq_weights"]), kw["run_refine"][i],
                g("run_realign", i, kw["realign"]), g("run_consistency_anchors", i, kw["consistency"]),
                g("run_consistency_weight", i, kw["consistency_weight"]), g("run_vsm_mode", i, kw["vsm_mode"]),
                kw["run_vsm_kappa"][i], kw["run_vsm_dzero"][i], kw["run_guide_dist"][i], g("run_gap_model", i, 0)] + \
            [g(k, i, -1) for k in ("run_tree_method", "run_tree_dist", "run_tree_alpha", "run_tree_maxlen",
                                   "run_tree_leaf", "run_tree_kmeans_dist", "run_realign_tree")]
        L.append("run " + " ".join(f"{float(v):.17g}" if isinstance(v, float) else str(int(v)) for v in vals))
    return "\n".join(L) + "\n"


def run_probe(probe, p, inp, tmp):
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        cf, out = Path(d) / "probe.cfg", Path(d) / "c.fa"
        cf.write_text(probe_config(p))
        r = subprocess.run([str(probe), str(cf), inp, str(out)], capture_output=True, text=True,
                           env={**os.environ, **E.SINGLE_THREAD_ENV})
        scores = [float(v) for _, v in re.findall(r"Run (\d+) score: ([-\d.]+)", r.stderr)]
        sel = re.findall(r"Selected run (\d+)", r.stderr)
        path = [l.split(":", 3)[-1].strip() for l in r.stderr.splitlines()
                if "consensus" in l.lower() or "Keeping" in l or "Refin" in l or "Filling" in l]
        return {"scores": scores, "selected": int(sel[0]) if sel else None, "sha_c": sha(out), "path": path,
                "rc": r.returncode, "err": r.stderr[-300:] if r.returncode else None}


def five_fast(fast, over_k=None, k=None, ens=None, acc=None):
    m1 = EF.fast_member(fast)
    ms = [dict(m1) for _ in range(EF.K)]
    if over_k:
        ms[k - 1].update(over_k)
    e = {kk: acc[kk] for kk in EF.ENS_KEYS} if ens is None else ens
    return EF.worker_dict(ms, e)


def run(ws):
    ws = Path(ws)
    out = ws / "runs" / "ensfit" / "identity"
    out.mkdir(parents=True, exist_ok=True)
    tmp = ws / "tmp"
    rows, meta, s20, f2 = cases(ws)
    ids = s20 + f2
    ev, b, cache = EF.make_evaluator(ws, rows, ids, threads=EF.SLOTS)
    slots = SlotPool(ws / "runs" / "slots", EF.SLOTS)
    fast, acc = EF.fast_setting(ws), EF.prerelease_accurate()
    probe = ws / PROBE_REL
    res = {"build": b.build_id, "cli_fast": str(CLI_FAST), "cli_fast_sha256": EF.sha256_file(CLI_FAST),
           "cli_acc": str(CLI_ACC), "cli_acc_sha256": EF.sha256_file(CLI_ACC), "probe": str(probe),
           "probe_sha256": EF.sha256_file(probe), "cases_small": s20, "cases_1k": f2,
           "started": time.strftime("%Y-%m-%dT%H:%M:%S")}
    if EF.sha256_file(CLI_ACC) != b.build_id:
        raise SystemExit("5b93950 CLI differs from the fit build")
    # aligner constants of the fit's _core
    q = subprocess.run([b.python, "-c", "import json; from kalign import _core; print(json.dumps({k: int(getattr("
                        "_core, k)) for k in %r if hasattr(_core, k)}))" % list(EXPECT)],
                       capture_output=True, text=True)
    const = json.loads(q.stdout)
    res["constants"] = {"found": const, "mismatch": {k: (v, EXPECT[k]) for k, v in const.items() if v != EXPECT[k]}}
    fh = open(out / "rows.jsonl", "w")

    def emit(r):
        fh.write(json.dumps(r) + "\n")
        fh.flush()

    def cli_sha(cli, mode, cid):
        with slots.slot():
            with tempfile.TemporaryDirectory(dir=tmp) as d:
                o = Path(d) / "cli.fa"
                wall, cpu, err = E.run_cli(str(cli), mode, rows[cid]["unaligned"], o, EF.TIMEOUT_S)
                return sha(o), err

    # ---- a, b
    kf, ka = pkey(fast), pkey(acc)
    got = ev.evaluate({kf: fast, ka: acc}, ids)
    with ThreadPoolExecutor(24) as ex:
        fa = {c: ex.submit(cli_sha, CLI_FAST, "fast", c) for c in ids}
        fb = {c: ex.submit(cli_sha, CLI_ACC, "accurate", c) for c in s20}
        for check, key, futs in (("a", kf, fa), ("b", ka, fb)):
            eqs = {}
            for c, f in futs.items():
                s_cli, err = f.result()
                s_ev = got[key][c].get("aln_sha256")
                eqs[c] = s_cli is not None and s_cli == s_ev
                emit({"check": check, "case_id": c, "sha_evaluator": s_ev, "sha_cli": s_cli, "identical": eqs[c],
                      "err": err or got[key][c].get("error")})
            res[check] = {"identical": sum(eqs.values()), "n": len(eqs),
                          "identical_small": sum(v for c, v in eqs.items() if c in s20),
                          "n_small": sum(1 for c in eqs if c in s20),
                          "identical_1k": sum(v for c, v in eqs.items() if c in f2),
                          "n_1k": sum(1 for c in eqs if c in f2)}
            print(f"{check}: {sum(eqs.values())}/{len(eqs)} byte-identical", flush=True)

    # ---- d: seed 1 (explicit shipped tree) == shipped tree at the engine default
    sp = EF.Space(fast, acc)
    seed1 = sp.seed_params()[0]
    dflt = copy.deepcopy(seed1)
    for kk in ("run_tree_method", "run_tree_dist", "run_tree_alpha", "run_tree_maxlen", "run_tree_leaf",
               "run_tree_kmeans_dist", "run_realign_tree"):
        dflt[kk] = [seed1[kk][0]] + [-1] * (EF.K - 1)
    g = ev.evaluate({pkey(seed1): seed1, pkey(dflt): dflt}, ids)
    n_eq = sum(g[pkey(seed1)][c].get("aln_sha256") == g[pkey(dflt)][c].get("aln_sha256")
               and g[pkey(seed1)][c].get("aln_sha256") is not None for c in ids)
    for c in ids:
        emit({"check": "d", "case_id": c, "sha_seed1": g[pkey(seed1)][c].get("aln_sha256"),
              "sha_default_tree": g[pkey(dflt)][c].get("aln_sha256")})
    res["d"] = {"identical": n_eq, "n": len(ids)}
    print(f"d: {n_eq}/{len(ids)} byte-identical", flush=True)

    # ---- c: per-member values
    specs, meta_c = {}, {}
    for es, ens in ENS_SETS.items():
        base = five_fast(fast, ens=ens, acc=acc)
        specs[pkey(base)] = base
        meta_c[pkey(base)] = (es, "equal", None, None)
        for s, over in SETTINGS.items():
            for k in range(2, EF.K + 1):
                p = five_fast(fast, over, k, ens=ens, acc=acc)
                specs[pkey(p)] = p
                meta_c[pkey(p)] = (es, s, k, over)
    got = ev.evaluate(specs, ids)
    jobs = [(key, c) for key in specs for c in ids]
    with ThreadPoolExecutor(EF.SLOTS) as ex:
        def pj(job):
            key, c = job
            with slots.slot():
                return job, run_probe(probe, specs[key], rows[c]["unaligned"], tmp)
        probes = dict(ex.map(pj, jobs))
    summ = {}
    for es in ENS_SETS:
        kb = next(k for k, m in meta_c.items() if m[0] == es and m[1] == "equal")
        base_sha = {c: got[kb][c].get("aln_sha256") for c in ids}
        for s in SETTINGS:
            d = {"n": 0, "output_differs": 0, "member_k_score_only": 0, "equal_scores_in_equal_run": 0,
                 "py_eq_c": 0, "winner": Counter(), "per_k": {}}
            for k in range(2, EF.K + 1):
                key = next(kk for kk, m in meta_c.items() if m[:3] == (es, s, k))
                pk = {"output_differs": 0, "n": 0}
                for c in ids:
                    pr, pb = probes[(key, c)], probes[(kb, c)]
                    sh = got[key][c].get("aln_sha256")
                    differs = sh != base_sha[c]
                    sc = pr["scores"]
                    others = [v for i, v in enumerate(sc) if i != k - 1]
                    konly = len(sc) == EF.K and len(set(others)) == 1 and sc[k - 1] != others[0]
                    d["n"] += 1
                    d["output_differs"] += differs
                    d["member_k_score_only"] += konly
                    d["equal_scores_in_equal_run"] += len(set(pb["scores"])) == 1
                    d["py_eq_c"] += (pr["sha_c"] == sh) + 0
                    d["winner"][str(pr["selected"])] += 1
                    pk["n"] += 1
                    pk["output_differs"] += differs
                    emit({"check": "c", "options": es, "setting": s, "member": k, "case_id": c,
                          "value": SETTINGS[s], "sha_evaluator": sh, "sha_equal_run": base_sha[c],
                          "output_differs": differs, "member_scores": sc, "equal_run_scores": pb["scores"],
                          "selected": pr["selected"], "equal_run_selected": pb["selected"],
                          "probe_eq_evaluator": pr["sha_c"] == sh, "path": pr["path"],
                          "err": pr["err"] or got[key][c].get("error")})
                d["per_k"][str(k)] = pk
            d["winner"] = dict(d["winner"])
            summ.setdefault(es, {})[s] = d
            print(f"c {es:11s} {s:12s} output differs {d['output_differs']}/{d['n']}, only member k's score differs "
                  f"{d['member_k_score_only']}/{d['n']}, probe == evaluator {d['py_eq_c']}/{d['n']}, winners "
                  f"{d['winner']}", flush=True)
    res["c"] = summ
    res["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    fh.close()
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    return run(Path(a.ws).expanduser())


if __name__ == "__main__":
    sys.exit(main())
