"""Align cases with one parameter dict (runs in the ALIGNER venv; imports only kalign).

Two modes:

1. one-shot (per-case worker):  <venv>/bin/python _align_worker.py JOB.json RESULT.json
   JOB.json: {"params": <decoded params>, "cases": [[input, output], ...], "seed": 42}
   RESULT.json: {"aligner": path, "rows": [{"wall": s, "cpu": s, "error": str|null}, ...]}

2. long-lived (batched worker):  <venv>/bin/python _align_worker.py --serve PARAMS.json
   PARAMS.json: {"params": ..., "seed": 42, "expect": {...}}. The process prints {"ready": true} once
   the aligner is imported, then answers one JSON line per request line on stdin:
     request  {"in": path, "out": path}      response {"wall": s, "cpu": s, "error": str|null}
   and exits at EOF of stdin. The protocol runs on a private copy of the original stdout; fd 1 is
   redirected to stderr so that anything the C library prints cannot corrupt it.

Ported from experiments/redesign/tools/refit/align_worker.py. The decoded-params layout is the
one of experiments/redesign/tools/refit/space.py decode(): run_gpo, run_gpe, run_tgpe, run_noise,
run_types, run_vsm_amax, run_refine, run_vsm_kappa, run_vsm_dzero, run_guide_dist (lists, one entry
per run) and min_support, realign, seq_weights, consistency, consistency_weight, consistency_merge,
consistency_merge_weight, vsm_mode, fill_in, fill_anchor_frac, fill_support (scalars).
CPU = user+sys of this process (RUSAGE_SELF delta around the single aligner call), n_threads=1.
"""
import json
import os
import resource
import sys
import time

for _k, _v in {"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1", "OMP_DYNAMIC": "FALSE",
               "OMP_WAIT_POLICY": "PASSIVE"}.items():
    os.environ[_k] = _v

from kalign import _core  # noqa: E402


TREE_KEYS = ("run_tree_method", "run_tree_dist", "run_tree_alpha", "run_tree_maxlen", "run_tree_leaf",
             "run_tree_kmeans_dist", "run_realign_tree")
MEMBER_KEYS = ("run_gap_model", "run_seq_weights", "run_realign", "run_consistency_anchors",
               "run_consistency_weight", "run_vsm_mode")


def cpu_self():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def build_kwargs(job):
    for name, val in job.get("expect", {}).items():
        if int(getattr(_core, name)) != int(val):
            raise SystemExit(f"constant mismatch {name}: aligner {getattr(_core, name)} != {val}")
    p = job["params"]
    kw = dict(
        run_gpo=p["run_gpo"], run_gpe=p["run_gpe"], run_tgpe=p["run_tgpe"], run_noise=p["run_noise"],
        run_types=p["run_types"], format="fasta", seq_type=p["run_types"][0], seed=job.get("seed", 42),
        min_support=p["min_support"], realign=p["realign"], seq_weights=p["seq_weights"],
        n_threads=1,
        consistency_anchors=p["consistency"], consistency_weight=p["consistency_weight"],
        run_vsm_amax=p["run_vsm_amax"], run_refine=p["run_refine"],
        consistency_merge=p["consistency_merge"], consistency_merge_weight=p["consistency_merge_weight"],
        vsm_mode=p["vsm_mode"], run_vsm_kappa=p["run_vsm_kappa"], run_vsm_dzero=p["run_vsm_dzero"],
        fill_in=p["fill_in"], fill_anchor_frac=p["fill_anchor_frac"], fill_support=p["fill_support"],
        run_guide_dist=p["run_guide_dist"],
    )
    # exploration (exp/fast-tree-fit): per-run guide-tree options, passed only when the dict has them, so
    # older dicts (and their cache keys) are unchanged
    kw.update({k: p[k] for k in TREE_KEYS if k in p})
    # ensemble fit (PLAN_ensemble.md, build 5b93950): per-member gap model, sequence weights, realign,
    # consistency and VSM mode, passed only when the dict has them (older dicts and keys unchanged)
    kw.update({k: p[k] for k in MEMBER_KEYS if k in p})
    return kw


def align_one(kw, inp, out):
    """One aligner call; CPU = RUSAGE_SELF delta (user+sys) around it. -> row dict."""
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    err = None
    c0, t0 = cpu_self(), time.perf_counter()
    try:
        _core.ensemble_custom_file_to_file(inp, out, **kw)
    except Exception as e:  # noqa: BLE001 - recorded in the cache as an error entry
        err = f"{type(e).__name__}: {e}"
    return {"wall": time.perf_counter() - t0, "cpu": cpu_self() - c0, "error": err}


def main():
    job = json.load(open(sys.argv[1]))
    kw = build_kwargs(job)
    rows = [align_one(kw, inp, out) for inp, out in job["cases"]]
    json.dump({"aligner": _core.__file__, "rows": rows}, open(sys.argv[2], "w"))


def serve(params_path):
    job = json.load(open(params_path))
    proto = os.fdopen(os.dup(1), "w", buffering=1)  # private protocol channel
    os.dup2(2, 1)                                    # library output -> stderr
    kw = build_kwargs(job)
    proto.write(json.dumps({"ready": True, "aligner": _core.__file__}) + "\n")
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        req = json.loads(line)
        proto.write(json.dumps(align_one(kw, req["in"], req["out"])) + "\n")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--serve":
        serve(sys.argv[2])
    else:
        main()
