"""Aligner worker of the confidence re-fit (docs/PLAN_ensconf.md): recipe/_align_worker.py plus two things, so the
frozen worker stays unchanged:
  - run_tree_sl_noise (single-linkage tree noise per member, build ea02ca9) is forwarded when the dict has it;
  - when the aligner returns the ensemble's confidence, the column confidence is written to OUT.conf.json
    (a JSON list, one value per output column), after the CPU of the call has been measured.
Same protocol and CPU measurement as recipe/_align_worker.py (run in the aligner venv, --serve PARAMS.json).
"""
import json
import os
import sys
import time

import _align_worker as W
from kalign import _core


def build_kwargs(job):
    kw = _base_build_kwargs(job)
    if "run_tree_sl_noise" in job["params"]:
        kw["run_tree_sl_noise"] = job["params"]["run_tree_sl_noise"]
    return kw


def align_one(kw, inp, out):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    err, conf = None, None
    c0, t0 = W.cpu_self(), time.perf_counter()
    try:
        r = _core.ensemble_custom_file_to_file(inp, out, **kw)
        if isinstance(r, dict):
            conf = r["column_confidence"]
    except Exception as e:  # noqa: BLE001 - recorded in the cache as an error entry
        err = f"{type(e).__name__}: {e}"
    row = {"wall": time.perf_counter() - t0, "cpu": W.cpu_self() - c0, "error": err}
    if conf is not None and err is None:
        with open(f"{out}.conf.json", "w") as fh:
            json.dump([float(x) for x in conf], fh)
    return row


_base_build_kwargs = W.build_kwargs
W.build_kwargs, W.align_one = build_kwargs, align_one

if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--serve":
        W.serve(sys.argv[2])
    else:
        W.main()
