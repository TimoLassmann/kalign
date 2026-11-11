"""Aligner worker of amendment 2 (docs/AMENDMENT_ensconf_2.md, selection mode): recipe/_align_worker_conf.py plus
  - select_only is forwarded when the dict has it (build 5d56876);
  - the selected member (or -1) is written to OUT.sel.json, next to OUT.conf.json.
Same protocol and CPU measurement as recipe/_align_worker.py (run in the aligner venv, --serve PARAMS.json).
"""
import json
import os
import sys
import time

import _align_worker_conf  # noqa: F401 - patches _align_worker (run_tree_sl_noise)
import _align_worker as W
from kalign import _core


def build_kwargs(job):
    kw = _conf_build_kwargs(job)
    if "select_only" in job["params"]:
        kw["select_only"] = int(job["params"]["select_only"])
    return kw


def align_one(kw, inp, out):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    err, r = None, None
    c0, t0 = W.cpu_self(), time.perf_counter()
    try:
        r = _core.ensemble_custom_file_to_file(inp, out, **kw)
    except Exception as e:  # noqa: BLE001 - recorded in the cache as an error entry
        err = f"{type(e).__name__}: {e}"
    row = {"wall": time.perf_counter() - t0, "cpu": W.cpu_self() - c0, "error": err}
    if err is None and isinstance(r, dict):
        with open(f"{out}.conf.json", "w") as fh:
            json.dump([float(x) for x in r["column_confidence"]], fh)
        with open(f"{out}.sel.json", "w") as fh:
            json.dump({"selected_member": int(r.get("selected_member", -1))}, fh)
    return row


_conf_build_kwargs = W.build_kwargs
W.build_kwargs, W.align_one = build_kwargs, align_one

if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--serve":
        W.serve(sys.argv[2])
    else:
        W.main()
