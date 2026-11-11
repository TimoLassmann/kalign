"""Cost pilot of the nucfit plan (PLAN_nucfit.md section 10): peak memory and CPU of single alignments of the largest
fit inputs under the fit's seed settings, before the freeze. Aligns only; nothing is scored.

    python -m recipe.nucfit_cost --ws WS [--cases ID,ID] [--jobs 1]   -> WS/runs/cost/cost.jsonl

Each (setting, case) runs the fit worker (recipe/_align_worker.py, one-shot mode) in its own process; peak RSS and
CPU come from os.wait4 of that process.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from . import ensfit as EF
from . import nucens as NE
from . import nucfast as NF
from .evaluate import Build, WORKER_SCRIPT
from .fit import _py, pkey

DEFAULT_CASES = ("sim_dna/dna_t10000_d1.0_ir0.02_r0", "sim_dna/dna_t1000_d1.0_ir0.02_r0")


def settings():
    out = {}
    for i, p in enumerate(NF.NucSpace().seeds()):
        out[f"fastseed{i}"] = p
    # a 5-member ensemble seed of the accurate fit, with the pre-release fast setting standing in for member 1
    stand_in = _py(NF.with_gap(NF.FF.with_tree(NF.v2_params(NF.NUC_FAST_KEY), NF.TREE_PFAST), NE.GAP_MODEL))
    for i, p in enumerate(NE.NucEnsSpace(stand_in).seed_params()):
        out[f"ensseed{i}"] = p
    return out


def one(build, name, params, row, tmp, timeout_s):
    d = Path(tempfile.mkdtemp(dir=tmp))
    job, res = d / "job.json", d / "res.json"
    job.write_text(json.dumps({"params": params, "cases": [[row["unaligned"], str(d / "out.fa")]], "seed": 42}))
    t0 = time.time()
    p = subprocess.Popen([build.python, str(WORKER_SCRIPT), str(job), str(res)], stdout=subprocess.DEVNULL,
                         stderr=subprocess.PIPE, env={**os.environ, "OMP_NUM_THREADS": "1"})
    try:
        _, status, ru = os.wait4(p.pid, 0)
    except Exception:
        p.kill()
        raise
    wall = time.time() - t0
    r = json.loads(res.read_text())["rows"][0] if res.exists() else {"error": f"status {status}"}
    for f in d.iterdir():
        f.unlink()
    d.rmdir()
    return {"setting": name, "key": pkey(params), "case_id": row["case_id"], "n_seq": int(row["n_seq"]),
            "peak_rss_gb": ru.ru_maxrss / 1e6, "cpu_s": ru.ru_utime + ru.ru_stime, "wall_s": wall,
            "error": r.get("error"), "worker_cpu": r.get("cpu")}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ws", default=str(NF.default_ws()))
    ap.add_argument("--cases", default=",".join(DEFAULT_CASES))
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=1800.0)
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    EF.use_ws_tmp(ws)
    rows = {}
    for c in a.cases.split(","):          # simulation rows straight from the data folder (no manifest needed)
        d = ws / "data" / "large" / "sim" / c.split("/", 1)[1]
        rows[c] = {"case_id": c, "unaligned": str(d / "unaligned.fasta"),
                   "n_seq": json.loads((d / "info.json").read_text())["n_seq"]}
    b = Build.load(ws / "kalign" / NF.BUILD)
    out = ws / "runs" / "cost"
    out.mkdir(parents=True, exist_ok=True)
    tasks = [(n, p, rows[c]) for c in a.cases.split(",") for n, p in settings().items()]
    with open(out / "cost.jsonl", "a") as fh, ThreadPoolExecutor(a.jobs) as ex:
        for r in ex.map(lambda t: one(b, t[0], t[1], t[2], ws / "tmp", a.timeout), tasks):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(f"{r['setting']:10s} {r['case_id']:42s} RSS {r['peak_rss_gb']:6.2f} GB  CPU {r['cpu_s']:7.1f}s "
                  f"wall {r['wall_s']:6.1f}s err {r['error']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
