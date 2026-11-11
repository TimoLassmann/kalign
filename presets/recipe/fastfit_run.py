"""Launch the fast-fit plan's fits and the CPU calibration (PLAN_fastfit.md sections 4, 5, 7).

    python -m recipe.fastfit_run start [--ws WS]     the 6 registered fits (folds 0-4, all), the seed-robustness
                                                     repeat (fold all, seed 4321, random start) and the calibration,
                                                     all at once, sharing the 60 aligner slots; returns when all
                                                     fits have exited
    python -m recipe.fastfit_run calib [--ws WS]     the calibration loop only (started by `start`)
    python -m recipe.fastfit_run status [--ws WS]    one line per fit from its progress.log

Output: WS/runs/fastfit/fold<k>/, final/, final_seed4321/ (recipe.fastfit output), WS/runs/fastfit/calib.jsonl.
Calibration (plan section 4): every 30 minutes, 4.0.0 fast parameters on the new build (member gap model) over a fixed
set (20 small cases of the fold-0 screening subset in case-id order, the first 3 training families of fold 0 at
1,000 and the first at 10,000), without the cache, each case in one aligner slot; one JSON line per case.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from . import fastfit as FF
from .evaluate import Build, ParamWorker
from .fit import SlotPool

CALIB_EVERY_S = 1800
PRESETS_DIR = Path(__file__).resolve().parents[1]
FITS = {**{f"fold{k}": ["--fold", str(k)] for k in range(5)}, "final": ["--fold", "all"],
        "final_seed4321": ["--fold", "all", "--seed", "4321", "--random-init"]}


def out_root(ws):
    return Path(ws) / "runs" / "fastfit"


def calib(ws):
    ws = Path(ws)
    rows, meta, core, screen, test, lockbox = FF.load_fold(ws, 0)
    ids = sorted(c for c in screen if meta[c]["kind"] == "small")[:20]
    ids += sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)[:3]
    ids += sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 10000)[:1]
    build = Build.load(ws / "kalign" / FF.NEW_BUILD)
    slots = SlotPool(ws / "runs" / "slots", FF.SLOTS)
    path = out_root(ws) / "calib.jsonl"
    stop = out_root(ws) / "CALIB_STOP"
    rnd = 0
    while not stop.exists():
        t0 = time.time()
        for cid in ids:
            with slots.slot():
                w = ParamWorker(build.python, FF.P400_FAST)
                try:
                    with tempfile.TemporaryDirectory(dir=os.environ.get("KALIGN_REFIT_TMP")) as d:
                        wall, cpu, err = w.align(rows[cid]["unaligned"], Path(d) / "a.fa", FF.TIMEOUT_S)
                finally:
                    w.close()
            with open(path, "a") as fh:
                fh.write(json.dumps({"round": rnd, "time": time.strftime("%Y-%m-%dT%H:%M:%S"), "case_id": cid,
                                     "cpu_s": cpu, "wall_s": wall, "error": err}) + "\n")
        rnd += 1
        while time.time() - t0 < CALIB_EVERY_S and not stop.exists():
            time.sleep(10)


def start(ws):
    ws = Path(ws)
    root = out_root(ws)
    root.mkdir(parents=True, exist_ok=True)
    (root / "CALIB_STOP").unlink(missing_ok=True)
    env = {**os.environ, "TMPDIR": str(ws / "tmp"), "KALIGN_REFIT_TMP": str(ws / "tmp")}
    procs = {}
    cal_fh = open(root / "calib.out", "a")
    cal = subprocess.Popen([sys.executable, "-m", "recipe.fastfit_run", "calib", "--ws", str(ws)], cwd=PRESETS_DIR,
                           stdout=cal_fh, stderr=subprocess.STDOUT, env=env)
    for name, args in FITS.items():
        d = root / name
        d.mkdir(parents=True, exist_ok=True)
        fh = open(d / "fit.out", "a")
        procs[name] = (subprocess.Popen([sys.executable, "-m", "recipe.fastfit", *args, "--out", str(d), "--ws", str(ws)],
                                        cwd=PRESETS_DIR, stdout=fh, stderr=subprocess.STDOUT, env=env), fh)
        print(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {name}: started pid {procs[name][0].pid}", flush=True)
    bad = []
    for name, (p, fh) in procs.items():
        rc = p.wait()
        fh.close()
        print(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {name}: exited {rc}", flush=True)
        if rc:
            bad.append(name)
    (root / "CALIB_STOP").touch()
    cal.wait()
    cal_fh.close()
    return 1 if bad else 0


def status(ws):
    root = out_root(ws)
    for name in FITS:
        log = root / name / "progress.log"
        last = log.read_text().strip().splitlines()[-1] if log.exists() and log.stat().st_size else "(no log)"
        print(f"{name:15s} {last}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["start", "calib", "status"])
    ap.add_argument("--ws", default=str(FF.default_ws()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    if a.cmd == "start":
        return start(ws)
    if a.cmd == "calib":
        calib(ws)
        return 0
    status(ws)
    return 0


if __name__ == "__main__":
    sys.exit(main())
