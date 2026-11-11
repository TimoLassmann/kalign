"""Launch the ensemble plan's fits (PLAN_ensemble.md): folds 0-4 and the final fit, all at once, sharing the 60
aligner slots (recipe.fit.SlotPool); returns when every fit has exited. Resumable: start again after a stop.

    python -m recipe.ensfit_run start [--ws WS]
    python -m recipe.ensfit_run status [--ws WS]     one line per fit from its progress.log

Output: WS/runs/ensfit/fold<k>/, final/ (recipe.ensfit output); launcher log WS/runs/ensfit/launch.out.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from . import ensfit as EF

PRESETS_DIR = Path(__file__).resolve().parents[1]
FITS = {"fold0": ["--fold", "0"], "fold1": ["--fold", "1"], "fold2": ["--fold", "2"], "fold3": ["--fold", "3"],
        "fold4": ["--fold", "4"], "final": ["--fold", "all"]}   # AMENDMENT_ensemble_3.md: all five folds again


def out_root(ws):
    return Path(ws) / "runs" / "ensfit"


def start(ws):
    ws = Path(ws)
    root = out_root(ws)
    if not EF.select10k_path(ws).exists():
        raise SystemExit("select10k.json missing: run recipe.ensfit_prep select10k first")
    env = {**os.environ, "TMPDIR": str(ws / "tmp"), "KALIGN_REFIT_TMP": str(ws / "tmp")}
    procs = {}
    for name, args in FITS.items():
        d = root / name
        d.mkdir(parents=True, exist_ok=True)
        fh = open(d / "fit.out", "a")
        procs[name] = (subprocess.Popen([sys.executable, "-m", "recipe.ensfit", *args, "--out", str(d), "--ws", str(ws)],
                                        cwd=PRESETS_DIR, stdout=fh, stderr=subprocess.STDOUT, env=env), fh)
        print(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {name}: started pid {procs[name][0].pid}", flush=True)
    bad = []
    for name, (p, fh) in procs.items():
        rc = p.wait()
        fh.close()
        print(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {name}: exited {rc}", flush=True)
        if rc:
            bad.append(name)
    return 1 if bad else 0


def status(ws):
    root = out_root(ws)
    for name in FITS:
        log = root / name / "progress.log"
        last = log.read_text().strip().splitlines()[-1] if log.exists() and log.stat().st_size else "(no log)"
        hv = root / name / "hv.json"
        g = len(json.loads(hv.read_text())["history"]) if hv.exists() else 0
        print(f"{name:7s} gen {g:3d} | {last}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["start", "status"])
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    return start(ws) if a.cmd == "start" else (status(ws) or 0)


if __name__ == "__main__":
    sys.exit(main())
