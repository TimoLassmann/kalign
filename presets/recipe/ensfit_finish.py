"""Finish the ensemble experiment without supervision (AMENDMENT_ensemble_4.md):

  as each fit's process logs "fit done": recipe.ensfit_select for that fit (selection_a4.json);
  when folds 0-4 all have selection_a4.json: recipe.ensfit_oof run, then recipe.ensfit_oof gate.

Idempotent: a fit with selection_a4.json is not selected again; the out-of-fold run resumes from its rows.

    python -m recipe.ensfit_finish [--ws WS]          (run detached; log: runs/ensfit/finish.log)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from . import ensfit as EF

FITS = ["final", "fold0", "fold1", "fold2", "fold3", "fold4"]
POLL_S = 300
PRESETS_DIR = Path(__file__).resolve().parents[1]


def log(fh, msg):
    fh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}\n")
    fh.flush()


def fit_done(d):
    p = d / "progress.log"
    if not p.exists():
        return False
    return any("] fit done" in line and "[a4]" not in line for line in open(p, errors="replace"))


def run_module(fh, args):
    log(fh, "start: python -m " + " ".join(args))
    rc = subprocess.call([sys.executable, "-m", *args], cwd=PRESETS_DIR, stdout=fh, stderr=subprocess.STDOUT)
    log(fh, f"exit {rc}: python -m {' '.join(args)}")
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ws", default=str(EF.default_ws()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    EF.use_ws_tmp(ws)
    root = ws / "runs" / "ensfit"
    fh = open(root / "finish.log", "a")
    log(fh, "finisher started")
    procs, failed = {}, set()
    while True:
        for name in FITS:
            d = root / name
            if name in procs or name in failed or (d / "selection_a4.json").exists() or not fit_done(d):
                continue
            log(fh, f"{name}: search and fit process done -> selection (amendment 4)")
            procs[name] = subprocess.Popen([sys.executable, "-m", "recipe.ensfit_select", "--fit", name,
                                            "--ws", str(ws)], cwd=PRESETS_DIR, stdout=fh, stderr=subprocess.STDOUT)
        for name, p in list(procs.items()):
            rc = p.poll()
            if rc is not None:
                log(fh, f"{name}: selection exit {rc}")
                del procs[name]
                if rc:
                    failed.add(name)
        if failed and not procs:
            log(fh, f"selection failed for {sorted(failed)}; stopping (fix and restart the finisher)")
            return 1
        if all((root / f"fold{k}" / "selection_a4.json").exists() for k in range(5)) and not procs:
            break
        time.sleep(POLL_S)
    if run_module(fh, ["recipe.ensfit_oof", "run", "--ws", str(ws)]):
        return 1
    if run_module(fh, ["recipe.ensfit_oof", "gate", "--ws", str(ws)]):
        return 1
    while not (root / "final" / "selection_a4.json").exists():
        time.sleep(POLL_S)
    log(fh, "finished: selection_a4.json for every fit, runs/ensfit/oof/gate.json written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
