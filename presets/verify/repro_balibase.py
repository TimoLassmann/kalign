"""Independent reproduction of the v2 presets on BAliBASE 3 (laptop).

Independent of the workstation pipeline in four ways: a fresh build of Kalign from source into a new venv (ARM,
not the workstation's x86), no cache, the original BAliBASE 3 release (lbgi.fr download in
benchmarks/data/downloads/bb3_release, not the drive5 copy the re-fit used), and the official bali_score program
(BAliBASE core blocks from the .xml files) instead of the re-fit scorer.

    <fresh venv>/bin/python presets/verify/repro_balibase.py --bb3 BB3_DIR --presets presets/results/v2/presets_v2.json \
        --bali-score BALI_SCORE --out OUT [--workers 8]

Older releases can be added as programs: --cli "Kalign 3.1=/path/kalign" (run as EXE -i IN -o OUT -f msf, CPU of the
child process); --only-cli runs just those.
Arms: the v2 presets (protein fast, accurate) through _core.ensemble_custom_file_to_file (the same call the re-fit
worker makes), and the 3.7.0 modes shipped in this build (fast, default, accurate) through the mode API. One thread
per alignment; CPU = user+sys of the worker process around the call. Output: OUT/rows.tsv (one row per case and
arm) and OUT/summary.tsv (mean SP and TC per arm and category).
"""

import argparse
import concurrent.futures as cf
import csv
import json
import os
import re
import resource
import subprocess
import time
from pathlib import Path

for _k, _v in {"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1"}.items():
    os.environ[_k] = _v

CATS = ("RV11", "RV12", "RV20", "RV30", "RV40", "RV50")


def cpu_self():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def kwargs(p):
    return dict(run_gpo=p["run_gpo"], run_gpe=p["run_gpe"], run_tgpe=p["run_tgpe"], run_noise=p["run_noise"],
                run_types=p["run_types"], format="msf", seq_type=p["run_types"][0], seed=42,
                min_support=p["min_support"], realign=p["realign"], seq_weights=p["seq_weights"], n_threads=1,
                consistency_anchors=p["consistency"], consistency_weight=p["consistency_weight"],
                run_vsm_amax=p["run_vsm_amax"], run_refine=p["run_refine"],
                consistency_merge=p["consistency_merge"], consistency_merge_weight=p["consistency_merge_weight"],
                vsm_mode=p["vsm_mode"], run_vsm_kappa=p["run_vsm_kappa"], run_vsm_dzero=p["run_vsm_dzero"],
                fill_in=p["fill_in"], fill_anchor_frac=p["fill_anchor_frac"], fill_support=p["fill_support"],
                run_guide_dist=p["run_guide_dist"])


def cpu_children():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def align(task):
    arm, spec, inp, out = task
    if spec["kind"] == "cli":   # an older Kalign release: EXE -i IN -o OUT -f msf; CPU of the child process
        c0, t0 = cpu_children(), time.perf_counter()
        p = subprocess.run([spec["exe"], "-i", inp, "-o", out, "-f", "msf"], capture_output=True, text=True)
        err = "" if p.returncode == 0 and os.path.exists(out) else f"exit {p.returncode}"
        return arm, inp, out, cpu_children() - c0, time.perf_counter() - t0, err
    from kalign import _core
    c0, t0 = cpu_self(), time.perf_counter()
    err = ""
    try:
        if spec["kind"] == "params":
            _core.ensemble_custom_file_to_file(inp, out, **kwargs(spec["params"]))
        else:
            _core.align_file_to_file_mode(inp, out, spec["mode"], format="msf", n_threads=1)
    except Exception as e:  # noqa: BLE001
        err = f"{type(e).__name__}: {e}"
    return arm, inp, out, cpu_self() - c0, time.perf_counter() - t0, err


def score(bali_score, xml, msf):
    p = subprocess.run([bali_score, xml, msf], capture_output=True, text=True)
    m = re.search(r"^auto \S+ ([\d.]+) ([\d.]+)\s*$", p.stdout, re.M)
    if p.returncode != 0 or not m:
        return None, None, f"bali_score exit {p.returncode}"
    return float(m.group(1)), float(m.group(2)), ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bb3", required=True, type=Path)
    ap.add_argument("--presets", required=True, type=Path)
    ap.add_argument("--bali-score", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cli", action="append", default=[], help="NAME=EXE: an older Kalign release run as a program")
    ap.add_argument("--only-cli", action="store_true", help="run only the --cli arms")
    a = ap.parse_args()
    pr = json.loads(a.presets.read_text())["presets"]["protein"]
    arms = {"new fast": {"kind": "params", "params": pr["fast"]["params"]},
            "new accurate": {"kind": "params", "params": pr["accurate"]["params"]},
            "3.7.0 fast": {"kind": "mode", "mode": "fast"},
            "3.7.0 default": {"kind": "mode", "mode": "default"},
            "3.7.0 accurate": {"kind": "mode", "mode": "accurate"}}
    if a.only_cli:
        arms = {}
    for c in a.cli:
        name, exe = c.split("=", 1)
        arms[name] = {"kind": "cli", "exe": exe}
    cases = sorted((c, f) for c in CATS for f in (a.bb3 / c).glob("*.tfa"))
    tasks = []
    for arm, spec in arms.items():
        d = a.out / "aln" / arm.replace(" ", "_")
        d.mkdir(parents=True, exist_ok=True)
        tasks += [(arm, spec, str(f), str(d / (f.stem + ".msf"))) for _, f in cases]
    cat_of = {str(f): c for c, f in cases}
    rows = []
    t0 = time.time()
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for i, (arm, inp, out, cpu, wall, err) in enumerate(ex.map(align, tasks, chunksize=4), 1):
            sp = tc = None
            if not err:
                sp, tc, err = score(a.bali_score, inp[:-4] + ".xml", out)
            rows.append({"arm": arm, "category": cat_of[inp], "case": Path(inp).stem, "sp": sp, "tc": tc,
                         "cpu_s": round(cpu, 4), "wall_s": round(wall, 4), "error": err})
            if i % 200 == 0:
                print(f"{i}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    with open(a.out / "rows.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    summ = []
    for arm in arms:
        for c in CATS + ("all",):
            rr = [r for r in rows if r["arm"] == arm and (c == "all" or r["category"] == c)]
            ok = [r for r in rr if not r["error"]]
            summ.append({"arm": arm, "category": c, "n": len(rr), "errors": len(rr) - len(ok),
                         "sp": sum(r["sp"] for r in ok) / max(len(ok), 1), "tc": sum(r["tc"] for r in ok) / max(len(ok), 1),
                         "cpu_s_per_case": sum(r["cpu_s"] for r in rr) / len(rr)})
    with open(a.out / "summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0]), delimiter="\t")
        w.writeheader()
        w.writerows(summ)
    print(f"done in {time.time() - t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
