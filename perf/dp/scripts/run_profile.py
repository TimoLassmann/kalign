#!/usr/bin/env python3
"""Run the instrumented driver (kalign_prof) over datasets x modes.

Each job is single threaded (-n 1, a clean per-phase breakdown needs one
thread); up to --jobs jobs run concurrently.  One JSON per job is written to
--results (<dataset>__<mode>.json); existing results are skipped, so the
script can be re-run to resume.

    uv run --no-project python perf/dp/scripts/run_profile.py --prof build/perf/dp/kalign_prof \
        --data $DATA --results $RES --modes fast default accurate --jobs 2

Optional DP dump: --dump DIR sets KALIGN_DP_DUMP (plus --dump-every,
--dump-max) so the same runs also write .kdp problems (timings of a dumping
run include file I/O; use separate runs for profiling and dumping).
"""
import argparse
import concurrent.futures as cf
import json
import os
import subprocess
import sys
import time


def jobs_for(manifest, modes, skip, only):
    out = []
    for d in manifest["datasets"]:
        if only and not any(s in d["name"] for s in only):
            continue
        for m in modes:
            key = f"{d['name']}__{m}"
            if key in skip:
                continue
            if d["n_seq"] > 500 or d.get("model_len", 0) > 500:
                sys.exit(f"refusing {d['name']}: exceeds N<=500, L<=500")
            out.append((d, m))
    # largest first for better packing
    out.sort(key=lambda x: -(x[0]["n_seq"] * x[0]["mean_len"] ** 2 * (5 if x[1] == "accurate" else 1)))
    return out


def run_one(args, d, mode):
    name = f"{d['name']}__{mode}"
    res = os.path.join(args.results, name + ".json")
    if os.path.exists(res) and not args.force:
        return name, "exists", 0.0
    fa = os.path.join(args.data, d["file"])
    out_fa = os.path.join(args.results, name + ".aln.fa")
    env = dict(os.environ)
    if args.dump:
        env["KALIGN_DP_DUMP"] = os.path.join(args.dump, name)
        os.makedirs(env["KALIGN_DP_DUMP"], exist_ok=True)
        env["KALIGN_DP_DUMP_TAG"] = name + "__"  # unique file names across runs (oracle files are keyed by basename)
        env["KALIGN_DP_DUMP_EVERY"] = str(args.dump_every)
        env["KALIGN_DP_DUMP_MAX"] = str(args.dump_max)
    tmp = res + ".tmp"
    cmd = [args.prof, "-i", fa, "-o", out_fa, "-m", mode, "-n", "1", "-j", tmp]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=args.timeout)
    except subprocess.TimeoutExpired:
        with open(os.path.join(args.results, name + ".timeout"), "w") as fh:
            fh.write(f"timeout after {args.timeout}s\n")
        return name, "timeout", time.time() - t0
    dt = time.time() - t0
    if p.returncode != 0:
        with open(os.path.join(args.results, name + ".err"), "w") as fh:
            fh.write(p.stderr)
        return name, f"fail rc={p.returncode}", dt
    with open(tmp) as fh:
        js = json.load(fh)
    js["dataset"] = d
    with open(res, "w") as fh:
        json.dump(js, fh)
    os.remove(tmp)
    if not args.keep_aln:
        os.remove(out_fa)
    return name, "ok", dt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prof", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--modes", nargs="+", default=["fast", "default", "accurate"])
    ap.add_argument("--jobs", type=int, default=1, help="concurrent jobs (hard cap 4; prefer 1-2)")
    ap.add_argument("--timeout", type=int, default=120, help="per-job timeout in s (keep runs small)")
    ap.add_argument("--skip", nargs="*", default=[],
                    help="dataset__mode keys to skip")
    ap.add_argument("--only", nargs="*", default=None, help="substring filter on dataset names")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--keep-aln", action="store_true")
    ap.add_argument("--dump", default=None)
    ap.add_argument("--dump-every", type=int, default=1)
    ap.add_argument("--dump-max", type=int, default=50)
    args = ap.parse_args()
    if args.jobs > 4:
        sys.exit("--jobs > 4 is not allowed (memory limits)")
    os.makedirs(args.results, exist_ok=True)
    with open(os.path.join(args.data, "manifest.json")) as fh:
        manifest = json.load(fh)
    jobs = jobs_for(manifest, args.modes, set(args.skip or []), args.only)
    print(f"{len(jobs)} jobs, {args.jobs} concurrent", flush=True)
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_one, args, d, m) for d, m in jobs]
        for f in cf.as_completed(futs):
            name, status, dt = f.result()
            print(f"{name:40s} {status:10s} {dt:9.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
