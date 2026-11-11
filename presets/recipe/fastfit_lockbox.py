"""One-shot lockbox of the fast-fit (PLAN_fastfit.md sections 3, 7, 9; DECISION_fastfit_selection.md).

  run --dry-run  checks only (case counts, leak guards, file hashes, arms and keys); runs nothing scored
  run            every lockbox input with every arm -> runs/fastfit/lockbox/rows.jsonl (+ cli_identity.jsonl)
  rss            peak RSS of the 50,000 / 200,000 runs re-measured through a fork shim (mechanics fix, not scored)
                 -> runs/fastfit/lockbox/rss.jsonl
  gate           rows -> per-arm, per-size table and lockbox criteria 3-4 -> runs/fastfit/lockbox/gate.json

Arms (fixed before any lockbox case was aligned):
  knee              the final fit's knee pick (runs/fastfit/final_pick.json; DECISION rule), new build
  registered_ultra  the final fit's registered ultra-fast pick (runs/fastfit/final/knee_selection.json), new build
  gm1               the gm1 fast preset with the 4.0.0 tree (knee_selection.json "gm1"; plan section 7 fallback)
  p400              Kalign 4.0.0 fast (P400_FAST on the 4.0.0 build): the reference
  famsa             FAMSA 2.5.2, default options (reported only)

Cases: every lockbox case of load_fold(ws, "all") (60 families at 1,000, 49 at 10,000, 10 at 50,000, the four
200,000-sequence inputs).
  1,000 and 10,000: one thread, timeout 1,800 s. Kalign arms through the fit's Evaluator (cache, ParamWorker, slot
     pool) exactly as recipe.fastfit_oof; FAMSA "-t 1" under one slot.
  50,000 and 200,000: 8 threads (8 slots held), timeout 4 h wall per input (criterion 4). Kalign arms: a one-shot
     process of the arm's build that builds its keyword arguments with recipe/_align_worker.py build_kwargs (the
     evaluator's code) and sets n_threads = 8; FAMSA "-t 8". Not written to the fit cache.
  A failed case (error or timeout) scores SP = TC = 0 and is counted.
  Knee only: the tagged release CLI (`kalign --mode fast -n 1`, f739560) on every 1,000 and 10,000 input; its output
     sha256 is compared with the evaluator's output sha256 (cache aln_sha256). Not scored.

Gate (plan section 7, lockbox part): criterion 3: Q at 1,000 and Q at 10,000 (Q = mean over families of
(SP + TC) / 2) above 4.0.0 fast; criterion 4: every 50,000 and 200,000 input completes. The same for gm1 (fallback);
registered_ultra with ">=" (plan section 7 ultra rule); the DECISION exception (knee fails 4, ultra passes 4).
Paired differences to 4.0.0 fast with 95 % bootstrap intervals over families (10,000 resamples, seed 1234) are
reported, not gated. Criteria 1-2 are copied from runs/fastfit/oof/gate.json.

No per-case score is printed by run; gate prints aggregates only.

    python -m recipe.fastfit_lockbox run [--dry-run] | rss | gate   [--ws WS]
"""
from __future__ import annotations

import argparse
import contextlib
import json
import os
import signal
import subprocess
import sys
import tempfile
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from . import cache as C
from . import fit as FIT
from .fastfit import NEW_BUILD, P400_FAST, REF_BUILD, SLOTS, load_fold
from .fit import Evaluator, LeakError, SlotPool, _atomic_json, log, pkey

FIT.MAX_SLOTS = SLOTS
ARMS = ("knee", "registered_ultra", "gm1", "p400", "famsa")
KALIGN_ARMS = ("knee", "registered_ultra", "gm1", "p400")
EXPECTED = {1000: 60, 10000: 49, 50000: 10, 200000: 4}
SMALL_K, LARGE_K = (1000, 10000), (50000, 200000)
SMALL_TIMEOUT_S = 1800.0
LARGE_TIMEOUT_S = 4 * 3600.0
LARGE_THREADS = 8
SMALL_WORKERS = 40            # concurrent single-thread aligners (the box also runs another 2 x 8-thread job)
LARGE_CONCURRENT = 5          # concurrent 8-thread aligners (40 slots)
N_BOOT, SEED = 10000, 1234
HOME = Path.home()
RELEASE_CLI = HOME / "tmp/kalign-rel/f739560/src/build/src/kalign"
FAMSA = HOME / "tmp/k4bench/tools/famsa/famsa"
WORKER = Path(__file__).with_name("_align_worker.py")
FIT_DIRS = ("fold0", "fold1", "fold2", "fold3", "fold4", "final", "final_seed4321")

# one-shot multi-threaded aligner process: thread environment before the library loads, then the evaluator's own
# build_kwargs (recipe/_align_worker.py; its single-thread environment settings come after the library is loaded and
# do not apply), n_threads replaced. argv: worker.py job.json input output n_threads
MT_SNIPPET = r"""
import importlib.util, json, os, sys
n = int(sys.argv[4])
os.environ.update({"OMP_NUM_THREADS": str(n), "OMP_THREAD_LIMIT": str(n), "OMP_DYNAMIC": "FALSE",
                   "OMP_WAIT_POLICY": "PASSIVE"})
from kalign import _core
spec = importlib.util.spec_from_file_location("_align_worker", sys.argv[1])
W = importlib.util.module_from_spec(spec)
spec.loader.exec_module(W)
kw = W.build_kwargs(json.load(open(sys.argv[2])))
kw["n_threads"] = n
_core.ensemble_custom_file_to_file(sys.argv[3], sys.argv[5], **kw)
"""


def lb_dir(ws):
    d = ws / "runs" / "fastfit" / "lockbox"
    d.mkdir(parents=True, exist_ok=True)
    return d


def tmp_root(ws):
    d = Path(os.environ.get("KALIGN_REFIT_TMP") or ws / "tmp") / "lockbox"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------------------------------ arms and cases
def load_arms(ws):
    from .evaluate import Build
    ff = ws / "runs" / "fastfit"
    pick = json.load(open(ff / "final_pick.json"))
    sel = json.load(open(ff / "final" / "knee_selection.json"))
    arms = {"knee": {"key": pick["key"], "params": pick["params"], "build": NEW_BUILD,
                     "source": "runs/fastfit/final_pick.json"},
            "registered_ultra": {"key": sel["registered_ultra"]["key"], "params": sel["registered_ultra"]["params"],
                                 "build": NEW_BUILD, "source": "final/knee_selection.json registered_ultra"},
            "gm1": {"key": sel["gm1"]["key"], "params": sel["gm1"]["params"], "build": NEW_BUILD,
                    "source": "final/knee_selection.json gm1"},
            "p400": {"key": pkey(P400_FAST), "params": P400_FAST, "build": REF_BUILD, "source": "fastfit.P400_FAST"}}
    if pick["key"] != sel["knee"]["key"] or pick["params_sha256"] != pick["key"]:
        raise SystemExit("final_pick.json and final/knee_selection.json disagree on the knee pick")
    for a, v in arms.items():
        if pkey(v["params"]) != v["key"]:
            raise SystemExit(f"arm {a}: params do not hash to their key {v['key'][:12]}")
    builds = {b: Build.load(ws / "kalign" / b) for b in (NEW_BUILD, REF_BUILD)}
    for v in arms.values():
        v["build_id"] = builds[v["build"]].build_id
    arms["famsa"] = {"key": "famsa:default:" + C.sha256_file(FAMSA)[:16], "params": None, "build": str(FAMSA),
                     "build_id": C.sha256_file(FAMSA), "source": "FAMSA default options"}
    return arms, builds


def check_cases(ws):
    """Lockbox ids per size + leak guards (plan section 9). Raises on any disagreement."""
    rows, meta, _, _, _, lockbox = load_fold(ws, "all")
    lock = sorted(lockbox)
    by_k = Counter(meta[c]["k"] for c in lock)
    if dict(by_k) != EXPECTED:
        raise SystemExit(f"lockbox counts {dict(by_k)} != expected {EXPECTED}")
    lfam = {meta[c]["family"] for c in lock}
    other_fam = {m["family"] for c, m in meta.items() if m["kind"] == "large" and c not in lockbox}
    if lfam & other_fam:
        raise LeakError(f"lockbox families with non-lockbox cases: {sorted(lfam & other_fam)[:5]}")
    seen = {}
    for f in (0, 1, 2, 3, 4, "all"):
        _, _, core, screen, test, lb = load_fold(ws, f)
        if set(lb) != set(lockbox):
            raise LeakError(f"fold {f}: lockbox differs")
        for name, ids in (("core", core), ("screen", screen), ("test", test)):
            bad = set(ids) & set(lockbox)
            fams = {meta[c]["family"] for c in ids if meta[c]["kind"] == "large"} & lfam
            if bad or fams:
                raise LeakError(f"fold {f} {name}: {len(bad)} lockbox cases, {len(fams)} lockbox families")
            seen[f"fold{f}.{name}"] = len(ids)
    ff = ws / "runs" / "fastfit"
    for d in FIT_DIRS:
        for name in ("cases_core.txt", "cases_screen.txt"):
            p = ff / d / name
            if not p.exists():
                raise SystemExit(f"{p} missing")
            ids = [l.strip() for l in open(p) if l.strip()]
            if set(ids) & set(lockbox):
                raise LeakError(f"{p}: lockbox case ids")
            seen[f"{d}/{name}"] = len(ids)
    oof = ff / "oof" / "rows.jsonl"
    n_oof = 0
    for l in open(oof):
        n_oof += 1
        if json.loads(l)["case_id"] in lockbox:
            raise LeakError(f"{oof}: a lockbox case")
    seen["oof_rows"] = n_oof
    bad_files = []
    for c in lock:
        r = rows[c]
        if C.sha256_file(r["unaligned"]) != r["sha256_unaligned"] or \
                C.sha256_file(r["reference"]) != r["sha256_reference"]:
            bad_files.append(c)
    if bad_files:
        raise SystemExit(f"lockbox files differ from the manifest hashes: {bad_files}")
    ids = {k: sorted(c for c in lock if meta[c]["k"] == k) for k in EXPECTED}
    return rows, meta, ids, {"counts": {str(k): len(v) for k, v in ids.items()}, "checked": seen,
                             "files_sha256_ok": len(lock)}


def cache_lockbox_rows(ws, lock):
    import sqlite3
    c = sqlite3.connect(f"file:{ws / 'cache' / 'cache.sqlite'}?mode=ro", uri=True)
    q = ",".join("?" * len(lock))
    out = c.execute(f"select build_id, param_key, count(*) from results where case_id in ({q}) group by 1, 2",
                    sorted(lock)).fetchall()
    c.close()
    return out


# ------------------------------------------------------------------------------------ processes
# Peak RSS: a child forked from this (large) process inherits its memory high-water mark, so ru_maxrss of a direct
# child is at least this process's peak (incident 2026-10-15: the run's large-input RSS values were the harness's
# VmHWM). The shim (a small exec'd python) forks the aligner and reports the aligner's own ru_maxrss and CPU.
# argv: result.json aligner argv...
RSS_SHIM = r"""
import json, os, sys
pid = os.fork()
if pid == 0:
    os.execv(sys.argv[2], sys.argv[2:])
_, status, ru = os.wait4(pid, 0)
json.dump({"maxrss_kb": ru.ru_maxrss, "cpu": ru.ru_utime + ru.ru_stime, "rc": os.waitstatus_to_exitcode(status)},
          open(sys.argv[1], "w"))
"""


def run_proc_rss(argv, timeout_s, env):
    """run_proc through RSS_SHIM -> (wall, cpu, maxrss_kb of the aligner, rc or None if killed, stderr tail)."""
    with tempfile.TemporaryDirectory(prefix="lb_rss_", dir=os.environ.get("KALIGN_REFIT_TMP")) as d:
        res = Path(d) / "rss.json"
        wall, cpu, _, rc, txt = run_proc([sys.executable, "-c", RSS_SHIM, str(res), *argv], timeout_s, env)
        if rc is None or not res.exists():
            return wall, cpu, None, rc, txt
        r = json.loads(res.read_text())
        return wall, r["cpu"], r["maxrss_kb"], r["rc"], txt


def run_proc(argv, timeout_s, env):
    """Own session, killed (group) after timeout_s. -> (wall, cpu, maxrss_kb, rc or None if killed, stderr tail).
    maxrss_kb is NOT the aligner's own peak when this process is large (see RSS_SHIM)."""
    t0 = time.perf_counter()
    with tempfile.TemporaryFile() as errf:
        p = subprocess.Popen(argv, stdout=subprocess.DEVNULL, stderr=errf, start_new_session=True, env=env)
        lock, st = threading.Lock(), {"done": False, "killed": False}

        def _kill():
            with lock:
                if not st["done"]:
                    st["killed"] = True
                    try:
                        os.killpg(p.pid, signal.SIGKILL)
                    except (ProcessLookupError, PermissionError):
                        pass
        timer = threading.Timer(timeout_s, _kill)
        timer.daemon = True
        timer.start()
        _, status, ru = os.wait4(p.pid, 0)
        with lock:
            st["done"] = True
        timer.cancel()
        rc = os.waitstatus_to_exitcode(status)
        errf.seek(0)
        txt = errf.read().decode(errors="replace")[-300:].replace("\n", " ")
    return (time.perf_counter() - t0, ru.ru_utime + ru.ru_stime, ru.ru_maxrss,
            None if st["killed"] else rc, txt)


def thread_env(n):
    return {**os.environ, "OMP_NUM_THREADS": str(n), "OMP_THREAD_LIMIT": str(n), "OMP_DYNAMIC": "FALSE",
            "OMP_WAIT_POLICY": "PASSIVE"}


def aligner_argv(arm, v, builds, row, threads, d, out):
    if arm == "famsa":
        return [str(FAMSA), "-t", str(threads), row["unaligned"], str(out)]
    job = Path(d) / "job.json"
    job.write_text(json.dumps({"params": v["params"], "seed": 42}))
    return [builds[v["build"]].python, "-c", MT_SNIPPET, str(WORKER), str(job), row["unaligned"], str(threads),
            str(out)]


def align_score(arm, v, builds, row, threads, timeout_s, tmp):
    """One aligner process (kalign arm via MT_SNIPPET, or FAMSA) + scoring. -> row dict (no cache)."""
    from .scoring import score_case
    with tempfile.TemporaryDirectory(prefix=f"lb_{arm}_", dir=tmp) as d:
        out = Path(d) / "aln.fa"
        argv = aligner_argv(arm, v, builds, row, threads, d, out)
        wall, cpu, rss, rc, txt = run_proc(argv, timeout_s, thread_env(threads))
        err, m, sha = None, {}, None
        if rc is None:
            err = "timeout"
        elif rc != 0 or not out.exists():
            err = f"exit {rc}: {txt}"
        else:
            try:
                sha = C.sha256_file(out)
                m = score_case(row, str(out))
            except Exception as e:  # noqa: BLE001 - recorded as a failed case
                err = f"score {type(e).__name__}: {e}"[:300]
    return {"recall": m.get("recall"), "tc": m.get("tc"), "precision": m.get("precision"), "cpu_s": cpu,
            "wall_s": wall, "peak_rss_kb": rss, "aln_sha256": sha, "error": err}


# ------------------------------------------------------------------------------------ run
def run(ws, dry):
    out = lb_dir(ws)
    rows, meta, ids, chk = check_cases(ws)
    arms, builds = load_arms(ws)
    small = ids[1000] + ids[10000]
    large = ids[200000] + ids[50000]
    lock = small + large
    rows_path, id_path = out / "rows.jsonl", out / "cli_identity.jsonl"
    prior = cache_lockbox_rows(ws, lock)
    allowed_prior = {(arms[a]["build_id"], arms[a]["key"]) for a in KALIGN_ARMS}
    if any((b, k) not in allowed_prior for b, k, _ in prior) or (prior and not rows_path.exists()):
        raise LeakError(f"fit cache holds lockbox rows not written by this lockbox run: {prior}")
    plan = {"cases": {str(k): v for k, v in ids.items()}, "checks": chk,
            "arms": {a: {kk: v[kk] for kk in ("key", "build", "build_id", "source")} for a, v in arms.items()},
            "small": {"threads": 1, "timeout_s": SMALL_TIMEOUT_S, "workers": SMALL_WORKERS},
            "large": {"threads": LARGE_THREADS, "timeout_s": LARGE_TIMEOUT_S, "concurrent": LARGE_CONCURRENT},
            "release_cli": {"path": str(RELEASE_CLI), "sha256": C.sha256_file(RELEASE_CLI)},
            "sha256": {"final_pick.json": C.sha256_file(ws / "runs/fastfit/final_pick.json"),
                       "final/knee_selection.json": C.sha256_file(ws / "runs/fastfit/final/knee_selection.json"),
                       "fastfit_lockbox.py": C.sha256_file(__file__),
                       "_align_worker.py": C.sha256_file(WORKER)},
            "cache_lockbox_rows_before": prior}
    if dry:
        for k, v in ids.items():
            print(f"{k:>6}: {len(v)} cases: {' '.join(v)}")
        for a, v in plan["arms"].items():
            print(f"arm {a:17s} key {v['key']}  build {v['build']} ({v['build_id'][:12]})  [{v['source']}]")
        print(json.dumps({k: plan[k] for k in ("checks", "small", "large", "release_cli", "sha256",
                                                "cache_lockbox_rows_before")}, indent=1))
        _atomic_json(out / "dryrun.json", {**plan, "written": time.strftime("%Y-%m-%dT%H:%M:%S")})
        print("dry run: nothing aligned")
        return 0
    logfh = open(out / "progress.log", "a")
    runs_meta = out / "run.json"
    rm = json.loads(runs_meta.read_text()) if runs_meta.exists() else {**plan, "started": time.strftime(
        "%Y-%m-%dT%H:%M:%S"), "resumed": []}
    if runs_meta.exists():
        if rm["arms"] != plan["arms"] or rm["cases"] != plan["cases"] or rm["large"] != plan["large"] \
                or rm["small"]["timeout_s"] != SMALL_TIMEOUT_S:
            raise SystemExit("arms, cases or limits differ from the recorded lockbox run; refusing to resume")
        rm["resumed"].append({"at": time.strftime("%Y-%m-%dT%H:%M:%S"), "sha256": plan["sha256"]})
    _atomic_json(runs_meta, rm)
    t_start = time.time()
    done = set()
    if rows_path.exists():
        for l in open(rows_path):
            r = json.loads(l)
            done.add((r["arm"], r["case_id"]))
    wlock = threading.Lock()

    def write(rec):
        with wlock, open(rows_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
            done.add((rec["arm"], rec["case_id"]))

    def base(arm, cid, threads, timeout_s, path):
        return {"arm": arm, "key": arms[arm]["key"], "case_id": cid, "family": meta[cid]["family"],
                "k": meta[cid]["k"], "threads": threads, "timeout_s": timeout_s, "path": path}

    slots = SlotPool(ws / "runs" / "slots", SLOTS)
    cache = ws / "cache" / "cache.sqlite"
    tmp = tmp_root(ws)
    log(f"[lockbox] start: {len(small)} inputs at 1,000/10,000, {len(large)} at 50,000/200,000; arms "
        f"{', '.join(f'{a}={v['key'][:10]}' for a, v in arms.items())}", logfh)

    # ---- phase 1: 1,000 and 10,000, Kalign arms through the evaluator (as recipe.fastfit_oof)
    for build, arm_list in ((NEW_BUILD, ("knee", "registered_ultra", "gm1")), (REF_BUILD, ("p400",))):
        todo = [a for a in arm_list if any((a, c) not in done for c in small)]
        if not todo:
            continue
        ev = Evaluator(builds[build], rows.values(), small, cache_path=cache, slots=slots, threads=SMALL_WORKERS,
                       timeout_s=SMALL_TIMEOUT_S, logfh=logfh, rerun_short_timeouts=True)
        specs = {arms[a]["key"]: arms[a]["params"] for a in todo}
        log(f"[lockbox] evaluator {build}: {', '.join(todo)} on {len(small)} inputs", logfh)
        res = ev.evaluate(specs, small)
        for a in todo:
            for cid in small:
                if (a, cid) in done:
                    continue
                r = res[arms[a]["key"]][cid]
                write({**base(a, cid, 1, SMALL_TIMEOUT_S, f"evaluator:{build}"), "recall": r.get("recall"),
                       "tc": r.get("tc"), "precision": r.get("precision"), "cpu_s": r.get("cpu_s"),
                       "wall_s": r.get("wall_s"), "peak_rss_kb": None, "aln_sha256": r.get("aln_sha256"),
                       "error": r.get("error"), "cached": r.get("cached")})
            n_err = sum(1 for cid in small if res[arms[a]["key"]][cid].get("error"))
            log(f"[lockbox] {a}: 1,000/10,000 done, {n_err} failed", logfh)

    # ---- phase 1b: FAMSA at 1,000 and 10,000 (one thread, one slot)
    todo = [c for c in sorted(small, key=lambda c: -meta[c]["k"]) if ("famsa", c) not in done]
    if todo:
        log(f"[lockbox] famsa: {len(todo)} inputs at 1,000/10,000", logfh)

        def fjob(cid):
            with slots.slot():
                r = align_score("famsa", arms["famsa"], builds, rows[cid], 1, SMALL_TIMEOUT_S, tmp)
            write({**base("famsa", cid, 1, SMALL_TIMEOUT_S, "famsa"), **r})
            return r
        with ThreadPoolExecutor(SMALL_WORKERS) as ex:
            errs = sum(bool(r["error"]) for r in ex.map(fjob, todo))
        log(f"[lockbox] famsa: 1,000/10,000 done, {errs} failed", logfh)

    # ---- phase 2: release CLI identity (knee only; sha256, not scored)
    idone = {json.loads(l)["case_id"] for l in open(id_path)} if id_path.exists() else set()
    ev = Evaluator(builds[NEW_BUILD], rows.values(), small, cache_path=cache)
    knee_sha = {c: r.get("aln_sha256") for c, r in ev.lookup(arms["knee"]["key"], small).items()}
    todo = [c for c in small if c not in idone]
    if todo:
        log(f"[lockbox] release CLI identity: {len(todo)} inputs", logfh)
        env1 = {**{k: v for k, v in os.environ.items() if not k.startswith("KALIGN_")}, "OMP_NUM_THREADS": "1",
                "OMP_THREAD_LIMIT": "1"}

        def ijob(cid):
            with slots.slot(), tempfile.TemporaryDirectory(prefix="lb_cli_", dir=tmp) as d:
                o = Path(d) / "aln.fa"
                wall, cpu, rss, rc, txt = run_proc([str(RELEASE_CLI), "--mode", "fast", "-n", "1", "-i",
                                                    rows[cid]["unaligned"], "-o", str(o)], SMALL_TIMEOUT_S, env1)
                sha = C.sha256_file(o) if rc == 0 and o.exists() else None
            rec = {"case_id": cid, "k": meta[cid]["k"], "sha_release": sha, "sha_evaluator": knee_sha.get(cid),
                   "identical": bool(sha and sha == knee_sha.get(cid)),
                   "error": None if rc == 0 else ("timeout" if rc is None else f"exit {rc}: {txt}"),
                   "cpu_s": cpu, "wall_s": wall}
            with wlock, open(id_path, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            return rec
        with ThreadPoolExecutor(SMALL_WORKERS) as ex:
            recs = list(ex.map(ijob, todo))
        log(f"[lockbox] release CLI identity: {sum(r['identical'] for r in recs)}/{len(recs)} identical "
            f"(this pass)", logfh)

    # ---- phase 3: 50,000 and 200,000, 8 threads, 4 h
    order = {a: i for i, a in enumerate(("p400", "famsa", "knee", "registered_ultra", "gm1"))}
    jobs = sorted(((a, c) for c in large for a in ARMS if (a, c) not in done),
                  key=lambda j: (-meta[j[1]]["k"], order[j[0]], j[1]))
    if jobs:
        log(f"[lockbox] {len(jobs)} runs at 50,000/200,000 ({LARGE_THREADS} threads, {LARGE_CONCURRENT} at once)",
            logfh)

        def ljob(arm, cid):
            with contextlib.ExitStack() as stack:
                for _ in range(LARGE_THREADS):
                    stack.enter_context(slots.slot())
                log(f"[lockbox] start {arm} {cid}", logfh)
                r = align_score(arm, arms[arm], builds, rows[cid], LARGE_THREADS, LARGE_TIMEOUT_S, tmp)
            write({**base(arm, cid, LARGE_THREADS, LARGE_TIMEOUT_S, "mt-worker" if arm != "famsa" else "famsa"), **r})
            log(f"[lockbox] done {arm} {cid}: {'FAILED ' + r['error'][:80] if r['error'] else 'completed'}, "
                f"wall {r['wall_s']:.0f}s cpu {r['cpu_s']:.0f}s rss {r['peak_rss_kb'] / 1e6:.1f}GB", logfh)
            return r
        errs = []
        with ThreadPoolExecutor(LARGE_CONCURRENT) as ex:
            futs = [ex.submit(ljob, a, c) for a, c in jobs]
            for f in as_completed(futs):
                try:
                    f.result()
                except Exception as e:  # noqa: BLE001 - harness failure: no row written, rerun on resume
                    errs.append(repr(e))
                    log(f"[lockbox] harness error (no row written): {e!r}", logfh)
        if errs:
            raise RuntimeError(f"{len(errs)} large runs failed in the harness; fix and resume")
    missing = [(a, c) for a in ARMS for c in lock if (a, c) not in done]
    log(f"[lockbox] run finished in {time.time() - t_start:.0f}s; missing rows: {len(missing)}", logfh)
    return 1 if missing else 0


# ------------------------------------------------------------------------------------ rss (mechanics re-measurement)
def measure_rss(ws):
    """Peak RSS of every 50,000 / 200,000 run, measured through RSS_SHIM (the run's own values were inherited from the
    harness; incident 2026-10-15). Same aligner command, threads and timeout; nothing is scored. The output sha256 is
    compared with the scored run's (determinism check). -> runs/fastfit/lockbox/rss.jsonl"""
    out = lb_dir(ws)
    rm = json.loads((out / "run.json").read_text())
    arms, builds = load_arms(ws)
    if {a: {kk: v[kk] for kk in ("key", "build", "build_id", "source")} for a, v in arms.items()} != rm["arms"]:
        raise SystemExit("arms differ from the recorded lockbox run")
    rows, _, _, _, _, _ = load_fold(ws, "all")
    by = {(r["arm"], r["case_id"]): r for r in map(json.loads, open(out / "rows.jsonl"))}
    large = rm["cases"]["200000"] + rm["cases"]["50000"]
    path = out / "rss.jsonl"
    done = {(r["arm"], r["case_id"]) for r in map(json.loads, open(path))} if path.exists() else set()
    jobs = [(a, c) for c in large for a in ARMS if (a, c) not in done]
    missing = [j for j in jobs if j not in by]
    if missing:
        raise SystemExit(f"{len(missing)} large runs not in rows.jsonl; finish the run first")
    slots, tmp, wlock = SlotPool(ws / "runs" / "slots", SLOTS), tmp_root(ws), threading.Lock()
    logfh = open(out / "progress.log", "a")
    log(f"[lockbox rss] {len(jobs)} runs (aligner peak RSS via shim; not scored)", logfh)

    def job(arm, cid):
        with contextlib.ExitStack() as stack:
            for _ in range(LARGE_THREADS):
                stack.enter_context(slots.slot())
            with tempfile.TemporaryDirectory(prefix=f"lb_rss_{arm}_", dir=tmp) as d:
                o = Path(d) / "aln.fa"
                argv = aligner_argv(arm, arms[arm], builds, rows[cid], LARGE_THREADS, d, o)
                wall, cpu, rss, rc, txt = run_proc_rss(argv, LARGE_TIMEOUT_S, thread_env(LARGE_THREADS))
                sha = C.sha256_file(o) if rc == 0 and o.exists() else None
        ref = by[(arm, cid)]
        rec = {"arm": arm, "case_id": cid, "peak_rss_kb": rss, "cpu_s": cpu, "wall_s": wall,
               "error": None if rc == 0 else ("timeout" if rc is None else f"exit {rc}: {txt}"),
               "aln_sha256": sha, "same_output_as_run": sha == ref.get("aln_sha256"),
               "same_status_as_run": bool(rc == 0) == (not ref.get("error"))}
        with wlock, open(path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        log(f"[lockbox rss] {arm} {cid}: rss {(rss or 0) / 1e6:.1f}GB, same status {rec['same_status_as_run']}, "
            f"same output {rec['same_output_as_run']}", logfh)
    with ThreadPoolExecutor(LARGE_CONCURRENT) as ex:
        for f in as_completed([ex.submit(job, a, c) for a, c in jobs]):
            f.result()
    return 0


# ------------------------------------------------------------------------------------ gate
def gate(ws):
    out = lb_dir(ws)
    rm = json.loads((out / "run.json").read_text())
    rows = [json.loads(l) for l in open(out / "rows.jsonl")]
    by = {(r["arm"], r["case_id"]): r for r in rows}
    cases = {int(k): v for k, v in rm["cases"].items()}
    for a in ARMS:
        miss = [c for v in cases.values() for c in v if (a, c) not in by]
        if miss:
            raise SystemExit(f"arm {a}: {len(miss)} lockbox rows missing; finish the run first")
    rp = out / "rss.jsonl"
    rssm = {(r["arm"], r["case_id"]): r for r in map(json.loads, open(rp))} if rp.exists() else {}
    rng = np.random.default_rng(SEED)
    boot = {k: rng.integers(0, len(v), (N_BOOT, len(v))) for k, v in cases.items()}
    table, vals = {}, {}
    for a in ARMS:
        table[a] = {}
        for k, ids in cases.items():
            rs = [by[(a, c)] for c in ids]
            fail = np.array([bool(r["error"]) for r in rs])
            sp = np.array([0.0 if r["error"] else float(r["recall"]) for r in rs])
            tc = np.array([0.0 if r["error"] else float(r["tc"]) for r in rs])
            q = 0.5 * (sp + tc)
            vals[(a, k)] = {"sp": sp, "tc": tc, "q": q}
            # peak RSS: only from the shim re-measurement (rss.jsonl; 50,000 / 200,000); the run's own values were
            # inherited from the harness process (incident 2026-10-15) and are not reported
            rss = [rssm[(a, c)]["peak_rss_kb"] for c in ids if (a, c) in rssm and rssm[(a, c)].get("peak_rss_kb")
                   and not rssm[(a, c)]["error"]]
            table[a][str(k)] = {
                "n": len(ids), "n_failed": int(fail.sum()),
                "failed": [{"case_id": r["case_id"], "error": r["error"][:120]} for r in rs if r["error"]],
                "sp": float(sp.mean()), "tc": float(tc.mean()), "q": float(q.mean()),
                "cpu_s": float(sum(float(r["cpu_s"] or 0.0) for r in rs)),
                "wall_s": float(sum(float(r["wall_s"] or 0.0) for r in rs)),
                "wall_max_s": float(max(float(r["wall_s"] or 0.0) for r in rs)),
                "peak_rss_gb_max": (max(rss) / 1e6) if rss else None, "threads": rs[0]["threads"],
                "timeout_s": rs[0]["timeout_s"]}
    for a in ARMS:
        if a == "p400":
            continue
        for k in cases:
            d = {}
            for m in ("sp", "tc", "q"):
                x = vals[(a, k)][m] - vals[("p400", k)][m]
                bm = x[boot[k]].mean(axis=1)
                d[m] = {"diff": float(x.mean()), "lo": float(np.percentile(bm, 2.5)),
                        "hi": float(np.percentile(bm, 97.5))}
            table[a][str(k)]["vs_p400"] = d
    ref = table["p400"]

    def crit(a, strict):
        t = table[a]
        above = (lambda x, y: x > y) if strict else (lambda x, y: x >= y)
        c3 = {f"q_{k}": {"arm": t[str(k)]["q"], "p400": ref[str(k)]["q"], "pass": above(t[str(k)]["q"],
                                                                                       ref[str(k)]["q"])}
              for k in SMALL_K}
        c4 = {f"completed_{k}": f"{t[str(k)]['n'] - t[str(k)]['n_failed']}/{t[str(k)]['n']}" for k in LARGE_K}
        return {"c3_lockbox_q_above_p400" if strict else "c3_lockbox_q_ge_p400": c3,
                "c3_pass": all(v["pass"] for v in c3.values()),
                "c4_detail": c4, "c4_pass": all(t[str(k)]["n_failed"] == 0 for k in LARGE_K)}
    crits = {"knee": crit("knee", True), "gm1": crit("gm1", True), "registered_ultra": crit("registered_ultra", False)}
    oof = json.loads((ws / "runs/fastfit/oof/gate.json").read_text())["criteria_oof"]
    k_oof = oof.get("knee", {})
    crits["knee"]["c1_c2_oof"] = k_oof
    crits["knee"]["all_four_pass"] = bool(k_oof and all(k_oof.values()) and crits["knee"]["c3_pass"]
                                          and crits["knee"]["c4_pass"])
    crits["registered_ultra"]["oof"] = oof.get("registered_ultra")
    crits["decision_exception"] = {"knee_fails_c4": not crits["knee"]["c4_pass"],
                                   "ultra_completes_all": crits["registered_ultra"]["c4_pass"],
                                   "ultra_reconsidered": (not crits["knee"]["c4_pass"])
                                   and crits["registered_ultra"]["c4_pass"]}
    ident = [json.loads(l) for l in open(out / "cli_identity.jsonl")]
    idn = {"n": len(ident), "identical": sum(r["identical"] for r in ident),
           "by_k": {str(k): f"{sum(r['identical'] for r in ident if r['k'] == k)}/"
                            f"{sum(1 for r in ident if r['k'] == k)}" for k in SMALL_K},
           "errors": sum(bool(r["error"]) for r in ident)}
    res = {"arms": rm["arms"], "cases": {k: len(v) for k, v in rm["cases"].items()}, "table": table,
           "criteria": crits, "release_cli_identity": idn, "n_boot": N_BOOT, "seed": SEED,
           "rss_pass": {"n": len(rssm), "same_status_as_run": sum(r["same_status_as_run"] for r in rssm.values()),
                        "same_output_as_run": sum(r["same_output_as_run"] for r in rssm.values()),
                        "errors": sum(bool(r["error"]) for r in rssm.values())},
           "notes": ["Q = mean over families of (SP + TC) / 2; failed case SP = TC = 0",
                     "1,000/10,000: 1 thread, CPU = per-case user+sys; 50,000/200,000: 8 threads, CPU = process "
                     "user+sys of the scored run; peak RSS = the aligner's own ru_maxrss from the shim re-measurement "
                     "(rss.jsonl), reported for 50,000/200,000 only (incident: the scored run's values were inherited "
                     "from the harness process)",
                     "bootstrap intervals of the differences are reported, not gated"],
           "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _atomic_json(out / "gate.json", res)
    print(f"{'arm':17s} {'k':>6s} {'n':>3s} {'fail':>4s} {'SP':>6s} {'TC':>6s} {'Q':>6s} {'dQ':>7s} "
          f"{'CPU s':>9s} {'wall s':>9s} {'RSS GB':>6s}")
    for a in ARMS:
        for k in cases:
            t = table[a][str(k)]
            dq = f"{t['vs_p400']['q']['diff']:+.4f}" if "vs_p400" in t else ""
            rss = f"{t['peak_rss_gb_max']:.1f}" if t["peak_rss_gb_max"] else "-"
            print(f"{a:17s} {k:6d} {t['n']:3d} {t['n_failed']:4d} {t['sp']:6.3f} {t['tc']:6.3f} {t['q']:6.3f} "
                  f"{dq:>7s} {t['cpu_s']:9.0f} {t['wall_s']:9.0f} {rss:>6s}")
    print(json.dumps({"criteria": crits, "release_cli_identity": idn, "rss_pass": res["rss_pass"]}, indent=1))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["run", "rss", "gate"])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--ws", default=str(Path("~/tmp/kalign-fastfit").expanduser()))
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    return run(ws, a.dry_run) if a.cmd == "run" else measure_rss(ws) if a.cmd == "rss" else gate(ws)


if __name__ == "__main__":
    sys.exit(main())
