"""Evaluate a spec (CLI preset or parameter dict) on a set of cases through the cache.

    evaluate(case_rows, spec, build, workers) -> list of per-case result rows

spec:  "fast" | "default" | "recall" | "accurate" | {"preset": mode}   -> CLI, ``kalign --mode M -n 1``
       a decoded parameter dict (see recipe/_align_worker.py)         -> worker process in the venv
Every case goes through the cache: hit -> returned as is; miss -> align into a temp file, score with
scoring.score_case(case_row, aln_path), store, delete the temp alignment (unless keep_dir is set).
CPU seconds: CLI = user+sys of the child via os.wait4; params = RUSAGE_SELF delta in the worker.
Everything is single-threaded (OMP_NUM_THREADS=1, -n 1) and seed 42.

Cache identity of a result: (build_id, param_key, case_id, sha256_unaligned, sha256_reference,
scorer_version); scorer_version is the SCORER_VERSION constant of the pool scorer module (see cache.py).

Timeouts: every aligner call (CLI, per-case worker, batched worker) is killed after ``timeout_s``
(default 1800 s wall of that call; wall >= CPU for a single-threaded run, so this is a slightly
conservative stand-in for 30 min CPU) and stored as error "timeout" (CPU = what the killed process used).

Batched parameter-set worker: ``evaluate(..., batch=True)`` runs ``workers`` long-lived worker processes
(one per slot, recycled every ``recycle`` cases); each aligns many cases, one request at a time over a
pipe, and measures per-case CPU as a RUSAGE_SELF delta around each ensemble_custom_file_to_file call.
Presets (CLI) have no worker and ignore ``batch``.

CLI:  python -m recipe.evaluate --cases TSV --spec preset:default --workers 30 --out TSV
      python -m recipe.evaluate --init-build --cli PATH --python PATH [--wheel PATH] --build-dir DIR
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import importlib
import json
import multiprocessing as mp
import os
import queue
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import traceback
from dataclasses import dataclass
from pathlib import Path

from . import cache as C

HARNESS = "evaluate-2"
DEFAULT_TIMEOUT_S = 1800.0
DEFAULT_RECYCLE = 500
SEED = 42
MODES = ("fast", "default", "recall", "accurate")
SINGLE_THREAD_ENV = {"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1", "OMP_DYNAMIC": "FALSE",
                     "OMP_WAIT_POLICY": "PASSIVE"}
MAX_WORKERS = 30
WORKER_SCRIPT = Path(__file__).with_name("_align_worker.py")
OUT_COLS = ("case_id", "pool", "benchmark", "subset", "recall", "precision", "f1", "tc", "width_ratio",
            "n_seq", "cpu_s", "wall_s", "aln_sha256", "error", "cached")


# ---------------------------------------------------------------- build description
@dataclass(frozen=True)
class Build:
    build_id: str          # sha256 of the CLI binary
    cli: str               # path to the kalign CLI
    python: str            # python of the aligner venv (kalign wheel installed)
    wheel_sha256: str = ""

    @classmethod
    def create(cls, cli, python, wheel=None):
        return cls(C.sha256_file(cli), str(Path(cli).resolve()), str(Path(python).absolute()),
                   C.sha256_file(wheel) if wheel else "")

    def save(self, build_dir):
        Path(build_dir, "build.json").write_text(json.dumps(self.__dict__, indent=1) + "\n")

    @classmethod
    def load(cls, build_dir):
        d = json.loads(Path(build_dir, "build.json").read_text())
        b = cls(**d)
        if C.sha256_file(b.cli) != b.build_id:
            raise RuntimeError(f"CLI {b.cli} does not match build_id {b.build_id[:12]} in build.json")
        return b


def normalise_spec(spec):
    """Return ("preset", mode) or ("params", dict)."""
    if isinstance(spec, str):
        spec = {"preset": spec}
    if set(spec) == {"preset"}:
        if spec["preset"] not in MODES:
            raise ValueError(f"unknown preset {spec['preset']!r}")
        return "preset", spec["preset"]
    return "params", spec


# ---------------------------------------------------------------- timeouts
def _kill_group(pid):
    try:
        os.killpg(pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def _run_killable(argv, timeout_s):
    """Run argv (own session) and kill it after timeout_s. -> (wall_s, cpu_s, rc_or_None, stderr_text).
    rc None = killed by the timeout. CPU = user+sys of the child via os.wait4."""
    t0 = time.perf_counter()
    with tempfile.TemporaryFile() as errf:
        p = subprocess.Popen(argv, stdout=subprocess.DEVNULL, stderr=errf, start_new_session=True,
                             env={**os.environ, **SINGLE_THREAD_ENV})
        lock, st = threading.Lock(), {"done": False, "killed": False}

        def _timeout():
            with lock:
                if not st["done"]:
                    st["killed"] = True
                    _kill_group(p.pid)

        timer = None
        if timeout_s:
            timer = threading.Timer(timeout_s, _timeout)
            timer.daemon = True
            timer.start()
        _, status, ru = os.wait4(p.pid, 0)
        with lock:
            st["done"] = True
        if timer:
            timer.cancel()
        p.returncode = os.waitstatus_to_exitcode(status)
        errf.seek(0)
        txt = errf.read().decode(errors="replace")
    wall = time.perf_counter() - t0
    return wall, ru.ru_utime + ru.ru_stime, (None if st["killed"] else p.returncode), txt


def _tail(txt, n=300):
    return txt[-n:].replace("\t", " ").replace("\n", " ")


# ---------------------------------------------------------------- aligners
def run_cli(cli, mode, inp, out, timeout_s=DEFAULT_TIMEOUT_S):
    """-> (wall_s, cpu_s, error). CPU from wait4 rusage of the child (user+sys). error "timeout" if killed."""
    wall, cpu, rc, txt = _run_killable([cli, "--mode", mode, "-n", "1", "-i", str(inp), "-o", str(out)],
                                       timeout_s)
    if rc is None:
        return wall, cpu, "timeout"
    return wall, cpu, (None if rc == 0 else f"exit {rc}: {_tail(txt)}")


def run_params(python, params, inp, out, tmp, timeout_s=DEFAULT_TIMEOUT_S):
    """-> (wall_s, cpu_s, error) via a one-shot worker process in the aligner venv (per-case path)."""
    job, res = Path(tmp, "job.json"), Path(tmp, "result.json")
    job.write_text(json.dumps({"params": params, "cases": [[str(inp), str(out)]], "seed": SEED}))
    wall, cpu, rc, txt = _run_killable([python, str(WORKER_SCRIPT), str(job), str(res)], timeout_s)
    if rc is None:
        return wall, cpu, "timeout"
    if rc != 0 or not res.exists():
        return wall, None, f"worker exit {rc}: {_tail(txt)}"
    r = json.loads(res.read_text())["rows"][0]
    return r["wall"], r["cpu"], r["error"]


class ParamWorker:
    """One long-lived _align_worker.py --serve process for one parameter dict (see module doc).

    align(inp, out, timeout_s) -> (wall_s, cpu_s, error). On a timeout the worker process group is killed
    (error "timeout", CPU = what the process used since its last reported case) and a fresh worker is
    started for the next call; the same on a crash (error "worker exit ..."). Recycled after ``recycle``
    cases to bound memory growth of the long-lived process."""

    def __init__(self, python, params, recycle=DEFAULT_RECYCLE):
        self.python, self.params, self.recycle = python, params, recycle
        self.proc = None
        self.tmp = None
        self.errf = None
        self.buf = b""
        self.n = 0
        self.cpu_reported = 0.0

    # -- lifecycle
    def _start(self, timeout_s):
        self.tmp = tempfile.mkdtemp(prefix="pw_", dir=os.environ.get("KALIGN_REFIT_TMP"))
        pj = Path(self.tmp, "params.json")
        pj.write_text(json.dumps({"params": self.params, "seed": SEED}))
        self.errf = tempfile.TemporaryFile()
        self.proc = subprocess.Popen([self.python, str(WORKER_SCRIPT), "--serve", str(pj)],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.errf,
                                     bufsize=0, start_new_session=True, env={**os.environ, **SINGLE_THREAD_ENV})
        self.buf, self.n, self.cpu_reported = b"", 0, 0.0
        line = self._readline(time.monotonic() + timeout_s)  # {"ready": true}
        if not json.loads(line).get("ready"):
            raise RuntimeError(f"worker did not become ready: {line!r}")

    def _stderr_tail(self):
        try:
            self.errf.seek(0)
            return _tail(self.errf.read().decode(errors="replace"))
        except Exception:  # noqa: BLE001
            return ""

    def _reap(self, kill):
        """Stop the worker; -> CPU (user+sys) it used in total."""
        p, cpu = self.proc, 0.0
        if p is not None:
            if kill:
                _kill_group(p.pid)
            else:
                try:
                    p.stdin.close()
                except OSError:
                    pass
            try:
                _, status, ru = os.wait4(p.pid, 0)
                p.returncode = os.waitstatus_to_exitcode(status)
                cpu = ru.ru_utime + ru.ru_stime
            except ChildProcessError:
                pass
            for f in (p.stdout, p.stdin):
                try:
                    f and f.close()
                except OSError:
                    pass
        self.proc = None
        if self.errf:
            self.errf.close()
            self.errf = None
        if self.tmp:
            shutil.rmtree(self.tmp, ignore_errors=True)
            self.tmp = None
        return cpu

    def close(self):
        self._reap(kill=False)

    # -- protocol
    def _readline(self, deadline):
        fd = self.proc.stdout.fileno()
        while b"\n" not in self.buf:
            rem = deadline - time.monotonic()
            if rem <= 0:
                raise TimeoutError
            if not select.select([fd], [], [], rem)[0]:
                raise TimeoutError
            chunk = os.read(fd, 65536)
            if not chunk:
                raise EOFError
            self.buf += chunk
        line, _, self.buf = self.buf.partition(b"\n")
        return line

    def align(self, inp, out, timeout_s=DEFAULT_TIMEOUT_S):
        t0 = time.perf_counter()
        try:
            if self.proc is None:
                self._start(timeout_s)
            os.write(self.proc.stdin.fileno(), (json.dumps({"in": str(inp), "out": str(out)}) + "\n").encode())
            r = json.loads(self._readline(time.monotonic() + timeout_s))
        except TimeoutError:
            used = self._reap(kill=True) - self.cpu_reported
            return time.perf_counter() - t0, max(used, 0.0), "timeout"
        except (EOFError, BrokenPipeError, OSError, ValueError, RuntimeError) as e:
            tail = self._stderr_tail()
            used = self._reap(kill=True) - self.cpu_reported
            return (time.perf_counter() - t0, max(used, 0.0),
                    f"worker exit ({type(e).__name__}): {tail}")
        self.n += 1
        self.cpu_reported += r["cpu"]
        if self.recycle and self.n >= self.recycle:
            self.close()
        return r["wall"], r["cpu"], r["error"]


# ---------------------------------------------------------------- one case
def _load_scorer(name):
    mod, fn = name.split(":")
    return getattr(importlib.import_module(mod), fn)


def scorer_version(scorer_name, case_row) -> str:
    """SCORER_VERSION of the module that scores this case (the pool module for recipe.scoring)."""
    mod = scorer_name.split(":")[0]
    if mod == "recipe.scoring":
        from .scoring import POOL_MODULES
        pool = case_row["pool"]
        if pool not in POOL_MODULES:
            raise ValueError(f"unknown pool {pool!r} for case {case_row.get('case_id')}")
        mod = POOL_MODULES[pool]
    v = getattr(importlib.import_module(mod), "SCORER_VERSION", None)
    if not v:
        raise ValueError(f"scorer module {mod} defines no SCORER_VERSION (needed for the cache key)")
    return str(v)


def _process_case(build, spec, case_row, scorer_name, keep_dir, cache, retry_errors, align):
    """Cache lookup, align (``align(inp, out, tmp) -> (wall, cpu, err)``), score, store. -> result row."""
    key = C.param_key(spec)
    cid = case_row["case_id"]
    idt = C.ident(case_row, scorer_version(scorer_name, case_row))
    hit = cache.get(build.build_id, key, idt)
    if hit is not None and not (retry_errors and hit["error"]):
        return {**hit, "cached": 1}
    tmp = tempfile.mkdtemp(prefix=f"ev_{cid[:40].replace('/', '_')}_", dir=os.environ.get("KALIGN_REFIT_TMP"))
    try:
        out = Path(tmp, "aln.fa")
        wall, cpu, err = align(case_row["unaligned"], out, tmp)
        metrics, sha = {}, None
        if err is None:
            try:
                sha = C.sha256_file(out)
                metrics = _load_scorer(scorer_name)(case_row, str(out)) or {}
            except Exception as e:  # noqa: BLE001 - scoring failure is recorded, not raised
                err = f"score {type(e).__name__}: {e}"[:400]
        if keep_dir and out.exists():
            d = Path(keep_dir, key[:16])
            d.mkdir(parents=True, exist_ok=True)
            shutil.copy(out, d / f"{cid.replace('/', '_')}.aln.fa")
        cache.put(build.build_id, key, idt, metrics, cpu_s=cpu, wall_s=wall, aln_sha256=sha, error=err,
                  spec=spec, harness=HARNESS, replace_error=bool(hit and retry_errors))
        row = cache.get(build.build_id, key, idt)  # what is stored is what is returned
        return {**row, "cached": 0}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def eval_one(build_d, spec, case_row, scorer_name, keep_dir, cache_path, retry_errors,
             timeout_s=DEFAULT_TIMEOUT_S):
    """Worker-process body (per-case path): one aligner process per case."""
    build = Build(**build_d)
    kind, val = normalise_spec(spec)

    def align(inp, out, tmp):
        if kind == "preset":
            return run_cli(build.cli, val, inp, out, timeout_s)
        return run_params(build.python, val, inp, out, tmp, timeout_s)

    return _process_case(build, spec, case_row, scorer_name, keep_dir, C.Cache(cache_path), retry_errors, align)


def _slot(build_d, spec, scorer_name, keep_dir, cache_path, retry_errors, timeout_s, recycle, q, done_q):
    """One pool slot of the batched path: a long-lived ParamWorker aligning the cases it pulls from q."""
    build = Build(**build_d)
    cache = C.Cache(cache_path)
    worker = ParamWorker(build.python, normalise_spec(spec)[1], recycle)
    try:
        while True:
            c = q.get()
            if c is None:
                break
            try:
                _process_case(build, spec, c, scorer_name, keep_dir, cache, retry_errors,
                              lambda inp, out, tmp: worker.align(inp, out, timeout_s))
                done_q.put((c["case_id"], None))
            except Exception:  # noqa: BLE001 - reported to the parent, which raises
                done_q.put((c["case_id"], traceback.format_exc()[-500:]))
    finally:
        worker.close()


def _cost(c):
    try:
        return int(c["n_seq"]) * int(c["ref_len"])
    except (KeyError, ValueError, TypeError):
        return 0


# ---------------------------------------------------------------- batch
def evaluate(case_rows, spec, build: Build, workers=1, *, cache_path=None, keep_dir=None,
             scorer="recipe.scoring:score_case", retry_errors=False, progress=False,
             timeout_s=DEFAULT_TIMEOUT_S, batch=False, recycle=DEFAULT_RECYCLE):
    """Evaluate ``spec`` on every case; returns rows in case order (dicts incl. pool/benchmark/subset).

    batch=True: parameter dicts are aligned by ``workers`` long-lived worker processes (cheapest-last
    scheduling); results are identical to the per-case path. Ignored for CLI presets."""
    if workers > MAX_WORKERS:
        raise ValueError(f"workers={workers} exceeds the limit of {MAX_WORKERS} concurrent aligners")
    kind, _ = normalise_spec(spec)
    cache = C.Cache(cache_path)
    key = C.param_key(spec)
    cases = list(case_rows)
    ids = [c["case_id"] for c in cases]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate case_id in case list")
    idents = {c["case_id"]: C.ident(c, scorer_version(scorer, c)) for c in cases}
    hits = cache.get_many(build.build_id, key, list(idents.values()))
    todo = [c for c in cases if c["case_id"] not in hits or (retry_errors and hits[c["case_id"]]["error"])]
    results = {cid: {**h, "cached": 1} for cid, h in hits.items()}
    args = (build.__dict__, spec)
    if todo and batch and kind == "params":
        _run_batched(args, todo, scorer, keep_dir, cache, retry_errors, timeout_s, recycle, workers, progress)
        got = cache.get_many(build.build_id, key, [idents[c["case_id"]] for c in todo])
        for c in todo:
            results[c["case_id"]] = {**got[c["case_id"]], "cached": 0}
    elif todo:
        tail = (scorer, str(keep_dir) if keep_dir else None, str(cache.path), retry_errors, timeout_s)
        if workers <= 1:
            for c in todo:
                results[c["case_id"]] = eval_one(*args, c, *tail)
        else:
            ctx = mp.get_context("fork")  # no __main__ guard needed by callers; workers open their own DB connections
            with cf.ProcessPoolExecutor(workers, mp_context=ctx) as ex:
                futs = {ex.submit(eval_one, *args, c, *tail): c["case_id"] for c in todo}
                for i, f in enumerate(cf.as_completed(futs), 1):
                    results[futs[f]] = f.result()
                    if progress and i % 25 == 0:
                        print(f"  {i}/{len(todo)} aligned", file=sys.stderr, flush=True)
    out = []
    for c in cases:
        r = results[c["case_id"]]
        out.append({**{k: c.get(k) for k in ("pool", "benchmark", "subset")}, **r})
    return out


def _run_batched(args, todo, scorer, keep_dir, cache, retry_errors, timeout_s, recycle, workers, progress):
    ctx = mp.get_context("fork")
    todo = sorted(todo, key=_cost, reverse=True)  # expensive first: better packing
    n = max(1, min(workers, len(todo)))
    q, done_q = ctx.Queue(), ctx.Queue()
    for c in todo:
        q.put(c)
    for _ in range(n):
        q.put(None)
    procs = [ctx.Process(target=_slot, args=(*args, scorer, str(keep_dir) if keep_dir else None,
                                             str(cache.path), retry_errors, timeout_s, recycle, q, done_q))
             for _ in range(n)]
    for p in procs:
        p.start()
    n_done, failures = 0, []
    try:
        while n_done < len(todo):
            try:
                cid, tb = done_q.get(timeout=2)
            except queue.Empty:
                if not any(p.is_alive() for p in procs):
                    break
                continue
            n_done += 1
            if tb:
                failures.append((cid, tb))
            if progress and n_done % 25 == 0:
                print(f"  {n_done}/{len(todo)} aligned", file=sys.stderr, flush=True)
    finally:
        for p in procs:
            if p.is_alive() and n_done < len(todo):
                p.terminate()
        for p in procs:
            p.join()
    if failures:
        raise RuntimeError(f"{len(failures)} cases failed in the batch slots, first: {failures[0]}")
    if n_done < len(todo):
        raise RuntimeError(f"batch workers died: {n_done}/{len(todo)} cases finished")


def write_rows(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(OUT_COLS)
        for r in rows:
            w.writerow(["" if r.get(c) is None else r.get(c) for c in OUT_COLS])


# ---------------------------------------------------------------- CLI
def _parse_spec(s):
    if s.startswith("preset:"):
        return s.split(":", 1)[1]
    if s.startswith("params:"):
        return json.loads(Path(s.split(":", 1)[1]).read_text())
    raise SystemExit("--spec must be preset:<mode> or params:<file.json>")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--init-build", action="store_true", help="write build.json for a built aligner")
    ap.add_argument("--cli"); ap.add_argument("--python"); ap.add_argument("--wheel")
    ap.add_argument("--build-dir", help="dir with build.json (aligner build)")
    ap.add_argument("--cases"); ap.add_argument("--spec"); ap.add_argument("--out")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--cache"); ap.add_argument("--keep-alignments")
    ap.add_argument("--retry-errors", action="store_true")
    ap.add_argument("--timeout-s", type=float, default=DEFAULT_TIMEOUT_S)
    ap.add_argument("--batch", action="store_true", help="long-lived parameter-set workers")
    a = ap.parse_args(argv)
    if a.init_build:
        b = Build.create(a.cli, a.python, a.wheel)
        b.save(a.build_dir)
        print(json.dumps(b.__dict__, indent=1))
        return 0
    if not (a.cases and a.spec and a.out and a.build_dir):
        ap.error("--cases, --spec, --out and --build-dir are required")
    from .scoring import read_manifest
    build = Build.load(a.build_dir)
    rows = evaluate(read_manifest(a.cases), _parse_spec(a.spec), build, a.workers, cache_path=a.cache,
                    keep_dir=a.keep_alignments, retry_errors=a.retry_errors, progress=True,
                    timeout_s=a.timeout_s, batch=a.batch)
    write_rows(rows, a.out)
    n_hit = sum(r["cached"] for r in rows)
    n_err = sum(1 for r in rows if r.get("error"))
    print(f"{len(rows)} cases, {n_hit} cached, {n_err} errors -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
