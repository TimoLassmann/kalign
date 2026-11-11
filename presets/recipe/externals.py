"""External aligners (MAFFT L-INS-i, MUSCLE 5, Clustal Omega) run on case rows through the cache.

    run_external(tool, case_row, ...)       one case, cache-aware, returns the cache row (+ "cached")
    evaluate_external(case_rows, tool, workers=...)   many cases, process pool

Comparison only: these tools are never tuned. Everything is single threaded; CPU seconds are
user+sys from ``os.wait4`` of the tool's process (descendants that the tool waited for are included,
so MAFFT's wrapper script and its helper binaries are counted). Timeout 3 h per case: the whole
process group is killed and the case is recorded as an error row (with the CPU spent).

Cache key: build_id = sha256 over the tool's executable files, param_key = param_key({"external":
tool, "version": version}); the version string also encodes the command line / policy, so a changed
command is a different key. The tool output is normalised before scoring (identical for all
tools): upper case, ``.``/``*`` -> ``-``, one line per sequence. Scoring is recipe.scoring.score_case,
so every integrity check of the scorers (all input sequences present and unchanged) applies.

Command lines (all single thread, protein or nucleotide auto-detected by the tool):
  mafft_linsi  mafft --localpair --maxiterate 1000 --thread 1 --quiet IN          (L-INS-i)
  muscle5      muscle -align IN -output OUT -threads 1      if n_seq <= 200
               muscle -super5 IN -output OUT -threads 1     if n_seq  > 200
  clustalo     clustalo -i IN -o OUT --outfmt=fa --threads=1 --force   (defaults otherwise)

  kalign361_<mode>  kalign -i IN -o OUT --mode <mode> -n 1     Kalign 3.6.1 release (commit c007ad0), CLI
  kalign370_<mode>  kalign -i IN -o OUT --mode <mode> -n 1     Kalign 3.7.0 CLI of the fitting build (960a4f6)
               (Addendum A1 of PLAN.md: descriptive comparison of releases, same harness for both; not part of
               the decision; mode in fast, default, recall, accurate)
  kalign31          kalign -i IN -o OUT                         Kalign 3.1 (tag v3.1, Oct 2019), defaults (Addendum A2)
  famsa             famsa -t 1 IN OUT                           FAMSA 2.5.2, defaults (Addendum A3, Kalign 4 paper)
  kalign335         kalign -i IN -o OUT -n 1                    Kalign 3.3.5 (tag v3.3.5, Nov 2022), defaults (Addendum A2)

Tools live in <WS>/tools (env KALIGN_REFIT_TOOLS overrides); see presets/docs/build/externals.md. Kalign builds live in
<WS>/kalign/c007ad0-3.6.1 and <WS>/kalign/960a4f6 (build/src/kalign).

CLI:
  python -m recipe.externals --info                      JSON: path, version, sha256 of every tool
  python -m recipe.externals --register INFO.json        merge that info into presets/resources/REGISTRY.yaml
  python -m recipe.externals --tool all --cases TSV --workers 4 --out TSV [--limit-per-pool N]
  python -m recipe.externals --smoke --protein TSV --nucleotide TSV --n 20 --out TSV
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import hashlib
import json
import multiprocessing as mp
import os
import random
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from . import cache as C

HARNESS = "externals-1"
TIMEOUT_S = 3 * 3600
MUSCLE_ALIGN_MAX_SEQ = 200
SINGLE_THREAD_ENV = {"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1", "OMP_DYNAMIC": "FALSE",
                     "OMP_WAIT_POLICY": "PASSIVE", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
TOOL_NAMES = ("mafft_linsi", "muscle5", "clustalo")  # "famsa" (Addendum A3) is run by name only
OLD_KALIGN = {"kalign31": ("v3.1", "kalign31"), "kalign335": ("v3.3.5", "build/src/kalign")}
KALIGN_MODES = ("fast", "default", "recall", "accurate")
KALIGN_BUILDS = {"kalign361": "c007ad0-3.6.1", "kalign370": "960a4f6"}
KALIGN_TOOLS = tuple(f"{b}_{m}" for b in KALIGN_BUILDS for m in KALIGN_MODES)


def _kalign_tool(tool: str):
    """'kalign361_accurate' -> ('kalign361', 'accurate'), or None for other tools."""
    b, _, m = tool.partition("_")
    return (b, m) if b in KALIGN_BUILDS and m in KALIGN_MODES else None
OUT_COLS = ("tool", "case_id", "pool", "benchmark", "subset", "recall", "precision", "f1", "tc",
            "width_ratio", "n_seq", "cpu_s", "wall_s", "aln_sha256", "error", "cached")


def tools_dir() -> Path:
    return Path(os.environ.get("KALIGN_REFIT_TOOLS") or (C.default_ws() / "tools")).expanduser()


# ---------------------------------------------------------------- tool descriptions
def _run_text(cmd, env=None) -> str:
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return (p.stdout + p.stderr).strip()


def tool_paths(tool: str, tdir: Path | None = None) -> dict:
    """Executables and hashed files of a tool."""
    t = tdir or tools_dir()
    if tool in OLD_KALIGN:
        exe = t.parent / "kalign" / OLD_KALIGN[tool][0] / OLD_KALIGN[tool][1]
        return {"exe": exe, "files": [exe], "env": {}}
    kt = _kalign_tool(tool)
    if kt:
        exe = t.parent / "kalign" / KALIGN_BUILDS[kt[0]] / "build" / "src" / "kalign"
        return {"exe": exe, "files": [exe], "env": {}}
    if tool == "mafft_linsi":
        root = t / "mafft"
        files = [root / "bin" / "mafft"] + sorted((root / "libexec" / "mafft").glob("*"))
        return {"exe": root / "bin" / "mafft", "files": [f for f in files if f.is_file()],
                "env": {"MAFFT_BINARIES": str(root / "libexec" / "mafft")}}
    if tool == "muscle5":
        return {"exe": t / "bin" / "muscle", "files": [t / "bin" / "muscle"], "env": {}}
    if tool == "famsa":
        return {"exe": t / "bin" / "famsa", "files": [t / "bin" / "famsa"], "env": {}}
    if tool == "clustalo":
        exe = t / "bin" / "clustalo"
        if not exe.exists():
            exe = Path(shutil.which("clustalo") or exe)
        return {"exe": exe, "files": [exe], "env": {}}
    raise ValueError(f"unknown tool {tool!r}")


def tool_version(tool: str, tdir: Path | None = None) -> str:
    """Version string; it also encodes the command line so that it is part of the cache key."""
    p = tool_paths(tool, tdir)
    exe = str(p["exe"])
    env = {**os.environ, **p["env"]}
    if tool in OLD_KALIGN:
        v = [ln for ln in _run_text([exe, "--version"], env).splitlines() if "align" in ln.lower()][:1]
        return f"{v[0].strip() if v else tool}; defaults" + ("; -n 1" if tool == "kalign335" else "")
    kt = _kalign_tool(tool)
    if kt:
        v = _run_text([exe, "--version"], env).split("\n")[0].strip()
        return f"{v}; --mode {kt[1]} -n 1"
    if tool == "mafft_linsi":
        v = _run_text([exe, "--version"], env).split("\n")[-1].strip()
        return f"{v}; L-INS-i: --localpair --maxiterate 1000 --thread 1"
    if tool == "muscle5":
        v = _run_text([exe, "-version"], env).split("\n")[0].strip()
        return f"{v}; -align if n_seq<={MUSCLE_ALIGN_MAX_SEQ} else -super5; -threads 1"
    if tool == "clustalo":
        v = _run_text([exe, "--version"], env).split("\n")[0].strip()
        return f"{v}; defaults; --threads=1"
    if tool == "famsa":
        v = [ln.strip() for ln in _run_text([exe], env).splitlines() if "version" in ln][:1]
        return f"FAMSA {v[0] if v else '?'}; defaults; -t 1"
    raise ValueError(tool)


def tool_build_id(tool: str, tdir: Path | None = None) -> str:
    """sha256 over (name, sha256) of every executable file of the tool."""
    h = hashlib.sha256()
    for f in tool_paths(tool, tdir)["files"]:
        h.update(f.name.encode() + b"\0" + C.sha256_file(f).encode() + b"\n")
    return h.hexdigest()


def tool_spec(tool: str, tdir: Path | None = None) -> dict:
    return {"external": tool, "version": tool_version(tool, tdir)}


def command(tool: str, inp, out, n_seq, tdir: Path | None = None):
    """-> (argv, stdout_path_or_None). MAFFT writes the alignment to stdout."""
    exe = str(tool_paths(tool, tdir)["exe"])
    if tool == "kalign31":
        return [exe, "-i", str(inp), "-o", str(out)], None
    if tool == "kalign335":
        return [exe, "-i", str(inp), "-o", str(out), "-n", "1"], None
    kt = _kalign_tool(tool)
    if kt:
        return [exe, "-i", str(inp), "-o", str(out), "--mode", kt[1], "-n", "1"], None
    if tool == "mafft_linsi":
        return [exe, "--localpair", "--maxiterate", "1000", "--thread", "1", "--quiet", str(inp)], out
    if tool == "muscle5":
        mode = "-align" if int(n_seq) <= MUSCLE_ALIGN_MAX_SEQ else "-super5"
        return [exe, mode, str(inp), "-output", str(out), "-threads", "1"], None
    if tool == "clustalo":
        return [exe, "-i", str(inp), "-o", str(out), "--outfmt=fa", "--threads=1", "--force"], None
    if tool == "famsa":
        return [exe, "-t", "1", str(inp), str(out)], None
    raise ValueError(tool)


# ---------------------------------------------------------------- running
def run_process(argv, stdout_path, env, timeout_s):
    """-> (wall_s, cpu_s, error). CPU = user+sys from wait4. Kills the process group on timeout."""
    t0 = time.perf_counter()
    sink = open(stdout_path, "wb") if stdout_path else subprocess.DEVNULL
    try:
        p = subprocess.Popen(argv, stdout=sink, stderr=subprocess.PIPE, env=env, start_new_session=True)
    finally:
        if stdout_path:
            sink.close()  # child holds its own copy
    timed_out = [False]

    def _kill():
        timed_out[0] = True
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    import threading
    timer = threading.Timer(timeout_s, _kill)
    timer.start()
    try:
        err_txt = p.stderr.read().decode(errors="replace")  # returns at exit (or after the kill)
        _, status, ru = os.wait4(p.pid, 0)
    finally:
        timer.cancel()
        p.stderr.close()
        try:  # orphaned helper processes of a killed tool
            os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    wall = time.perf_counter() - t0
    cpu = ru.ru_utime + ru.ru_stime
    rc = os.waitstatus_to_exitcode(status)
    if timed_out[0]:
        return wall, cpu, f"timeout after {timeout_s:g} s"
    if rc != 0:
        return wall, cpu, f"exit {rc}: {err_txt[-300:]}".replace("\t", " ").replace("\n", " ")
    return wall, cpu, None


def normalise_alignment(path) -> int:
    """Rewrite ``path`` in place: upper case, non-letters -> '-', unwrapped. Returns n sequences."""
    recs, name, buf = [], None, []
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n{''.join(c.upper() if c.isalpha() else '-' for c in s)}\n")
    return len(recs)


def scorer_version(scorer: str, case_row) -> str:
    """SCORER_VERSION of the module that scores this case (cache key part). The pool scorer modules
    define it; "<pool>-v1" is the step-0 label used until they do."""
    import importlib
    mod = importlib.import_module(scorer.split(":")[0])
    if mod.__name__ == "recipe.scoring":
        mod = importlib.import_module(mod.POOL_MODULES[case_row["pool"]])
    return getattr(mod, "SCORER_VERSION", None) or f"{case_row.get('pool', 'x')}-v1"


def run_external(tool, case_row, *, cache_path=None, scorer="recipe.scoring:score_case",
                 timeout_s=TIMEOUT_S, retry_errors=False, tdir=None, keep_dir=None):
    """One case through the cache. Returns the stored row (dict) with ``cached`` 0/1."""
    import importlib
    cache = C.Cache(cache_path)
    spec = tool_spec(tool, tdir)
    key, build_id, cid = C.param_key(spec), tool_build_id(tool, tdir), case_row["case_id"]
    ident = C.ident(case_row, scorer_version(scorer, case_row))
    hit = cache.get(build_id, key, ident)
    if hit is not None and not (retry_errors and hit["error"]):
        return {**hit, "cached": 1}
    tmp = tempfile.mkdtemp(prefix=f"ext_{tool}_{cid[:30].replace('/', '_')}_",
                           dir=os.environ.get("KALIGN_REFIT_TMP"))
    try:
        out = Path(tmp, "aln.fa")
        argv, stdout_path = command(tool, case_row["unaligned"], out, case_row["n_seq"], tdir)
        env = {**os.environ, **SINGLE_THREAD_ENV, **tool_paths(tool, tdir)["env"], "TMPDIR": tmp}
        wall, cpu, err = run_process(argv, stdout_path, env, timeout_s)
        metrics, sha = {}, None
        if err is None:
            try:
                if not out.exists() or out.stat().st_size == 0:
                    raise ValueError("no alignment written")
                normalise_alignment(out)
                sha = C.sha256_file(out)
                mod, fn = scorer.split(":")
                metrics = getattr(importlib.import_module(mod), fn)(case_row, str(out)) or {}
                if metrics.get("integrity_bad"):  # nucleotide scorer records instead of raising
                    raise ValueError(f"integrity_bad={metrics['integrity_bad']}")
            except Exception as e:  # noqa: BLE001 - recorded, not raised
                err = f"score {type(e).__name__}: {e}"[:400]
        if keep_dir and out.exists():
            d = Path(keep_dir, tool)
            d.mkdir(parents=True, exist_ok=True)
            shutil.copy(out, d / f"{cid.replace('/', '_')}.aln.fa")
        cache.put(build_id, key, ident, metrics, cpu_s=cpu, wall_s=wall, aln_sha256=sha, error=err,
                  spec=spec, harness=HARNESS, replace_error=bool(hit and retry_errors))
        return {**cache.get(build_id, key, ident), "cached": 0}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def evaluate_external(case_rows, tool, workers=1, *, cache_path=None, scorer="recipe.scoring:score_case",
                      timeout_s=TIMEOUT_S, retry_errors=False, tdir=None, keep_dir=None, progress=False):
    """Rows in case order; each carries pool/benchmark/subset/tool."""
    cases = list(case_rows)
    ids = [c["case_id"] for c in cases]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate case_id in case list")
    kw = dict(cache_path=cache_path, scorer=scorer, timeout_s=timeout_s, retry_errors=retry_errors,
              tdir=tdir, keep_dir=keep_dir)
    results = {}
    if workers <= 1:
        for c in cases:
            results[c["case_id"]] = run_external(tool, c, **kw)
    else:
        ctx = mp.get_context("fork")
        with cf.ProcessPoolExecutor(workers, mp_context=ctx) as ex:
            futs = {ex.submit(run_external, tool, c, **kw): c["case_id"] for c in cases}
            for i, f in enumerate(cf.as_completed(futs), 1):
                results[futs[f]] = f.result()
                if progress and i % 25 == 0:
                    print(f"  {tool}: {i}/{len(cases)}", file=sys.stderr, flush=True)
    return [{"tool": tool, **{k: c.get(k) for k in ("pool", "benchmark", "subset")}, **results[c["case_id"]]}
            for c in cases]


def write_rows(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(OUT_COLS)
        for r in rows:
            w.writerow(["" if r.get(c) is None else r.get(c) for c in OUT_COLS])


# ---------------------------------------------------------------- info / registry
def tools_info(tdir=None) -> dict:
    out = {}
    for t in TOOL_NAMES:
        p = tool_paths(t, tdir)
        out[t] = {"exe": str(p["exe"]), "version": tool_version(t, tdir),
                  "exe_sha256": C.sha256_file(p["exe"]), "build_id": tool_build_id(t, tdir),
                  "n_hashed_files": len(p["files"])}
    return out


def register(info: dict, registry_path=None):
    """Merge tool entries (kind: tool) into the data registry (a YAML list; entries are kept)."""
    import yaml
    from . import download as D
    if registry_path:
        D.REGISTRY = Path(registry_path)
    reg = D.load_registry()
    for name, e in info.items():
        reg[f"tool_{name}"] = {"name": f"tool_{name}", "kind": "tool", **e}
    D.save_registry(reg)
    return reg


# ---------------------------------------------------------------- smoke test
def pick(rows, n, seed=1234):
    return sorted(random.Random(seed).sample(rows, min(n, len(rows))), key=lambda r: r["case_id"])


def main(argv=None):
    from .scoring import read_manifest
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--info", action="store_true")
    ap.add_argument("--register")
    ap.add_argument("--tool", default="all")
    ap.add_argument("--cases")
    ap.add_argument("--limit-per-pool", type=int)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--protein"); ap.add_argument("--nucleotide")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=TIMEOUT_S)
    ap.add_argument("--cache"); ap.add_argument("--out")
    ap.add_argument("--keep-alignments"); ap.add_argument("--retry-errors", action="store_true")
    a = ap.parse_args(argv)
    if a.workers > 30:
        ap.error("at most 30 concurrent aligner processes")
    if a.info:
        print(json.dumps(tools_info(), indent=1))
        return 0
    if a.register:
        reg = register(json.loads(Path(a.register).read_text()))
        print(f"registry updated ({len(reg)} entries)")
        return 0
    if a.smoke:
        rows = pick(read_manifest(a.protein), a.n) + pick(read_manifest(a.nucleotide), a.n)
    elif a.cases:
        rows = read_manifest(a.cases)
        if a.limit_per_pool:
            rows = pick(rows, a.limit_per_pool)
    else:
        ap.error("need --info, --register, --smoke or --cases")
    tools = TOOL_NAMES if a.tool == "all" else a.tool.split(",")
    allrows = []
    for t in tools:
        res = evaluate_external(rows, t, a.workers, cache_path=a.cache, timeout_s=a.timeout,
                                retry_errors=a.retry_errors, keep_dir=a.keep_alignments, progress=True)
        allrows += res
        n_err = sum(1 for r in res if r.get("error"))
        print(f"{t}: {len(res)} cases, {sum(r['cached'] for r in res)} cached, {n_err} errors",
              file=sys.stderr)
    if a.out:
        write_rows(allrows, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
