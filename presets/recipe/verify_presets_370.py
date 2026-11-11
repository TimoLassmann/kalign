"""Byte-for-byte check: 3.7.0 preset via the CLI (--mode M -n 1) vs via the parameter-set worker.

    python -m recipe.verify_presets_370 --build-dir BUILD --workers 28 --out DIR [--n 40]

Protein: n cases x fast/default/accurate (+ recall on a subset, alias of accurate);
nucleotide: n cases x fast/default/recall/accurate.  Cases are spread over benchmarks (round robin
over benchmark, seeded shuffle within).  Cases with > max_cpu_s CPU for the CLI accurate preset in
the step 0 cache are skipped (keeps the check cheap).  Output alignments are compared as bytes
(sha256).  Writes DIR/verify_presets_370.tsv and prints a summary; exit 1 on any mismatch.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import hashlib
import json
import multiprocessing as mp
import random
import sys
import tempfile
from pathlib import Path

from . import evaluate as E
from .presets_370 import MODES, params_370
from .scoring import read_manifest

HERE = Path(__file__).resolve().parent


def pick_cases(rows, n, seed=1234, skip=()):
    rows = [r for r in rows if r["case_id"] not in skip]
    rng = random.Random(seed)
    by = {}
    for r in rows:
        by.setdefault(r["benchmark"], []).append(r)
    for v in by.values():
        rng.shuffle(v)
    out = []
    while len(out) < n and any(by.values()):
        for b in sorted(by):
            if by[b] and len(out) < n:
                out.append(by[b].pop())
    return out


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def one(args):
    build_d, pool, mode, case = args
    b = E.Build(**build_d)
    with tempfile.TemporaryDirectory(prefix="v370_") as tmp:
        o1, o2 = Path(tmp, "cli.fa"), Path(tmp, "wrk.fa")
        w1, c1, e1 = E.run_cli(b.cli, mode, case["unaligned"], o1)
        w2, c2, e2 = E.run_params(b.python, params_370(pool, mode), case["unaligned"], o2, tmp)
        s1 = _sha(o1) if e1 is None and o1.exists() else None
        s2 = _sha(o2) if e2 is None and o2.exists() else None
    return {"pool": pool, "mode": mode, "case_id": case["case_id"], "benchmark": case["benchmark"],
            "cli_sha": s1, "worker_sha": s2, "identical": int(s1 is not None and s1 == s2),
            "cli_cpu": c1, "worker_cpu": c2, "cli_err": e1, "worker_err": e2}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--workers", type=int, default=28)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--recall-protein-n", type=int, default=10)
    a = ap.parse_args(argv)
    build = E.Build.load(a.build_dir)
    jobs = []
    for pool in ("protein", "nucleotide"):
        rows = read_manifest(Path(a.data_dir, f"cases_{pool}.tsv"))
        cases = pick_cases(rows, a.n, a.seed)
        for mode in MODES[pool]:
            if pool == "protein" and mode == "recall":
                use = cases[:a.recall_protein_n]
            else:
                use = cases
            jobs += [(build.__dict__, pool, mode, c) for c in use]
    Path(a.out).mkdir(parents=True, exist_ok=True)
    res = []
    with cf.ProcessPoolExecutor(a.workers, mp_context=mp.get_context("fork")) as ex:
        for i, r in enumerate(ex.map(one, jobs, chunksize=1), 1):
            res.append(r)
            if i % 20 == 0:
                print(f"{i}/{len(jobs)}", file=sys.stderr, flush=True)
    cols = list(res[0])
    with open(Path(a.out, "verify_presets_370.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in res:
            w.writerow(["" if r[c] is None else r[c] for c in cols])
    summ = {}
    for r in res:
        s = summ.setdefault((r["pool"], r["mode"]), {"n": 0, "identical": 0, "cli_cpu": 0.0, "wrk_cpu": 0.0})
        s["n"] += 1
        s["identical"] += r["identical"]
        s["cli_cpu"] += r["cli_cpu"] or 0
        s["wrk_cpu"] += r["worker_cpu"] or 0
    bad = 0
    for (pool, mode), s in sorted(summ.items()):
        print(f"{pool:10s} {mode:9s} identical {s['identical']}/{s['n']}  cpu cli {s['cli_cpu']:.1f}s worker {s['wrk_cpu']:.1f}s")
        bad += s["n"] - s["identical"]
    for r in res:
        if not r["identical"]:
            print("MISMATCH", r["pool"], r["mode"], r["case_id"], r["cli_err"], r["worker_err"])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
