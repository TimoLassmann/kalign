"""Verify that the batched parameter-set worker gives the same alignments as the per-case worker.

    python -m recipe.verify_batch_worker [--n 50] [--workers 8] [--out DIR]

For each (pool, 3.7.0 mode) spec: n random cases (seed 1234) are aligned twice with throw-away caches,
once per case (``batch=False``) and once through long-lived workers (``batch=True``); the kept
alignments are compared byte for byte, and the stored metrics and error entries must be identical.
CPU sums of both paths are reported; they agree to a few percent (both are RUSAGE deltas around the
aligner call; the difference is run-to-run noise at 8 concurrent workers, not a systematic bias claim).
"""
from __future__ import annotations

import argparse
import filecmp
import json
import random
import sys
import tempfile
from pathlib import Path

from . import cache as C
from . import evaluate as E
from .presets_370 import params_370
from .scoring import read_manifest
from .step0 import BUILD_DIR, MANIFESTS

SPECS = (("protein", "default"), ("protein", "accurate"), ("nucleotide", "default"), ("nucleotide", "accurate"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    build = E.Build.load(BUILD_DIR)
    out = Path(a.out or tempfile.mkdtemp(prefix="verify_batch_"))
    out.mkdir(parents=True, exist_ok=True)
    report, ok_all = [], True
    for pool, mode in SPECS:
        cases = read_manifest(MANIFESTS[pool])
        cases = random.Random(1234).sample(cases, a.n)
        spec = params_370(pool, mode)
        res = {}
        for name, batch in (("per_case", False), ("batch", True)):
            kd = out / f"{pool}_{mode}_{name}"
            rows = E.evaluate(cases, spec, build, a.workers, cache_path=out / f"{pool}_{mode}_{name}.sqlite",
                              keep_dir=kd, batch=batch)
            res[name] = (rows, kd / C.param_key(spec)[:16])
        n_same = n_diff = 0
        for r1, r2 in zip(res["per_case"][0], res["batch"][0]):
            fn = f"{r1['case_id'].replace('/', '_')}.aln.fa"
            p1, p2 = res["per_case"][1] / fn, res["batch"][1] / fn
            same = (r1["error"] == r2["error"] and r1["aln_sha256"] == r2["aln_sha256"]
                    and all(r1[k] == r2[k] for k in ("recall", "precision", "f1", "tc", "width_ratio"))
                    and ((not p1.exists() and not p2.exists()) or filecmp.cmp(p1, p2, shallow=False)))
            n_same += same
            n_diff += not same
        cpu = {k: sum(r["cpu_s"] or 0 for r in v[0]) for k, v in res.items()}
        n_err = sum(1 for r in res["batch"][0] if r["error"])
        line = {"pool": pool, "mode": mode, "n": a.n, "identical": n_same, "different": n_diff,
                "errors": n_err, "cpu_per_case_s": round(cpu["per_case"], 2), "cpu_batch_s": round(cpu["batch"], 2),
                "wall_s_per_case_sum": round(sum(r["wall_s"] or 0 for r in res["per_case"][0]), 2),
                "wall_s_batch_sum": round(sum(r["wall_s"] or 0 for r in res["batch"][0]), 2)}
        print(json.dumps(line), flush=True)
        report.append(line)
        ok_all &= n_diff == 0
    (out / "verify_batch_worker.json").write_text(json.dumps(report, indent=1) + "\n")
    print("ALL IDENTICAL" if ok_all else "MISMATCH", f"(details in {out})")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
