"""Identity checks of the fast-fit plan (PLAN_fastfit.md section 11), run before the freeze.

    python -m recipe.fastfit_identity --ws WS [--workers 60]

I1  4.0.0 build (daff955): the parameter-set worker with fastfit.P400_FAST gives the same alignment (sha256) as the
    CLI preset ``--mode fast`` on every check case.
I2  new build (99aecb2): P400_FAST and GM1_FAST with the 4.0.0 tree values written out (fastfit.with_tree(.,
    TREE_400)) give the same alignments as the same dicts without tree keys.
Check cases: the first 60 cases (case-id order) of the small screening subset of fold 0, and the first 3 / 2
training families of fold 0 at 1,000 / 10,000 (never a lockbox case). Results go through the fit cache.
Writes WS/runs/identity/identity.json; exit 1 on any difference or failed case.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import fastfit as FF
from .evaluate import Build, evaluate


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ws", default=str(FF.default_ws()))
    ap.add_argument("--workers", type=int, default=30)
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    rows, meta, core, screen, test, lockbox = FF.load_fold(ws, 0)
    small = sorted(c for c in screen if meta[c]["kind"] == "small")[:60]
    l1 = sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 1000)[:3]
    l10 = sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == 10000)[:2]
    ids = small + l1 + l10
    assert not set(ids) & set(lockbox)
    cases = [rows[c] for c in ids]
    cache = ws / "cache" / "cache.sqlite"
    b_ref, b_new = Build.load(ws / "kalign" / FF.REF_BUILD), Build.load(ws / "kalign" / FF.NEW_BUILD)

    def shas(spec, build):
        out = evaluate(cases, spec, build, workers=a.workers, cache_path=cache, timeout_s=FF.TIMEOUT_S)
        return {r["case_id"]: (r.get("aln_sha256"), r.get("error")) for r in out}

    res, ok = {}, True
    pairs = {"I1_400_fast_cli_vs_worker": ((("fast"), b_ref), (FF.P400_FAST, b_ref)),
             "I2_400_fast_explicit_tree": ((FF.P400_FAST, b_new), (FF.with_tree(FF.P400_FAST, FF.TREE_400), b_new)),
             "I2_gm1_fast_explicit_tree": ((FF.GM1_FAST, b_new), (FF.with_tree(FF.GM1_FAST, FF.TREE_400), b_new))}
    for name, ((s1, bb1), (s2, bb2)) in pairs.items():
        x, y = shas(s1, bb1), shas(s2, bb2)
        err = [c for c in ids if x[c][1] or y[c][1]]
        diff = [c for c in ids if x[c][0] != y[c][0]]
        res[name] = {"n": len(ids), "identical": len(ids) - len(diff), "differ": diff, "errors": err}
        ok &= not diff and not err
        print(f"{name}: identical {len(ids) - len(diff)}/{len(ids)}, errors {len(err)}", flush=True)
    out = ws / "runs" / "identity"
    out.mkdir(parents=True, exist_ok=True)
    (out / "identity.json").write_text(json.dumps({"cases": ids, "results": res, "pass": ok}, indent=1) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
