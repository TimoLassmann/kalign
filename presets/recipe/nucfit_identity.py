"""Identity checks of the nucfit plan (PLAN_nucfit.md section 10), run before the freeze.

    python -m recipe.nucfit_identity --ws WS [--workers 60]

N1  build 9137658: the parameter-set worker with the pre-release nucleotide fast dict (4ff7e4d1) gives the same
    alignment (sha256) as the CLI preset ``--mode fast``; the same for accurate (bd03cdb7) and ``--mode accurate``.
N2  build 9137658: the fast dict with the 4.0.0 guide tree and gap model 0 written out (the fit's reference) gives the
    same alignments as the plain dict.
N3  every reference alignment of the manifest (fit and lockbox; no alignment program runs) scores SP = TC = 1 (projected cases: the full true alignment against the subset)
    against itself.
Check cases (N1, N2): the first 60 cases (case-id order) of the small screening subset of fold 0, the first 3 Rfam
large training cases of fold 0, the first 2 simulations at 1,000 and the first 1 at 10,000 (never a lockbox case).
Writes WS/runs/identity/identity.json; exit 1 on any difference or failed case.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import fastfit as FF
from . import nucfast as NF
from .evaluate import Build, evaluate
from .scoring import score_case


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ws", default=str(NF.default_ws()))
    ap.add_argument("--workers", type=int, default=30)
    a = ap.parse_args(argv)
    ws = Path(a.ws).expanduser()
    rows, meta, core, screen, test, lockbox = NF.load_fold(ws, 0)
    small = sorted(c for c in screen if meta[c]["kind"] == "small")[:60]
    large = []
    for k, n in (("rfam", 3), ("1000", 2), ("10000", 1)):
        large += sorted(c for c in core if meta[c]["kind"] == "large" and meta[c]["k"] == k)[:n]
    ids = small + large
    assert not set(ids) & set(lockbox)
    cases = [rows[c] for c in ids]
    cache = ws / "cache" / "cache.sqlite"
    b = Build.load(ws / "kalign" / NF.BUILD)

    def shas(spec):
        out = evaluate(cases, spec, b, workers=a.workers, cache_path=cache, timeout_s=FF.TIMEOUT_S)
        return {r["case_id"]: (r.get("aln_sha256"), r.get("error")) for r in out}

    fast, acc = NF.v2_params(NF.NUC_FAST_KEY), NF.v2_params(NF.NUC_ACC_KEY)
    res, ok = {}, True
    pairs = {"N1_fast_cli_vs_worker": ("fast", fast),
             "N1_accurate_cli_vs_worker": ("accurate", acc),
             "N2_fast_explicit_tree_gap": (fast, NF.with_gap(FF.with_tree(fast, FF.TREE_400), 0))}
    for name, (s1, s2) in pairs.items():
        x, y = shas(s1), shas(s2)
        err = [c for c in ids if x[c][1] or y[c][1]]
        diff = [c for c in ids if x[c][0] != y[c][0]]
        res[name] = {"n": len(ids), "identical": len(ids) - len(diff), "differ": diff, "errors": err}
        ok &= not diff and not err
        print(f"{name}: identical {len(ids) - len(diff)}/{len(ids)}, errors {len(err)}", flush=True)
    bad = []
    for cid, r in sorted(rows.items()):
        # projected cases (simulations at 10,000): the reference is a 1,000-row subset, so the known answer is the
        # full true alignment scored against it
        test = (str(Path(r["reference"]).with_name("true.fasta")) if r.get("score_rule") == "projected"
                else r["reference"])
        s = score_case(r, test)
        if abs(float(s["recall"]) - 1.0) > 1e-12 or abs(float(s["tc"]) - 1.0) > 1e-12:
            bad.append((cid, s["recall"], s["tc"]))
    res["N3_reference_self_score"] = {"n": len(rows), "not_one": bad[:50], "n_not_one": len(bad)}
    ok &= not bad
    print(f"N3_reference_self_score: {len(rows) - len(bad)}/{len(rows)} score 1", flush=True)
    out = ws / "runs" / "identity"
    out.mkdir(parents=True, exist_ok=True)
    (out / "identity.json").write_text(json.dumps({"cases": ids, "results": res, "pass": ok}, indent=1) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
