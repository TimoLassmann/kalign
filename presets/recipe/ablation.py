"""Component ablation of the accurate presets on held-out cases (Kalign 4 paper, block F).

For every fold k, the fold's own accurate preset (selection_fold{k}.json of the v2 picks, chosen on that fold's
training core only) is run with one component switched off, on that fold's test cases; pooled over the five
folds this gives out-of-fold scores for each variant. Nothing is selected or tuned here: the variants are fixed
in VARIANTS before any score is seen, and the full preset is one of them (cache hits from the v2 test run).

    python -m recipe.ablation --pool protein --ws WS --sel-dir ../presets/results/v2/protein --out OUT [--threads 12]

Writes OUT/ablation_<pool>_fold<k>.tsv (one row per case and variant); a finished fold is skipped on rerun, an
unfinished one resumes from the cache.
"""

import argparse
import copy
import json
from pathlib import Path

from . import fit as FIT
from .nested import ROW_COLS, _make_evaluator, write_rows


def _single_run(p):
    q = copy.deepcopy(p)
    for k, v in p.items():
        if k.startswith("run_") and isinstance(v, list):
            q[k] = v[:1]
    q["n_runs"] = 1
    q["min_support"] = 0
    q["fill_in"] = 0
    q["fill_anchor_frac"] = 0.0
    return q


def _set(**kw):
    def f(p):
        q = copy.deepcopy(p)
        for k, v in kw.items():
            q[k] = [v] * len(q[k]) if isinstance(q[k], list) else v
        return q
    return f


VARIANTS = {
    "protein": {
        "full": lambda p: copy.deepcopy(p),
        "single_run": _single_run,                 # ensemble off: run 1 of the preset alone
        "no_consistency": _set(consistency=0),
        "no_realign": _set(realign=0),
        "no_seq_weights": _set(seq_weights=0.0),
        "no_vsm": _set(run_vsm_amax=0.0),
        "no_refine": _set(run_refine=0),
        "no_fill_in": _set(fill_in=0, fill_anchor_frac=0.0),
    },
    "nucleotide": {
        "full": lambda p: copy.deepcopy(p),
        "no_consistency": _set(consistency=0),
        "no_realign": _set(realign=0),
        "no_seq_weights": _set(seq_weights=0.0),
        "no_vsm": _set(run_vsm_amax=0.0),
        "no_refine": _set(run_refine=0),
    },
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--ws", required=True)
    ap.add_argument("--sel-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--build-dir")
    ap.add_argument("--slots-dir")
    ap.add_argument("--cache")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    for k in range(5):
        dst = a.out / f"ablation_{a.pool}_fold{k}.tsv"
        if dst.exists():
            continue
        sel = json.loads((a.sel_dir / f"selection_fold{k}.json").read_text())
        base = sel["selection"]["accurate"]["params"]
        arms = {name: f(base) for name, f in VARIANTS[a.pool].items()}
        keys = {name: FIT.pkey(p) for name, p in arms.items()}
        rows, meta, core, _, test = FIT.load_fold_inputs(a.pool, k, a.ws)
        if set(test) & set(core):
            raise FIT.LeakError(f"{a.pool} fold {k}: test cases in the training core")
        test = sorted(test)
        ev = _make_evaluator(a.ws, rows, frozenset(test), a)
        res = ev.evaluate({keys[n]: arms[n] for n in arms}, test)
        out = []
        for name in arms:
            for cid in test:
                r = res[keys[name]][cid]
                m = meta[cid]
                out.append({"pool": a.pool, "fold": k, "scope": "oof", "mode": "accurate", "arm": name,
                            "key": keys[name], "is_370": 0, "fallback_370": 0, "smoke": 0, "case_id": cid,
                            "benchmark": m["benchmark"], "difficulty": m["difficulty"],
                            "homolog_in_train": m["homolog_in_train"],
                            **{c: r.get(c) for c in ROW_COLS[13:]}})
        write_rows(dst, out)
        print(f"{a.pool} fold {k}: {len(out)} rows -> {dst}", flush=True)


if __name__ == "__main__":
    main()
