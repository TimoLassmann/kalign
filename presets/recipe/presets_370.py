"""Kalign 3.7.0 presets as decoded parameter dicts for the parameter-set worker.

Source of truth: lib/src/aln_wrap.c (kalign_get_mode_preset, preset_protein, preset_nucleotide,
preset_protein_ensemble, preset_run, preset_consistency, kalign_run_config_defaults,
kalign_ensemble_config_defaults).  NOT experiments/redesign space.py SHIPPED (that table is the 3.5
protein presets and is outdated).  The numbers below were transcribed from aln_wrap.c; the test
tests/test_presets_370.py parses aln_wrap.c again and compares every number, so a drift in the C
source is caught.

Layout = decoded params of experiments/redesign/tools/refit/space.py decode() (see _align_worker.py).

    params_370(pool, mode)      -> deep copy of the decoded dict, pool "protein" | "nucleotide"
    MODES[pool]                 -> modes that exist for the pool (protein "recall" is an alias of "accurate")
    check_in_space(pool, p)     -> list of problems (empty = valid point of the search space)
    protein_stage1_members(p)   -> the single-run (stage-1) individuals of a protein ensemble preset

Why the worker reproduces the CLI (verified by tests/verify_presets_370.py, results in
docs/build/presets_370.md):
  * the worker calls _core.ensemble_custom_file_to_file with exactly the per-run values of the
    preset and the shared (seq_weights, realign, consistency, vsm_mode, ...) values; the C side
    builds kalign_run_config_defaults() + overrides, as preset_run() does;
  * the seed only enters through tree_seed, which is used only when tree_noise > 0 (aln_wrap.c
    single_run: `cfg->tree_seed != 0 && cfg->tree_noise > 0`).  Every 3.7.0 preset has noise 0 in
    every run, so the CLI seeds (0, or 42+k) and the worker seed 42 (tree_seed 42+k) give the same
    alignments.  The seed would matter only for individuals with noise > 0, which the stage-1/-2
    search space masks to 0 whenever realign > 0 or n_runs = 1 (space.mask) and the stage-2 space
    fixes to 0 (STAGE2_FIXED).
"""

from __future__ import annotations

import copy
import os
import sys
from pathlib import Path

# aligner constants (kalign.h); the worker re-checks the ones it is given via job["expect"]
PFASUM43, PFASUM60, CORBLOSUM66 = 1, 2, 3
K1, K20, K200 = 10, 11, 12
REF_NONE, REF_ALL, REF_CONF, REF_INLINE = 0, 1, 2, 3
GUIDE_BPM, GUIDE_PAIRALN = 0, 1
FILL_OFF = 0

MODES = {"protein": ("fast", "default", "accurate", "recall"),
         "nucleotide": ("fast", "default", "recall", "accurate")}
ALIAS = {("protein", "recall"): "accurate"}   # preset_protein(): recall == accurate


def _params(runs, *, seq_weights, realign, consistency=0, consistency_weight=1.0, min_support=0,
            vsm_mode=0, kappa=None, dzero=None, guide=None):
    """runs: list of (matrix, gpo, gpe, tgpe, vsm_amax, refine) as the arguments of preset_run().

    Unset run fields take kalign_run_config_defaults(): vsm_kappa 1, vsm_dzero 1, guide_dist BPM,
    tree_noise 0.  Ensemble defaults (kalign_ensemble_config_defaults + preset overrides):
    consistency_merge 0, consistency_merge_weight 2, fill_in OFF, fill_anchor_frac 0, fill_support 1.
    """
    n = len(runs)
    return {
        "n_runs": n,
        "run_gpo": [r[1] for r in runs], "run_gpe": [r[2] for r in runs], "run_tgpe": [r[3] for r in runs],
        "run_noise": [0.0] * n, "run_types": [r[0] for r in runs],
        "run_vsm_amax": [r[4] for r in runs], "run_refine": [r[5] for r in runs],
        "run_vsm_kappa": list(kappa) if kappa else [1.0] * n,
        "run_vsm_dzero": list(dzero) if dzero else [1.0] * n,
        "run_guide_dist": list(guide) if guide else [GUIDE_BPM] * n,
        "seq_weights": seq_weights, "consistency": consistency, "consistency_weight": consistency_weight,
        "realign": realign, "min_support": min_support, "consistency_merge": 0,
        "consistency_merge_weight": 2.0, "vsm_mode": vsm_mode, "fill_in": FILL_OFF,
        "fill_anchor_frac": 0.0, "fill_support": 1,
    }


def _build():
    # ------------------------------------------------------------------ protein (preset_protein)
    protein = {}
    # fast: preset_run(PFASUM60, 8.626, 0.843, 0.433, vsm_amax 0.592, seq_weights 1.534, realign 0, NONE)
    protein["fast"] = _params([(PFASUM60, 8.626, 0.843, 0.433, 0.592, REF_NONE)],
                              seq_weights=1.534, realign=0)
    # default: single run, realign 2, refine ALL, consistency 1 round with weight 1.5850479405163567
    protein["default"] = _params(
        [(PFASUM60, 5.621009959254649, 0.8981372818573407, 0.4737930605697324, 0.4863285321570886,
          REF_ALL)],
        seq_weights=1.2332522985904468, realign=2, consistency=1,
        consistency_weight=1.5850479405163567, vsm_mode=0)
    # accurate (= recall): preset_protein_ensemble(), 5 runs, vsm_mode 1, shared seq_weights,
    # consistency_anchors 0 (weight 1.0), realign 0, min_support 2, fill_in OFF.
    # vsm_amax is 0.0 in every run (unused with vsm_mode 1).
    matrix = [PFASUM60, PFASUM60, CORBLOSUM66, PFASUM43, PFASUM43]
    gpo = [6.123414948157974, 3.9422242450518326, 3.951413346944608, 5.951508465248285, 6.129249240073377]
    gpe = [0.8905358238624422, 0.8451893377347646, 0.7953340627319447, 0.9097516251997103,
           0.8977023562247453]
    tgpe = [0.49010011751845955, 0.4332371708704297, 0.35817139949635746, 0.4911277927807959,
            0.4898265278919298]
    refine = [REF_INLINE, REF_NONE, REF_NONE, REF_ALL, REF_ALL]
    kappa = [0.9266191602810352, 0.9988193293912846, 0.9954160682480878, 0.42732812652138796,
             0.8908869809114087]
    dzero = [7.675591889662474, 0.9937625844005973, 1.0026464536496102, 7.671180058242367,
             7.718281725286436]
    guide = [GUIDE_BPM, GUIDE_PAIRALN, GUIDE_PAIRALN, GUIDE_PAIRALN, GUIDE_BPM]
    protein["accurate"] = _params(
        [(matrix[k], gpo[k], gpe[k], tgpe[k], 0.0, refine[k]) for k in range(5)],
        seq_weights=0.930328121699592, realign=0, consistency=0, consistency_weight=1.0,
        min_support=2, vsm_mode=1, kappa=kappa, dzero=dzero, guide=guide)

    # ------------------------------------------------------------------ nucleotide (preset_nucleotide)
    nuc = {}
    nuc["fast"] = _params([(K200, 19.873, 0.642, 2.582, 0.874, REF_NONE)], seq_weights=0.199, realign=1)
    nuc["default"] = _params(
        [(K20, 19.275, 0.628, 2.854, 0.445, REF_NONE),
         (K200, 16.581, 0.274, 2.491, 0.706, REF_NONE),
         (K200, 16.499, 3.754, 2.896, 0.491, REF_CONF)],
        seq_weights=0.448, realign=1, min_support=2)
    nuc["recall"] = _params([(K200, 19.882, 0.646, 2.854, 0.874, REF_CONF)], seq_weights=0.449, realign=1,
                            consistency=1, consistency_weight=1.616)
    nuc["accurate"] = _params(
        [(K200, 19.475, 0.728, 2.856, 0.492, REF_CONF),
         (K1, 19.377, 1.379, 2.486, 1.211, REF_CONF),
         (K200, 15.869, 0.353, 2.430, 0.721, REF_NONE),
         (K200, 15.989, 2.331, 2.764, 0.188, REF_CONF),
         (K1, 5.346, 3.483, 0.610, 1.406, REF_NONE)],
        seq_weights=0.590, realign=1, consistency=2, consistency_weight=0.800, min_support=3)
    return {"protein": protein, "nucleotide": nuc}


_PRESETS = _build()


def params_370(pool: str, mode: str) -> dict:
    """Decoded parameter dict of the 3.7.0 preset ``mode`` for ``pool`` (a fresh deep copy)."""
    if pool not in _PRESETS:
        raise ValueError(f"unknown pool {pool!r}")
    mode = ALIAS.get((pool, mode), mode)
    if mode not in _PRESETS[pool]:
        raise ValueError(f"no 3.7.0 preset {mode!r} for pool {pool!r}")
    return copy.deepcopy(_PRESETS[pool][mode])


# ---------------------------------------------------------------------- search-space access
def _tool_dirs():
    here = Path(__file__).resolve()
    dirs = []
    if os.environ.get("KALIGN_REFIT_TOOLS"):
        dirs.append(Path(os.environ["KALIGN_REFIT_TOOLS"]))
    dirs.append(here.parents[2] / "experiments" / "redesign" / "tools")      # repo checkout
    dirs.append(Path.home() / "tmp" / "kalign-refit" / "redesign_tools")      # workstation copy
    return dirs


# sha256 of the search-space modules the fit was frozen with (experiments/redesign at 960a4f6, unchanged).
# The freeze manifest covers this file, so pinning the hashes here puts the spaces under Freeze 1 as well;
# recipe.fit refuses to start with any other space file.
# Nucleotide: the frozen file had sha256 869c22f5d10f24ef7fea6be4e9b5cb110fc0a2215c5f28bf043f795cd05cc1a3; for
# publication (2026-11-11) a local path in its docstring was replaced, nothing else (see presets/README.md).
SPACE_SHA256 = {"protein": "c1ce7be04ec90094621c4b83f196b16a2420dc4133d3d6bb5d0bff2381326e13",
                "nucleotide": "a3c7d53ba534c8dd7e9b7c4fe411b2da02c784f1b64e9c0af89ab9e1f87776a4"}


def load_space(pool: str):
    """Import the (read-only) protein ``refit/space.py`` or nucleotide ``nested_nuc/nuc_space.py``."""
    import importlib.util
    rel = "refit/space.py" if pool == "protein" else "nested_nuc/nuc_space.py"
    for d in _tool_dirs():
        f = d / rel
        if f.exists():
            spec = importlib.util.spec_from_file_location(f"_space_{pool}", f)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = mod
            spec.loader.exec_module(mod)
            return mod
    raise FileNotFoundError(f"search space {rel} not found in {[str(d) for d in _tool_dirs()]}")


def _bounds(space, pool):
    """(per-run bounds, shared bounds, choice sets) read from the space modules."""
    if pool == "protein":
        per = {n: b for n, kind, b in space.PER_RUN if kind != "choice"}
        sh = {n: b for n, kind, b in space.SHARED if kind in ("real", "int")}
    else:
        per = {"gpo": (space.PER_RUN_LO[0], space.PER_RUN_HI[0]), "gpe": (space.PER_RUN_LO[1], space.PER_RUN_HI[1]),
               "tgpe": (space.PER_RUN_LO[2], space.PER_RUN_HI[2]), "noise": (space.PER_RUN_LO[3], space.PER_RUN_HI[3]),
               "vsm_amax": (space.SHARED_LO[0], space.SHARED_HI[0])}
        sh = {"seq_weights": (space.SHARED_LO[1], space.SHARED_HI[1]),
              "consistency_weight": (space.SHARED_LO[2], space.SHARED_HI[2]),
              "realign": (0, 2), "min_support": (0, 5), "consistency_merge_weight": (0.5, 10.0)}
    return per, sh


def check_in_space(pool: str, p: dict) -> list[str]:
    """Problems that keep decoded params ``p`` from being a point of the full search space.

    Checks (1) every variable inside its bounds / choice set, (2) space.encode -> space.decode
    returns ``p`` unchanged (so the individual has the same cache key after a round trip).
    """
    sp = load_space(pool)
    bad = []
    per, sh = _bounds(sp, pool)
    key = {"gpo": "run_gpo", "gpe": "run_gpe", "tgpe": "run_tgpe", "noise": "run_noise",
           "vsm_amax": "run_vsm_amax", "vsm_kappa": "run_vsm_kappa", "vsm_dzero": "run_vsm_dzero"}
    for name, (lo, hi) in per.items():
        vals = p.get(key[name], [])
        if pool == "protein" and name == "vsm_amax" and p["vsm_mode"] == 1:
            continue            # masked to the canonical 0.0 in vsm_mode 1
        if pool == "protein" and name in ("vsm_kappa", "vsm_dzero") and p["vsm_mode"] == 0:
            continue
        for k, v in enumerate(vals):
            if not lo <= v <= hi:
                bad.append(f"{name}_{k}={v} outside [{lo}, {hi}]")
    for name, (lo, hi) in sh.items():
        if name in p and not lo <= p[name] <= hi:
            bad.append(f"{name}={p[name]} outside [{lo}, {hi}]")
    if p["n_runs"] not in sp.N_RUNS_MAP:
        bad.append(f"n_runs={p['n_runs']} not in {sp.N_RUNS_MAP}")
    if p["consistency"] not in sp.CONSISTENCY_MAP:
        bad.append(f"consistency={p['consistency']} not in {sp.CONSISTENCY_MAP}")
    for t in p["run_types"]:
        if t not in sp.MATRIX_MAP:
            bad.append(f"matrix {t} not in {sp.MATRIX_MAP}")
    for r in p["run_refine"]:
        if r not in sp.REFINE_MAP:
            bad.append(f"refine {r} not in {sp.REFINE_MAP}")
    if pool == "protein":
        for g in p["run_guide_dist"]:
            if g not in (0, 1):
                bad.append(f"guide_dist {g}")
    # round trip
    full = copy.deepcopy(p)
    if pool == "nucleotide":
        for k in ("vsm_mode", "run_vsm_kappa", "run_vsm_dzero", "run_guide_dist", "fill_in",
                  "fill_anchor_frac", "fill_support"):
            full.pop(k, None)        # not part of the nucleotide space (3.6+ mechanisms off)
    back = sp.decode(sp.encode(full))
    if back != full:
        diff = {k: (full[k], back.get(k)) for k in full if back.get(k) != full[k]}
        bad.append(f"encode/decode round trip changes {diff}")
    if pool == "nucleotide":
        w = sp.to_worker(back)
        for k in ("vsm_mode", "run_vsm_kappa", "run_vsm_dzero", "run_guide_dist", "fill_in",
                  "fill_anchor_frac", "fill_support"):
            if w[k] != p[k]:
                bad.append(f"nucleotide worker default {k}={w[k]} differs from preset {p[k]}")
    return bad


def protein_stage1_members(p: dict) -> list[dict]:
    """Single-run (stage-1) individuals of a protein ensemble preset: run k with the ensemble-wide
    single-run parameters (seq_weights, consistency, realign, vsm_mode) of the lead.  A stage-2
    individual built from these members (stage2_decode) reproduces the ensemble when the
    ensemble parameters (min_support, fill_in, ...) are set as in the preset."""
    sp = load_space("protein")
    out = []
    for k in range(p["n_runs"]):
        q = {"n_runs": 1,
             "run_gpo": [p["run_gpo"][k]], "run_gpe": [p["run_gpe"][k]], "run_tgpe": [p["run_tgpe"][k]],
             "run_noise": [0.0], "run_types": [p["run_types"][k]], "run_vsm_amax": [p["run_vsm_amax"][k]],
             "run_refine": [p["run_refine"][k]], "run_vsm_kappa": [p["run_vsm_kappa"][k]],
             "run_vsm_dzero": [p["run_vsm_dzero"][k]], "run_guide_dist": [p["run_guide_dist"][k]],
             "seq_weights": p["seq_weights"], "consistency": p["consistency"],
             "consistency_weight": p["consistency_weight"], "realign": p["realign"],
             "min_support": 0, "consistency_merge": 0, "consistency_merge_weight": 2.0,
             "vsm_mode": p["vsm_mode"], "fill_in": FILL_OFF, "fill_anchor_frac": 0.0, "fill_support": 1}
        out.append(sp.mask(q))
    return out


def all_presets():
    """[(pool, mode, params)] for every distinct preset (recall alias of protein accurate skipped)."""
    out = []
    for pool, modes in MODES.items():
        for m in modes:
            if (pool, m) in ALIAS:
                continue
            out.append((pool, m, params_370(pool, m)))
    return out
