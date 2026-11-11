"""Regenerate the golden hashes of tests/preset_golden.cmake (Kalign 4 modes).

The golden outputs are produced through the Python path the re-fit evaluated
(kalign._core.ensemble_custom_file_to_file with the parameters of
presets/results/v2/presets_v2.json, and for protein fast those of
presets/results/fastfit/final_pick.json with the member gap model, and for
protein accurate those of presets/results/ensconf/final_pick.json, for nucleotide
fast and accurate those of presets/results/nucfit/{fast,accurate}_final_pick.json, seed 42,
n_threads 1), NOT through the C
presets, so the ctest checks that lib/src/mode_presets.h reproduces the
evaluated presets byte for byte. "default" is an alias of "fast".

Usage (from the repository root, after `uv sync`):

    uv run python tests/preset_golden_regen.py

and paste the printed lines into the GOLDEN list of tests/preset_golden.cmake.
"""

import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

from kalign import _core

ROOT = Path(__file__).resolve().parents[1]
P = json.loads((ROOT / "presets/results/v2/presets_v2.json").read_text())["presets"]
P["protein"]["fast"] = json.loads((ROOT / "presets/results/fastfit/final_pick.json").read_text())
P["protein"]["fast"]["gap_model"] = 1  # the fast fit ran with the member gap model
# the ensemble fit's pick carries its gap model per member (run_gap_model)
P["protein"]["accurate"] = json.loads((ROOT / "presets/results/ensconf/final_pick.json").read_text())
# nucleotide fast and accurate: the nucfit final picks, as presets/release/emit_c_presets.py takes them
for _m in ("fast", "accurate"):
    _path = ROOT / f"presets/results/nucfit/{_m}_final_pick.json"
    if _path.exists():
        P["nucleotide"][_m] = json.loads(_path.read_text())
TREE_KEYS = ("run_tree_method", "run_tree_dist", "run_tree_alpha", "run_tree_maxlen", "run_tree_leaf",
             "run_tree_kmeans_dist", "run_realign_tree")
# per-member settings; when present they override the scalar (seq_weights, realign, consistency, ...)
MEMBER_KEYS = ("run_seq_weights", "run_realign", "run_consistency_anchors", "run_consistency_weight", "run_vsm_mode",
               "run_gap_model", "run_tree_sl_noise")
# (mode name given to the CLI, pool, preset mode, input)
CASES = [(m, "protein", "fast" if m == "default" else m, f)
         for m in ("fast", "accurate", "default") for f in ("BB11001.tfa", "BB12006.tfa", "BB30014.tfa")]
CASES += [(m, "nucleotide", m, "small.fa") for m in ("fast", "accurate")]


def kwargs(preset):
    p = preset["params"]
    kw = dict(
        run_gpo=p["run_gpo"], run_gpe=p["run_gpe"], run_tgpe=p["run_tgpe"], run_noise=p["run_noise"],
        run_types=p["run_types"], format="fasta", seq_type=p["run_types"][0], seed=42,
        min_support=p["min_support"], realign=p["realign"], seq_weights=p["seq_weights"],
        n_threads=1, consistency_anchors=p["consistency"], consistency_weight=p["consistency_weight"],
        run_vsm_amax=p["run_vsm_amax"], run_refine=p["run_refine"],
        consistency_merge=p["consistency_merge"], consistency_merge_weight=p["consistency_merge_weight"],
        vsm_mode=p["vsm_mode"], run_vsm_kappa=p["run_vsm_kappa"], run_vsm_dzero=p["run_vsm_dzero"],
        fill_in=p["fill_in"], fill_anchor_frac=p["fill_anchor_frac"], fill_support=p["fill_support"],
        run_guide_dist=p["run_guide_dist"],
    )
    for k in TREE_KEYS + MEMBER_KEYS:
        if k in p:
            kw[k] = p[k]
    if preset.get("gap_model"):
        kw["run_gap_model"] = [preset["gap_model"]] * p["n_runs"]
    return kw


def main():
    data = ROOT / "tests" / "data"
    with tempfile.TemporaryDirectory() as tmp:
        for mode, pool, preset, name in CASES:
            out = os.path.join(tmp, f"{mode}_{name}.fasta")
            _core.ensemble_custom_file_to_file(str(data / name), out, **kwargs(P[pool][preset]))
            h = hashlib.sha256(open(out, "rb").read()).hexdigest()
            sys.stdout.write(f"  \"{mode}:{name}:{h}\"\n")


if __name__ == "__main__":
    main()
