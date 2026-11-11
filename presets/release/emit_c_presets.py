"""Write lib/src/mode_presets.h from presets/results/v2/presets_v2.json (the Kalign 4 modes), for protein fast
presets/results/fastfit/final_pick.json (the fast-fit knee setting: member gap model and guide-tree options) and, for
protein accurate, presets/results/ensconf/final_pick.json (the confidence re-fit's pick: the ensemble fit's 5 members,
member 1 = fast, every member on member 1's guide tree with single-linkage tree noise per member) and, for nucleotide
fast and accurate, presets/results/nucfit/{fast,accurate}_final_pick.json (the final-fit picks of the nucfit plan,
presets/docs/nucfit/PLAN_nucfit.md: the nucfast knee pick and the nucconf pick) when they exist.

    uv run --no-project python presets/release/emit_c_presets.py [--check]

Every number is the float32 rounding of the JSON value (struct.pack('f'), round-to-nearest-even), printed with
9 significant digits, so the C literal is exactly the float the Python path (_core.ensemble_custom_file_to_file,
std::vector<float>) receives. --check exits 1 if the header on disk differs from what would be written.
"""
import hashlib
import json
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "presets/results/v2/presets_v2.json"
FASTFIT = ROOT / "presets/results/fastfit/final_pick.json"
ENSCONF = ROOT / "presets/results/ensconf/final_pick.json"
NUCFIT = {m: ROOT / f"presets/results/nucfit/{m}_final_pick.json" for m in ("fast", "accurate")}
# The fast fit ran on the build where the member gap model was the default (1196865, exp/fast-tree-fit), so its
# setting is reproduced with gap_model = 1.  The v2 settings ran with Kalign's gap rule (gap_model = 0).
FASTFIT_GAP_MODEL = 1
TREE_KEYS = ("run_tree_method", "run_tree_dist", "run_tree_alpha", "run_tree_maxlen", "run_tree_leaf",
             "run_tree_kmeans_dist", "run_realign_tree")
# Settings that act inside each member run (kalign_single_run): a per-run list under the run_* key (the
# _core.ensemble_custom_file_to_file argument name, which overrides the scalar as it does there), else the scalar
# broadcast to every run; the third value pads unused slots.
MEMBER_KEYS = (("run_seq_weights", "seq_weights", 0.0), ("run_realign", "realign", 0),
               ("run_consistency_anchors", "consistency", 0), ("run_consistency_weight", "consistency_weight", 2.0),
               ("run_vsm_mode", "vsm_mode", 0))
OUT = ROOT / "lib/src/mode_presets.h"
MAX = 5


def f32(v):
    x = struct.unpack("f", struct.pack("f", float(v)))[0]
    s = f"{x:.9g}"
    if not any(c in s for c in ".en"):
        s += ".0"
    return s + "f"


def row(vals, fmt, pad):
    vals = list(vals) + [pad] * (MAX - len(vals))
    return "{ " + ", ".join(fmt(v) for v in vals) + " }"


def entry(name, key, p, gap_model):
    n = p["n_runs"]
    for k in ("run_gpo", "run_gpe", "run_tgpe", "run_types", "run_vsm_amax", "run_refine", "run_vsm_kappa",
              "run_vsm_dzero", "run_guide_dist", "run_noise"):
        assert len(p[k]) == n, (name, k)
    # gap model: one value per preset (the ensemble fit's setting carries it per member, all equal)
    if "run_gap_model" in p:
        assert len(set(p["run_gap_model"])) == 1 and p["run_gap_model"][0] == gap_model, name
    # guide-tree options: -1 (the shipped tree) when the setting has none, as the fit's Python call left them unset
    tree = {k: p.get(k, [-1] * n) for k in TREE_KEYS}
    for k in TREE_KEYS:
        assert len(tree[k]) == n, (name, k)
    mem = {rk: p[rk] if rk in p else [p[sk]] * n for rk, sk, _ in MEMBER_KEYS}
    pad = {rk: d for rk, _, d in MEMBER_KEYS}
    for rk, sk, _ in MEMBER_KEYS:
        assert len(mem[rk]) == n, (name, rk)
    i = str
    return f"""        /* {name}: parameter key {key} */
        {{ {n},
          {row(p['run_types'], i, 0)},
          {row(p['run_gpo'], f32, 0)},
          {row(p['run_gpe'], f32, 0)},
          {row(p['run_tgpe'], f32, 0)},
          {row(p['run_vsm_amax'], f32, 0)},
          {row(p['run_noise'], f32, 0)},
          {row(p['run_refine'], i, 0)},
          {row(p['run_vsm_kappa'], f32, 1)},
          {row(p['run_vsm_dzero'], f32, 1)},
          {row(p['run_guide_dist'], i, 0)},
          {gap_model},
          {row(tree['run_tree_method'], i, -1)},
          {row(tree['run_tree_dist'], i, -1)},
          {row(tree['run_tree_alpha'], i, -1)},
          {row(tree['run_tree_maxlen'], i, -1)},
          {row(tree['run_tree_leaf'], i, -1)},
          {row(tree['run_tree_kmeans_dist'], i, -1)},
          {row(tree['run_realign_tree'], i, -1)},
          {row(p.get('run_tree_sl_noise', [0.0] * n), f32, 0)},
          {row(mem['run_seq_weights'], f32, pad['run_seq_weights'])},
          {row(mem['run_realign'], i, pad['run_realign'])},
          {row(mem['run_consistency_anchors'], i, pad['run_consistency_anchors'])},
          {row(mem['run_consistency_weight'], f32, pad['run_consistency_weight'])},
          {row(mem['run_vsm_mode'], i, pad['run_vsm_mode'])},
          {p['min_support']}, {p['consistency_merge']}, {f32(p['consistency_merge_weight'])},
          {p['fill_in']}, {f32(p['fill_anchor_frac'])}, {p['fill_support']} }}"""


def render():
    raw = SRC.read_bytes()
    P = json.loads(raw)["presets"]
    sha = hashlib.sha256(raw).hexdigest()
    ff_raw = FASTFIT.read_bytes()
    ff = json.loads(ff_raw)
    ff_sha = hashlib.sha256(ff_raw).hexdigest()
    ef_raw = ENSCONF.read_bytes()
    ef = json.loads(ef_raw)
    ef_sha = hashlib.sha256(ef_raw).hexdigest()
    sel = {(pool, m): (P[pool][m]["key"], P[pool][m]["params"], 0)
           for pool in ("protein", "nucleotide") for m in ("fast", "accurate")}
    sel[("protein", "fast")] = (ff["key"], ff["params"], FASTFIT_GAP_MODEL)
    sel[("protein", "accurate")] = (ef["key"], ef["params"], ef["params"]["run_gap_model"][0])
    nuc_src = []
    for m, path in NUCFIT.items():
        if path.exists():
            nraw = path.read_bytes()
            n = json.loads(nraw)
            sel[("nucleotide", m)] = (n["key"], n["params"], n["params"]["run_gap_model"][0])
            nuc_src.append(f" *         {path.relative_to(ROOT)} (sha256 {hashlib.sha256(nraw).hexdigest()}), nucleotide {m}\n")
    body = ",\n".join(entry(f"{pool} {m}", *sel[(pool, m)])
                      for pool in ("protein", "nucleotide") for m in ("fast", "accurate"))
    return f"""/* GENERATED by presets/release/emit_c_presets.py -- do not edit.
 * Source: presets/results/v2/presets_v2.json (sha256 {sha})
 *         presets/results/fastfit/final_pick.json (sha256 {ff_sha}), protein fast
 *         presets/results/ensconf/final_pick.json (sha256 {ef_sha}), protein accurate
{"".join(nuc_src)} * Kalign 4 modes: each the final fit on all training cases after grouped nested cross-validation;
 * presets_v2.json is used only for a mode whose pick file is absent.
 * protein fast: the fast-fit knee setting (member gap model, hybrid single-linkage guide tree);
 * protein accurate: the confidence re-fit's pick (the ensemble fit's 5 members, member 1 = protein fast, all on member
 * 1's guide tree with per-member single-linkage tree noise; consensus min_support 2; presets/results/ensconf/oof_gate.json).
 * nucleotide: the final-fit picks of presets/docs/nucfit/PLAN_nucfit.md where listed above (else the v2 settings).
 * Order: protein fast, protein accurate, nucleotide fast, nucleotide accurate. */
#ifndef KALIGN_MODE_PRESETS_H
#define KALIGN_MODE_PRESETS_H

#define KALIGN_PRESET_MAX_RUNS {MAX}

struct kalign_mode_preset {{
        int n_runs;
        int matrix[KALIGN_PRESET_MAX_RUNS];
        float gpo[KALIGN_PRESET_MAX_RUNS];
        float gpe[KALIGN_PRESET_MAX_RUNS];
        float tgpe[KALIGN_PRESET_MAX_RUNS];
        float vsm_amax[KALIGN_PRESET_MAX_RUNS];
        float tree_noise[KALIGN_PRESET_MAX_RUNS];
        int refine[KALIGN_PRESET_MAX_RUNS];
        float vsm_kappa[KALIGN_PRESET_MAX_RUNS];
        float vsm_dzero[KALIGN_PRESET_MAX_RUNS];
        int guide_dist[KALIGN_PRESET_MAX_RUNS];
        int gap_model;
        int tree_method[KALIGN_PRESET_MAX_RUNS];
        int tree_dist[KALIGN_PRESET_MAX_RUNS];
        int tree_alpha[KALIGN_PRESET_MAX_RUNS];
        int tree_maxlen[KALIGN_PRESET_MAX_RUNS];
        int tree_leaf[KALIGN_PRESET_MAX_RUNS];
        int tree_kmeans_dist[KALIGN_PRESET_MAX_RUNS];
        int realign_tree[KALIGN_PRESET_MAX_RUNS];
        float tree_sl_noise[KALIGN_PRESET_MAX_RUNS];
        float seq_weights[KALIGN_PRESET_MAX_RUNS];
        int realign[KALIGN_PRESET_MAX_RUNS];
        int consistency_anchors[KALIGN_PRESET_MAX_RUNS];
        float consistency_weight[KALIGN_PRESET_MAX_RUNS];
        int vsm_mode[KALIGN_PRESET_MAX_RUNS];
        int min_support;
        int consistency_merge;
        float consistency_merge_weight;
        int fill_in;
        float fill_anchor_frac;
        int fill_support;
}};

#define KALIGN_PRESET_PROTEIN_FAST        0
#define KALIGN_PRESET_PROTEIN_ACCURATE    1
#define KALIGN_PRESET_NUCLEOTIDE_FAST     2
#define KALIGN_PRESET_NUCLEOTIDE_ACCURATE 3

static const struct kalign_mode_preset kalign_mode_presets[4] = {{
{body}
}};

#endif
"""


if __name__ == "__main__":
    text = render()
    if "--check" in sys.argv:
        sys.exit(0 if OUT.exists() and OUT.read_text() == text else 1)
    OUT.write_text(text)
    print("written", OUT)
