"""Search space of the Kalign 3.5 nucleotide preset fit, reconstructed from optimize_unified.py.

Source (copied, never edited there): scripts/optimizers/optimize_unified.py of the manuscript
repository (identical to kalign commit 32e555d benchmarks/optimize_unified.py).
The profiles "rna" (BRAliBASE fit), "dna" (MDSA fit) and "nucleotide" (combined fit) have the same
bounds and matrices, so one space serves all three stages:

  per run k (k < 5): gpo [1,20], gpe [0.2,5], tgpe [0.05,3], noise [0,0.5],
                     matrix {K1, K20, K200} (Kimura 1/20/200 PAM), vsm_amax [0,3],
                     refine {NONE, ALL, CONFIDENT, INLINE}
  shared           : seq_weights [0,3], consistency_weight [0.5,5], n_runs {1,3,5},
                     consistency {0,1,2,3,5,8,10}, realign int [0,2], min_support int [0,5],
                     consistency_merge {0,1}, consistency_merge_weight [0.5,10]

decode() applies the masking rules of optimize_unified.decode_unified_params verbatim.
to_worker() adds the 3.6+ mechanism parameters at their off/default values (they did not exist in
3.5 and are not searched): vsm_mode 0, vsm_kappa 1, vsm_dzero 1, guide_dist 0, fill_in 0.
"""

import hashlib
import json

# Aligner constants (checked against kalign._core of the aligner at every job)
MATRIX_NUC_1PAM, MATRIX_NUC_20PAM, MATRIX_NUC_200PAM = 10, 11, 12
REFINE_NONE, REFINE_ALL, REFINE_CONFIDENT, REFINE_INLINE = 0, 1, 2, 3
EXPECT = {"MATRIX_NUC_1PAM": MATRIX_NUC_1PAM, "MATRIX_NUC_20PAM": MATRIX_NUC_20PAM,
          "MATRIX_NUC_200PAM": MATRIX_NUC_200PAM, "REFINE_NONE": REFINE_NONE, "REFINE_ALL": REFINE_ALL,
          "REFINE_CONFIDENT": REFINE_CONFIDENT, "REFINE_INLINE": REFINE_INLINE}
SEQ_TYPE_INT = MATRIX_NUC_200PAM  # optimize_unified: seq_type=_active_profile["seq_type_int"]

MAX_RUNS = 5
N_RUNS_MAP = [1, 3, 5]
CONSISTENCY_MAP = [0, 1, 2, 3, 5, 8, 10]
MATRIX_MAP = [MATRIX_NUC_1PAM, MATRIX_NUC_20PAM, MATRIX_NUC_200PAM]
MATRIX_NAMES = {MATRIX_NUC_1PAM: "K1", MATRIX_NUC_20PAM: "K20", MATRIX_NUC_200PAM: "K200"}
REFINE_MAP = [REFINE_NONE, REFINE_ALL, REFINE_CONFIDENT, REFINE_INLINE]
REFINE_NAMES = {REFINE_NONE: "N", REFINE_ALL: "A", REFINE_CONFIDENT: "C", REFINE_INLINE: "I"}

PER_RUN_LO = (1.0, 0.2, 0.05, 0.0)      # gpo, gpe, tgpe, noise
PER_RUN_HI = (20.0, 5.0, 3.0, 0.5)
SHARED_LO = (0.0, 0.0, 0.5)             # vsm_amax (per run), seq_weights, consistency_weight
SHARED_HI = (3.0, 3.0, 5.0)


def pymoo_vars(max_runs=MAX_RUNS):
    """optimize_unified.get_vars(max_runs) for the nucleotide profile (same names and order)."""
    from pymoo.core.variable import Choice, Integer, Real
    v = {}
    for k in range(max_runs):
        v[f"gpo_{k}"] = Real(bounds=(PER_RUN_LO[0], PER_RUN_HI[0]))
        v[f"gpe_{k}"] = Real(bounds=(PER_RUN_LO[1], PER_RUN_HI[1]))
        v[f"tgpe_{k}"] = Real(bounds=(PER_RUN_LO[2], PER_RUN_HI[2]))
        v[f"noise_{k}"] = Real(bounds=(PER_RUN_LO[3], PER_RUN_HI[3]))
        v[f"matrix_{k}"] = Choice(options=list(range(len(MATRIX_MAP))))
        v[f"vsm_amax_{k}"] = Real(bounds=(SHARED_LO[0], SHARED_HI[0]))
        v[f"refine_{k}"] = Choice(options=REFINE_MAP)
    v["seq_weights"] = Real(bounds=(SHARED_LO[1], SHARED_HI[1]))
    v["consistency_weight"] = Real(bounds=(SHARED_LO[2], SHARED_HI[2]))
    v["n_runs"] = Choice(options=N_RUNS_MAP)
    v["consistency"] = Choice(options=CONSISTENCY_MAP)
    v["realign"] = Integer(bounds=(0, 2))
    v["min_support"] = Integer(bounds=(0, max_runs))
    v["consistency_merge"] = Choice(options=[0, 1])
    v["consistency_merge_weight"] = Real(bounds=(0.5, 10.0))
    return v


def decode(x):
    """optimize_unified.decode_unified_params (masking rules verbatim)."""
    n = int(x["n_runs"])
    consistency = int(x["consistency"])
    realign = int(x["realign"])
    cm = int(x.get("consistency_merge", 0))
    cmw = float(x.get("consistency_merge_weight", 2.0))
    cw = float(x["consistency_weight"])
    p = {
        "n_runs": n,
        "run_gpo": [float(x[f"gpo_{k}"]) for k in range(n)],
        "run_gpe": [float(x[f"gpe_{k}"]) for k in range(n)],
        "run_tgpe": [float(x[f"tgpe_{k}"]) for k in range(n)],
        "run_noise": [float(x[f"noise_{k}"]) for k in range(n)],
        "run_types": [MATRIX_MAP[int(x[f"matrix_{k}"])] for k in range(n)],
        "run_vsm_amax": [float(x[f"vsm_amax_{k}"]) for k in range(n)],
        "run_refine": [int(x[f"refine_{k}"]) for k in range(n)],
        "seq_weights": float(x["seq_weights"]),
        "consistency": consistency, "consistency_weight": cw, "realign": realign,
        "min_support": int(x["min_support"]), "consistency_merge": cm, "consistency_merge_weight": cmw,
    }
    return mask(p)


def mask(p):
    n = p["n_runs"]
    if n == 1:
        p["run_noise"] = [0.0]
        p["min_support"] = 0
        p["consistency_merge"] = 0
    else:
        p["min_support"] = min(p["min_support"], n)
    if p["realign"] > 0:
        p["run_noise"] = [0.0] * n
    if p["consistency"] == 0:
        p["consistency_weight"] = 1.0
    if p["consistency_merge"] == 0:
        p["consistency_merge_weight"] = 2.0
    return p


def encode(p, max_runs=MAX_RUNS):
    """optimize_unified.encode_unified_params (unused slots copy run 0)."""
    x = {}
    for k in range(max_runs):
        j = k if k < p["n_runs"] else 0
        x[f"gpo_{k}"] = p["run_gpo"][j]
        x[f"gpe_{k}"] = p["run_gpe"][j]
        x[f"tgpe_{k}"] = p["run_tgpe"][j]
        x[f"noise_{k}"] = p["run_noise"][j]
        x[f"matrix_{k}"] = MATRIX_MAP.index(p["run_types"][j])
        x[f"vsm_amax_{k}"] = p["run_vsm_amax"][j]
        x[f"refine_{k}"] = p["run_refine"][j]
    for k in ("seq_weights", "consistency_weight", "n_runs", "consistency", "realign", "min_support",
              "consistency_merge", "consistency_merge_weight"):
        x[k] = p[k]
    return x


def to_worker(p):
    """Decoded params -> the keyword set of tools/nested_nuc/nuc_align_worker.py (3.6+ mechanisms off)."""
    n = p["n_runs"]
    return {**p, "vsm_mode": 0, "run_vsm_kappa": [1.0] * n, "run_vsm_dzero": [1.0] * n,
            "run_guide_dist": [0] * n, "fill_in": 0, "fill_anchor_frac": 0.0, "fill_support": 1}


def key(p) -> str:
    return hashlib.sha1(json.dumps(p, sort_keys=True).encode()).hexdigest()[:16]


def short(p) -> str:
    n = p["n_runs"]
    mats = "/".join(MATRIX_NAMES[t] for t in p["run_types"])
    ref = "/".join(REFINE_NAMES[r] for r in p["run_refine"])
    s = (f"{'single' if n == 1 else f'ens{n}'} {mats} gpo={'/'.join(f'{v:.2f}' for v in p['run_gpo'])} "
         f"vsm={'/'.join(f'{v:.2f}' for v in p['run_vsm_amax'])} sw={p['seq_weights']:.2f} "
         f"c={p['consistency']} re={p['realign']} ref={ref}")
    if n > 1:
        s += f" ms={p['min_support']} cm={p['consistency_merge']}"
    return s


def _preset(runs, seq_weights, realign, consistency=0, consistency_weight=2.0, min_support=0):
    """runs: (matrix, gpo, gpe, tgpe, vsm_amax, refine) as in preset_run() of aln_wrap.c."""
    n = len(runs)
    return mask({
        "n_runs": n, "run_gpo": [r[1] for r in runs], "run_gpe": [r[2] for r in runs],
        "run_tgpe": [r[3] for r in runs], "run_noise": [0.0] * n, "run_types": [r[0] for r in runs],
        "run_vsm_amax": [r[4] for r in runs], "run_refine": [r[5] for r in runs],
        "seq_weights": seq_weights, "consistency": consistency, "consistency_weight": consistency_weight,
        "realign": realign, "min_support": min_support, "consistency_merge": 0, "consistency_merge_weight": 2.0,
    })


K1, K20, K200 = MATRIX_NUC_1PAM, MATRIX_NUC_20PAM, MATRIX_NUC_200PAM
# Shipped nucleotide presets: lib/src/aln_wrap.c preset_nucleotide(), unchanged from kalign 9208173
# (3.5) to redesign/experiment. Tree seeds 42+k with noise 0 = ensemble_custom_file_to_file(seed=42).
SHIPPED = {
    "fast": _preset([(K200, 19.873, 0.642, 2.582, 0.874, REFINE_NONE)], seq_weights=0.199, realign=1),
    "default": _preset([
        (K20, 19.275, 0.628, 2.854, 0.445, REFINE_NONE),
        (K200, 16.581, 0.274, 2.491, 0.706, REFINE_NONE),
        (K200, 16.499, 3.754, 2.896, 0.491, REFINE_CONFIDENT),
    ], seq_weights=0.448, realign=1, min_support=2),
    "recall": _preset([(K200, 19.882, 0.646, 2.854, 0.874, REFINE_CONFIDENT)], seq_weights=0.449, realign=1,
                      consistency=1, consistency_weight=1.616),
    "accurate": _preset([
        (K200, 19.475, 0.728, 2.856, 0.492, REFINE_CONFIDENT),
        (K1, 19.377, 1.379, 2.486, 1.211, REFINE_CONFIDENT),
        (K200, 15.869, 0.353, 2.430, 0.721, REFINE_NONE),
        (K200, 15.989, 2.331, 2.764, 0.188, REFINE_CONFIDENT),
        (K1, 5.346, 3.483, 0.610, 1.406, REFINE_NONE),
    ], seq_weights=0.590, realign=1, consistency=2, consistency_weight=0.800, min_support=3),
}
MODES = ("fast", "default", "recall", "accurate")
