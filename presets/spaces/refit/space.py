"""Re-fit search space (pre-registration section 6): existing protein preset parameters of
optimize_unified.py (manuscript repo, copied and adapted) + the new mechanism parameters of
redesign/integration.

Existing (as in optimize_unified.py, same bounds):
  per run k : gpo_k [2,15], gpe_k [0.5,5], tgpe_k [0.1,3], noise_k [0,0.5],
              matrix_k {PFASUM43, PFASUM60, CorBLOSUM66}, vsm_amax_k [0,5],
              refine_k {NONE, ALL, CONFIDENT, INLINE}
  shared    : n_runs {1,3,5}, seq_weights [0,5], consistency {0,1,2,3,5,8,10},
              consistency_weight [0.5,5], realign {0,1,2}, min_support {0..5},
              consistency_merge {0,1}, consistency_merge_weight [0.5,10]
New mechanism parameters:
  M2 VSM    : vsm_mode {0 legacy, 1 per-merge} (shared; the API takes one value),
              vsm_kappa_k [0,1], vsm_dzero_k [0.5,10] (per run, as vsm_amax)
  M3 fill-in: fill_in {0 off, 1 realign, 2 project, 3 attach}, fill_anchor_frac [0,1],
              fill_support {1..5} (shared; ensemble-level)
  M4 guide  : guide_dist_k {0 BPM, 1 pairwise alignment} (per run)

Masking (decode): parameters that cannot affect the output are set to a canonical value so
equivalent individuals share one cache key. Rules of optimize_unified.py are kept; new rules:
vsm_mode 0 -> kappa/dzero canonical; vsm_mode 1 -> vsm_amax canonical (ignored by the C code);
fill-in only for ensembles without consistency_merge; fill_anchor_frac only for fill_in 1/2;
fill_support only for fill_in 3 (capped at n_runs).
"""

import hashlib
import json

# Aligner constants (verified against the aligner's kalign._core at every job by align_worker.py)
MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66 = 1, 2, 3
REFINE_NONE, REFINE_ALL, REFINE_CONFIDENT, REFINE_INLINE = 0, 1, 2, 3
FILL_OFF, FILL_REALIGN, FILL_PROJECT, FILL_ATTACH = 0, 1, 2, 3
EXPECT = {
    "MATRIX_PFASUM43": MATRIX_PFASUM43, "MATRIX_PFASUM60": MATRIX_PFASUM60,
    "MATRIX_CORBLOSUM66": MATRIX_CORBLOSUM66, "REFINE_NONE": REFINE_NONE, "REFINE_ALL": REFINE_ALL,
    "REFINE_CONFIDENT": REFINE_CONFIDENT, "REFINE_INLINE": REFINE_INLINE, "FILL_OFF": FILL_OFF,
    "FILL_REALIGN": FILL_REALIGN, "FILL_PROJECT": FILL_PROJECT, "FILL_ATTACH": FILL_ATTACH,
}

MAX_RUNS = 5
N_RUNS_MAP = [1, 3, 5]
CONSISTENCY_MAP = [0, 1, 2, 3, 5, 8, 10]
MATRIX_MAP = [MATRIX_PFASUM43, MATRIX_PFASUM60, MATRIX_CORBLOSUM66]
MATRIX_NAMES = {MATRIX_PFASUM43: "P43", MATRIX_PFASUM60: "P60", MATRIX_CORBLOSUM66: "CB66"}
REFINE_MAP = [REFINE_NONE, REFINE_ALL, REFINE_CONFIDENT, REFINE_INLINE]

# (name, kind, bounds/options) ; kind in real/int/choice
PER_RUN = [
    ("gpo", "real", (2.0, 15.0)),
    ("gpe", "real", (0.5, 5.0)),
    ("tgpe", "real", (0.1, 3.0)),
    ("noise", "real", (0.0, 0.5)),
    ("matrix", "choice", list(range(len(MATRIX_MAP)))),
    ("vsm_amax", "real", (0.0, 5.0)),
    ("refine", "choice", REFINE_MAP),
    # new
    ("vsm_kappa", "real", (0.0, 1.0)),
    ("vsm_dzero", "real", (0.5, 10.0)),
    ("guide_dist", "choice", [0, 1]),
]
SHARED = [
    ("n_runs", "choice", N_RUNS_MAP),
    ("seq_weights", "real", (0.0, 5.0)),
    ("consistency", "choice", CONSISTENCY_MAP),
    ("consistency_weight", "real", (0.5, 5.0)),
    ("realign", "int", (0, 2)),
    ("min_support", "int", (0, MAX_RUNS)),
    ("consistency_merge", "choice", [0, 1]),
    ("consistency_merge_weight", "real", (0.5, 10.0)),
    # new
    ("vsm_mode", "choice", [0, 1]),
    ("fill_in", "choice", [FILL_OFF, FILL_REALIGN, FILL_PROJECT, FILL_ATTACH]),
    ("fill_anchor_frac", "real", (0.0, 1.0)),
    ("fill_support", "int", (1, MAX_RUNS)),
]
NEW_PARAMS = ["vsm_mode", "vsm_kappa_k", "vsm_dzero_k", "fill_in", "fill_anchor_frac", "fill_support",
              "guide_dist_k"]

# Canonical values for masked parameters (= kalign_run_config / ensemble defaults)
CANON = {"vsm_kappa": 1.0, "vsm_dzero": 1.0, "vsm_amax": 0.0, "fill_anchor_frac": 0.0, "fill_support": 1,
         "consistency_weight_off": 1.0, "consistency_merge_weight": 2.0}


def pymoo_vars():
    from pymoo.core.variable import Choice, Integer, Real

    def mk(kind, b):
        return Real(bounds=b) if kind == "real" else Integer(bounds=b) if kind == "int" else Choice(options=b)

    v = {}
    for k in range(MAX_RUNS):
        for name, kind, b in PER_RUN:
            v[f"{name}_{k}"] = mk(kind, b)
    for name, kind, b in SHARED:
        v[name] = mk(kind, b)
    return v


def decode(x):
    """Mixed-variable dict -> decoded params (what align_worker.py consumes). Applies masking."""
    n = int(x["n_runs"])
    runs = {f"run_{name}": [x[f"{name}_{k}"] for k in range(n)] for name, _, _ in PER_RUN}
    p = {
        "n_runs": n,
        "run_gpo": [float(v) for v in runs["run_gpo"]],
        "run_gpe": [float(v) for v in runs["run_gpe"]],
        "run_tgpe": [float(v) for v in runs["run_tgpe"]],
        "run_noise": [float(v) for v in runs["run_noise"]],
        "run_types": [MATRIX_MAP[int(v)] for v in runs["run_matrix"]],
        "run_vsm_amax": [float(v) for v in runs["run_vsm_amax"]],
        "run_refine": [int(v) for v in runs["run_refine"]],
        "run_vsm_kappa": [float(v) for v in runs["run_vsm_kappa"]],
        "run_vsm_dzero": [float(v) for v in runs["run_vsm_dzero"]],
        "run_guide_dist": [int(v) for v in runs["run_guide_dist"]],
        "seq_weights": float(x["seq_weights"]),
        "consistency": int(x["consistency"]),
        "consistency_weight": float(x["consistency_weight"]),
        "realign": int(x["realign"]),
        "min_support": int(x["min_support"]),
        "consistency_merge": int(x["consistency_merge"]),
        "consistency_merge_weight": float(x["consistency_merge_weight"]),
        "vsm_mode": int(x["vsm_mode"]),
        "fill_in": int(x["fill_in"]),
        "fill_anchor_frac": float(x["fill_anchor_frac"]),
        "fill_support": int(x["fill_support"]),
    }
    return mask(p)


def mask(p):
    """Canonicalise parameters that cannot change the output (rules of optimize_unified.py + new)."""
    n = p["n_runs"]
    # optimize_unified.py rules
    if n == 1:
        p["run_noise"] = [0.0]
        p["min_support"] = 0
        p["consistency_merge"] = 0
    else:
        p["min_support"] = min(p["min_support"], n)
    if p["realign"] > 0:
        p["run_noise"] = [0.0] * n
    if p["consistency"] == 0:
        p["consistency_weight"] = CANON["consistency_weight_off"]
    if p["consistency_merge"] == 0:
        p["consistency_merge_weight"] = CANON["consistency_merge_weight"]
    # new: VSM
    if p["vsm_mode"] == 0:
        p["run_vsm_kappa"] = [CANON["vsm_kappa"]] * n
        p["run_vsm_dzero"] = [CANON["vsm_dzero"]] * n
    else:
        p["run_vsm_amax"] = [CANON["vsm_amax"]] * n
    # new: fill-in acts on the POAR consensus of an ensemble only
    if n == 1 or p["consistency_merge"] == 1:
        p["fill_in"] = FILL_OFF
    if p["fill_in"] not in (FILL_REALIGN, FILL_PROJECT):
        p["fill_anchor_frac"] = CANON["fill_anchor_frac"]
    if p["fill_in"] != FILL_ATTACH:
        p["fill_support"] = CANON["fill_support"]
    else:
        p["fill_support"] = max(1, min(p["fill_support"], n))
    return p


def encode(p):
    """Decoded params -> mixed-variable dict (unused run slots copy run 0)."""
    x = {}
    per = {"gpo": "run_gpo", "gpe": "run_gpe", "tgpe": "run_tgpe", "noise": "run_noise",
           "vsm_amax": "run_vsm_amax", "refine": "run_refine", "vsm_kappa": "run_vsm_kappa",
           "vsm_dzero": "run_vsm_dzero", "guide_dist": "run_guide_dist"}
    for k in range(MAX_RUNS):
        j = k if k < p["n_runs"] else 0
        for name, key in per.items():
            x[f"{name}_{k}"] = p[key][j]
        x[f"matrix_{k}"] = MATRIX_MAP.index(p["run_types"][j])
    for name, _, _ in SHARED:
        x[name] = p[name]
    return x


def key(p) -> str:
    return hashlib.sha1(json.dumps(p, sort_keys=True).encode()).hexdigest()[:16]


def short(p) -> str:
    n = p["n_runs"]
    mats = "/".join(MATRIX_NAMES[t] for t in p["run_types"])
    gd = "".join(str(g) for g in p["run_guide_dist"])
    vsm = (f"vsm0 amax={'/'.join(f'{v:.2f}' for v in p['run_vsm_amax'])}" if p["vsm_mode"] == 0 else
           f"vsm1 k={'/'.join(f'{v:.2f}' for v in p['run_vsm_kappa'])} "
           f"d0={'/'.join(f'{v:.1f}' for v in p['run_vsm_dzero'])}")
    s = (f"n={n} {mats} gpo={'/'.join(f'{v:.2f}' for v in p['run_gpo'])} {vsm} gd={gd} "
         f"c={p['consistency']} re={p['realign']} sw={p['seq_weights']:.2f}")
    if n > 1:
        s += f" ms={p['min_support']} cm={p['consistency_merge']} fill={p['fill_in']}"
    return s


def _preset(runs, seq_weights, realign, consistency, consistency_weight, min_support):
    """Build decoded params for a shipped protein preset (lib/src/aln_wrap.c preset_protein)."""
    n = len(runs)
    return mask({
        "n_runs": n,
        "run_gpo": [r[1] for r in runs], "run_gpe": [r[2] for r in runs], "run_tgpe": [r[3] for r in runs],
        "run_noise": [0.0] * n, "run_types": [r[0] for r in runs], "run_vsm_amax": [r[4] for r in runs],
        "run_refine": [r[5] for r in runs],
        "run_vsm_kappa": [1.0] * n, "run_vsm_dzero": [1.0] * n, "run_guide_dist": [0] * n,
        "seq_weights": seq_weights, "consistency": consistency, "consistency_weight": consistency_weight,
        "realign": realign, "min_support": min_support, "consistency_merge": 0,
        "consistency_merge_weight": 2.0, "vsm_mode": 0, "fill_in": FILL_OFF, "fill_anchor_frac": 0.0,
        "fill_support": 1,
    })


# Shipped protein presets at 3caadc7 (lib/src/aln_wrap.c, preset_protein). Tree seeds 42+k are
# produced by ensemble_custom_file_to_file(seed=42); with noise 0 the seed has no effect.
# (matrix, gpo, gpe, tgpe, vsm_amax, refine)
SHIPPED = {
    "fast": _preset([(MATRIX_PFASUM60, 8.626, 0.843, 0.433, 0.592, REFINE_NONE)],
                    seq_weights=1.534, realign=0, consistency=0, consistency_weight=2.0, min_support=0),
    "default": _preset([(MATRIX_PFASUM60, 8.121, 0.684, 0.560, 1.661, REFINE_INLINE)],
                       seq_weights=2.356, realign=0, consistency=1, consistency_weight=1.640, min_support=0),
    "recall": _preset([
        (MATRIX_CORBLOSUM66, 5.416, 1.071, 1.620, 2.536, REFINE_ALL),
        (MATRIX_CORBLOSUM66, 12.091, 1.024, 1.284, 2.078, REFINE_INLINE),
        (MATRIX_PFASUM60, 6.970, 2.919, 0.632, 0.877, REFINE_NONE),
        (MATRIX_PFASUM43, 6.173, 1.102, 0.510, 1.276, REFINE_ALL),
        (MATRIX_PFASUM43, 5.278, 1.764, 1.088, 2.315, REFINE_CONFIDENT),
    ], seq_weights=0.282, realign=1, consistency=3, consistency_weight=2.120, min_support=0),
    "accurate": _preset([
        (MATRIX_PFASUM43, 8.682, 0.650, 1.465, 2.532, REFINE_INLINE),
        (MATRIX_PFASUM43, 6.148, 1.174, 1.297, 2.011, REFINE_ALL),
        (MATRIX_CORBLOSUM66, 3.622, 1.004, 0.631, 0.872, REFINE_CONFIDENT),
        (MATRIX_CORBLOSUM66, 5.988, 0.552, 0.495, 1.914, REFINE_CONFIDENT),
        (MATRIX_PFASUM43, 13.939, 2.629, 1.406, 1.830, REFINE_NONE),
    ], seq_weights=1.468, realign=1, consistency=8, consistency_weight=2.457, min_support=3),
}


# =============================================================================================
# Two-stage re-fit (PREREGISTRATION.md, Amendment A1 to section 6)
# =============================================================================================
# Stage 1: single-run presets (n_runs = 1). Variables: the per-run parameters of run slot 0
# (matrix, gpo, gpe, tgpe, vsm_amax, vsm_kappa, vsm_dzero, refine, guide_dist) and the single-run
# shared parameters below. noise is masked to 0 for n_runs = 1 (optimize_unified.py rule).
STAGE1_PER_RUN = [v for v in PER_RUN if v[0] != "noise"]
STAGE1_SHARED = [v for v in SHARED if v[0] in
                 ("seq_weights", "consistency", "consistency_weight", "realign", "vsm_mode")]
# Stage 2: 5-run ensembles. Run k = stage-1 Pareto member m_k (all its per-run parameters);
# the ensemble-wide single-run parameters (seq_weights, consistency, consistency_weight, realign,
# vsm_mode) are those of the lead member m_0. Ensemble parameters are searched.
STAGE2_N_RUNS = 5
STAGE2_ENSEMBLE = [v for v in SHARED if v[0] in ("min_support", "fill_in", "fill_anchor_frac", "fill_support")]
STAGE2_FIXED = {"consistency_merge": 0, "consistency_merge_weight": CANON["consistency_merge_weight"],
                "run_noise": 0.0}


def _mk(kind, b):
    from pymoo.core.variable import Choice, Integer, Real
    return Real(bounds=b) if kind == "real" else Integer(bounds=b) if kind == "int" else Choice(options=b)


def stage1_vars():
    v = {name: _mk(kind, b) for name, kind, b in STAGE1_PER_RUN}
    v.update({name: _mk(kind, b) for name, kind, b in STAGE1_SHARED})
    return v


def stage1_decode(x):
    p = {
        "n_runs": 1,
        "run_gpo": [float(x["gpo"])], "run_gpe": [float(x["gpe"])], "run_tgpe": [float(x["tgpe"])],
        "run_noise": [0.0], "run_types": [MATRIX_MAP[int(x["matrix"])]],
        "run_vsm_amax": [float(x["vsm_amax"])], "run_refine": [int(x["refine"])],
        "run_vsm_kappa": [float(x["vsm_kappa"])], "run_vsm_dzero": [float(x["vsm_dzero"])],
        "run_guide_dist": [int(x["guide_dist"])],
        "seq_weights": float(x["seq_weights"]), "consistency": int(x["consistency"]),
        "consistency_weight": float(x["consistency_weight"]), "realign": int(x["realign"]),
        "min_support": 0, "consistency_merge": 0,
        "consistency_merge_weight": CANON["consistency_merge_weight"],
        "vsm_mode": int(x["vsm_mode"]), "fill_in": FILL_OFF, "fill_anchor_frac": 0.0, "fill_support": 1,
    }
    return mask(p)


def stage1_encode(p):
    assert p["n_runs"] == 1
    x = {"gpo": p["run_gpo"][0], "gpe": p["run_gpe"][0], "tgpe": p["run_tgpe"][0],
         "matrix": MATRIX_MAP.index(p["run_types"][0]), "vsm_amax": p["run_vsm_amax"][0],
         "refine": p["run_refine"][0], "vsm_kappa": p["run_vsm_kappa"][0],
         "vsm_dzero": p["run_vsm_dzero"][0], "guide_dist": p["run_guide_dist"][0]}
    for name, _, _ in STAGE1_SHARED:
        x[name] = p[name]
    return x


def stage2_vars(pool_size):
    v = {f"m_{k}": _mk("choice", list(range(pool_size))) for k in range(STAGE2_N_RUNS)}
    v.update({name: _mk(kind, b) for name, kind, b in STAGE2_ENSEMBLE})
    return v


def stage2_members(x, pool):
    """Member indices actually used. m_0 is the lead. Because vsm_mode is ensemble-wide, a member
    m_k (k >= 1) whose vsm_mode differs from the lead's is replaced by compat[m_k mod len(compat)],
    compat = pool indices (ascending) with the lead's vsm_mode (the lead itself is always in it)."""
    lead = int(x["m_0"])
    mode = pool[lead]["vsm_mode"]
    compat = [i for i, q in enumerate(pool) if q["vsm_mode"] == mode]
    out = [lead]
    for k in range(1, STAGE2_N_RUNS):
        m = int(x[f"m_{k}"])
        out.append(m if pool[m]["vsm_mode"] == mode else compat[m % len(compat)])
    return out


def stage2_decode(x, pool):
    """pool: list of decoded stage-1 params (n_runs = 1), in the frozen pool order."""
    idx = stage2_members(x, pool)
    lead = pool[idx[0]]
    per = ["run_gpo", "run_gpe", "run_tgpe", "run_types", "run_vsm_amax", "run_refine", "run_vsm_kappa",
           "run_vsm_dzero", "run_guide_dist"]
    p = {k: [pool[i][k][0] for i in idx] for k in per}
    p.update({
        "n_runs": STAGE2_N_RUNS, "run_noise": [STAGE2_FIXED["run_noise"]] * STAGE2_N_RUNS,
        "seq_weights": lead["seq_weights"], "consistency": lead["consistency"],
        "consistency_weight": lead["consistency_weight"], "realign": lead["realign"],
        "vsm_mode": lead["vsm_mode"],
        "consistency_merge": STAGE2_FIXED["consistency_merge"],
        "consistency_merge_weight": STAGE2_FIXED["consistency_merge_weight"],
        "min_support": int(x["min_support"]), "fill_in": int(x["fill_in"]),
        "fill_anchor_frac": float(x["fill_anchor_frac"]), "fill_support": int(x["fill_support"]),
    })
    return mask(p)
