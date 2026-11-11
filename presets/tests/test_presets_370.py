"""3.7.0 presets as parameter dicts: transcription check against lib/src/aln_wrap.c, validity in
the search spaces, stage-1/-2 reconstruction.  The byte-for-byte CLI comparison needs the
workstation build: recipe/verify_presets_370.py (results in docs/build/presets_370.md)."""
import re
from pathlib import Path

import pytest

from recipe import presets_370 as P

C_SRC = Path(__file__).resolve().parents[2] / "lib" / "src" / "aln_wrap.c"
# Kalign 4 replaced the 3.7.0 presets; the 3.7.0 source is read from the 3.7.0 release commit.
REV_370 = "e2c0605"
NUM = r"-?\d+\.\d+f?|-?\d+"
MATRIX = {"KALIGN_MATRIX_PFASUM43": 1, "KALIGN_MATRIX_PFASUM60": 2, "KALIGN_MATRIX_CORBLOSUM66": 3,
          "KALIGN_MATRIX_NUC_1PAM": 10, "KALIGN_MATRIX_NUC_20PAM": 11, "KALIGN_MATRIX_NUC_200PAM": 12}
REFINE = {"KALIGN_REFINE_NONE": 0, "KALIGN_REFINE_ALL": 1, "KALIGN_REFINE_CONFIDENT": 2,
          "KALIGN_REFINE_INLINE": 3}
GUIDE = {"KALIGN_GUIDE_DIST_BPM": 0, "KALIGN_GUIDE_DIST_PAIRALN": 1}


def _f(tok):
    return float(tok.rstrip("f"))


def _func(src, name):
    i = src.index(f"static int {name}(")
    j = src.index("\n}\n", i)
    return src[i:j]


def _calls(body):
    """[(matrix, gpo, gpe, tgpe, amax, sw, realign, refine, seed, noise)] of each preset_run( call."""
    out = []
    for m in re.finditer(r"preset_run\(&runs\[\d\],\s*(\w+),\s*(%s),\s*(%s),\s*(%s),\s*(%s),\s*(%s),\s*(\d+),\s*(\w+),\s*(\d+),\s*(%s)\)"
                         % ((NUM,) * 5 + (NUM,)), body):
        g = m.groups()
        out.append((MATRIX[g[0]], _f(g[1]), _f(g[2]), _f(g[3]), _f(g[4]), _f(g[5]), int(g[6]),
                    REFINE[g[7]], int(g[8]), _f(g[9])))
    return out


@pytest.fixture(scope="module")
def csrc():
    import subprocess
    r = subprocess.run(["git", "-C", str(C_SRC.parents[2]), "show", f"{REV_370}:lib/src/aln_wrap.c"],
                       capture_output=True, text=True)
    if r.returncode:
        pytest.skip("3.7.0 source not available (git show failed)")
    return r.stdout


def _check_runs(p, calls):
    assert p["n_runs"] == len(calls)
    for k, (mat, gpo, gpe, tgpe, amax, sw, re_, ref, seed, noise) in enumerate(calls):
        assert p["run_types"][k] == mat
        assert p["run_gpo"][k] == pytest.approx(gpo, rel=0, abs=0)
        assert p["run_gpe"][k] == gpe and p["run_tgpe"][k] == tgpe
        assert p["run_vsm_amax"][k] == amax and p["run_refine"][k] == ref
        assert p["seq_weights"] == sw and p["realign"] == re_
        assert noise == 0.0 and p["run_noise"][k] == 0.0


def test_nucleotide_matches_c(csrc):
    body = _func(csrc, "preset_nucleotide")
    parts = {m: body[body.index(f'strcasecmp(m, "{m}")'):] for m in P.MODES["nucleotide"]}
    order = sorted(parts, key=lambda m: -len(parts[m]))   # as they appear in the source
    for i, m in enumerate(order):
        seg = parts[m] if i == len(order) - 1 else parts[m][:len(parts[m]) - len(parts[order[i + 1]])]
        calls = _calls(seg)
        p = P.params_370("nucleotide", m)
        _check_runs(p, calls)
        mc = re.search(r"preset_consistency\(runs,\s*(\d+),\s*(\d+),\s*(%s)\)" % NUM, seg)
        ms = re.search(r"ens->min_support = (\d+)", seg)
        assert p["consistency"] == (int(mc.group(2)) if mc else 0)
        assert p["consistency_weight"] == (_f(mc.group(3)) if mc else 1.0)
        assert p["min_support"] == (int(ms.group(1)) if ms else 0)
        assert p["vsm_mode"] == 0 and p["fill_in"] == 0 and p["consistency_merge"] == 0


def test_protein_fast_default_match_c(csrc):
    body = _func(csrc, "preset_protein")
    fast = body[body.index('"fast"'):body.index('"default"')]
    dflt = body[body.index('"default"'):body.index('"accurate"')]
    _check_runs(P.params_370("protein", "fast"), _calls(fast))
    p = P.params_370("protein", "default")
    _check_runs(p, _calls(dflt))
    mc = re.search(r"preset_consistency\(runs,\s*1,\s*(\d+),\s*(%s)\)" % NUM, dflt)
    assert p["consistency"] == int(mc.group(1)) and p["consistency_weight"] == _f(mc.group(2))
    assert re.search(r"vsm_mode = 0", dflt) and p["vsm_mode"] == 0 and p["run_guide_dist"] == [0]
    assert p["run_vsm_kappa"] == [1.0] and p["run_vsm_dzero"] == [1.0]


def test_protein_ensemble_matches_c(csrc):
    i = csrc.index("static void preset_protein_ensemble(")
    body = csrc[i:csrc.index("\n}\n", i)]

    def arr(name, conv):
        m = re.search(r"static const \w+\s+%s\[5\]\s*=\s*\{([^}]*)\}" % name, body)
        return [conv(t.strip()) for t in m.group(1).split(",")]

    p = P.params_370("protein", "accurate")
    assert p["run_types"] == arr("matrix", lambda t: MATRIX[t])
    assert p["run_gpo"] == arr("gpo", _f) and p["run_gpe"] == arr("gpe", _f) and p["run_tgpe"] == arr("tgpe", _f)
    assert p["run_refine"] == arr("refine", lambda t: REFINE[t])
    assert p["run_vsm_kappa"] == arr("kappa", _f) and p["run_vsm_dzero"] == arr("dzero", _f)
    assert p["run_guide_dist"] == arr("guide", lambda t: GUIDE[t])
    sw = re.search(r"(%s)\s*,\s*/\* seq_weights \*/" % NUM, body)
    assert p["seq_weights"] == _f(sw.group(1))
    assert "0.0f,                 /* vsm_amax" in body and p["run_vsm_amax"] == [0.0] * 5
    assert "runs[k].vsm_mode = 1" in body and p["vsm_mode"] == 1
    assert "42 + (uint64_t)k, 0.0f" in body and p["run_noise"] == [0.0] * 5
    assert "runs[k].consistency_anchors = 0" in body and p["consistency"] == 0
    assert "runs[k].consistency_weight = 1.0f" in body and p["consistency_weight"] == 1.0
    assert "ens->min_support = 2" in body and p["min_support"] == 2
    assert "ens->fill_in = KALIGN_FILL_OFF" in body and p["fill_in"] == 0
    assert re.search(r"preset_run\(&runs\[k\].*?\n.*?\n.*?\n\s*0, refine\[k\]", body, re.S)  # realign 0
    assert p["realign"] == 0
    # recall is an alias of accurate
    assert re.search(r'strcasecmp\(m, "accurate"\) == 0 \|\| strcasecmp\(m, "recall"\) == 0', csrc)
    assert P.params_370("protein", "recall") == p


def test_defaults_used_for_unset_fields(csrc):
    i = csrc.index("struct kalign_run_config kalign_run_config_defaults")
    body = csrc[i:csrc.index("\n}\n", i)]
    for frag in ("cfg.tree_noise = 0.0f", "cfg.vsm_kappa = 1.0f", "cfg.vsm_dzero = 1.0f",
                 "cfg.guide_dist = KALIGN_GUIDE_DIST_BPM", "cfg.consistency_anchors = 0"):
        assert frag in body
    j = csrc.index("struct kalign_ensemble_config kalign_ensemble_config_defaults")
    eb = csrc[j:csrc.index("\n}\n", j)]
    for frag in ("ens.consistency_merge = 0", "ens.consistency_merge_weight = 2.0f", "ens.fill_in = 0",
                 "ens.fill_anchor_frac = 0.0f", "ens.fill_support = 1"):
        assert frag in eb


def test_worker_keys_complete():
    need = {"run_gpo", "run_gpe", "run_tgpe", "run_noise", "run_types", "run_vsm_amax", "run_refine",
            "run_vsm_kappa", "run_vsm_dzero", "run_guide_dist", "min_support", "realign", "seq_weights",
            "consistency", "consistency_weight", "consistency_merge", "consistency_merge_weight",
            "vsm_mode", "fill_in", "fill_anchor_frac", "fill_support"}
    for pool, mode, p in P.all_presets():
        assert need <= set(p), (pool, mode)
        n = p["n_runs"]
        for k in need:
            if k.startswith("run_"):
                assert len(p[k]) == n


def test_params_are_copies():
    a = P.params_370("protein", "fast")
    a["run_gpo"][0] = -1
    assert P.params_370("protein", "fast")["run_gpo"][0] == 8.626


@pytest.mark.parametrize("pool,mode", [(p, m) for p, ms in P.MODES.items() for m in ms])
def test_valid_point_of_search_space(pool, mode):
    try:
        P.load_space(pool)
    except FileNotFoundError:
        pytest.skip("search space module not available")
    assert P.check_in_space(pool, P.params_370(pool, mode)) == []


def test_out_of_bounds_is_detected():
    try:
        P.load_space("protein")
    except FileNotFoundError:
        pytest.skip("search space module not available")
    p = P.params_370("protein", "fast")
    p["run_gpo"][0] = 99.0
    assert any("gpo_0" in s for s in P.check_in_space("protein", p))


def test_protein_stage1_members_roundtrip_and_stage2_reproduces():
    try:
        sp = P.load_space("protein")
    except FileNotFoundError:
        pytest.skip("search space module not available")
    acc = P.params_370("protein", "accurate")
    members = P.protein_stage1_members(acc)
    for m in members:
        assert sp.stage1_decode(sp.stage1_encode(m)) == m
    x = {f"m_{k}": k for k in range(5)}
    x.update(min_support=acc["min_support"], fill_in=0, fill_anchor_frac=0.0, fill_support=1)
    assert sp.stage2_decode(x, members) == acc
    d = P.params_370("protein", "default")
    assert sp.stage1_decode(sp.stage1_encode(d)) == d
    f = P.params_370("protein", "fast")
    assert sp.stage1_decode(sp.stage1_encode(f)) == f
