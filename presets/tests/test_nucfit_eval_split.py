"""Reporting rule A2 (PLAN_nucfit.md section 7): real and simulated large-tier values, known answers."""
import numpy as np
import pytest

from recipe.nucfit_eval import Pool, summarise


def _meta_rows():
    meta, rows = {}, []
    # small pool: two groups of one benchmark; large: 2 Rfam families, 2 families per simulation size
    for g in ("g1", "g2"):
        meta[f"s_{g}"] = {"kind": "small", "group": g, "benchmark": "bralibase"}
    for k, fams in (("rfam", ("R1", "R2")), ("1000", ("S1", "S2")), ("10000", ("T1", "T2"))):
        for f in fams:
            meta[f"l_{f}"] = {"kind": "large", "family": f, "k": k}
    # label a: SP = TC = value below (so Q = value); label b: everything 0.5
    val_a = {"s_g1": 0.8, "s_g2": 0.6, "l_R1": 0.9, "l_R2": 0.7, "l_S1": 0.4, "l_S2": 0.2, "l_T1": 0.3, "l_T2": 0.1}
    for lab in ("a", "b"):
        for c in meta:
            v = val_a[c] if lab == "a" else 0.5
            rows.append({"label": lab, "case_id": c, "error": "", "recall": v, "tc": v, "precision": v,
                         "brier": None, "cpu_s": 1.0})
    return meta, rows


def test_split_known_answers():
    meta, rows = _meta_rows()
    P = Pool(rows, meta, ["a", "b"])
    out = summarise(P, ["a", "b"], [("a", "b")])
    pa = out["point"]["a"]
    # category means: rfam 0.8, 1000 0.3, 10000 0.2 -> real large 0.8, simulated 0.25, Q_large 1.3/3
    assert pa["q_rfam"] == pytest.approx(0.8)
    assert pa["q_1000"] == pytest.approx(0.3)
    assert pa["q_10000"] == pytest.approx(0.2)
    assert pa["q_sim"] == pytest.approx(0.25)
    assert pa["q_large"] == pytest.approx(1.3 / 3)
    assert pa["q_small"] == pytest.approx(0.7)
    d = out["diff"]["a - b"]
    assert d["q_rfam"]["diff"] == pytest.approx(0.3)
    assert d["q_sim"]["diff"] == pytest.approx(-0.25)
    # interval of the simulated difference: every value of a's simulations is below 0.5, so the whole interval is < 0
    assert d["q_sim"]["hi"] < 0
    assert d["q_rfam"]["lo"] <= d["q_rfam"]["diff"] <= d["q_rfam"]["hi"]


def test_no_sim_key_without_simulations():
    meta, rows = _meta_rows()
    keep = {c for c, m in meta.items() if m.get("k") not in ("1000", "10000")}
    meta = {c: m for c, m in meta.items() if c in keep}
    rows = [r for r in rows if r["case_id"] in keep]
    out = summarise(Pool(rows, meta, ["a", "b"]), ["a", "b"], [("a", "b")])
    assert "q_sim" not in out["point"]["a"]
    assert out["point"]["a"]["q_large"] == pytest.approx(0.8)
