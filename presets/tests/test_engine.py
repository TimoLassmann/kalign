"""Tests of the fitting engine (recipe.fit / nested / pool_results / decide) with fake evaluations."""
import json
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pymoo")

from recipe import cache as C  # noqa: E402
from recipe import decide as D  # noqa: E402
from recipe import evaluate as E  # noqa: E402
from recipe import fit as F  # noqa: E402
from recipe import pool_results as PR  # noqa: E402
from recipe import presets_370 as P370  # noqa: E402

try:
    P370.load_space("protein")
    HAVE_SPACE = True
except FileNotFoundError:
    HAVE_SPACE = False
needs_space = pytest.mark.skipif(not HAVE_SPACE, reason="search-space modules not found")


# ------------------------------------------------------------------------------------ HV stop
def test_hv_stop_min_max_and_tolerance():
    flat = [1.0] * 30
    assert F.hv_stop(flat[:19], 20, 100) == (False, "")           # below the minimum
    stop, why = F.hv_stop(flat[:20], 20, 100)
    assert stop and "hypervolume" in why                           # flat at the minimum -> stop
    rising = [1.00021 ** g for g in range(100)]                    # +0.105 % per 5 generations
    for g in range(20, 100):
        assert not F.hv_stop(rising[:g], 20, 100)[0]
    assert F.hv_stop(rising[:100], 20, 100) == (True, "max generations (100)")
    slow = [1.00019 ** g for g in range(100)]                      # +0.095 % per 5 generations
    assert not F.hv_stop(slow[:19], 20, 100)[0] and F.hv_stop(slow[:20], 20, 100)[0]
    # a late plateau: stops 5 generations after the improvement ends
    late = [1.01 ** g for g in range(30)] + [1.01 ** 29] * 10
    stops = [g for g in range(20, 41) if F.hv_stop(late[:g], 20, 100)[0]]
    assert stops[0] == 35
    assert not F.hv_stop([0.0] * 30, 20, 100)[0]                   # empty box: run on to the maximum


def test_hypervolume_box_and_reference_filter():
    hv = F.hypervolume([[0.1, 0.1, 0.1, 0.1]], 1.1)
    assert hv == pytest.approx(1.0)
    # a point outside the reference box contributes nothing; a dominated point neither
    assert F.hypervolume([[0.1] * 4, [1.2, 0, 0, 0], [0.5] * 4], 1.1) == pytest.approx(1.0)
    assert F.hypervolume([], 1.1) == 0.0


def test_normaliser():
    refs = {"fast": {"sp": 0.7, "tc": 0.5, "precision": 0.6, "cpu": 10.0},
            "default": {"sp": 0.7, "tc": 0.5, "precision": 0.6, "cpu": 40.0},
            "accurate": {"sp": 0.75, "tc": 0.55, "precision": 0.65, "cpu": 200.0}}
    n1, n2 = F.Normaliser(refs, 1), F.Normaliser(refs, 2)
    f = n1.F({"sp": 0.6, "tc": 0.4, "precision": 0.5, "cpu": 5.0})
    assert f[:3] == pytest.approx([1.0, 1.0, 1.0])                 # fast - 0.10 -> 1
    assert n1.F({"sp": 1, "tc": 1, "precision": 1, "cpu": 80.0}) == pytest.approx([0, 0, 0, 1])
    assert f[3] == pytest.approx(0.0)                              # 0.5 x fast -> 0
    assert n2.F({"sp": 1, "tc": 1, "precision": 1, "cpu": 400.0})[3] == pytest.approx(1.0)


def test_summarise_benchmark_weighting_and_errors():
    meta = {"a1": {"benchmark": "A"}, "a2": {"benchmark": "A"}, "b1": {"benchmark": "B"}}
    rows = {"a1": {"recall": 1.0, "tc": 1.0, "precision": 1.0, "f1": 1.0, "cpu_s": 1.0, "error": None},
            "a2": {"recall": 0.0, "tc": 0.0, "precision": 0.0, "f1": 0.0, "cpu_s": 1.0, "error": None},
            "b1": {"recall": 1.0, "tc": 0.5, "precision": 0.5, "f1": 0.5, "cpu_s": 2.0, "error": None}}
    s = F.summarise(rows, meta)
    assert s["sp"] == pytest.approx(0.75) and s["case_mean"]["sp"] == pytest.approx(2 / 3)
    assert s["cpu"] == 4.0 and s["valid"]
    rows["b1"] = {**rows["b1"], "error": "timeout", "recall": None}
    s = F.summarise(rows, meta)
    assert not s["valid"] and s["n_err"] == 1


# ------------------------------------------------------------------------------------ selection / budgets
def _e(key, sp, tc, p, cpu, valid=True):
    return {"key": key, "params": {"k": key}, "sp": sp, "tc": tc, "precision": p, "f1": sp, "cpu": cpu,
            "valid": valid, "n_err": 0}


def test_select_modes_budget_z_and_recall_rule():
    refs = {"fast": _e("rf", 0.70, 0.50, 0.60, 10), "default": _e("rd", 0.71, 0.51, 0.61, 50),
            "accurate": _e("ra", 0.74, 0.55, 0.63, 200)}
    refs["recall"] = refs["accurate"]
    cands = [_e("c_cheap", 0.72, 0.52, 0.60, 11.9),     # within fast budget (12)
             _e("c_over", 0.80, 0.60, 0.70, 12.1),      # just over the fast budget
             _e("c_def", 0.75, 0.56, 0.64, 59.0),
             _e("c_acc", 0.78, 0.57, 0.66, 230.0),
             _e("c_rec", 0.79, 0.50, 0.55, 239.0),      # higher SP, lower z
             _e("c_bad", 0.99, 0.99, 0.99, 1.0, valid=False),
             refs["fast"], refs["default"], refs["accurate"]]
    res = F.select_modes(cands, refs, "protein")
    sel = {m: (r["entry"] or {}).get("key") for m, r in res["selection"].items()}
    assert res["budgets"]["fast"] == pytest.approx(12.0) and res["budgets"]["recall"] == pytest.approx(240.0)
    assert sel["fast"] == "c_cheap"                     # c_over is over budget, invalid is ignored
    assert sel["default"] == "c_over"                   # dominates c_def, within the default budget
    for m, r in res["selection"].items():
        assert r["entry"]["cpu"] <= res["budgets"][m]
    acc = res["selection"]["accurate"]["entry"]
    rec = res["selection"]["recall"]["entry"]
    assert rec["sp"] >= acc["sp"] and res["selection"]["recall"]["need_sp"] == acc["sp"]
    assert "c_bad" not in [e["key"] for e in res["front"]]


def test_select_modes_none_within_budget():
    refs = {m: _e("r" + m, 0.7, 0.5, 0.6, 10) for m in ("fast", "default", "recall", "accurate")}
    res = F.select_modes([_e("x", 0.9, 0.9, 0.9, 100)], refs, "nucleotide")
    assert all(r["entry"] is None for r in res["selection"].values())


# ------------------------------------------------------------------------------------ fake fits
def _cases(n=60, benchmarks=("A", "B", "C")):
    rows = [{"case_id": f"{benchmarks[i % len(benchmarks)]}{i:03d}", "benchmark": benchmarks[i % len(benchmarks)],
             "n_seq": "5", "ref_len": "100"} for i in range(n)]
    meta = {r["case_id"]: {"benchmark": r["benchmark"]} for r in rows}
    return rows, meta


SMOKE = dict(pop_size=10, ref_partitions=2, stage1_gens=(3, 4), stage2_gens=(2, 3), note="test")


def _fit(tmp_path, pool, rows, meta, core, screen, test=frozenset(), name="f", fail_if=None):
    ev = F.FakeEvaluator(rows, core, fail_if=fail_if)
    f = F.Fit(pool, 0, tmp_path / name, core, screen, ev, meta, F.Settings(**SMOKE), test_ids=test)
    return f, ev


@needs_space
@pytest.mark.parametrize("pool", ["protein", "nucleotide"])
def test_fake_fit_end_to_end_resume_and_no_leak(tmp_path, pool):
    rows, meta = _cases()
    ids = sorted(meta)
    test = frozenset(ids[:12])
    core = ids[12:]
    screen = core[::2]
    f, ev = _fit(tmp_path, pool, rows, meta, core, screen, test)
    sel = f.run()
    # every evaluated case is a training case, screens are screen-only during the search
    for _, cids, _ in ev.calls:
        assert set(cids) <= set(core) and not set(cids) & test
    assert any(set(c) == set(screen) for _, c, _ in ev.calls) and any(set(c) == set(core) for _, c, _ in ev.calls)
    for m in F.MODES[pool]:
        s = sel["selection"][m]
        if s["key"] is not None:
            assert s["core"]["cpu"] <= s["budget_cpu"] + 1e-9
            ref = sel["refs_core"][F.BUDGET_REF[pool][m]]["cpu"]
            assert s["budget_cpu"] == pytest.approx(1.2 * ref)
    out = tmp_path / "f"
    for p in ("run.json", "refs.json", "selection.json", "stage1/hv.json", "stage2/hv.json", "stage1/pool.json",
              "stage1/front_core.json", "stage2/front_core.json", "stage1/evaluated.jsonl"):
        assert (out / p).exists(), p
    hv = json.loads((out / "stage1/hv.json").read_text())
    assert hv["stop"] and 3 <= hv["n_gen"] <= 4
    # resume: same selection, no re-evaluation needed beyond cache-free fake recomputation
    f2, _ = _fit(tmp_path, pool, rows, meta, core, screen, test)
    sel2 = f2.run()
    assert {m: s["key"] for m, s in sel2["selection"].items()} == {m: s["key"] for m, s in sel["selection"].items()}
    assert json.loads((out / "stage1/hv.json").read_text())["n_gen"] == hv["n_gen"]
    # the 3.7.0 presets are seeds of stage 1 / 2 and appear in the archives
    keys1 = {json.loads(l)["key"] for l in open(out / "stage1/evaluated.jsonl")}
    keys2 = {json.loads(l)["key"] for l in open(out / "stage2/evaluated.jsonl")}
    sp = F.get_space(pool)
    assert {F.pkey(p) for p in sp.stage1_seeds()} <= keys1
    assert {F.pkey(p) for p in sp.stage2_seed_ensembles()} <= keys2


@needs_space
def test_leak_assertions(tmp_path):
    rows, meta = _cases(30)
    ids = sorted(meta)
    with pytest.raises(F.LeakError):
        _fit(tmp_path, "protein", rows, meta, ids, ids[:10], test=frozenset(ids[:3]))
    ev = F.FakeEvaluator(rows, ids[5:])
    with pytest.raises(F.LeakError):
        ev.evaluate({"k": {}}, ids[:6])
    with pytest.raises(ValueError):
        _fit(tmp_path, "protein", rows, meta, ids[5:], ids[:10])  # screen outside core


@needs_space
def test_failed_individuals_are_infeasible(tmp_path):
    rows, meta = _cases(30)
    ids = sorted(meta)
    # every ensemble with consistency > 0 times out on one case
    f, ev = _fit(tmp_path, "nucleotide", rows, meta, ids, ids[::2],
                 fail_if=lambda p, cid: p["consistency"] > 0 and p["n_runs"] == 1 and cid == ids[0]
                 and p["run_gpo"][0] < 10)
    sel = f.run()
    arch = [json.loads(l) for l in open(tmp_path / "f/stage1/evaluated.jsonl")]
    bad = [e for e in arch if not e["screen"]["valid"]]
    fr = json.loads((tmp_path / "f/stage1/front_core.json").read_text())["front"]
    assert not {e["key"] for e in bad} & {e["key"] for e in fr}
    for s in sel["selection"].values():
        assert s["key"] not in {e["key"] for e in bad}


# ------------------------------------------------------------------------------------ spaces
@needs_space
def test_nucleotide_stage1_roundtrip_and_stage2_seeds():
    sp = F.get_space("nucleotide")
    for p in sp.stage1_seeds():
        assert sp.stage1_decode(sp.stage1_encode(p)) == p
    pool = sp.stage1_seeds()
    for ens in sp.stage2_seed_ensembles():
        x = F.stage2_seed_x(sp, ens, pool)
        assert F.pkey(sp.stage2_decode(x, pool)) == F.pkey(ens)
    assert F.pkey(pool[0]) == F.pkey(F._py(P370.params_370("nucleotide", "fast")))


@needs_space
def test_protein_stage2_matches_space_py_for_five_runs():
    sp = F.get_space("protein")
    pool = sp.stage1_seeds()
    for p in pool:
        assert F.pkey(sp.stage1_decode(sp.stage1_encode(p))) == F.pkey(p)
    acc = F._py(P370.params_370("protein", "accurate"))
    x = F.stage2_seed_x(sp, acc, pool)
    assert F.pkey(sp.stage2_decode(x, pool)) == F.pkey(acc)
    assert F.pkey(sp.stage2_decode(x, pool)) == F.pkey(F._py(sp.sp.stage2_decode(x, pool)))
    x3 = dict(x, n_runs=3)
    p3 = sp.stage2_decode(x3, pool)
    assert p3["n_runs"] == 3 and p3["run_gpo"] == acc["run_gpo"][:3]


# ------------------------------------------------------------------------------------ slots / evaluator
def test_slot_pool_bounds_concurrency(tmp_path):
    pools = [F.SlotPool(tmp_path / "slots", 3) for _ in range(2)]   # two "processes" share 3 slots
    cur, peak, lock = [0], [0], threading.Lock()

    def job(sp):
        with sp.slot():
            with lock:
                cur[0] += 1
                peak[0] = max(peak[0], cur[0])
            time.sleep(0.02)
            with lock:
                cur[0] -= 1

    ts = [threading.Thread(target=job, args=(pools[i % 2],)) for i in range(24)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    assert peak[0] <= 3
    with pytest.raises(ValueError):
        F.SlotPool(tmp_path / "s2", F.MAX_SLOTS + 1)   # fit modules raise MAX_SLOTS on import


def test_slot_pool_hold_others(tmp_path):
    a = F.SlotPool(tmp_path / "slots", 30)
    a.hold_others(4)
    b = F.SlotPool(tmp_path / "slots", 30)
    got = []
    for _ in range(4):
        cm = b.slot()
        got.append((cm, cm.__enter__()))
    assert sorted(i for _, i in got) == [0, 1, 2, 3]   # only the 4 free slots exist
    for cm, _ in got:
        cm.__exit__(None, None, None)
    a.release_others()


FAKE_P = {"n_runs": 1, "consistency": 0, "realign": 0, "run_refine": [0], "run_gpo": [5.0]}


def _real_cases(tmp_path, ids):
    rows = []
    for cid in ids:
        f = tmp_path / f"{cid}.fa"
        f.write_text(">a\nACGT\n>b\nACGA\n")
        rows.append({"case_id": cid, "pool": "protein", "benchmark": "T", "subset": "s", "unaligned": str(f),
                     "reference": str(f), "score_rule": "x", "n_seq": "2", "ref_len": "4",
                     "sha256_unaligned": "", "sha256_reference": ""})
    return rows


def test_real_evaluator_with_fake_worker(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "WORKER_SCRIPT", Path(__file__).with_name("fake_worker.py"))
    (tmp_path / "unused").write_text("x")
    build = E.Build.create(str(tmp_path / "unused"), sys.executable)
    ids = [f"c{i}" for i in range(25)] + ["hang_c"]
    rows = _real_cases(tmp_path, ids)
    ev = F.Evaluator(build, rows, ids[:-1] + ["hang_c"], cache_path=tmp_path / "c.sqlite",
                     slots=F.SlotPool(tmp_path / "slots", 3), threads=3, timeout_s=1.0, chunk=7,
                     scorer="fake_scorer:score_case")
    p2 = dict(FAKE_P, run_gpo=[6.0])
    specs = {F.pkey(FAKE_P): FAKE_P, F.pkey(p2): p2}
    res = ev.evaluate(specs, ids)
    for k in specs:
        assert set(res[k]) == set(ids)
        assert res[k]["hang_c"]["error"] == "timeout"
        assert all(res[k][c]["error"] is None for c in ids[:-1])
    assert C.Cache(tmp_path / "c.sqlite").count() == 2 * len(ids)
    with pytest.raises(F.LeakError):
        F.Evaluator(build, rows, ids[:3], cache_path=tmp_path / "c.sqlite").evaluate(specs, ids[:4])


def test_real_evaluator_stop_on_error_and_cpu_cap(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "WORKER_SCRIPT", Path(__file__).with_name("fake_worker.py"))
    (tmp_path / "unused").write_text("x")
    build = E.Build.create(str(tmp_path / "unused"), sys.executable)
    ids = ["crash_c"] + [f"c{i}" for i in range(10)]
    cpu_ids = [f"cpu_{i}" for i in range(10)]
    rows = _real_cases(tmp_path, ids + cpu_ids)
    mk = lambda: F.Evaluator(build, rows, ids + cpu_ids, cache_path=tmp_path / "c.sqlite", threads=1, chunk=100,
                             timeout_s=5.0, scorer="fake_scorer:score_case")
    for r in rows:                       # crash first (most expensive), then the rest
        r["n_seq"] = "9" if r["case_id"] == "crash_c" else "2"
    res = mk().evaluate({F.pkey(FAKE_P): FAKE_P}, ids, stop_on_error=True)[F.pkey(FAKE_P)]
    assert res["crash_c"]["error"].startswith("worker exit")
    assert all(res[c]["error"] == "skipped: error" for c in ids[1:])
    assert C.Cache(tmp_path / "c.sqlite").count() == 1          # skipped rows are never cached
    res = mk().evaluate({F.pkey(FAKE_P): FAKE_P}, cpu_ids, cpu_cap=0.7)[F.pkey(FAKE_P)]
    n_done = sum(1 for c in cpu_ids if res[c]["error"] is None)
    assert 3 <= n_done <= 4 and sum(1 for c in cpu_ids if res[c]["error"] == "skipped: cpu_cap") == 10 - n_done
    # a second call counts the cached CPU and stops at once
    res2 = mk().evaluate({F.pkey(FAKE_P): FAKE_P}, cpu_ids, cpu_cap=0.7)[F.pkey(FAKE_P)]
    assert sum(1 for c in cpu_ids if res2[c]["error"] is None) == n_done



def test_short_cached_timeouts_are_rerun_under_the_longer_test_limit(tmp_path, monkeypatch):
    """Review 2026-10-05: a "timeout" cached under the 600 s fitting limit must not be reused by an out-of-fold
    test (1800 s): it is re-run and replaced; a timeout under the full limit is kept."""
    monkeypatch.setattr(E, "WORKER_SCRIPT", Path(__file__).with_name("fake_worker.py"))
    (tmp_path / "unused").write_text("x")
    build = E.Build.create(str(tmp_path / "unused"), sys.executable)
    ids = ["hang_c", "c0"]
    rows = _real_cases(tmp_path, ids)
    k = F.pkey(FAKE_P)
    fitting = F.Evaluator(build, rows, ids, cache_path=tmp_path / "c.sqlite", threads=1, timeout_s=0.5,
                          scorer="fake_scorer:score_case")
    assert fitting.evaluate({k: FAKE_P}, ids)[k]["hang_c"]["error"] == "timeout"
    # rename the hanging input so that the re-run succeeds: proves the cached row was not reused
    hang = next(r for r in rows if r["case_id"] == "hang_c")
    newp = tmp_path / "now_ok.fa"
    Path(hang["unaligned"]).rename(newp)
    hang["unaligned"] = hang["reference"] = str(newp)
    plain = F.Evaluator(build, rows, ids, cache_path=tmp_path / "c.sqlite", threads=1, timeout_s=2.0,
                        scorer="fake_scorer:score_case")
    assert plain.evaluate({k: FAKE_P}, ids)[k]["hang_c"]["error"] == "timeout"     # default: cached row reused
    test = F.Evaluator(build, rows, ids, cache_path=tmp_path / "c.sqlite", threads=1, timeout_s=2.0,
                       scorer="fake_scorer:score_case", rerun_short_timeouts=True)
    assert test.evaluate({k: FAKE_P}, ids)[k]["hang_c"]["error"] is None

@needs_space
def test_fit_caps_search_candidates_not_references(tmp_path):
    rows, meta = _cases(30)
    ids = sorted(meta)
    f, ev = _fit(tmp_path, "protein", rows, meta, ids, ids[::2])
    f.run()
    acc_key = F.pkey(F._py(P370.params_370("protein", "accurate")))
    caps = [cap for keys, _, cap in ev.calls if acc_key in keys and cap is None]
    assert caps, "references are evaluated without a cap"
    capped = [cap for _, _, cap in ev.calls if cap is not None]
    assert capped and all(c == pytest.approx(3.0 * f.refs_screen["accurate"]["cpu"]) or
                          c == pytest.approx(3.0 * f.refs_core["accurate"]["cpu"]) for c in capped)


# ------------------------------------------------------------------------------------ bootstrap / pooling
def test_bca_matches_normal_theory_and_strata():
    rng = np.random.default_rng(0)
    x = rng.normal(0.5, 0.1, 2000)
    est, lo, hi = PR.bca(x, B=4000)
    se = x.std(ddof=1) / np.sqrt(len(x))
    assert est[0] == pytest.approx(x.mean())
    assert lo[0] == pytest.approx(x.mean() - 1.96 * se, abs=0.3 * se)
    assert hi[0] == pytest.approx(x.mean() + 1.96 * se, abs=0.3 * se)
    strata = ["a"] * 1900 + ["b"] * 100
    y = np.r_[np.zeros(1900), np.ones(100)]
    est, lo, hi = PR.bca(y, strata, B=2000)
    assert est[0] == pytest.approx(0.5) and lo[0] == hi[0] == pytest.approx(0.5)   # no within-stratum spread
    e2 = PR.bca(np.c_[x, 2 * x], B=1000)
    assert e2[0][1] == pytest.approx(2 * e2[0][0]) and e2[1][1] == pytest.approx(2 * e2[1][0])


def test_decision_rule():
    assert D.rule(0.005, 0.0) and D.rule(-0.005, 0.005) and D.rule(0.01, 0.01)
    assert not D.rule(0.004, 0.004) and not D.rule(0.02, -0.0051) and not D.rule(-0.006, 0.02)


def _write_oof(root, pool, fold, cases, d_sp=0.01, sel_key="kn"):
    import csv
    from recipe import nested as N
    d = root / pool / f"fold{fold}"
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    rng = np.random.default_rng(fold)
    for m in F.MODES[pool]:
        qs = {cid: 0.7 + 0.1 * rng.random() for cid, _, _ in cases}
        for arm in ("new", "v370"):
            for cid, b, h in cases:
                q = qs[cid]
                rows.append({"pool": pool, "fold": fold, "scope": "oof", "mode": m, "arm": arm,
                             "key": sel_key if arm == "new" else "k370", "is_370": int(arm == "v370"),
                             "fallback_370": 0, "smoke": 0, "case_id": cid, "benchmark": b, "difficulty": "x",
                             "homolog_in_train": h, "recall": q + (d_sp if arm == "new" else 0), "tc": q,
                             "precision": q, "f1": q, "width_ratio": 1.0, "cpu_s": 1.0, "wall_s": 1.0, "error": ""})
    N.write_rows(d / "oof_rows.tsv", rows)
    sel = {"fold": fold, "refs_core": {m: {"cpu": 1.0} for m in ("fast", "default", "recall", "accurate")},
           "selection": {m: {"key": sel_key, "is_370": False, "fallback_370": False, "budget_cpu": 1.2,
                             "core": {"sp": 0.8, "tc": 0.7, "precision": 0.7, "cpu": 1.0}}
                         for m in F.MODES[pool]}}
    (d / "selection.json").write_text(json.dumps(sel))


def test_pool_results_tables(tmp_path):
    pool = "nucleotide"
    for k in range(5):
        cases = [(f"f{k}_{i}", "A" if i % 3 else "B", i % 2) for i in range(30)]
        _write_oof(tmp_path, pool, k, cases)
    res = PR.pool_results(pool, tmp_path, log=lambda *a: None)
    di = res["decision_inputs"]["default"]
    assert di["benchmark"]["d_sp"] == pytest.approx(0.01) and di["case"]["d_tc"] == pytest.approx(0.0)
    assert di["n"] == 150
    rd = tmp_path / pool / "results"
    for f in ("oof_summary.tsv", "oof_bw.tsv", "oof_diff.tsv", "per_fold.tsv", "selected.tsv", "results.json"):
        assert (rd / f).exists()
    # optimism: final-fit selection in-sample (all cases) minus pooled OOF, paired by case
    import csv
    from recipe import nested as N
    ins = []
    for k in range(5):
        for r in csv.DictReader(open(tmp_path / pool / f"fold{k}" / "oof_rows.tsv"), delimiter="\t"):
            if r["arm"] == "new":
                ins.append({**r, "scope": "insample", "recall": float(r["recall"]) + 0.02})
    (tmp_path / pool / "final").mkdir()
    N.write_rows(tmp_path / pool / "final" / "insample_rows.tsv", ins)
    res = PR.pool_results(pool, tmp_path, log=lambda *a: None)
    opt = [o for o in res["optimism"] if o["mode"] == "fast" and o["weighting"] == "case"][0]
    assert opt["opt_sp"] == pytest.approx(0.02) and opt["opt_tc"] == pytest.approx(0.0) and opt["n"] == 150
    with pytest.raises(SystemExit):     # partial pooling is refused by default
        (tmp_path / pool / "fold4" / "oof_rows.tsv").unlink()
        PR.pool_results(pool, tmp_path, log=lambda *a: None)


def test_budget_check_uses_final_core_retime():
    pool = "protein"
    ref = {m: F.pkey(F._py(P370.params_370(pool, m))) for m in ("fast", "default", "accurate")}
    sel = {"selection": {m: {"key": f"new_{m}"} for m in F.MODES[pool]}}
    rows = []
    cpu = {ref["fast"]: 10, ref["default"]: 50, ref["accurate"]: 200, "new_fast": 12.0, "new_default": 61.0,
           "new_recall": 239.0, "new_accurate": 100.0}
    for k, c in cpu.items():
        rows.append({"scope": "final_core", "key": k, "cpu_s": str(c), "error": ""})
        rows.append({"scope": "oof", "key": k, "cpu_s": "9999", "error": ""})        # ignored
    b = D.budget_check(pool, rows, sel)
    assert b["fast"]["within"] and not b["default"]["within"] and b["recall"]["within"]
    assert b["recall"]["ref_mode"] == "accurate" and b["accurate"]["ratio"] == pytest.approx(0.5)
