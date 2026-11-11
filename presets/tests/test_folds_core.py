import random
from collections import Counter

from recipe import core, folds
from recipe import grouping as g


def synth(pool, n_groups=120, seed=3):
    rng = random.Random(seed)
    benches = ["bali3", "prefab4", "ox", "sabre", "balifam"] if pool == "protein" else ["bralibase", "mdsa", "rfam"]
    cases = []
    for gi in range(n_groups):
        size = rng.choice([1, 1, 1, 2, 3, 8, 30]) if gi else 60
        for k in range(size):
            b = rng.choice(benches)
            pid = rng.random()
            cases.append(dict(case_id=f"{b}.{gi}_{k}", benchmark=b, subset="x", group_id=f"G{gi:04d}", pid=pid,
                              bin=folds.difficulty_bin(pool, pid), n_seq=5, ref_len=100, row={}))
    return cases


def fold_table(cases, pool):
    gf = folds.assign_folds(cases, pool)
    return {c["case_id"]: dict(fold=gf[c["group_id"]]) for c in cases}, gf


def test_waterfill():
    a = core.waterfill(10, {"a": 100, "b": 100}, {"a": 1, "b": 1})
    assert a == {"a": 5, "b": 5}
    a = core.waterfill(10, {"a": 2, "b": 100}, {"a": 1, "b": 1})
    assert a == {"a": 2, "b": 8}
    a = core.waterfill(10, {"a": 100, "b": 100}, {"a": 1, "b": 0.5})
    assert sum(a.values()) == 10 and a["a"] == 7
    assert sum(core.waterfill(50, {"a": 3, "b": 4}, {"a": 1, "b": 1}).values()) == 7


def test_difficulty_bin():
    assert [folds.difficulty_bin("protein", x) for x in (0.0, 0.19, 0.2, 0.35, 0.4, 0.9)] == [0, 0, 1, 2, 3, 3]
    assert folds.difficulty_bin("protein", float("nan")) == 0


def test_groups_never_split_and_balanced():
    cases = synth("protein")
    ft, gf = fold_table(cases, "protein")
    by_group = {}
    for c in cases:
        by_group.setdefault(c["group_id"], set()).add(ft[c["case_id"]]["fold"])
    assert all(len(v) == 1 for v in by_group.values())
    sizes = Counter(ft[c["case_id"]]["fold"] for c in cases)
    assert len(sizes) == 5 and max(sizes.values()) - min(sizes.values()) <= 0.15 * len(cases) / 5 + 30
    assert folds.assign_folds(cases, "protein") == gf            # deterministic


def test_family_groups_one_per_fold():
    cases = []
    for fam in range(5):
        cases += [dict(case_id=f"bralibase.f{fam}_{i}", benchmark="bralibase", subset="s", group_id=f"F{fam}", pid=0.6,
                       bin=1, n_seq=5, ref_len=100, row={}) for i in range(40)]
    cases += [dict(case_id=f"rfam.{i}", benchmark="rfam", subset="s", group_id=f"R{i}", pid=0.8, bin=2, n_seq=5,
                   ref_len=100, row={}) for i in range(100)]
    gf = folds.assign_folds(cases, "nucleotide")
    assert sorted(gf[f"F{f}"] for f in range(5)) == [0, 1, 2, 3, 4]


def test_protein_core_rules():
    cases = synth("protein", n_groups=400, seed=5)
    ft, _ = fold_table(cases, "protein")
    cpu = {c["case_id"]: (50.0 if i % 17 == 0 else 1.0) for i, c in enumerate(cases)}
    res = core.build_core("protein", 2, data=cases, cpu=cpu, folds=ft)
    ids = set(res["core"])
    assert ids and all(ft[i]["fold"] != 2 for i in ids)
    assert all(cpu[i] <= 10 for i in ids)
    assert set(res["screen"]) <= ids and not (set(res["guard"]) & set(res["screen"]))
    assert len(res["core"]) <= 840 and len(res["screen"]) <= 300
    grp = Counter(res["cases"][i]["group_id"] for i in ids)
    assert max(grp.values()) <= 5
    assert all(res["cases"][i]["bin"] == 3 for i in res["guard"])
    # determinism and the "all" variant
    assert core.build_core("protein", 2, data=cases, cpu=cpu, folds=ft)["core"] == res["core"]
    allres = core.build_core("protein", "all", data=cases, cpu=cpu, folds=ft)
    assert len(allres["core"]) >= len(res["core"])


def test_nucleotide_core_rules():
    cases = synth("nucleotide", n_groups=300, seed=9)
    ft, _ = fold_table(cases, "nucleotide")
    res = core.build_core("nucleotide", 0, data=cases, folds=ft)
    assert all(ft[i]["fold"] != 0 for i in res["core"])
    assert max(Counter(res["cases"][i]["group_id"] for i in res["core"]).values()) <= 25
    sb = Counter(res["cases"][i]["benchmark"] for i in res["screen"])
    assert len(res["screen"]) <= 300 and max(sb.values()) - min(sb.values()) <= 2


def test_clan_links_and_extra_merge():
    import numpy as np
    rows = [dict(benchmark="rfam", family=f) for f in ("RF1", "RF2", "RF3", "RF4")] + [dict(benchmark="mdsa", family="")]
    pairs = g.clan_pairs(rows, {"RF1": "CL1", "RF2": "CL1", "RF3": "CL2"})
    roots, merged = g.merge_extra(np.arange(5), pairs, None)
    assert roots.tolist()[0] == roots.tolist()[1] and len(set(roots.tolist())) == 4 and len(merged) == 1


def test_homolog_flags(tmp_path):
    from test_grouping import make_rows, mutate, rnd
    rng = random.Random(7)
    a, b, c = rnd(120, rng), rnd(130, rng), rnd(140, rng)
    rows = make_rows(tmp_path, {"A": [a], "B": [mutate(a, 0.3, rng)], "C": [b], "D": [c]})
    flags, _ = g.homolog_flags(rows, "protein", [0, 1, 1, 1], threads=2, log=lambda *_: None)
    assert flags.tolist() == [True, True, False, False]     # A (fold 0) ~ B (fold 1); C, D unrelated
    flags, _ = g.homolog_flags(rows, "protein", [0, 0, 1, 1], threads=2, log=lambda *_: None)
    assert not flags.any()                                   # the related pair is in one fold
