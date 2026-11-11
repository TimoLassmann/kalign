import multiprocessing as mp
import os
import stat
from pathlib import Path

import pytest

from recipe import cache as C
from recipe import evaluate as E

B = "b" * 64
V = "fake-v1"


def I(cid, su="su", sr="sr", v=V):  # noqa: E743
    return C.Ident(cid, su, sr, v)


# ---------------------------------------------------------------- param_key
def test_param_key_stable_golden():
    # golden values: a change here means every cached entry is orphaned
    assert C.canonical_json("fast") == '{"preset":"fast"}'
    assert C.param_key("fast") == C.param_key({"preset": "fast"})


def test_param_key_golden_hash():
    import hashlib
    assert C.param_key({"preset": "fast"}) == hashlib.sha256(b'{"preset":"fast"}').hexdigest()


def test_param_key_order_and_float_rounding():
    a = {"run_gpo": [5.5, 6.0], "realign": 1, "seq_weights": 0.1 + 0.2}
    b = {"seq_weights": 0.3, "realign": 1.0, "run_gpo": [5.5000000000004, 6]}
    assert C.param_key(a) == C.param_key(b)
    assert C.param_key({"x": 0.3}) != C.param_key({"x": 0.3000001})
    assert C.param_key({"x": -0.0}) == C.param_key({"x": 0})


def test_param_key_numpy_and_tuple():
    np = pytest.importorskip("numpy")
    assert C.param_key({"a": np.float64(1.5), "b": (1, 2), "c": np.array([1, 2])}) == \
        C.param_key({"a": 1.5, "b": [1, 2], "c": [1, 2]})


def test_param_key_rejects_nan():
    with pytest.raises(ValueError):
        C.param_key({"x": float("nan")})


# ---------------------------------------------------------------- put/get
def test_put_get_never_overwrite(tmp_path):
    c = C.Cache(tmp_path / "c.sqlite")
    k = C.param_key("fast")
    assert c.get(B, k, I("case1")) is None
    m = {"recall": 0.9, "precision": 0.8, "f1": 0.85, "tc": 0.5, "width_ratio": 1.1, "n_seq": 12, "foo": 3}
    assert c.put(B, k, I("case1"), m, cpu_s=1.5, wall_s=2.0, aln_sha256="ab", spec="fast")
    assert not c.put(B, k, I("case1"), {"recall": 0.1}, cpu_s=9)  # never overwritten
    r = c.get(B, k, I("case1"))
    assert r["recall"] == 0.9 and r["n_seq"] == 12 and r["cpu_s"] == 1.5 and r["extra"] == {"foo": 3}
    assert r["error"] is None and r["host"] and r["ts"]
    assert c.spec_json(k) == '{"preset":"fast"}'
    assert c.get("c" * 64, k, I("case1")) is None  # other build -> miss
    assert c.get(B, k, I("case1", su="other")) is None  # changed input file -> miss
    assert c.get(B, k, I("case1", sr="other")) is None  # changed reference file -> miss
    assert c.get(B, k, I("case1", v="fake-v2")) is None  # new scorer version -> miss
    got = c.get_many(B, k, [I("case1"), I("case2")])
    assert c.get_many(B, k, [I("case1", v="fake-v2")]) == {}
    assert set(got) == {"case1"}


def test_error_entry_replace_only_if_error(tmp_path):
    c = C.Cache(tmp_path / "c.sqlite")
    k = C.param_key("fast")
    c.put(B, k, I("x"), error="boom")
    assert not c.put(B, k, I("x"), {"recall": 1.0})
    assert c.put(B, k, I("x"), {"recall": 1.0}, replace_error=True)
    assert not c.put(B, k, I("x"), {"recall": 0.0}, replace_error=True)  # ok entry is kept
    assert c.get(B, k, I("x"))["recall"] == 1.0
    # a new scorer version / file hash is a separate entry, the old one is untouched
    assert c.put(B, k, I("x", v="fake-v2"), {"recall": 0.2})
    assert c.get(B, k, I("x"))["recall"] == 1.0 and c.get(B, k, I("x", v="fake-v2"))["recall"] == 0.2


def test_get_many_chunking(tmp_path):
    c = C.Cache(tmp_path / "c.sqlite")
    k = C.param_key("fast")
    for i in range(1200):
        c.put(B, k, I(f"c{i}"), {"recall": i})
    assert len(c.get_many(B, k, [I(f"c{i}") for i in range(1300)])) == 1200
    assert c.count(B, k) == 1200


# ---------------------------------------------------------------- concurrency
def _writer(args):
    path, wid, n = args
    c = C.Cache(path)
    k = C.param_key("fast")
    wrote = 0
    for i in range(n):
        # every writer inserts its own cases plus a shared set (duplicates must be ignored)
        wrote += c.put(B, k, I(f"w{wid}_{i}"), {"recall": float(i)}, cpu_s=0.1)
        wrote += c.put(B, k, I(f"shared_{i}"), {"recall": float(wid)}, cpu_s=0.1)
    return wrote


def test_concurrent_writers(tmp_path):
    path, workers, n = str(tmp_path / "c.sqlite"), 12, 60
    C.Cache(path)
    with mp.get_context("spawn").Pool(workers) as p:
        wrote = p.map(_writer, [(path, w, n) for w in range(workers)])
    c = C.Cache(path)
    assert c.count(B) == workers * n + n
    assert sum(wrote) == workers * n + n  # each shared entry written exactly once
    assert c._conn().execute("PRAGMA journal_mode").fetchone()[0] == "wal"
    assert c._conn().execute("PRAGMA integrity_check").fetchone()[0] == "ok"


# ---------------------------------------------------------------- evaluate (fake aligner)
@pytest.fixture
def fake_build(tmp_path):
    cli = tmp_path / "fakekalign"
    cli.write_text('#!/bin/sh\n[ "$2" = "boom" ] && { echo oops >&2; exit 3; }\ncp "$6" "$8"\n')
    cli.chmod(cli.stat().st_mode | stat.S_IEXEC)
    return E.Build.create(str(cli), "python")


def _cases(tmp_path, ids):
    rows = []
    for cid in ids:
        f = tmp_path / f"{cid}.fa"
        f.write_text(">a\nACGT\n>b\nACGA\n")
        rows.append({"case_id": cid, "pool": "protein", "benchmark": "T", "subset": "s",
                     "unaligned": str(f), "reference": str(f), "score_rule": "x", "n_seq": "2",
                     "ref_len": "4", "sha256_unaligned": "", "sha256_reference": ""})
    return rows


@pytest.mark.parametrize("workers", [1, 3])
def test_evaluate_cli_cache_roundtrip(tmp_path, fake_build, workers):
    cases = _cases(tmp_path, ["c1", "c2", "c3", "c4bad"])
    kw = dict(cache_path=tmp_path / "c.sqlite", scorer="fake_scorer:score_case")
    r1 = E.evaluate(cases, "fast", fake_build, workers, **kw)
    assert [r["cached"] for r in r1] == [0, 0, 0, 0]
    assert r1[0]["recall"] == 0.5 and r1[0]["cpu_s"] >= 0 and r1[0]["aln_sha256"]
    assert r1[3]["error"].startswith("score ValueError")  # scoring failure recorded
    r2 = E.evaluate(cases, "fast", fake_build, workers, **kw)
    assert [r["cached"] for r in r2] == [1, 1, 1, 1]
    assert [r["aln_sha256"] for r in r2] == [r["aln_sha256"] for r in r1]
    # temp alignments are gone
    assert not list(tmp_path.glob("ev_*"))


def test_evaluate_aligner_failure_and_workers_cap(tmp_path, fake_build):
    cases = _cases(tmp_path, ["c1"])
    kw = dict(cache_path=tmp_path / "c.sqlite", scorer="fake_scorer:score_case")
    # the fake CLI fails for mode "boom"; real presets are validated, so call run_cli directly
    wall, cpu, err = E.run_cli(fake_build.cli, "boom", cases[0]["unaligned"], tmp_path / "o.fa")
    assert err.startswith("exit 3") and "oops" in err
    with pytest.raises(ValueError):
        E.evaluate(cases, "fast", fake_build, 31, **kw)
    with pytest.raises(ValueError):
        E.evaluate(cases, "nonsense", fake_build, 1, **kw)


def test_evaluate_keep_alignment(tmp_path, fake_build):
    cases = _cases(tmp_path, ["c1"])
    E.evaluate(cases, "fast", fake_build, 1, cache_path=tmp_path / "c.sqlite",
               scorer="fake_scorer:score_case", keep_dir=tmp_path / "keep")
    assert len(list((tmp_path / "keep").rglob("*.aln.fa"))) == 1


# ---------------------------------------------------------------- key hygiene in evaluate
def test_evaluate_changed_file_or_scorer_is_a_miss(tmp_path, fake_build, monkeypatch):
    import fake_scorer
    cases = _cases(tmp_path, ["c1", "c2"])
    kw = dict(cache_path=tmp_path / "c.sqlite", scorer="fake_scorer:score_case")
    E.evaluate(cases, "fast", fake_build, 1, **kw)
    assert [r["cached"] for r in E.evaluate(cases, "fast", fake_build, 1, **kw)] == [1, 1]
    changed = [{**cases[0], "sha256_reference": "newhash"}, cases[1]]  # same id, other reference file
    assert [r["cached"] for r in E.evaluate(changed, "fast", fake_build, 1, **kw)] == [0, 1]
    monkeypatch.setattr(fake_scorer, "SCORER_VERSION", "fake-v2")
    assert [r["cached"] for r in E.evaluate(cases, "fast", fake_build, 1, **kw)] == [0, 0]


def test_scorer_version_required(tmp_path):
    with pytest.raises(ValueError):
        E.scorer_version("os.path:join", {"case_id": "x"})
    assert E.scorer_version("recipe.scoring:score_case", {"pool": "protein"}) == "protein-v1"
    assert E.scorer_version("recipe.scoring:score_case", {"pool": "nucleotide"}) == "nucleotide-v2"


# ---------------------------------------------------------------- schema-1 upgrade + migration
def _make_v1(path, rows):
    import sqlite3
    c = sqlite3.connect(path)
    c.executescript("""
    CREATE TABLE results (build_id TEXT NOT NULL, param_key TEXT NOT NULL, case_id TEXT NOT NULL,
      recall REAL, precision REAL, f1 REAL, tc REAL, width_ratio REAL, n_seq INTEGER, extra TEXT,
      cpu_s REAL, wall_s REAL, aln_sha256 TEXT, host TEXT, ts TEXT, error TEXT, harness TEXT,
      PRIMARY KEY (build_id, param_key, case_id)) WITHOUT ROWID;
    CREATE TABLE specs (param_key TEXT PRIMARY KEY, spec_json TEXT NOT NULL) WITHOUT ROWID;
    CREATE TABLE meta (k TEXT PRIMARY KEY, v TEXT) WITHOUT ROWID;
    INSERT INTO meta VALUES ('schema','1');""")
    for cid in rows:
        c.execute("INSERT INTO results VALUES (?,?,?,0.5,0.5,0.5,0.5,1.0,3,NULL,1.0,1.0,'aa','h','t',NULL,'old')",
                  (B, "k1", cid))
    c.commit()
    c.close()


def _manifest(tmp_path, ids, pool="protein"):
    rows = []
    for cid in ids:
        u, r = tmp_path / f"{cid}.u", tmp_path / f"{cid}.r"
        u.write_text("u" + cid)
        r.write_text("r" + cid)
        rows.append((cid, pool, str(u), str(r), C.sha256_file(u), C.sha256_file(r)))
    p = tmp_path / "m.tsv"
    p.write_text("case_id\tpool\tunaligned\treference\tsha256_unaligned\tsha256_reference\n"
                 + "".join("\t".join(r) + "\n" for r in rows))
    return p, rows


def test_v1_upgrade_hides_old_rows_and_migration_fills_hashes(tmp_path):
    db = str(tmp_path / "old.sqlite")
    _make_v1(db, ["a", "b", "gone"])
    c = C.Cache(db)  # upgrade on open
    assert c.count() == 0  # nothing from schema 1 is readable
    assert c._conn().execute("SELECT COUNT(*) FROM results_v1").fetchone()[0] == 3
    man, rows = _manifest(tmp_path, ["a", "b"])
    res = C.migrate_v1(c, [man], {"protein": "protein-v1"})
    assert res["legacy_rows"] == 3 and res["migrated"] == 2 and res["skipped_not_in_manifest"] == 1
    a = next(r for r in rows if r[0] == "a")
    got = c.get(B, "k1", C.Ident("a", a[4], a[5], "protein-v1"))
    assert got and got["recall"] == 0.5 and got["harness"] == "old"
    assert c.get(B, "k1", C.Ident("a", a[4], a[5], "protein-v2")) is None
    assert c.get(B, "k1", C.Ident("gone", "x", "y", "protein-v1")) is None
    assert C.migrate_v1(c, [man], {"protein": "protein-v1"})["migrated"] == 0  # idempotent
    c.close()
    assert C.Cache(db).count() == 2  # re-open: no second upgrade


def test_migration_skips_files_that_changed(tmp_path):
    db = str(tmp_path / "old.sqlite")
    _make_v1(db, ["a", "b"])
    c = C.Cache(db)
    man, rows = _manifest(tmp_path, ["a", "b"])
    (tmp_path / "a.r").write_text("tampered")  # file no longer matches its manifest hash
    res = C.migrate_v1(c, [man], {"protein": "protein-v1"})
    assert res["migrated"] == 1 and res["skipped_file_mismatch_cases"] == 1


# ---------------------------------------------------------------- timeouts
def test_cli_timeout_kills_and_records(tmp_path):
    import time
    cli = tmp_path / "sleepy"
    cli.write_text("#!/bin/sh\nsleep 60\n")
    cli.chmod(0o755)
    t0 = time.time()
    wall, cpu, err = E.run_cli(str(cli), "fast", tmp_path / "i", tmp_path / "o", timeout_s=0.5)
    assert err == "timeout" and 0.4 < wall < 5 and time.time() - t0 < 10
    assert cpu is not None


def test_cli_timeout_kills_process_group(tmp_path):
    import time
    marker = tmp_path / "child_alive"
    cli = tmp_path / "forky"
    cli.write_text(f"#!/bin/sh\n(sleep 2; touch {marker}) &\nsleep 60\n")
    cli.chmod(0o755)
    E.run_cli(str(cli), "fast", "i", "o", timeout_s=0.3)
    time.sleep(2.5)
    assert not marker.exists()  # grandchild died with the group


def test_evaluate_timeout_stored_as_error(tmp_path):
    cli = tmp_path / "sleepy"
    cli.write_text("#!/bin/sh\nsleep 60\n")
    cli.chmod(0o755)
    b = E.Build.create(str(cli), "python")
    cases = _cases(tmp_path, ["c1"])
    kw = dict(cache_path=tmp_path / "c.sqlite", scorer="fake_scorer:score_case", timeout_s=0.5)
    r = E.evaluate(cases, "fast", b, 1, **kw)
    assert r[0]["error"] == "timeout" and r[0]["recall"] is None
    assert E.evaluate(cases, "fast", b, 1, **kw)[0]["cached"] == 1


@pytest.fixture
def fake_params_build(tmp_path, monkeypatch):
    import sys
    monkeypatch.setattr(E, "WORKER_SCRIPT", Path(__file__).with_name("fake_worker.py"))
    (tmp_path / "unused").write_text("x")
    return E.Build.create(str(tmp_path / "unused"), sys.executable)


PARAMS = {"run_gpo": [5.0], "x": 1}


def test_params_oneshot_timeout(tmp_path, fake_params_build):
    cases = _cases(tmp_path, ["c1", "hang_c2"])
    kw = dict(cache_path=tmp_path / "c.sqlite", scorer="fake_scorer:score_case", timeout_s=1.0)
    r = E.evaluate(cases, PARAMS, fake_params_build, 2, **kw)
    assert r[0]["error"] is None and r[1]["error"] == "timeout"


# ---------------------------------------------------------------- batched worker
def test_batch_worker_matches_per_case_and_survives_timeout_crash(tmp_path, fake_params_build):
    names = [f"c{i}" for i in range(12)] + ["hang_c", "crash_c", "cpu_c", "c_after"]
    cases = _cases(tmp_path, names)
    kwb = dict(cache_path=tmp_path / "b.sqlite", scorer="fake_scorer:score_case", timeout_s=1.5)
    kwp = dict(cache_path=tmp_path / "p.sqlite", scorer="fake_scorer:score_case", timeout_s=1.5)
    rb = E.evaluate(cases, PARAMS, fake_params_build, 3, batch=True, **kwb)
    rp = E.evaluate(cases, PARAMS, fake_params_build, 3, batch=False, **kwp)
    assert [r["case_id"] for r in rb] == names
    for b, p in zip(rb, rp):
        assert b["error"] == p["error"] or (b["error"].startswith("worker exit") and p["error"].startswith("worker exit"))
        assert b["aln_sha256"] == p["aln_sha256"] and b["recall"] == p["recall"]
    by = {r["case_id"]: r for r in rb}
    assert by["hang_c"]["error"] == "timeout"
    assert by["crash_c"]["error"].startswith("worker exit")
    assert by["cpu_c"]["error"] is None and by["cpu_c"]["cpu_s"] >= 0.25  # per-case CPU delta
    assert by["c_after"]["error"] is None and by["c0"]["cpu_s"] < 0.2
    assert all(r["error"] is None for r in rb if r["case_id"].startswith("c") and r["case_id"][1:].isdigit())
    assert [r["cached"] for r in E.evaluate(cases, PARAMS, fake_params_build, 3, batch=True, **kwb)] == [1] * len(names)


def test_param_worker_recycles_and_cpu_is_per_case(tmp_path, fake_params_build):
    f = tmp_path / "cpu1.fa"
    f.write_text(">a\nAC\n")
    g = tmp_path / "g.fa"
    g.write_text(">a\nAC\n")
    w = E.ParamWorker(fake_params_build.python, PARAMS, recycle=2)
    try:
        r1 = w.align(f, tmp_path / "o1", 10)
        r2 = w.align(g, tmp_path / "o2", 10)
        assert r1[1] >= 0.25 and r2[1] < 0.2 and r1[2] is None
        assert w.proc is None  # recycled after 2 cases
        r3 = w.align(g, tmp_path / "o3", 10)
        assert r3[2] is None and w.n == 1
    finally:
        w.close()
