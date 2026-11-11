import os
import stat
import sys
import textwrap
from pathlib import Path

import pytest

from recipe import externals as X

SCORER = "fake_scorer:score_case"


@pytest.fixture(autouse=True)
def _path():
    sys.path.insert(0, str(Path(__file__).parent))
    yield
    sys.path.remove(str(Path(__file__).parent))


def _script(path: Path, body: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n" + textwrap.dedent(body))
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


@pytest.fixture
def tdir(tmp_path):
    t = tmp_path / "tools"
    # fake muscle: -align|-super5 IN -output OUT ; writes a lower-case wrapped alignment, logs the mode
    _script(t / "bin" / "muscle", """\
        if [ "$1" = "-version" ]; then echo "muscle 5.3.fake"; exit 0; fi
        echo "$1" >> "$MODELOG"
        printf ">a\\nac.gt\\n*\\n>b\\nacgt-\\n" > "$4"
        """)
    _script(t / "bin" / "clustalo", """\
        if [ "$1" = "--version" ]; then echo 1.2.4; exit 0; fi
        if grep -q fail "$2"; then echo boom >&2; exit 3; fi
        cp "$2" "$4"
        """)
    (t / "mafft" / "bin").mkdir(parents=True)
    _script(t / "mafft" / "bin" / "mafft", """\
        if [ "$1" = "--version" ]; then echo "v7.fake (2024)" >&2; exit 0; fi
        eval last=\\${$#}; cat "$last"
        """)
    (t / "mafft" / "libexec" / "mafft").mkdir(parents=True)
    (t / "mafft" / "libexec" / "mafft" / "tbfast").write_text("x")
    return t


def _case(tmp_path, cid="c1", n_seq=2, body=">a\nacgt\n>b\nacg\n"):
    inp = tmp_path / f"{cid.replace('/', '_')}.fa"
    inp.write_text(body)
    return {"case_id": cid, "pool": "protein", "benchmark": "t", "subset": "s", "unaligned": str(inp),
            "n_seq": str(n_seq)}


def test_normalise(tmp_path):
    p = tmp_path / "a.fa"
    p.write_text(">x y\nac.g\nt*\n>z\nAC-GT-\n")
    assert X.normalise_alignment(p) == 2
    assert p.read_text() == ">x y\nAC-GT-\n>z\nAC-GT-\n"


def test_run_process_timeout_and_error(tmp_path):
    w, cpu, err = X.run_process(["sh", "-c", "sleep 30"], None, os.environ, 0.5)
    assert w < 10 and err.startswith("timeout")
    w, cpu, err = X.run_process(["sh", "-c", "echo bad >&2; exit 2"], None, os.environ, 5)
    assert err.startswith("exit 2") and "bad" in err
    w, cpu, err = X.run_process(["sh", "-c", "i=0; while [ $i -lt 200000 ]; do i=$((i+1)); done"], None,
                                os.environ, 60)
    assert err is None and cpu > 0.01  # CPU of the child is measured


def test_versions_and_command_policy(tdir):
    assert "5.3.fake" in X.tool_version("muscle5", tdir) and "-super5" in X.tool_version("muscle5", tdir)
    assert "v7.fake" in X.tool_version("mafft_linsi", tdir)
    a, _ = X.command("muscle5", "in", "out", 200, tdir)
    b, _ = X.command("muscle5", "in", "out", 201, tdir)
    assert a[1] == "-align" and b[1] == "-super5" and "-threads" in a
    m, so = X.command("mafft_linsi", "in", "out", 5, tdir)
    assert "--localpair" in m and "1000" in m and so == "out"
    assert X.tool_build_id("mafft_linsi", tdir) != X.tool_build_id("clustalo", tdir)


@pytest.mark.parametrize("tool", ["mafft_linsi", "muscle5", "clustalo"])
def test_run_external_cache(tmp_path, tdir, tool, monkeypatch):
    monkeypatch.setenv("MODELOG", str(tmp_path / "mode.log"))
    case = _case(tmp_path, body=">a\nacgt\n>b\nacg\n")
    kw = dict(cache_path=tmp_path / "c.sqlite", scorer=SCORER, tdir=tdir)
    r = X.run_external(tool, case, **kw)
    if tool == "muscle5":  # fake output has other residues -> still stored (fake scorer does not check)
        assert (tmp_path / "mode.log").read_text().strip() == "-align"
    assert r["error"] is None and r["cached"] == 0 and r["cpu_s"] is not None and r["aln_sha256"]
    r2 = X.run_external(tool, case, **kw)
    assert r2["cached"] == 1 and r2["aln_sha256"] == r["aln_sha256"]


def test_failure_recorded_and_retry(tmp_path, tdir):
    case = _case(tmp_path, cid="fail1", body=">a\nfail\n>b\nacg\n")
    kw = dict(cache_path=tmp_path / "c.sqlite", scorer=SCORER, tdir=tdir)
    r = X.run_external("clustalo", case, **kw)
    assert r["error"].startswith("exit 3") and "boom" in r["error"]
    assert X.run_external("clustalo", case, **kw)["cached"] == 1


def test_score_error_recorded(tmp_path, tdir):
    case = _case(tmp_path, cid="x/bad")  # fake scorer raises on ids ending with "bad"
    r = X.run_external("clustalo", case, cache_path=tmp_path / "c.sqlite", scorer=SCORER, tdir=tdir)
    assert r["error"].startswith("score ValueError")


def test_evaluate_parallel_and_register(tmp_path, tdir):
    cases = [_case(tmp_path, cid=f"p/{i}") for i in range(6)]
    rows = X.evaluate_external(cases, "clustalo", 3, cache_path=tmp_path / "c.sqlite", scorer=SCORER, tdir=tdir)
    assert [r["case_id"] for r in rows] == [c["case_id"] for c in cases]
    assert all(r["tool"] == "clustalo" and r["error"] is None for r in rows)
    yaml = pytest.importorskip("yaml")
    reg = tmp_path / "REG.yaml"
    reg.write_text(yaml.safe_dump([{"name": "keep", "sha256": "ab"}]))
    X.register(X.tools_info(tdir), reg)
    got = {e["name"]: e for e in yaml.safe_load(reg.read_text())}
    assert set(got) == {"keep", "tool_mafft_linsi", "tool_muscle5", "tool_clustalo"}
    assert got["tool_clustalo"]["kind"] == "tool" and len(got["tool_clustalo"]["exe_sha256"]) == 64
