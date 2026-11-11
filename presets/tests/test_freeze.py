"""Freeze manifest (recipe.freeze), smoke-test subset helper and run_all.sh syntax."""
import json
import shutil
import subprocess
from pathlib import Path

from recipe import freeze as F
from recipe import nested as N

PRESETS = Path(__file__).resolve().parents[1]


def make_tree(root):
    (root / "recipe").mkdir(parents=True)
    (root / "tests").mkdir()
    (root / "resources").mkdir()
    (root / "recipe" / "a.py").write_text("x = 1\n")
    (root / "tests" / "test_a.py").write_text("def test(): pass\n")
    (root / "resources" / "REGISTRY.yaml").write_text("- id: x\n")
    (root / "run_all.sh").write_text("#!/bin/sh\n")
    (root / "docs").mkdir()
    (root / "docs" / "PLAN.md").write_text("plan\n")
    (root / "docs" / "other.md").write_text("not frozen\n")
    (root / "recipe" / "__pycache__").mkdir()
    (root / "recipe" / "__pycache__" / "a.pyc").write_text("junk")


def manifest(root):
    return {"commit": "abc", "files": F.tree_files(root), "n_files": len(F.tree_files(root))}


def test_tree_files_selects_frozen_files_only(tmp_path):
    make_tree(tmp_path)
    assert set(F.tree_files(tmp_path)) == {
        "presets/recipe/a.py", "presets/tests/test_a.py", "presets/resources/REGISTRY.yaml",
        "presets/run_all.sh", "presets/docs/PLAN.md"}


def test_check_detects_change_missing_and_extra(tmp_path):
    make_tree(tmp_path)
    man = manifest(tmp_path)
    assert F.check_tree(man, tmp_path) == []
    (tmp_path / "docs" / "other.md").write_text("edited\n")            # not frozen: fine
    assert F.check_tree(man, tmp_path) == []
    (tmp_path / "recipe" / "a.py").write_text("x = 2\n")
    assert F.check_tree(man, tmp_path) == ["differs: presets/recipe/a.py"]
    (tmp_path / "recipe" / "a.py").write_text("x = 1\n")
    (tmp_path / "recipe" / "b.py").write_text("")
    assert F.check_tree(man, tmp_path) == ["not in the freeze: presets/recipe/b.py"]
    (tmp_path / "recipe" / "b.py").unlink()
    (tmp_path / "tests" / "test_a.py").unlink()
    assert F.check_tree(man, tmp_path) == ["missing: presets/tests/test_a.py"]


def test_cli_check_exit_codes(tmp_path):
    make_tree(tmp_path)
    mf = tmp_path / "m.json"
    mf.write_text(json.dumps(manifest(tmp_path)))
    assert F.main(["check", "--manifest", str(mf), "--root", str(tmp_path)]) == 0
    (tmp_path / "recipe" / "a.py").write_text("changed")
    assert F.main(["check", "--manifest", str(mf), "--root", str(tmp_path)]) == 1


def test_data_hash_roundtrip(tmp_path):
    d = tmp_path / "data" / "cores"
    d.mkdir(parents=True)
    (tmp_path / "data" / "folds_protein.tsv").write_text("a\n")
    (d / "protein_fold0_core.txt").write_text("c1\n")
    out = tmp_path / "h.json"
    assert F.main(["hash-data", "--ws", str(tmp_path), "--out", str(out)]) == 0
    assert F.main(["hash-data", "--ws", str(tmp_path), "--verify", str(out)]) == 0
    (d / "protein_fold0_core.txt").write_text("c2\n")
    assert F.main(["hash-data", "--ws", str(tmp_path), "--verify", str(out)]) == 1


def test_limit_ids_even_spacing_and_noop():
    ids = [f"c{i:03d}" for i in range(100)]
    sub = N.limit_ids(ids, 10)
    assert len(sub) == 10 and sub[0] == "c000" and sub[-1] == "c090" and sub == sorted(sub)
    assert N.limit_ids(ids, 0) == ids and N.limit_ids(ids, 500) == ids
    assert N.limit_ids(list(reversed(ids)), 10) == sub


def test_run_all_syntax_and_refusal(tmp_path):
    sh = PRESETS / "run_all.sh"
    assert subprocess.run(["bash", "-n", str(sh)]).returncode == 0
    if shutil.which("flock") is None:       # run_all.sh is a workstation (Linux) script
        return
    r = subprocess.run(["bash", str(sh), "run", "--ws", str(tmp_path), "--out", str(tmp_path / "o")],
                       capture_output=True, text=True)
    assert r.returncode == 4 and "needs --freeze" in r.stderr   # real run without a freeze manifest is refused
