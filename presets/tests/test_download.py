import hashlib
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
pytest.importorskip("yaml")
from recipe import download as d  # noqa: E402


def test_sha256_and_tree_hash(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "x").write_bytes(b"hello")
    assert d.sha256(tmp_path / "a" / "x") == hashlib.sha256(b"hello").hexdigest()
    h1 = d.tree_hash(tmp_path)
    (tmp_path / "a" / "x").write_bytes(b"hellp")
    assert d.tree_hash(tmp_path) != h1


def test_registry_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(d, "REGISTRY", tmp_path / "r" / "REGISTRY.yaml")
    d.save_registry({"x": {"name": "x", "sha256": "ab"}})
    assert d.load_registry() == {"x": {"name": "x", "sha256": "ab"}}


def test_mismatch_is_error(tmp_path):
    entry = ("b", "n", "http://invalid/", "f.txt", "0" * 64, "v", "l", "")
    (tmp_path / "b").mkdir()
    (tmp_path / "b" / "f.txt").write_text("data")
    with pytest.raises(SystemExit):
        d.do_file(entry, tmp_path, {}, True)


def test_pinned_rfam_not_current():
    assert all("CURRENT" not in e[2] for e in d.FILES)
