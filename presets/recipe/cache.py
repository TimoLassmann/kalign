"""One SQLite cache for the whole project:
(build_id, param_key, case_id, sha256_unaligned, sha256_reference, scorer_version) -> per-case result.

Alignments are deterministic, so a cached entry is exact. Rules:
  * an existing entry is never overwritten (INSERT OR IGNORE); the only exception is
    ``put(..., replace_error=True)`` which may replace an entry whose ``error`` is set;
  * safe for ~30 concurrent writer processes: WAL mode, ``busy_timeout``, one short transaction per
    insert, retry with jittered back-off on "database is locked";
  * connections are per process (re-opened after fork) and never shared between threads.

Key parts:
  build_id   SHA-256 of the aligner CLI binary (the "aligner build id")
  param_key  sha256 of the canonical JSON of the spec: {"preset": mode} or a parameter dict
  case_id    from the case manifest
  sha256_unaligned / sha256_reference
             content hashes of the case's input and reference files (manifest columns): a rebuilt
             manifest with different files under the same case_id is a cache miss, never a stale hit
  scorer_version
             constant ``SCORER_VERSION`` of the pool scorer module; bump it on ANY change to scoring

Schema 2 (this file). A schema-1 database (key without the hashes and scorer version) is upgraded on
first open: its table is renamed ``results_v1`` and is NEVER read by ``get``/``get_many``; rows are
carried over only by the explicit ``migrate_v1`` (fills the hashes from the manifests, scorer v1).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import random
import socket
import sqlite3
import threading
import time
from typing import NamedTuple
from datetime import datetime, timezone
from pathlib import Path

SCHEMA_VERSION = 2
METRIC_COLS = ("recall", "precision", "f1", "tc", "width_ratio")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS results (
  build_id   TEXT NOT NULL,
  param_key  TEXT NOT NULL,
  case_id    TEXT NOT NULL,
  sha256_unaligned  TEXT NOT NULL,
  sha256_reference  TEXT NOT NULL,
  scorer_version    TEXT NOT NULL,
  recall     REAL, precision REAL, f1 REAL, tc REAL, width_ratio REAL,
  n_seq      INTEGER,
  extra      TEXT,            -- JSON: any further scorer outputs
  cpu_s      REAL,            -- user+sys of the aligner process (single thread)
  wall_s     REAL,
  aln_sha256 TEXT,
  host       TEXT,
  ts         TEXT,
  error      TEXT,            -- NULL = ok
  harness    TEXT,
  PRIMARY KEY (build_id, param_key, case_id, sha256_unaligned, sha256_reference, scorer_version)
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS specs (
  param_key TEXT PRIMARY KEY,
  spec_json TEXT NOT NULL
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT) WITHOUT ROWID;
"""


# ---------------------------------------------------------------- canonical keys
def _canon(x):
    if hasattr(x, "tolist") and not isinstance(x, (list, tuple, dict, str, bytes)):
        x = x.tolist()  # numpy array -> list; numpy scalar -> python scalar
    if isinstance(x, bool) or x is None or isinstance(x, str):
        return x
    if isinstance(x, int):
        return x
    if isinstance(x, float):
        if not math.isfinite(x):
            raise ValueError(f"non-finite value in spec: {x}")
        r = round(x, 9)
        if r == 0:
            return 0  # also folds -0.0
        if r == int(r) and abs(r) < 1e15:
            return int(r)  # 1.0 and 1 are the same parameter
        return r
    if isinstance(x, (list, tuple)):
        return [_canon(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _canon(v) for k, v in sorted(x.items(), key=lambda kv: str(kv[0]))}
    raise TypeError(f"cannot canonicalise {type(x).__name__}")


def canonical_json(spec) -> str:
    """Canonical JSON of a spec (a preset name, {"preset": m} or a parameter dict)."""
    if isinstance(spec, str):
        spec = {"preset": spec}
    return json.dumps(_canon(spec), sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def param_key(spec) -> str:
    """Stable key of a spec: sha256 hex of the canonical JSON."""
    return hashlib.sha256(canonical_json(spec).encode()).hexdigest()


def sha256_file(path, bufsize=1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(bufsize):
            h.update(chunk)
    return h.hexdigest()


def default_ws() -> Path:
    return Path(os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit")).expanduser()


def default_cache_path() -> Path:
    return default_ws() / "cache" / "cache.sqlite"


class Ident(NamedTuple):
    """Everything besides (build, spec) that a cached result depends on."""
    case_id: str
    sha256_unaligned: str
    sha256_reference: str
    scorer_version: str


def ident(case_row, scorer_version) -> Ident:
    """Ident of a manifest row. Hashes may be "" only in tests; real manifests always carry them."""
    return Ident(case_row["case_id"], case_row.get("sha256_unaligned") or "",
                 case_row.get("sha256_reference") or "", scorer_version)


# ---------------------------------------------------------------- the cache
class Cache:
    def __init__(self, path=None, timeout=60.0):
        self.path = Path(path) if path else default_cache_path()
        self.timeout = timeout
        self._local = threading.local()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn()  # create schema now

    # one connection per (process, thread)
    def _conn(self) -> sqlite3.Connection:
        c = getattr(self._local, "conn", None)
        if c is not None and getattr(self._local, "pid", None) == os.getpid():
            return c
        c = sqlite3.connect(str(self.path), timeout=self.timeout, isolation_level=None)
        c.row_factory = sqlite3.Row
        c.execute(f"PRAGMA busy_timeout={int(self.timeout * 1000)}")
        self._retry(lambda: c.execute("PRAGMA journal_mode=WAL").fetchall())
        c.execute("PRAGMA synchronous=NORMAL")
        self._retry(lambda: self._upgrade(c))
        self._retry(lambda: c.executescript(_SCHEMA))
        self._retry(lambda: c.execute("INSERT OR REPLACE INTO meta VALUES ('schema', ?)",
                                      (str(SCHEMA_VERSION),)))
        self._local.conn, self._local.pid = c, os.getpid()
        return c

    @staticmethod
    def _upgrade(c):
        """Schema 1 -> 2: park the old table as results_v1 (unreadable by get/get_many)."""
        def legacy():
            cols = [r[1] for r in c.execute("PRAGMA table_info(results)").fetchall()]
            return bool(cols) and "scorer_version" not in cols
        if not legacy():
            return
        c.execute("BEGIN IMMEDIATE")
        try:
            if legacy():  # re-check under the write lock: another process may have done it
                c.execute("ALTER TABLE results RENAME TO results_v1")
            c.execute("COMMIT")
        except BaseException:
            c.execute("ROLLBACK")
            raise

    @staticmethod
    def _retry(fn, tries=12):
        for i in range(tries):
            try:
                return fn()
            except sqlite3.OperationalError as e:
                msg = str(e).lower()
                if ("locked" not in msg and "busy" not in msg) or i == tries - 1:
                    raise
                time.sleep(min(5.0, 0.05 * 2 ** i) * (0.5 + random.random()))

    def close(self):
        c = getattr(self._local, "conn", None)
        if c is not None and getattr(self._local, "pid", None) == os.getpid():
            c.close()
        self._local.conn = None

    # ------------------------------------------------------------ reads
    def get(self, build_id, key, ident: Ident):
        """Return the cached row as a dict, or None. ``ident`` = cache.ident(case_row, scorer_version)."""
        r = self._retry(lambda: self._conn().execute(
            "SELECT * FROM results WHERE build_id=? AND param_key=? AND case_id=? AND sha256_unaligned=? "
            "AND sha256_reference=? AND scorer_version=?", (build_id, key, *ident)).fetchone())
        return _row(r) if r else None

    def get_many(self, build_id, key, idents):
        """Bulk lookup: {case_id: row} for the idents that are cached (all six key parts must match)."""
        out, want = {}, {i.case_id: i for i in idents}
        if len(want) != len(list(idents)):
            raise ValueError("duplicate case_id in idents")
        ids = list(want)
        for i in range(0, len(ids), 500):
            chunk = ids[i:i + 500]
            q = ("SELECT * FROM results WHERE build_id=? AND param_key=? AND case_id IN (%s)"
                 % ",".join("?" * len(chunk)))
            for r in self._retry(lambda: self._conn().execute(q, (build_id, key, *chunk)).fetchall()):
                w = want[r["case_id"]]
                if (r["sha256_unaligned"], r["sha256_reference"], r["scorer_version"]) == w[1:]:
                    out[r["case_id"]] = _row(r)
        return out

    def count(self, build_id=None, key=None) -> int:
        q, a = "SELECT COUNT(*) FROM results WHERE 1=1", []
        if build_id:
            q += " AND build_id=?"; a.append(build_id)
        if key:
            q += " AND param_key=?"; a.append(key)
        return self._retry(lambda: self._conn().execute(q, a).fetchone()[0])

    def spec_json(self, key):
        r = self._retry(lambda: self._conn().execute(
            "SELECT spec_json FROM specs WHERE param_key=?", (key,)).fetchone())
        return r[0] if r else None

    # ------------------------------------------------------------ writes
    def put(self, build_id, key, ident: Ident, metrics=None, *, cpu_s=None, wall_s=None, aln_sha256=None,
            error=None, spec=None, host=None, harness="", replace_error=False) -> bool:
        """Insert one result in one short transaction. Returns True if written, False if an entry
        already existed (never overwritten, except an errored one when replace_error=True)."""
        m = dict(metrics or {})
        vals = [m.pop(c, None) for c in METRIC_COLS]
        n_seq = m.pop("n_seq", None)
        row = (build_id, key, *ident, *vals, n_seq, json.dumps(m, sort_keys=True) if m else None,
               cpu_s, wall_s, aln_sha256, host or socket.gethostname(),
               datetime.now(timezone.utc).isoformat(timespec="seconds"), error, harness)

        def _do():
            c = self._conn()
            c.execute("BEGIN IMMEDIATE")
            try:
                if spec is not None:
                    c.execute("INSERT OR IGNORE INTO specs VALUES (?,?)", (key, canonical_json(spec)))
                if replace_error:
                    c.execute("DELETE FROM results WHERE build_id=? AND param_key=? AND case_id=? "
                              "AND sha256_unaligned=? AND sha256_reference=? AND scorer_version=? "
                              "AND error IS NOT NULL", (build_id, key, *ident))
                cur = c.execute("INSERT OR IGNORE INTO results VALUES (%s)" % ",".join("?" * 20), row)
                n = cur.rowcount
                c.execute("COMMIT")
                return n == 1
            except BaseException:
                c.execute("ROLLBACK")
                raise
        return self._retry(_do)


def _row(r) -> dict:
    d = dict(r)
    extra = d.pop("extra", None)
    if extra:
        d["extra"] = json.loads(extra)
    return d


# ---------------------------------------------------------------- schema-1 migration
def migrate_v1(cache: Cache, manifests, scorer_versions, *, verify_files=True) -> dict:
    """Carry the schema-1 rows (``results_v1``) over into ``results``.

    manifests        paths of the case manifests (TSV with case_id, pool, sha256_unaligned,
                     sha256_reference, unaligned, reference) that were current when the rows were made
    scorer_versions  {pool: scorer_version} the rows were scored with (step 0: "protein-v1", "nucleotide-v1")
    verify_files     recompute the sha256 of every input/reference file and skip rows whose file no
                     longer matches the manifest hash (so a manifest hash can never be attached to
                     stale content)
    Rows whose case_id is not in the manifests, or whose files do not match, are left in ``results_v1``
    (never readable). Existing ``results`` rows win (INSERT OR IGNORE). Returns counts.
    """
    import csv
    man = {}
    for p in manifests:
        with open(p, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                man[r["case_id"]] = r
    ok = {}
    n_bad_file = 0
    for cid, r in man.items():
        if r["pool"] not in scorer_versions:
            continue
        if verify_files and (sha256_file(r["unaligned"]) != r["sha256_unaligned"]
                             or sha256_file(r["reference"]) != r["sha256_reference"]):
            n_bad_file += 1
            continue
        ok[cid] = r
    c = cache._conn()
    if not c.execute("SELECT 1 FROM sqlite_master WHERE name='results_v1'").fetchone():
        return {"legacy_rows": 0, "migrated": 0, "skipped_not_in_manifest": 0, "skipped_file_mismatch": 0}
    legacy = cache._retry(lambda: c.execute("SELECT * FROM results_v1").fetchall())
    migrated = not_in = 0
    cols = [d[1] for d in c.execute("PRAGMA table_info(results_v1)").fetchall()]
    rest = [k for k in cols if k not in ("build_id", "param_key", "case_id")]

    def _do():
        nonlocal migrated, not_in
        migrated = not_in = 0
        c.execute("BEGIN IMMEDIATE")
        try:
            for row in legacy:
                r = ok.get(row["case_id"])
                if r is None:
                    not_in += 1
                    continue
                vals = [row[k] for k in rest]
                cur = c.execute(
                    "INSERT OR IGNORE INTO results (build_id,param_key,case_id,sha256_unaligned,"
                    "sha256_reference,scorer_version,%s) VALUES (?,?,?,?,?,?,%s)"
                    % (",".join(rest), ",".join("?" * len(rest))),
                    (row["build_id"], row["param_key"], row["case_id"], r["sha256_unaligned"],
                     r["sha256_reference"], scorer_versions[r["pool"]], *vals))
                migrated += cur.rowcount
            c.execute("INSERT OR REPLACE INTO meta VALUES ('migrated_v1', ?)",
                      (datetime.now(timezone.utc).isoformat(timespec="seconds"),))
            c.execute("COMMIT")
        except BaseException:
            c.execute("ROLLBACK")
            raise
    cache._retry(_do)
    return {"legacy_rows": len(legacy), "migrated": migrated, "skipped_not_in_manifest": not_in,
            "skipped_file_mismatch_cases": n_bad_file}


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="cache maintenance")
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("migrate-v1", help="carry schema-1 rows over (hashes from the manifests)")
    m.add_argument("--manifests", nargs="+", required=True)
    m.add_argument("--cache")
    m.add_argument("--protein-scorer-version", default="protein-v1")
    m.add_argument("--nucleotide-scorer-version", default="nucleotide-v1")
    m.add_argument("--no-verify", action="store_true")
    a = ap.parse_args(argv)
    cache = Cache(a.cache)
    print(json.dumps(migrate_v1(cache, a.manifests, {"protein": a.protein_scorer_version,
                                                     "nucleotide": a.nucleotide_scorer_version},
                                verify_files=not a.no_verify)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
