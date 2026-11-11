"""Download and register every benchmark of the plan (section 2a).

Run on the workstation (data never goes into git):

    uv run --with pyyaml python -m recipe.download [--ws ~/tmp/kalign-refit] [--only NAME ...] [--verify-only]

Layout: <ws>/data/raw/<benchmark>/{archive file, unpacked/ or repo/}.
Registry: presets/resources/REGISTRY.yaml (one entry per downloaded file; the git clone is
registered by commit SHA and a tree hash of its used subdirectory).

Re-runnable: a file whose sha256 matches the registry (or the pinned expectation below) is skipped
and re-verified; a mismatch is an error (nothing is overwritten silently).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

import yaml

PRESETS = Path(__file__).resolve().parents[1]
REGISTRY = PRESETS / "resources" / "REGISTRY.yaml"
RFAM_RELEASE = "15.1"  # pinned (latest numbered release dir on 2026-10-05); never CURRENT
WAYBACK1 = "https://web.archive.org/web/20160408045442id_/http://projects.binf.ku.dk/pgardner/bralibase/data-set1.tar.gz"
WAYBACK2 = "https://web.archive.org/web/20160408045453id_/http://projects.binf.ku.dk/pgardner/bralibase/data-set2.tar.gz"
MDSA = "https://dna.cs.byu.edu/mdsas/tarred/{}_mdsa_all.tar.gz"

# (benchmark dir, entry name, url, filename, expected sha256 or None, version, licence statement, note)
FILES = [
    ("drive5_bench", "drive5_bench_tar", "https://www.drive5.com/bench/bench.tar.gz", "bench.tar.gz", None,
     "drive5 bench.tar.gz (Last-Modified 2025-09-07): BAliBASE 3, PREFAB 4, OXBench, SABRE in FASTA",
     "none found (no licence file in the archive or on https://www.drive5.com/bench/); freely distributed for benchmarking",
     "Edgar, R.C."),
    ("rfam", "rfam_seed", f"https://ftp.ebi.ac.uk/pub/databases/Rfam/{RFAM_RELEASE}/Rfam.seed.gz", "Rfam.seed.gz", None,
     f"Rfam {RFAM_RELEASE}",
     "Rfam README: 'Rfam is freely available and in the public domain under the Creative Commons Zero licence.' (see COPYING)",
     "pinned release directory, not CURRENT"),
    ("rfam", "rfam_readme", f"https://ftp.ebi.ac.uk/pub/databases/Rfam/{RFAM_RELEASE}/README", "README", None,
     f"Rfam {RFAM_RELEASE}", "see rfam_seed", "licence statement source"),
    ("rfam", "rfam_copying", f"https://ftp.ebi.ac.uk/pub/databases/Rfam/{RFAM_RELEASE}/COPYING", "COPYING", None,
     f"Rfam {RFAM_RELEASE}", "see rfam_seed", "licence text"),
    ("bralibase", "bralibase_data_set1", WAYBACK1, "data-set1.tar.gz", None, "BRAliBASE II data-set1 (Wayback 2016-04-08)",
     "none found in archive; article (Gardner et al. 2005, NAR) CC BY; derived from Rfam", "Wayback Machine copy of projects.binf.ku.dk"),
    ("bralibase", "bralibase_data_set2", WAYBACK2, "data-set2.tar.gz", None, "BRAliBASE II data-set2 (Wayback 2016-04-08)",
     "none found in archive; article (Gardner et al. 2005, NAR) CC BY; derived from Rfam", "Wayback Machine copy of projects.binf.ku.dk"),
    ("mdsa", "mdsa_balibase", MDSA.format("balibase"), "balibase_mdsa_all.tar.gz",
     "499cee4a2a4942b9c17db3ba2b61a549b2abf5e9539963f71223ec30b6089c29", "MDSA 'all' (Carroll et al. 2007)",
     "none found (index page gone, files live)", "expected sha256 from experiments/redesign/resources/REGISTRY.yaml"),
    ("mdsa", "mdsa_oxbench", MDSA.format("oxbench"), "oxbench_mdsa_all.tar.gz",
     "be8b064ba800827a17de1fdbabf242556bb61a8895d04bd02a54b5236dcbb0aa", "MDSA 'all' (Carroll et al. 2007)",
     "none found (index page gone, files live)", "expected sha256 from experiments/redesign/resources/REGISTRY.yaml"),
    ("mdsa", "mdsa_smart", MDSA.format("smart"), "smart_mdsa_all.tar.gz",
     "e4e1fc4adbeaea69be101fa64886d5ae54d25fad0be459456b1c2f93145301b4", "MDSA 'all' (Carroll et al. 2007)",
     "none found (index page gone, files live)", "expected sha256 from experiments/redesign/resources/REGISTRY.yaml"),
]
BALIFAM_URL = "https://github.com/rcedgar/balifam"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_registry() -> dict:
    if REGISTRY.exists():
        return {e["name"]: e for e in yaml.safe_load(REGISTRY.read_text()) or []}
    return {}


def save_registry(reg: dict) -> None:
    REGISTRY.parent.mkdir(parents=True, exist_ok=True)
    head = ("# Data registry for the pooled re-fit (presets/). One entry per downloaded file.\n"
            "# Written by presets/recipe/download.py; do not edit sha256 values by hand.\n")
    REGISTRY.write_text(head + yaml.safe_dump(list(reg.values()), sort_keys=False, width=100))


def fetch(url: str, dest: Path, tries: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "kalign-refit-download/1"})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            tmp.rename(dest)
            return
        except Exception as e:  # noqa: BLE001
            print(f"  attempt {i + 1} failed: {e}", file=sys.stderr)
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"download failed: {url}")


def safe_extract(tar: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar) as t:
        t.extractall(out, filter="data")


def do_file(entry, raw: Path, reg: dict, verify_only: bool) -> None:
    bench, name, url, fname, expected, version, lic, note = entry
    dest = raw / bench / fname
    old = reg.get(name)
    pinned = expected or (old or {}).get("sha256")
    if dest.exists():
        got = sha256(dest)
        if pinned and got != pinned:
            raise SystemExit(f"SHA256 MISMATCH {name}: expected {pinned}, got {got} ({dest})")
        status = "verified"
        retrieved = dt.datetime.fromtimestamp(dest.stat().st_mtime, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    elif verify_only:
        raise SystemExit(f"missing {dest}")
    else:
        print(f"downloading {name} <- {url}")
        fetch(url, dest)
        got = sha256(dest)
        if pinned and got != pinned:
            raise SystemExit(f"SHA256 MISMATCH {name}: expected {pinned}, got {got}")
        status = "downloaded"
        retrieved = now()
    print(f"{status:10s} {name} {got[:12]} {dest.stat().st_size} bytes")
    if old is None or old["sha256"] != got:
        reg[name] = dict(name=name, benchmark=bench, url=url, version=version, license=lic,
                         retrieved=retrieved, sha256=got, size_bytes=dest.stat().st_size,
                         local_path=str(dest.relative_to(raw.parent)), note=note,
                         sha256_pinned_before_download=bool(expected))
    else:
        old["verified"] = now()
    if fname.endswith(".tar.gz") and not verify_only:
        out = raw / bench / "unpacked" / fname[: -len(".tar.gz")]
        if not out.exists():
            safe_extract(dest, out)
    if fname.endswith(".gz") and not fname.endswith(".tar.gz"):
        out = dest.with_suffix("")
        if not out.exists() and not verify_only:
            import gzip, shutil
            with gzip.open(dest) as a, open(out, "wb") as b:
                shutil.copyfileobj(a, b)


def tree_hash(root: Path) -> str:
    h = hashlib.sha256()
    for p in sorted(x for x in root.rglob("*") if x.is_file()):
        h.update(f"{sha256(p)}  {p.relative_to(root)}\n".encode())
    return h.hexdigest()


def do_balifam(raw: Path, reg: dict, verify_only: bool) -> None:
    repo = raw / "balifam" / "repo"
    if not repo.exists():
        if verify_only:
            raise SystemExit(f"missing {repo}")
        print(f"cloning {BALIFAM_URL}")
        subprocess.run(["git", "clone", "--quiet", BALIFAM_URL, str(repo)], check=True)
    sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
    old = reg.get("balifam_repo")
    if old and old["commit"] != sha:
        # The recorded commit is authoritative: check it out again rather than silently drift.
        subprocess.run(["git", "-C", str(repo), "checkout", "--quiet", old["commit"]], check=True)
        sha = old["commit"]
    lic = ""
    for cand in ("LICENSE", "LICENSE.md", "COPYING", "UNLICENSE"):
        if (repo / cand).exists():
            lic = (repo / cand).read_text()[:300].strip()
            break
    th = tree_hash(repo / "balifam100") if (repo / "balifam100").exists() else None
    date = subprocess.run(["git", "-C", str(repo), "log", "-1", "--format=%cI"], capture_output=True, text=True).stdout.strip()
    print(f"balifam commit {sha} ({date}) balifam100 tree sha256 {th}")
    if old and old.get("balifam100_tree_sha256") not in (None, th):
        raise SystemExit("balifam100 tree hash MISMATCH")
    if not old:
        reg["balifam_repo"] = dict(name="balifam_repo", benchmark="balifam", url=BALIFAM_URL, version=f"git {sha}",
                                   commit=sha, commit_date=date, license=lic or "none found", retrieved=now(),
                                   sha256=th, balifam100_tree_sha256=th, size_bytes=sum(
                                       p.stat().st_size for p in (repo / "balifam100").rglob("*") if p.is_file()) if th else None,
                                   local_path="raw/balifam/repo",
                                   note="git clone; only balifam100/ is used; sha256 = tree hash over (file sha256, relative path) of balifam100/")
    else:
        old["verified"] = now()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ws", default=os.path.expanduser("~/tmp/kalign-refit"))
    ap.add_argument("--only", nargs="*", help="benchmark dirs to process")
    ap.add_argument("--verify-only", action="store_true")
    a = ap.parse_args()
    raw = Path(a.ws).expanduser() / "data" / "raw"
    reg = load_registry()
    for e in FILES:
        if not a.only or e[0] in a.only:
            do_file(e, raw, reg, a.verify_only)
    if not a.only or "balifam" in a.only:
        do_balifam(raw, reg, a.verify_only)
    save_registry(reg)
    print(f"registry: {REGISTRY} ({len(reg)} entries)")


if __name__ == "__main__":
    main()
