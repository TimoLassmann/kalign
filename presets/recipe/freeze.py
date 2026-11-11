"""Freeze manifest: proves that the code that runs is the code that was committed and tagged (Freeze 1).

    python -m recipe.freeze record --rev HEAD --out FILE [--tag NAME]   (on the laptop, inside the git repo)
    python -m recipe.freeze record --from-tree --out FILE               (smoke tests: hash the working tree)
    python -m recipe.freeze check  --manifest FILE [--root PRESETS_DIR]  (on the workstation, no git needed)
    python -m recipe.freeze hash-data --ws WS --out FILE | --verify FILE

The manifest holds the commit id and the sha256 of every frozen file: presets/recipe/*.py, presets/tests/*.py,
presets/resources/*, presets/run_all.sh and presets/docs/PLAN.md. `check` fails if any frozen file differs, is
missing, or if an extra .py file sits in recipe/ or tests/ (so nothing can be added unnoticed).
`hash-data` records the sha256 of the data tables the fit reads (case manifests, folds, cores, features), so a
resumed run can show that its inputs did not change either.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

PATTERNS = ("recipe/*.py", "tests/*.py", "resources/*", "run_all.sh", "run_v2.sh", "docs/PLAN.md")
DATA_PATTERNS = ("cases_*.final.tsv", "folds_*.tsv", "features_*.final.tsv", "groups_*.final.tsv",
                 "cores/*.txt")


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def tree_files(root: Path):
    out = {}
    for pat in PATTERNS:
        for p in sorted(root.glob(pat)):
            if p.is_file() and "__pycache__" not in p.parts:
                out[f"presets/{p.relative_to(root).as_posix()}"] = sha(p.read_bytes())
    return out


def git_files(rev: str, repo: Path):
    names = subprocess.run(["git", "ls-tree", "-r", "--name-only", rev, "presets/"], cwd=repo, check=True,
                           capture_output=True, text=True).stdout.split()
    out = {}
    for n in names:
        rel = n[len("presets/"):]
        if any(Path(rel).match(pat) and Path(rel).parent == Path(pat).parent for pat in PATTERNS):
            blob = subprocess.run(["git", "show", f"{rev}:{n}"], cwd=repo, check=True, capture_output=True).stdout
            out[n] = sha(blob)
    return out


def record(args):
    root = Path(__file__).resolve().parents[1]
    repo = root.parent
    if args.from_tree:
        files, commit, dirty = tree_files(root), "WORKTREE", True
    else:
        commit = subprocess.run(["git", "rev-parse", args.rev], cwd=repo, check=True, capture_output=True,
                                text=True).stdout.strip()
        files = git_files(commit, repo)
        dirty = False
        if tree_files(root) != files:
            raise SystemExit("the working tree differs from the recorded revision; commit first or use --from-tree")
    man = {"commit": commit, "tag": args.tag, "worktree_snapshot": dirty,
           "recorded": time.strftime("%Y-%m-%dT%H:%M:%S"), "n_files": len(files), "files": files}
    Path(args.out).write_text(json.dumps(man, indent=1) + "\n")
    print(f"freeze manifest: {commit[:12]} {args.tag or ''} {len(files)} files -> {args.out}")


def check_tree(man, root: Path):
    """-> list of problems (empty when the tree is exactly the frozen one)."""
    have = tree_files(root)
    want = man["files"]
    bad = [f"differs: {k}" for k in want if k in have and have[k] != want[k]]
    bad += [f"missing: {k}" for k in want if k not in have]
    bad += [f"not in the freeze: {k}" for k in have if k not in want]
    return bad


def check(args):
    man = json.loads(Path(args.manifest).read_text())
    root = Path(args.root) if args.root else Path(__file__).resolve().parents[1]
    bad = check_tree(man, root)
    if bad:
        print(f"FREEZE CHECK FAILED against commit {man['commit']}:", file=sys.stderr)
        for b in bad:
            print("  " + b, file=sys.stderr)
        return 1
    print(f"freeze check ok: code equals commit {man['commit'][:12]} {man.get('tag') or ''} "
          f"({man['n_files']} files)")
    return 0


def data_hashes(ws: Path):
    d = Path(ws) / "data"
    out = {}
    for pat in DATA_PATTERNS:
        for p in sorted(d.glob(pat)):
            out[p.relative_to(d).as_posix()] = sha(p.read_bytes())
    if not out:
        raise SystemExit(f"{d}: no data tables found")
    return out


def hash_data(args):
    ws = Path(args.ws).expanduser()
    cur = data_hashes(ws)
    if args.verify:
        old = json.loads(Path(args.verify).read_text())["files"]
        bad = [k for k in sorted(set(old) | set(cur)) if old.get(k) != cur.get(k)]
        if bad:
            print("DATA CHANGED since the run started: " + ", ".join(bad), file=sys.stderr)
            return 1
        print(f"data tables unchanged ({len(cur)} files)")
        return 0
    Path(args.out).write_text(json.dumps({"recorded": time.strftime("%Y-%m-%dT%H:%M:%S"), "files": cur},
                                         indent=1) + "\n")
    print(f"data hashes: {len(cur)} files -> {args.out}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("record")
    r.add_argument("--rev", default="HEAD")
    r.add_argument("--from-tree", action="store_true")
    r.add_argument("--tag")
    r.add_argument("--out", required=True)
    c = sub.add_parser("check")
    c.add_argument("--manifest", required=True)
    c.add_argument("--root")
    h = sub.add_parser("hash-data")
    h.add_argument("--ws", required=True)
    h.add_argument("--out")
    h.add_argument("--verify")
    a = ap.parse_args(argv)
    if a.cmd == "record":
        record(a)
        return 0
    if a.cmd == "check":
        return check(a)
    if not (a.out or a.verify):
        ap.error("--out or --verify")
    return hash_data(a)


if __name__ == "__main__":
    sys.exit(main())
