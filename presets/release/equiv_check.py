"""Equivalence gate: the built-in modes (CLI and Python mode API) reproduce the evaluated parameter path.

    <venv>/bin/python presets/release/equiv_check.py --kalign BIN --cases LIST --out DIR [--workers N] [--modes fast,accurate]

LIST: text file, one input file per line, optionally "<pool>\t<path>"; the pool is decided by Kalign itself
(biotype detection), the reference preset is chosen to match that biotype.
For every case and mode: (a) _core.ensemble_custom_file_to_file with the presets_v2.json parameters (the path the
re-fit worker used), (b) the CLI `BIN -i IN -o OUT --mode M -n 1`, (c) _core.align_file_to_file_mode(IN, OUT, M).
All three in FASTA; (b) and (c) must equal (a) byte for byte. Resumable: a case whose three outputs exist is only
compared. Writes DIR/equiv.tsv and prints a summary; exit 1 on any mismatch.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import subprocess
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
ROOT = Path(__file__).resolve().parents[2]
P = json.loads((ROOT / "presets/results/v2/presets_v2.json").read_text())["presets"]


def kw(p):
    return dict(run_gpo=p["run_gpo"], run_gpe=p["run_gpe"], run_tgpe=p["run_tgpe"], run_noise=p["run_noise"],
                run_types=p["run_types"], format="fasta", seq_type=p["run_types"][0], seed=42,
                min_support=p["min_support"], realign=p["realign"], seq_weights=p["seq_weights"], n_threads=1,
                consistency_anchors=p["consistency"], consistency_weight=p["consistency_weight"],
                run_vsm_amax=p["run_vsm_amax"], run_refine=p["run_refine"],
                consistency_merge=p["consistency_merge"], consistency_merge_weight=p["consistency_merge_weight"],
                vsm_mode=p["vsm_mode"], run_vsm_kappa=p["run_vsm_kappa"], run_vsm_dzero=p["run_vsm_dzero"],
                fill_in=p["fill_in"], fill_anchor_frac=p["fill_anchor_frac"], fill_support=p["fill_support"],
                run_guide_dist=p["run_guide_dist"])


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest() if Path(p).exists() else "missing"


def one(task):
    from kalign import _core
    kalign, pool_hint, inp, mode, d = task
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    ref, cli, api = d / "ref.fa", d / "cli.fa", d / "api.fa"
    err = []
    nuc = None
    if not cli.exists():
        r = subprocess.run([kalign, "-i", inp, "-o", str(cli), "--mode", mode, "-n", "1"], capture_output=True, text=True)
        if r.returncode:
            err.append(f"cli exit {r.returncode}")
    if pool_hint:
        nuc = pool_hint == "nucleotide"
    if not api.exists():
        try:
            _core.align_file_to_file_mode(inp, str(api), mode, format="fasta", n_threads=1)
        except Exception as e:  # noqa: BLE001
            err.append(f"api {e}")
    if not ref.exists():
        pools = ["nucleotide", "protein"] if nuc is None else (["nucleotide"] if nuc else ["protein"])
        for pool in pools:   # unknown biotype: the reference that matches the CLI output is the one Kalign chose
            tmp = d / f"ref_{pool}.fa"
            try:
                _core.ensemble_custom_file_to_file(inp, str(tmp), **kw(P[pool][mode]["params"]))
            except Exception as e:  # noqa: BLE001
                err.append(f"ref {pool} {e}")
        cands = [d / f"ref_{p}.fa" for p in pools if (d / f"ref_{p}.fa").exists()]
        hit = [x for x in cands if sha(x) == sha(cli)]
        (hit[0] if hit else cands[0]).rename(ref) if cands else None
        for x in cands:
            if x.exists():
                x.unlink()
    a, b, c2 = sha(ref), sha(cli), sha(api)
    return inp, mode, a == b, a == c2, ";".join(err)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kalign", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--modes", default="fast,accurate")
    a = ap.parse_args()
    cases = [(ln.split("\t")[0] if "\t" in ln else "", ln.split("\t")[-1].strip()) for ln in open(a.cases) if ln.strip()]
    tasks = [(a.kalign, pool, c, m, str(a.out / "aln" / m / hashlib.sha1(c.encode()).hexdigest()[:16]))
             for m in a.modes.split(",") for pool, c in cases]
    rows, bad = [], 0
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for i, r in enumerate(ex.map(one, tasks, chunksize=2), 1):
            rows.append(r)
            bad += not (r[2] and r[3])
            if i % 500 == 0:
                print(f"{i}/{len(tasks)} mismatches {bad}", flush=True)
    a.out.mkdir(parents=True, exist_ok=True)
    with open(a.out / "equiv.tsv", "w") as fh:
        fh.write("input\tmode\tcli_equal\tapi_equal\terror\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    for m in a.modes.split(","):
        rr = [r for r in rows if r[1] == m]
        print(f"{m}: {len(rr)} cases, CLI identical {sum(r[2] for r in rr)}, API identical {sum(r[3] for r in rr)}, "
              f"errors {sum(bool(r[4]) for r in rr)}")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
