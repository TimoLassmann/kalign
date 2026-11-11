#!/usr/bin/env python3
"""Generate the simulated perf datasets with Kalign's DSSim simulator.

Grid (protein): N in {50, 200, 500} x L in {150, 500} x
divergence {lo: n_obs=20, hi: n_obs=2}.  DNA: N=200, L=500, n_obs=8.

HARD LIMIT: at most 500 sequences and 500 residues per sequence (a laptop
went into heavy swap with larger sets).  main() refuses anything larger.
    uv run --no-project python perf/dp/scripts/make_datasets.py ...

L is DSSim's model length; realised mean sequence lengths are within a few
percent of L and are recorded in manifest.json together with every generator
parameter and seed, so the files can be regenerated bit for bit:

    uv run --no-project python perf/dp/scripts/make_datasets.py \
        --dssim-gen build/perf/dp/dssim_gen --out $DATA

Datasets are not committed; put them outside the repository.
"""
import argparse
import json
import os
import subprocess
import sys

PROT_N = [50, 200, 500]
PROT_L = [150, 500]
DIV = {"lo": 20, "hi": 2}  # DSSim n_obs: larger = less divergent
DNA = [(200, 500, 8)]
MAX_N = 500
MAX_L = 500
BASE_SEED = 20261002


def grid():
    rows = []
    for di, (dname, nobs) in enumerate(DIV.items()):
        for ni, n in enumerate(PROT_N):
            for li, length in enumerate(PROT_L):
                seed = BASE_SEED + 100 * ni + 10 * li + di
                rows.append(dict(name=f"prot_N{n}_L{length}_{dname}", biotype="protein",
                                 n_seq=n, n_obs=nobs, dna=0, model_len=length,
                                 divergence=dname, seed=seed))
    for k, (n, length, nobs) in enumerate(DNA):
        rows.append(dict(name=f"dna_N{n}_L{length}", biotype="dna", n_seq=n, n_obs=nobs,
                         dna=1, model_len=length, divergence="mid",
                         seed=BASE_SEED + 5000 + k))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dssim-gen", required=True, help="path to the dssim_gen binary")
    ap.add_argument("--out", required=True, help="output directory for FASTA files + manifest.json")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    manifest = dict(generator="tests/dssim.c via perf/dp/tools/dssim_gen",
                    command="dssim_gen <n_seq> <n_obs> <dna> <model_len> <seed> <out.fa> <max_len>",
                    datasets=[])
    for row in grid():
        if row["n_seq"] > MAX_N or row["model_len"] > MAX_L:
            sys.exit(f"refusing {row['name']}: exceeds N<={MAX_N}, L<={MAX_L}")
        fa = os.path.join(args.out, row["name"] + ".fa")
        cmd = [args.dssim_gen, str(row["n_seq"]), str(row["n_obs"]), str(row["dna"]),
               str(row["model_len"]), str(row["seed"]), fa, str(MAX_L)]
        out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.split()
        row["file"] = os.path.basename(fa)
        row["numseq"] = int(out[0])
        row["mean_len"] = float(out[1])
        row["max_len"] = int(out[2])
        row["n_truncated_to_max"] = int(out[3])
        manifest["datasets"].append(row)
        print(f"{row['name']:24s} n={row['numseq']:5d} mean_len={row['mean_len']:8.1f} "
              f"max_len={row['max_len']} truncated={row['n_truncated_to_max']} seed={row['seed']}")
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
