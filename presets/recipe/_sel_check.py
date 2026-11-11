"""Known-answer tests S2-S4 of docs/AMENDMENT_ensconf_2.md (runs in the ALIGNER venv of build 5d56876).

    <venv>/bin/python _sel_check.py JOB.json RESULT.json
JOB.json: {"params": <select_only worker dict>, "cases": [[case_id, input, outdir], ...]}

Per case:
  S2 the select output equals, byte for byte, the selected member aligned alone (n_runs = 1, the member's parameters,
     seed 42 + k);
  S3 the returned column confidence equals the confidence computed here from the K member alignments (|d| <= 1e-5);
  S4 the selected member is the argmax (ties: lower index) of the POAR scores computed here with integers:
     N_k = sum over member pairs l != k of the residue pairs aligned by both k and l
     (the C score is N_k / (K - 1)).
"""
import hashlib
import json
import os
import sys

import numpy as np

import _align_worker_sel  # noqa: F401 - patches _align_worker (select_only, run_tree_sl_noise)
import _align_worker as W
from kalign import _core


def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def residue_cols(recs, names):
    """Flattened residue -> column, sequences in `names` order."""
    d = dict(recs)
    out = []
    for n in names:
        s = d[n]
        out.extend(j for j, c in enumerate(s) if c.isalpha())
    return np.asarray(out, dtype=np.int64)


def group_pairs(*keys):
    """Number of residue pairs that share every key."""
    k = np.stack(keys, axis=1)
    _, cnt = np.unique(k, axis=0, return_counts=True)
    return int((cnt * (cnt - 1) // 2).sum())


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def member_kwargs(kw, k):
    m = {}
    for a, v in kw.items():
        m[a] = [v[k]] if a.startswith("run_") and isinstance(v, list) else v
    m["seed"] = kw["seed"] + k
    m.pop("select_only", None)
    return m


def main(job_path, out_path):
    job = json.load(open(job_path))
    kw = W.build_kwargs({"params": job["params"], "seed": 42})
    K = len(kw["run_gpo"])
    res = []
    for cid, inp, od in job["cases"]:
        os.makedirs(od, exist_ok=True)
        out = os.path.join(od, "select.fa")
        r = _core.ensemble_custom_file_to_file(inp, out, **kw)
        sel = int(r["selected_member"])
        names = list(r["names"])
        conf_c = np.asarray(r["column_confidence"], dtype=np.float64)
        mem = []
        for k in range(K):
            mo = os.path.join(od, f"member{k}.fa")
            _core.ensemble_custom_file_to_file(inp, mo, **member_kwargs(kw, k))
            mem.append(mo)
        cols = [residue_cols(read_fasta(m), names) for m in mem]
        o = residue_cols(read_fasta(out), names)
        # S2
        s2 = sha(out) == sha(mem[sel]) if sel >= 0 else False
        # S3: per residue, the number of other residues of its output column that member l puts in its column
        n_o = np.bincount(o)
        isum = np.zeros(len(o), dtype=np.int64)
        for c in cols:
            pair = o * (int(c.max()) + 1) + c
            _, inv, cnt = np.unique(pair, return_inverse=True, return_counts=True)
            isum += cnt[inv] - 1
        m = n_o[o]
        rconf = np.where(m > 1, isum / np.maximum((m - 1) * K, 1), 0.0)
        width = len(read_fasta(out)[0][1])
        csum = np.bincount(o, weights=rconf, minlength=width)
        conf_p = np.divide(csum, np.bincount(o, minlength=width), out=np.zeros(width),
                           where=np.bincount(o, minlength=width) > 0)
        s3_max = float(np.max(np.abs(conf_p - conf_c))) if len(conf_c) == width else float("inf")
        # S4
        N = [sum(group_pairs(cols[k], cols[l]) for l in range(K) if l != k) for k in range(K)]
        best = max(range(K), key=lambda k: (N[k], -k))
        res.append({"case_id": cid, "selected": sel, "N": N, "argmax": best, "S2": bool(s2),
                    "S3_max_abs": s3_max, "S3": s3_max <= 1e-5, "S4": best == sel})
        print(json.dumps(res[-1]), flush=True)
    json.dump(res, open(out_path, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
