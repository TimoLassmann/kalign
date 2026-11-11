"""Held-out comparison tables, one program per row (Addendum A2 pooling).

    python presets/verify/stacked_tables.py --a2 A2_DIR --externals V11/externals_protein.tsv

    (--oof REFIT_DIR: read the Kalign rows from REFIT_DIR/protein/fold*/oof_rows.tsv instead, e.g. the final
    rerun's merged rows: Kalign 4 on the shipped build, Kalign 3.7.0 on its own build)

A2_DIR: WS/runs/addendum_a2 (refit/<pool>/results/oof_bw.tsv, refit/protein/fold*/oof_rows.tsv via the v2 links,
<pool>_kalign31.tsv, <pool>_kalign335.tsv). Overall numbers are benchmark-weighted out-of-fold means; BAliBASE
category numbers are case means (BB and BBS together).
"""

import argparse
import csv
from pathlib import Path

ORDER = [("kalign31", "-", "Kalign 3.1"), ("kalign335", "-", "Kalign 3.3.5"),
         ("kalign361_fast", "-", "Kalign 3.6.1 fast"), ("kalign361_default", "-", "Kalign 3.6.1 default"),
         ("kalign361_recall", "-", "Kalign 3.6.1 recall"), ("kalign361_accurate", "-", "Kalign 3.6.1 accurate"),
         ("v370", "fast", "Kalign 3.7.0 fast"), ("v370", "default", "Kalign 3.7.0 default"),
         ("v370", "accurate", "Kalign 3.7.0 accurate"),
         ("new", "fast", "Kalign new fast"), ("new", "accurate", "Kalign new accurate"),
         ("famsa", "-", "FAMSA"), ("clustalo", "-", "Clustal Omega"), ("mafft_linsi", "-", "MAFFT L-INS-i"),
         ("muscle5", "-", "MUSCLE 5")]
CATS = ("RV11", "RV12", "RV20", "RV30", "RV40", "RV50")


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def overall(a2, pool, pooled=None):
    rows = {(r["arm"], r["mode"]): r for r in tsv(Path(pooled or a2) / "refit" / pool / "results" / "oof_bw.tsv")
            if r["subset_type"] == "all"}
    print(f"\n### {pool}: held out, benchmark-weighted (n = {next(iter(rows.values()))['n']})\n")
    print("| program | SP | TC | precision | CPU per case (s) |\n|---|---|---|---|---|")
    for arm, mode, name in ORDER:
        r = rows.get((arm, mode))
        if r:
            print(f"| {name} | {float(r['sp']):.3f} | {float(r['tc']):.3f} | {float(r['precision']):.3f} | "
                  f"{float(r['cpu']):.3f} |")


def bali(a2, externals, a3=None, oof=None):
    per, cat = {}, {}
    name = {a: n for a, _, n in ORDER}
    for r in tsv(externals):
        if r["benchmark"] == "bali3":
            cat[r["case_id"]] = r["subset"][-4:]
            per[(name[r["tool"]], r["case_id"])] = r
    for t in ("kalign31", "kalign335"):
        for r in tsv(Path(a2) / f"protein_{t}.tsv"):
            if r["benchmark"] == "bali3":
                per[(name[t], r["case_id"])] = r
    if a3:
        for r in tsv(Path(a3) / "protein_famsa.tsv"):
            if r["benchmark"] == "bali3":
                per[(name["famsa"], r["case_id"])] = r
    a1 = Path(a2).parent / "addendum_a1"
    for m in ("fast", "default", "recall", "accurate"):
        for r in tsv(a1 / f"protein_kalign361_{m}.tsv"):
            if r["benchmark"] == "bali3":
                per[(f"Kalign 3.6.1 {m}", r["case_id"])] = r
    for k in range(5):
        for r in tsv((Path(oof) if oof else Path(a2) / "refit") / "protein" / f"fold{k}" / "oof_rows.tsv"):
            if r["benchmark"] == "bali3":
                key = {("new", "fast"): "Kalign new fast", ("new", "accurate"): "Kalign new accurate",
                       ("v370", "fast"): "Kalign 3.7.0 fast", ("v370", "default"): "Kalign 3.7.0 default",
                       ("v370", "accurate"): "Kalign 3.7.0 accurate"}.get((r["arm"], r["mode"]))
                if key:
                    per[(key, r["case_id"].strip())] = r
    names = [n for _, _, n in ORDER]
    n = {c: sum(1 for v in cat.values() if v == c) for c in CATS}
    for met, col, fmt in (("SP", "recall", "{:.3f}"), ("TC", "tc", "{:.3f}"), ("CPU per case (s)", "cpu_s", "{:.3f}")):
        print(f"\n### BAliBASE 3, held out: {met}\n")
        print("| program | " + " | ".join(f"{c} (n={n[c]})" for c in CATS) + " | all (n=385) |")
        print("|---|" + "---|" * (len(CATS) + 1))
        for nm in names:
            vals, allv = [], []
            ok = True
            for c in CATS:
                ids = [i for i, v in cat.items() if v == c]
                if not all((nm, i) in per for i in ids):
                    ok = False
                    break
                x = [0.0 if per[(nm, i)]["error"] and col != "cpu_s" else float(per[(nm, i)][col]) for i in ids]
                vals.append(sum(x) / len(x))
                allv += x
            if ok:
                print(f"| {nm} | " + " | ".join(fmt.format(v) for v in vals) + f" | {fmt.format(sum(allv)/len(allv))} |")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a2", required=True)
    ap.add_argument("--externals", required=True)
    ap.add_argument("--a3", help="Addendum A3 dir (FAMSA rows and pooled tables incl. FAMSA)")
    ap.add_argument("--pooled", help="dir with refit/<pool>/results/oof_bw.tsv (default: --a3, else --a2)")
    ap.add_argument("--oof", help="refit dir with protein/fold*/oof_rows.tsv for the Kalign rows (default: A2_DIR/refit)")
    a = ap.parse_args()
    for pool in ("protein", "nucleotide"):
        overall(a.a2, pool, a.pooled or a.a3)
    bali(a.a2, a.externals, a.a3, a.oof)


if __name__ == "__main__":
    main()
