"""Compare the laptop BAliBASE reproduction (repro_balibase.py) with the workstation re-fit rows, case by case.

    python presets/verify/compare_repro.py --repro OUT/rows.tsv --v2 V2_REFIT_DIR

Workstation rows: the final v2 presets on every case (final/insample_rows.tsv, arm "new", modes fast and
accurate) and the 3.7.0 presets (arm "v370", same file). Laptop: fresh build, original BAliBASE 3 release,
official bali_score. Reports per arm: cases matched, mean SP and TC on both sides, mean and max absolute
difference, and the share of cases within 0.001.
"""

import argparse
import csv
from pathlib import Path

ARM_MAP = {("new", "fast"): "new fast", ("new", "accurate"): "new accurate", ("v370", "fast"): "3.7.0 fast",
           ("v370", "default"): "3.7.0 default", ("v370", "accurate"): "3.7.0 accurate"}


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repro", required=True)
    ap.add_argument("--v2", required=True)
    a = ap.parse_args()
    lap = {(r["arm"], r["case"]): r for r in tsv(a.repro) if not r["error"]}
    ws = {}
    for r in tsv(Path(a.v2) / "protein" / "final" / "insample_rows.tsv"):
        arm = ARM_MAP.get((r["arm"], r["mode"]))
        if arm and r["benchmark"] == "bali3":
            ws[(arm, r["case_id"].split(".", 1)[1])] = r
    print(f"{'arm':16s} {'n':>4s} {'SP laptop':>9s} {'SP ws':>7s} {'TC laptop':>9s} {'TC ws':>7s} "
          f"{'mean|dSP|':>9s} {'max|dSP|':>8s} {'SP within .001':>14s}")
    for arm in ARM_MAP.values():
        keys = sorted(k for k in ws if k[0] == arm and k in lap)
        if not keys:
            continue
        n = len(keys)
        spl = [float(lap[k]["sp"]) for k in keys]
        spw = [0.0 if ws[k]["error"] else float(ws[k]["recall"]) for k in keys]
        tcl = [float(lap[k]["tc"]) for k in keys]
        tcw = [0.0 if ws[k]["error"] else float(ws[k]["tc"]) for k in keys]
        d = [abs(x - y) for x, y in zip(spl, spw)]
        print(f"{arm:16s} {n:4d} {sum(spl)/n:9.4f} {sum(spw)/n:7.4f} {sum(tcl)/n:9.4f} {sum(tcw)/n:7.4f} "
              f"{sum(d)/n:9.4f} {max(d):8.4f} {sum(x < 0.001 for x in d)/n:14.1%}")
    big = sorted(((abs(float(lap[k]["sp"]) - float(ws[k]["recall"] or 0)), k) for k in ws if k in lap), reverse=True)[:8]
    print("largest SP differences (laptop vs workstation):")
    for dd, k in big:
        print(f"  {k[0]:16s} {k[1]:10s} laptop {float(lap[k]['sp']):.3f}  workstation {float(ws[k]['recall']):.3f}")


if __name__ == "__main__":
    main()
