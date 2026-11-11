"""Report tables of a finished run_all.sh run.

    python -m recipe.report --out WS/runs/full      -> <out>/report/report.md and <out>/report/tables/*.tsv

Reads only files that earlier steps wrote (results.json / TSVs of recipe.pool_results, decision.json of
recipe.decide, fit.py's run.json / selection.json / hv.json, the progress log, the freeze records) and does no
computation beyond formatting, so every number in the report is traceable to one of those files.
Tables (per pool):
  1  pooled out-of-fold quality, benchmark-weighted (the fit objective) and case-weighted, 95% BCa CI:
     new presets, 3.7.0 presets, MAFFT L-INS-i, MUSCLE 5, Clustal Omega
  2  paired differences new - 3.7.0 with the decision rule, budget check and decision
  3  the same differences split by homolog_in_train (has a >=50% identity link into another fold), by benchmark
     and by difficulty
  4  CPU per case: fitting-concurrency measurement and 4-worker re-timing
  5  selected preset per fold (core scores) and search effort (generations, stop reason)
  6  optimism: final-fit in-sample minus pooled out-of-fold
  7  provenance: freeze commit/tag, build id, data hashes, step wall times
A table or column that is missing in the inputs is reported as "n/a", never filled in.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

POOLS = ("protein", "nucleotide")
MODES = ("fast", "default", "recall", "accurate")
EXT_NAMES = {"mafft_linsi": "MAFFT L-INS-i", "muscle5": "MUSCLE 5", "clustalo": "Clustal Omega"}


def f(x, nd=3, sign=False):
    if x in (None, ""):
        return "n/a"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"


def ci(r, k, nd=3, sign=False):
    if r is None or r.get(k) in (None, ""):
        return "n/a"
    lo, hi = r.get(f"{k}_lo"), r.get(f"{k}_hi")
    if lo in (None, "") or hi in (None, ""):
        return f(r[k], nd, sign)
    return f"{f(r[k], nd, sign)} [{f(lo, nd, sign)}, {f(hi, nd, sign)}]"


def md(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


def write_tsv(path, headers, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(headers)
        w.writerows(rows)


def jload(p):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else None


def arm_name(r):
    a = r["arm"]
    if a == "new":
        return "new (refit)"
    if a == "v370":
        return "3.7.0"
    return EXT_NAMES.get(a, a)


def arm_order(r):
    return {"new": 0, "v370": 1}.get(r["arm"], 2), MODES.index(r["mode"]) if r["mode"] in MODES else 9, r["arm"]


class Report:
    def __init__(self, out):
        self.out = Path(out)
        self.refit = self.out / "refit"
        self.dir = self.out / "report"
        self.tdir = self.dir / "tables"
        self.tdir.mkdir(parents=True, exist_ok=True)
        self.parts = []

    def add(self, text):
        self.parts.append(text)

    def table(self, name, title, headers, rows, note=""):
        write_tsv(self.tdir / f"{name}.tsv", headers, rows)
        self.add(f"### {title}\n\n" + (md(headers, rows) if rows else "n/a (no data)\n") +
                 (f"\n{note}\n" if note else "") + "\n")

    # ---------------------------------------------------------------------------------- provenance
    def provenance(self):
        fz1, fz2 = jload(self.out / "freeze1.json"), jload(self.out / "freeze2.json")
        dh = jload(self.out / "data_hashes.json")
        self.add("## Provenance\n")
        lines = []
        if fz1:
            lines.append(f"- Freeze 1: commit `{fz1['commit']}`" + (f", tag `{fz1['tag']}`" if fz1.get("tag") else "") +
                         f", {fz1['n_files']} frozen files, manifest recorded {fz1['recorded']}"
                         + (" (WORKTREE SNAPSHOT, not a commit)" if fz1.get("worktree_snapshot") else ""))
        else:
            lines.append("- Freeze 1: NOT RECORDED (this run was started without a freeze manifest: smoke only)")
        if fz2:
            lines.append(f"- Freeze 2: out-of-fold rows of {len(fz2['oof_rows_sha256'])} fold files hashed "
                         f"{fz2['recorded']}, before the first final-fit step; verified unchanged at the pooling step")
        if dh:
            lines.append(f"- Data tables hashed at the start: {len(dh['files'])} files ({dh['recorded']}); a resumed run "
                         "refuses to continue if any changed")
        for pool in POOLS:
            dec = jload(self.refit / pool / "decision.json")
            if dec:
                lines.append(f"- {pool}: aligner build id `{dec.get('build_id')}`, settings "
                             f"{'SMOKE (non-registered)' if dec.get('smoke') else 'registered (full)'}")
        self.add("\n".join(lines) + "\n\n")
        log = self.out / "progress.log"
        rows = []
        if log.exists():
            for l in log.read_text().splitlines():
                p = l.split(" ", 2)
                if len(p) == 3 and p[1] in ("DONE", "FAIL", "FREEZE2", "STOP"):
                    rows.append((p[0], p[1], p[2]))
        self.table("steps", "Step log", ["time", "event", "step"], rows)

    # ---------------------------------------------------------------------------------- per pool
    def pool(self, pool):
        res = jload(self.refit / pool / "results" / "results.json")
        dec = jload(self.refit / pool / "decision.json")
        self.add(f"## {pool.capitalize()} pool\n")
        if not res:
            self.add("n/a: no pooled results found.\n\n")
            return
        tag = ("PARTIAL " if res.get("partial") else "") + ("SMOKE " if res.get("smoke") else "")
        self.add(f"{tag.strip() + ': ' if tag else ''}{res['n_cases']} cases in the pool; folds pooled "
                 f"{res['folds']}; bootstrap {res['bootstrap']['method']} {res['bootstrap']['B']} resamples, seed "
                 f"{res['bootstrap']['seed']}, {int((1 - res['bootstrap']['alpha']) * 100)}% intervals.\n\n")
        # -- 1 quality
        for key, title in (("benchmark_weighted", "benchmark-weighted (fit objective)"), ("summary", "case-weighted")):
            rows = [r for r in res[key] if r["subset_type"] == "all"]
            rows.sort(key=arm_order)
            table = [(r["mode"], arm_name(r), r["n"], r.get("n_err", ""), ci(r, "sp"), ci(r, "tc"),
                      ci(r, "precision"), ci(r, "f1"), f(r.get("cpu4"), 2), f(r.get("cpu"), 2)) for r in rows]
            self.table(f"{pool}_oof_{'bw' if key != 'summary' else 'case'}",
                       f"Table 1 ({pool}). Pooled out-of-fold quality, {title}",
                       ["mode", "arm", "n", "errors", "SP", "TC", "precision", "F1",
                        "CPU s/case, 4 workers", "CPU s/case, fitting/run concurrency"], table,
                       "SP = fraction of reference residue pairs recovered; precision = fraction of aligned pairs in "
                       "the reference; per-case means. Failed alignments count as 0. External aligners are not tuned: "
                       "they have no mode (mode '-') and their CPU was measured once at the run concurrency (30).")
        # -- 2 differences + decision
        rows = []
        diff_all = {(d["mode"], d["weighting"]): d for d in res["diff"] if d["subset_type"] == "all"}
        for m in MODES:
            d = diff_all.get((m, "benchmark"))
            c = diff_all.get((m, "case"))
            dm = (dec or {}).get("modes", {}).get(m, {})
            b = dm.get("budget", {})
            rows.append((m, ci(d, "d_sp", 4, True), ci(d, "d_tc", 4, True), ci(d, "d_precision", 4, True),
                         ci(c, "d_sp", 4, True), ci(c, "d_tc", 4, True),
                         f(b.get("ratio"), 2), ("yes" if dm.get("ships_new") else "no") if dm else "n/a",
                         dm.get("reason", "n/a")))
        self.table(f"{pool}_decision", f"Table 2 ({pool}). New - 3.7.0 out of fold and the decision",
                   ["mode", "dSP (bench-weighted)", "dTC (bench-weighted)", "dPrecision (bench-weighted)",
                    "dSP (case)", "dTC (case)", "CPU ratio new/3.7.0 (4 workers, final core)", "ships", "reason"], rows,
                   "Rule (fixed before any result): ships iff dSP >= -0.005 and dTC >= -0.005 and at least one "
                   ">= +0.005 (benchmark-weighted, pooled), final selection <= 1.2x the 3.7.0 CPU at 4 workers, and "
                   "the final selection is not simply 3.7.0.")
        # -- 3 splits
        for st, label in (("homolog_in_train", "homolog_in_train (1 = >=50% identity link into another fold)"),
                          ("benchmark", "benchmark"), ("difficulty", "difficulty bin")):
            rows = []
            for d in res["diff"]:
                if d["subset_type"] == st and d["weighting"] == "case":
                    rows.append((d["mode"], d["subset"], d["n"], ci(d, "d_sp", 4, True), ci(d, "d_tc", 4, True),
                                 ci(d, "d_precision", 4, True)))
            rows.sort(key=lambda r: (MODES.index(r[0]) if r[0] in MODES else 9, r[1]))
            self.table(f"{pool}_diff_by_{st}", f"Table 3 ({pool}). New - 3.7.0 by {label}, case-weighted",
                       ["mode", st, "n", "dSP", "dTC", "dPrecision"], rows)
        rows = [(r["mode"], arm_name(r), r["subset"], r["n"], ci(r, "sp"), ci(r, "tc"), ci(r, "precision"))
                for r in sorted((r for r in res["benchmark_weighted"] if r["subset_type"] == "homolog_in_train"),
                                key=arm_order)]
        self.table(f"{pool}_oof_by_homolog", f"Table 3b ({pool}). Out-of-fold quality by homolog_in_train, "
                   "benchmark-weighted", ["mode", "arm", "homolog_in_train", "n", "SP", "TC", "precision"], rows)
        # -- 4 CPU
        rows = []
        for r in sorted((r for r in res["summary"] if r["subset_type"] == "all"), key=arm_order):
            rows.append((r["mode"], arm_name(r), f(r.get("cpu"), 3), f(r.get("cpu4"), 3), f(r.get("cpu_total"), 1)))
        self.table(f"{pool}_cpu", f"Table 4 ({pool}). Single-thread CPU per case (s)",
                   ["mode", "arm", "mean, fitting concurrency", "mean, re-timed at 4 workers", "total (s)"], rows,
                   "Reported CPU is the 4-worker re-timing where available; the fitting-concurrency values are "
                   "5-30% higher (step-0 finding) and are used only inside the fits.")
        # -- 5 per fold selection and effort
        rows = []
        for s in res["selected"]:
            rows.append((s["fold"], s["mode"], "3.7.0" if s["is_370"] else "new",
                         "fallback" if s["fallback_370"] else "", f(s.get("core_sp")), f(s.get("core_tc")),
                         f(s.get("core_precision")), f(s.get("core_cpu"), 2), f(s.get("ref_core_cpu"), 2),
                         s.get("short", "")))
        self.table(f"{pool}_selected", f"Table 5a ({pool}). Preset selected inside each fold (training-core scores)",
                   ["fold", "mode", "source", "flag", "core SP", "core TC", "core precision", "core CPU s",
                    "3.7.0 core CPU s (budget basis)", "preset"], rows)
        rows = []
        for k in [*res["folds"], "all"]:
            d = self.refit / pool / ("final" if k == "all" else f"fold{k}")
            run = jload(d / "run.json")
            row = [k, run["n_core"] if run else "n/a", run["n_screen"] if run else "n/a"]
            for st in (1, 2):
                hv = jload(d / f"stage{st}" / "hv.json")
                row += [hv["n_gen"] if hv else "n/a", hv["reason"] if hv else "n/a",
                        f(hv["history"][-1], 4) if hv and hv.get("history") else "n/a"]
            rows.append(tuple(row))
        self.table(f"{pool}_effort", f"Table 5b ({pool}). Search effort per fit",
                   ["fold", "core cases", "screen cases", "stage 1 generations", "stage 1 stop", "stage 1 HV",
                    "stage 2 generations", "stage 2 stop", "stage 2 HV"], rows)
        # -- 6 optimism
        rows = [(o["mode"], o["weighting"], o["n"], ci(o, "opt_sp", 4, True), ci(o, "opt_tc", 4, True),
                 ci(o, "opt_precision", 4, True)) for o in res.get("optimism", [])]
        self.table(f"{pool}_optimism", f"Table 6 ({pool}). Optimism: final-fit in-sample minus pooled out-of-fold",
                   ["mode", "weighting", "n", "dSP", "dTC", "dPrecision"], rows,
                   "Positive = the final preset scores higher on the cases it was fitted on than the fold-wise "
                   "procedure does on unseen cases. Descriptive; no decision uses it.")

    def finish(self):
        head = (f"# Pooled nested-CV re-fit of the Kalign presets: report tables\n\nGenerated "
                f"{time.strftime('%Y-%m-%dT%H:%M:%S')} by `recipe.report` from `{self.out}`.\n\n")
        (self.dir / "report.md").write_text(head + "\n".join(self.parts))
        print(f"report -> {self.dir / 'report.md'} ({len(list(self.tdir.glob('*.tsv')))} tables)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    r = Report(Path(a.out))
    r.provenance()
    for p in POOLS:
        r.pool(p)
    r.finish()
    return 0


if __name__ == "__main__":
    sys.exit(main())
