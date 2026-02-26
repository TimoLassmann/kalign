"""Comprehensive BAliBASE benchmark comparing all kalign modes.

Modes:
  1. baseline        - no VSM, no refinement
  2. +vsm            - VSM only (vsm_amax=2.0)
  3. +vsm+ref        - VSM + refinement
  4. ens3            - ensemble(3), no VSM, no refinement in runs
  5. ens3+vsm        - ensemble(3), VSM in each run
  6. ens3+vsm+ref    - ensemble(3), VSM + refinement in each run

Usage:
  uv run python -m benchmarks.vsm_ensemble_experiment
  uv run python -m benchmarks.vsm_ensemble_experiment --max-cases 10
  uv run python -m benchmarks.vsm_ensemble_experiment --categories RV11
"""

import argparse
import json
import statistics
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import kalign
from .datasets import get_cases
from .scoring import parse_balibase_xml


def _read_fasta_seqs(path):
    seqs = []
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def _alignment_stats(path):
    seqs = _read_fasta_seqs(path)
    if not seqs:
        return {}
    alnlen = len(seqs[0][1])
    nseq = len(seqs)
    total_chars = 0
    total_gaps = 0
    total_gap_opens = 0
    for _, s in seqs:
        in_gap = False
        for ch in s:
            if ch == "-":
                total_gaps += 1
                if not in_gap:
                    total_gap_opens += 1
                    in_gap = True
            else:
                total_chars += 1
                in_gap = False
    gap_frac = total_gaps / (nseq * alnlen) if alnlen > 0 else 0
    mean_seqlen = total_chars / nseq if nseq > 0 else 0
    alnlen_ratio = alnlen / mean_seqlen if mean_seqlen > 0 else 0
    return {
        "alnlen": alnlen,
        "nseq": nseq,
        "gap_frac": gap_frac,
        "gap_opens_per_seq": total_gap_opens / nseq if nseq > 0 else 0,
        "alnlen_ratio": alnlen_ratio,
        "mean_seqlen": mean_seqlen,
    }


def _score_case(case, output_path):
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(str(case.reference), str(output_path), column_mask=mask)
    return kalign.compare_detailed(str(case.reference), str(output_path))


CONFIGS = [
    {
        "label": "baseline",
        "refine": "none",
        "vsm_amax": 0.0,  # explicitly disable VSM
    },
    {
        "label": "+vsm",
        "refine": "none",
        "vsm_amax": 2.0,
    },
    {
        "label": "+vsm+ref",
        "refine": "confident",
        "vsm_amax": 2.0,
    },
    {
        "label": "+vsm+iref",
        "refine": "inline",
        "vsm_amax": 2.0,
    },
    {
        "label": "ens3",
        "ensemble": 3,
        "refine": "none",
        "vsm_amax": 0.0,  # no VSM in individual runs
    },
    {
        "label": "ens3+vsm",
        "ensemble": 3,
        "refine": "none",
        "vsm_amax": 2.0,
    },
    {
        "label": "ens3+vsm+ref",
        "ensemble": 3,
        "refine": "confident",
        "vsm_amax": 2.0,
    },
    {
        "label": "ens3+vsm+ref+ra1",
        "ensemble": 3,
        "refine": "confident",
        "vsm_amax": 2.0,
        "realign": 1,
    },
    {
        "label": "ens3+vsm+ref+ra2",
        "ensemble": 3,
        "refine": "confident",
        "vsm_amax": 2.0,
        "realign": 2,
    },
    {
        "label": "ens3+vsm+iref",
        "ensemble": 3,
        "refine": "inline",
        "vsm_amax": 2.0,
    },
    {
        "label": "ens3+vsm+iref+ra1",
        "ensemble": 3,
        "refine": "inline",
        "vsm_amax": 2.0,
        "realign": 1,
    },
]


def _run_one(case, config):
    label = config["label"]
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        t0 = time.perf_counter()
        kwargs = dict(format="fasta", seq_type=case.seq_type)

        # Alignment mode
        if config.get("ensemble"):
            kwargs["ensemble"] = config["ensemble"]
            kwargs["refine"] = config.get("refine", "none")
            kwargs["vsm_amax"] = config.get("vsm_amax", -1.0)
            if config.get("realign"):
                kwargs["realign"] = config["realign"]
        else:
            # Standard kalign
            kwargs["refine"] = config.get("refine", "none")
            kwargs["vsm_amax"] = config.get("vsm_amax", -1.0)

        kalign.align_file_to_file(str(case.unaligned), tmp_path, **kwargs)
        wall = time.perf_counter() - t0
        scores = _score_case(case, tmp_path)
        stats = _alignment_stats(tmp_path)
        return {
            "family": case.family, "dataset": case.dataset,
            "method": label,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
            "wall_time": wall,
            **stats,
        }
    except Exception as e:
        return {
            "family": case.family, "dataset": case.dataset,
            "method": label,
            "recall": 0, "precision": 0, "f1": 0, "tc": 0,
            "wall_time": 0,
            "error": str(e),
        }
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _print_summary(results, method_names):
    print(f"\n{'Method':>20} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}"
          f" {'GapFrac':>8} {'AlnRatio':>9} {'Time':>7}")
    print("-" * 90)
    groups = defaultdict(list)
    for r in results:
        groups[r["method"]].append(r)
    for method in method_names:
        entries = groups.get(method, [])
        valid = [r for r in entries if "error" not in r]
        errs = len(entries) - len(valid)
        if not valid:
            print(f"{method:>20}  (no results)")
            continue
        rec = statistics.mean(r["recall"] for r in valid)
        prec = statistics.mean(r["precision"] for r in valid)
        f1 = statistics.mean(r["f1"] for r in valid)
        tc = statistics.mean(r["tc"] for r in valid)
        gf = statistics.mean(r.get("gap_frac", 0) for r in valid)
        ar = statistics.mean(r.get("alnlen_ratio", 0) for r in valid)
        wt = sum(r.get("wall_time", 0) for r in valid)
        suffix = f" ({errs} err)" if errs else ""
        print(f"{method:>20} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}"
              f" {gf:>8.3f} {ar:>9.2f} {wt:>6.1f}s{suffix}")


def _print_per_category(results, method_names):
    all_cats = sorted({r["dataset"].replace("balibase_", "") for r in results})
    for cat in all_cats:
        cat_results = [r for r in results if cat in r["dataset"]]
        cat_groups = defaultdict(list)
        for r in cat_results:
            cat_groups[r["method"]].append(r)

        n = len(cat_groups.get(method_names[0], []))
        print(f"\n=== {cat} ({n} cases) ===")
        print(f"{'Method':>20} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}")
        print("-" * 55)
        for method in method_names:
            entries = [r for r in cat_groups.get(method, []) if "error" not in r]
            if not entries:
                continue
            rec = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            print(f"{method:>20} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}")


def main():
    parser = argparse.ArgumentParser(description="Comprehensive BAliBASE benchmark")
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--categories", nargs="*", default=None)
    parser.add_argument("--configs", nargs="*", default=None,
                        help="Only run specific configs by label")
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)
    if args.categories:
        cats = [c.upper() for c in args.categories]
        cases = [c for c in cases if any(cat in c.dataset.upper() for cat in cats)]

    configs = CONFIGS
    if args.configs:
        configs = [c for c in CONFIGS if c["label"] in args.configs]

    print(f"{len(cases)} BAliBASE cases, {len(configs)} configs", flush=True)

    method_names = [c["label"] for c in configs]
    n_tasks = len(configs) * len(cases)
    print(f"{n_tasks} total tasks (sequential)", flush=True)

    t0 = time.perf_counter()
    results = []
    done = 0

    for cfg in configs:
        print(f"\n  Config: {cfg['label']}", flush=True)
        for case in cases:
            r = _run_one(case, cfg)
            results.append(r)
            done += 1
            if done % 20 == 0:
                elapsed = time.perf_counter() - t0
                eta = elapsed / done * (n_tasks - done)
                print(f"    {done}/{n_tasks} ({elapsed:.0f}s, ETA {eta:.0f}s)", flush=True)

    elapsed = time.perf_counter() - t0
    print(f"\nAll done in {elapsed:.0f}s")

    _print_summary(results, method_names)
    _print_per_category(results, method_names)

    # Save
    out = Path("benchmarks/data/vsm_ensemble_experiment.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
