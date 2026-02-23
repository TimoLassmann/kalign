"""Compare ProbMSA vs baseline kalign on BAliBASE.

Runs both methods on every case, scores against references, and collects
alignment diagnostics (alignment length, gap fraction, gap-open count)
to help identify systematic issues.

WARNING: ProbMSA allocates O(N^2) matrices per pair. The largest BAliBASE
cases can use ~4.5 GB RAM each. All runs are SEQUENTIAL to avoid OOM.
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
    """Read FASTA file, return list of (name, sequence) tuples."""
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
    """Compute diagnostics for an alignment file."""
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
    # Mean ungapped length
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


def _run_one(case, config):
    label = config["label"]
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        t0 = time.perf_counter()
        kwargs = dict(
            format="fasta", seq_type=case.seq_type, refine="none",
        )
        if config.get("probmsa_5state"):
            kwargs["probmsa_5state"] = True
            for k in ("pm_delta_s", "pm_epsilon_s", "pm_delta_l", "pm_epsilon_l"):
                if k in config:
                    kwargs[k] = config[k]
        elif config.get("probmsa"):
            kwargs["probmsa"] = True
            for k in ("pm_delta", "pm_epsilon", "pm_threshold", "pm_gpo", "pm_gpe"):
                if k in config:
                    kwargs[k] = config[k]
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
    """Print overall summary table."""
    groups = defaultdict(list)
    for r in results:
        groups[r["method"]].append(r)

    print(f"\n{'Method':>20} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}"
          f" {'GapFrac':>8} {'AlnRatio':>9} {'GapOpens':>9} {'Time':>7}")
    print("-" * 100)
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
        go = statistics.mean(r.get("gap_opens_per_seq", 0) for r in valid)
        wt = sum(r.get("wall_time", 0) for r in valid)
        suffix = f" ({errs} err)" if errs else ""
        print(f"{method:>20} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}"
              f" {gf:>8.3f} {ar:>9.2f} {go:>9.1f} {wt:>6.1f}s{suffix}")


def _print_per_category(results, method_names):
    """Print per-category breakdown."""
    all_cats = sorted({r["dataset"].replace("balibase_", "") for r in results})
    for cat in all_cats:
        cat_results = [r for r in results if cat in r["dataset"]]
        cat_groups = defaultdict(list)
        for r in cat_results:
            cat_groups[r["method"]].append(r)

        n = len(cat_groups.get(method_names[0], []))
        print(f"\n=== {cat} ({n} cases) ===")
        print(f"{'Method':>20} {'Recall':>8} {'Prec':>8} {'F1':>8} {'TC':>8}"
              f" {'GapFrac':>8} {'AlnRatio':>9}")
        print("-" * 75)
        for method in method_names:
            entries = [r for r in cat_groups.get(method, []) if "error" not in r]
            if not entries:
                continue
            rec = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            gf = statistics.mean(r.get("gap_frac", 0) for r in entries)
            ar = statistics.mean(r.get("alnlen_ratio", 0) for r in entries)
            print(f"{method:>20} {rec:>8.3f} {prec:>8.3f} {f1:>8.3f} {tc:>8.3f}"
                  f" {gf:>8.3f} {ar:>9.2f}")


def _print_deltas(results, baseline_label, compare_label):
    """Print per-case F1 deltas between two methods."""
    baseline_by_fam = {r["family"]: r for r in results
                       if r["method"] == baseline_label and "error" not in r}
    compare_by_fam = {r["family"]: r for r in results
                      if r["method"] == compare_label and "error" not in r}

    deltas = []
    for fam in sorted(set(baseline_by_fam) & set(compare_by_fam)):
        b = baseline_by_fam[fam]
        p = compare_by_fam[fam]
        df1 = p["f1"] - b["f1"]
        dgf = p.get("gap_frac", 0) - b.get("gap_frac", 0)
        dar = p.get("alnlen_ratio", 0) - b.get("alnlen_ratio", 0)
        deltas.append({"family": fam, "dataset": b["dataset"],
                        "df1": df1, "dgap_frac": dgf, "dalnlen_ratio": dar,
                        "b_f1": b["f1"], "p_f1": p["f1"],
                        "b_gap_frac": b.get("gap_frac", 0),
                        "p_gap_frac": p.get("gap_frac", 0),
                        "b_alnlen_ratio": b.get("alnlen_ratio", 0),
                        "p_alnlen_ratio": p.get("alnlen_ratio", 0)})

    improved = [d for d in deltas if d["df1"] > 0.01]
    worsened = [d for d in deltas if d["df1"] < -0.01]
    print(f"\n── Per-case F1 deltas ({compare_label} - {baseline_label}) ──")
    print(f"Improved (>0.01): {len(improved)}")
    print(f"Worsened (<-0.01): {len(worsened)}")
    print(f"Neutral:  {len(deltas) - len(improved) - len(worsened)}")
    if deltas:
        df1s = [d["df1"] for d in deltas]
        print(f"Mean dF1: {statistics.mean(df1s):+.4f}")
        print(f"Median dF1: {statistics.median(df1s):+.4f}")

    deltas.sort(key=lambda d: d["df1"])
    print(f"\n  5 worst cases:")
    hdr = f"  {'Family':>15} {'Cat':>6} {'Base F1':>8} {'Test F1':>8} {'dF1':>7}"
    print(hdr)
    for d in deltas[:5]:
        cat = d["dataset"].replace("balibase_", "")
        print(f"  {d['family']:>15} {cat:>6} {d['b_f1']:>8.3f} {d['p_f1']:>8.3f} {d['df1']:>+7.3f}")

    print(f"\n  5 best cases:")
    print(hdr)
    for d in deltas[-5:][::-1]:
        cat = d["dataset"].replace("balibase_", "")
        print(f"  {d['family']:>15} {cat:>6} {d['b_f1']:>8.3f} {d['p_f1']:>8.3f} {d['df1']:>+7.3f}")


def main():
    parser = argparse.ArgumentParser(description="ProbMSA vs baseline on BAliBASE")
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--categories", nargs="*", default=None)
    parser.add_argument("--sweep", action="store_true",
                        help="Sweep MEA gap parameters instead of single comparison")
    parser.add_argument("--sweep5", action="store_true",
                        help="Sweep 5-state HMM parameters (δs, εs, δl, εl)")
    args = parser.parse_args()

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)
    if args.categories:
        cats = [c.upper() for c in args.categories]
        cases = [c for c in cases if any(cat in c.dataset.upper() for cat in cats)]

    print(f"{len(cases)} BAliBASE cases", flush=True)

    if args.sweep5:
        configs = [
            {"label": "baseline"},
            {"label": "3state",         "probmsa": True},
            # Default 5-state: δs=0.10, εs=0.50, δl=0.01, εl=0.95
            {"label": "5s_default",     "probmsa_5state": True},
            # Lower short gap-open (less frequent short gaps)
            {"label": "5s_ds03",        "probmsa_5state": True, "pm_delta_s": 0.03},
            {"label": "5s_ds05",        "probmsa_5state": True, "pm_delta_s": 0.05},
            # Higher short gap-extend (longer short gaps)
            {"label": "5s_es70",        "probmsa_5state": True, "pm_epsilon_s": 0.70},
            # Lower short gap-open + higher short gap-extend
            {"label": "5s_ds05_es70",   "probmsa_5state": True, "pm_delta_s": 0.05, "pm_epsilon_s": 0.70},
            # Lower long gap-open (rarer long gaps)
            {"label": "5s_dl005",       "probmsa_5state": True, "pm_delta_l": 0.005},
        ]
    elif args.sweep:
        configs = [
            {"label": "baseline", "probmsa": False},
            {"label": "probmsa_3s", "probmsa": True},
            {"label": "probmsa_5s", "probmsa_5state": True},
        ]

    method_names = [c["label"] for c in configs]
    n_tasks = len(configs) * len(cases)
    print(f"{len(configs)} configs x {len(cases)} cases = {n_tasks} tasks (sequential)", flush=True)

    t0 = time.perf_counter()
    results = []
    done = 0

    # Run SEQUENTIALLY — ProbMSA can use 4+ GB per case
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

    # ── Summary ──
    _print_summary(results, method_names)

    # ── Per-category ──
    _print_per_category(results, method_names)

    # ── Deltas for each probmsa config vs baseline ──
    for label in method_names[1:]:
        _print_deltas(results, "baseline", label)

    # Save
    out = Path("benchmarks/results/probmsa_experiment.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
