"""Verify MUMSA claim: consensus alignments at higher support thresholds
have higher precision (fraction of aligned residue pairs that are correct).

Runs kalign ensemble_consensus at min_support = 1..N for each BAliBASE case,
scores against references, and compares precision with baseline/external tools.
"""

import argparse
import json
import statistics
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kalign
from .datasets import get_cases
from .scoring import parse_balibase_xml


def score_case(case, output_path):
    """Score a test alignment against the BAliBASE reference.
    Returns dict with recall, precision, f1, tc."""
    xml_path = case.reference.with_suffix(".xml")
    if xml_path.exists():
        mask = parse_balibase_xml(xml_path)
        return kalign.compare_detailed(
            str(case.reference), str(output_path), column_mask=mask
        )
    return kalign.compare_detailed(str(case.reference), str(output_path))


def _run_one_case(case, n_runs, max_support):
    """Process a single BAliBASE case: run all support thresholds + baseline + ensemble.
    Returns list of result dicts."""
    results = []

    for ms in range(1, max_support + 1):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            kalign.align_file_to_file(
                str(case.unaligned), tmp_path, format="fasta",
                seq_type=case.seq_type, ensemble=n_runs, min_support=ms,
            )
            scores = score_case(case, tmp_path)
            results.append({
                "family": case.family, "dataset": case.dataset,
                "method": f"consensus_ms{ms}", "min_support": ms,
                "n_runs": n_runs,
                "recall": scores["recall"], "precision": scores["precision"],
                "f1": scores["f1"], "tc": scores["tc"],
            })
        except Exception as e:
            results.append({
                "family": case.family, "dataset": case.dataset,
                "method": f"consensus_ms{ms}", "min_support": ms,
                "n_runs": n_runs,
                "recall": 0, "precision": 0, "f1": 0, "tc": 0,
                "error": str(e),
            })
        finally:
            Path(tmp_path).unlink(missing_ok=True)

    # Normal ensemble (no consensus)
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type, ensemble=n_runs,
        )
        scores = score_case(case, tmp_path)
        results.append({
            "family": case.family, "dataset": case.dataset,
            "method": "ensemble", "min_support": 0, "n_runs": n_runs,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
        })
    except Exception:
        pass
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    # Baseline kalign
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        kalign.align_file_to_file(
            str(case.unaligned), tmp_path, format="fasta",
            seq_type=case.seq_type,
        )
        scores = score_case(case, tmp_path)
        results.append({
            "family": case.family, "dataset": case.dataset,
            "method": "kalign_baseline", "min_support": 0, "n_runs": 0,
            "recall": scores["recall"], "precision": scores["precision"],
            "f1": scores["f1"], "tc": scores["tc"],
        })
    except Exception:
        pass
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    return results


def _worker(args):
    """Picklable worker for ProcessPoolExecutor."""
    case, n_runs, max_support = args
    return _run_one_case(case, n_runs, max_support)


def run_consensus_sweep(cases, n_runs=8, max_support=None, categories=None, parallel=1):
    """Run ensemble consensus at each support threshold for each case."""
    if max_support is None:
        max_support = n_runs

    # Filter cases
    filtered = []
    for case in cases:
        if categories and not any(cat in case.dataset.upper() for cat in categories):
            continue
        filtered.append(case)

    all_results = []

    if parallel > 1:
        tasks = [(case, n_runs, max_support) for case in filtered]
        done = 0
        with ProcessPoolExecutor(max_workers=parallel) as pool:
            futures = {pool.submit(_worker, t): t[0].family for t in tasks}
            for future in as_completed(futures):
                done += 1
                fam = futures[future]
                try:
                    results = future.result()
                    all_results.extend(results)
                    print(f"  [{done}/{len(filtered)}] {fam} done")
                except Exception as e:
                    print(f"  [{done}/{len(filtered)}] {fam} FAILED: {e}")
    else:
        for ci, case in enumerate(filtered):
            print(f"  [{ci+1}/{len(filtered)}] {case.family} ({case.dataset})", end="", flush=True)
            results = _run_one_case(case, n_runs, max_support)
            all_results.extend(results)
            print(" done")

    return all_results


def load_external_scores(json_path="benchmarks/results/full_comparison.json"):
    """Load external tool scores (mafft, muscle, clustalo) from saved results."""
    p = Path(json_path)
    if not p.exists():
        return {}
    data = json.load(open(p))
    ext = {}
    for r in data["results"]:
        if r["method"] in ("mafft", "muscle", "clustalo"):
            key = (r["family"], r["method"])
            ext[key] = {
                "recall": r.get("recall", 0),
                "precision": r.get("precision", 0),
                "f1": r.get("f1", 0),
                "tc": r.get("tc", 0),
            }
    return ext


def summarize(results, external_scores=None):
    """Print summary table of precision across support thresholds."""
    by_method = defaultdict(list)
    by_method_cat = defaultdict(lambda: defaultdict(list))

    for r in results:
        by_method[r["method"]].append(r)
        cat = r["dataset"].replace("balibase_", "")
        by_method_cat[r["method"]][cat].append(r)

    if external_scores:
        families = {r["family"] for r in results}
        family_datasets = {r["family"]: r["dataset"] for r in results}
        for tool in ("mafft", "muscle", "clustalo"):
            for fam in families:
                key = (fam, tool)
                if key in external_scores:
                    s = external_scores[key]
                    ds = family_datasets.get(fam, "")
                    cat = ds.replace("balibase_", "")
                    entry = {"family": fam, "dataset": ds, **s}
                    by_method[tool].append(entry)
                    by_method_cat[tool][cat].append(entry)

    def method_sort_key(m):
        if m == "kalign_baseline":
            return (0, 0)
        if m.startswith("consensus_ms"):
            return (1, int(m.replace("consensus_ms", "")))
        if m == "ensemble":
            return (2, 0)
        return (3, {"mafft": 0, "muscle": 1, "clustalo": 2}.get(m, 9))

    methods = sorted(by_method.keys(), key=method_sort_key)

    print("\n=== Overall (all categories) ===")
    print(f"{'Method':<22} {'Recall':>8} {'Precision':>10} {'F1':>8} {'TC':>8}  N")
    print("-" * 70)
    for m in methods:
        entries = by_method[m]
        n = len(entries)
        rec = statistics.mean(r["recall"] for r in entries)
        prec = statistics.mean(r["precision"] for r in entries)
        f1 = statistics.mean(r["f1"] for r in entries)
        tc = statistics.mean(r["tc"] for r in entries)
        print(f"{m:<22} {rec:>8.3f} {prec:>10.3f} {f1:>8.3f} {tc:>8.3f}  {n}")

    all_cats = sorted({r["dataset"].replace("balibase_", "") for r in results})
    for cat in all_cats:
        print(f"\n=== {cat} ===")
        print(f"{'Method':<22} {'Recall':>8} {'Precision':>10} {'F1':>8} {'TC':>8}  N")
        print("-" * 70)
        for m in methods:
            entries = by_method_cat[m].get(cat, [])
            if not entries:
                continue
            n = len(entries)
            rec = statistics.mean(r["recall"] for r in entries)
            prec = statistics.mean(r["precision"] for r in entries)
            f1 = statistics.mean(r["f1"] for r in entries)
            tc = statistics.mean(r["tc"] for r in entries)
            print(f"{m:<22} {rec:>8.3f} {prec:>10.3f} {f1:>8.3f} {tc:>8.3f}  {n}")


def main():
    parser = argparse.ArgumentParser(description="MUMSA precision verification")
    parser.add_argument("--n-runs", type=int, default=8,
                        help="Number of ensemble runs (default: 8)")
    parser.add_argument("--max-support", type=int, default=None,
                        help="Maximum min_support threshold (default: n_runs)")
    parser.add_argument("--categories", nargs="*", default=None,
                        help="BAliBASE categories to test (e.g. RV11 RV12)")
    parser.add_argument("--max-cases", type=int, default=0,
                        help="Limit number of cases (0 = all)")
    parser.add_argument("--parallel", "-j", type=int, default=4,
                        help="Number of parallel workers (default: 4)")
    args = parser.parse_args()

    categories = [c.upper() for c in args.categories] if args.categories else None

    cases = get_cases("balibase", max_cases=args.max_cases if args.max_cases else None)
    print(f"Running MUMSA precision analysis on {len(cases)} BAliBASE cases")
    print(f"Ensemble runs: {args.n_runs}, max support: {args.max_support or args.n_runs}")
    print(f"Parallel workers: {args.parallel}")

    results = run_consensus_sweep(
        cases, n_runs=args.n_runs,
        max_support=args.max_support,
        categories=categories,
        parallel=args.parallel,
    )

    external = load_external_scores()
    summarize(results, external)

    # Save results for plotting
    out_dir = Path("benchmarks/results")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "mumsa_precision.json"
    # Merge external scores into a flat list alongside our results
    ext_list = []
    if external:
        families = {r["family"] for r in results}
        family_datasets = {r["family"]: r["dataset"] for r in results}
        for tool in ("mafft", "muscle", "clustalo"):
            for fam in families:
                key = (fam, tool)
                if key in external:
                    ext_list.append({
                        "family": fam, "dataset": family_datasets.get(fam, ""),
                        "method": tool, "min_support": 0, "n_runs": 0,
                        **external[key],
                    })
    json.dump({"results": results + ext_list, "n_runs": args.n_runs},
              open(out_path, "w"), indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
