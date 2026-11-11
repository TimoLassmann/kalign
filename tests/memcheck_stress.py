"""Python-level memory stress test for kalign.

Exercises the Python API paths (fast and accurate modes, sequence types,
thread counts, output formats, comparison) repeatedly to surface memory bugs
in the C library or pybind11 bindings.

Usage:
    python tests/memcheck_stress.py [n_iters]
"""

import gc
import sys
import tempfile
import traceback
from pathlib import Path

DATA = Path(__file__).parent / "data"
UNALIGNED = DATA / "BB11001.tfa"
REFERENCE = DATA / "BB11001.msf"
UNALIGNED2 = DATA / "BB30014.tfa"
REFERENCE2 = DATA / "BB30014.msf"

n_passed = 0
n_failed = 0


def run_test(fn, *args, **kwargs):
    global n_passed, n_failed
    name = fn.__name__
    print(f"--- {name} ---", flush=True)
    try:
        fn(*args, **kwargs)
        n_passed += 1
        print(f"  PASSED", flush=True)
    except Exception as e:
        n_failed += 1
        print(f"  *** FAILED: {e} ***", flush=True)
        traceback.print_exc()


def test_align_file_to_file(n):
    """Basic align_file_to_file loop."""
    import kalign
    for i in range(n):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=True) as tmp:
            kalign.align_file_to_file(str(UNALIGNED), tmp.name)
        gc.collect()


def test_align_file_to_file_with_params(n):
    """align_file_to_file with parameter combinations allowed for each mode."""
    import kalign
    params_list = [
        {"mode": "fast"},
        {"mode": "fast", "gap_open": 5.5, "gap_extend": 2.0,
         "terminal_gap_extend": 1.0},
        {"mode": "fast", "seq_type": "protein", "n_threads": 2},
        {"mode": "accurate"},
        {"mode": "accurate", "n_threads": 2},
        {"mode": "accurate", "fill_in": "realign", "fill_anchor_frac": 0.5},
        {"mode": "accurate", "fill_in": "attach", "fill_support": 2},
    ]
    formats = ["fasta", "msf", "clu"]
    for i in range(n):
        params = params_list[i % len(params_list)]
        fmt = formats[i % len(formats)]
        with tempfile.NamedTemporaryFile(suffix=".out", delete=True) as tmp:
            kalign.align_file_to_file(str(UNALIGNED), tmp.name, format=fmt,
                                      **params)
        gc.collect()


def test_align_inmem(n):
    """In-memory alignment via align()."""
    import kalign
    seqs = [
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQD",
        "MKQAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV",
        "MKTVYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEV",
    ]
    for i in range(n):
        result = kalign.align(seqs, mode="fast")
        assert len(result) == len(seqs)
        gc.collect()


def test_align_inmem_accurate(n):
    """In-memory alignment with mode="accurate" (ensemble, confidence)."""
    import kalign
    seqs = [
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQD",
        "MKQAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV",
        "MKTVYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEV",
    ]
    for i in range(n):
        result = kalign.align(seqs, mode="accurate", n_threads=2)
        assert len(result) == len(seqs)
        gc.collect()


def test_align_from_file(n):
    """align_from_file loop."""
    import kalign
    for i in range(n):
        result = kalign.align_from_file(str(UNALIGNED), mode="fast")
        assert len(result.names) > 0
        gc.collect()


def test_align_from_file_accurate(n):
    """align_from_file with mode="accurate"."""
    import kalign
    for i in range(n):
        result = kalign.align_from_file(str(UNALIGNED), mode="accurate")
        assert result.column_confidence is not None
        gc.collect()


def test_align_dna_types(n):
    """DNA and RNA input with explicit and detected sequence types."""
    import kalign
    seqs = [
        "ACGTACGTTAGCTAGCTAGGATCGATCGATTAGC",
        "ACGTACGTTAGCTAGGATCGATCGATTAGC",
        "ACGTTCGTTAGCTAGCTAGGATCGATCGTTAGC",
        "ACGTACGTTAGCTAGCTAGGATCGATTAGC",
    ]
    rna = [x.replace("T", "U") for x in seqs]
    for i in range(n):
        for mode in ("fast", "accurate"):
            assert len(kalign.align(seqs, seq_type="auto", mode=mode)) == 4
            assert len(kalign.align(seqs, seq_type="dna", mode=mode)) == 4
            assert len(kalign.align(rna, seq_type="rna", mode=mode)) == 4
        gc.collect()


def test_compare(n):
    """Repeated MSA comparison."""
    import kalign
    for i in range(n):
        score = kalign.compare(str(REFERENCE), str(REFERENCE))
        assert score > 0
        gc.collect()


def test_compare_detailed(n):
    """Repeated detailed MSA comparison."""
    import kalign
    for i in range(n):
        result = kalign.compare_detailed(str(REFERENCE), str(REFERENCE))
        assert result["recall"] > 0
        gc.collect()


def test_compare_detailed_all_cols(n):
    """Repeated detailed comparison with max_gap_frac=-1.0."""
    import kalign
    for i in range(n):
        result = kalign.compare_detailed(str(REFERENCE), str(REFERENCE),
                                          max_gap_frac=-1.0)
        assert result["recall"] > 0
        gc.collect()


def test_full_benchmark_loop(n):
    """Simulates what the benchmark does: align + compare + compare_detailed."""
    import kalign
    for i in range(n):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=True) as tmp:
            kalign.align_file_to_file(str(UNALIGNED), tmp.name, mode="fast")
            sp = kalign.compare(str(REFERENCE), tmp.name)
            detailed = kalign.compare_detailed(str(REFERENCE), tmp.name)
        gc.collect()


def test_full_benchmark_loop_accurate(n):
    """Full loop with mode="accurate" (the expensive path)."""
    import kalign
    for i in range(n):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=True) as tmp:
            kalign.align_file_to_file(str(UNALIGNED), tmp.name,
                                       mode="accurate", n_threads=2)
            sp = kalign.compare(str(REFERENCE), tmp.name)
            detailed = kalign.compare_detailed(str(REFERENCE), tmp.name)
        gc.collect()


def test_generate_sequences(n):
    """Repeated sequence generation."""
    import kalign
    for i in range(n):
        seqs = kalign.generate_test_sequences(20, 10, False, 100, seed=i)
        assert len(seqs) == 20
        gc.collect()


def test_multiple_datasets(n):
    """Alternate between different datasets (like a real benchmark)."""
    import kalign
    datasets = [
        (UNALIGNED, REFERENCE),
        (UNALIGNED2, REFERENCE2),
    ]
    for i in range(n):
        unaln, ref = datasets[i % len(datasets)]
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=True) as tmp:
            kalign.align_file_to_file(str(unaln), tmp.name, mode="fast")
            sp = kalign.compare(str(ref), tmp.name)
            detailed = kalign.compare_detailed(str(ref), tmp.name)
        gc.collect()


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    print(f"=== Python Memory Stress Test ({n} iters per test) ===\n", flush=True)

    run_test(test_align_file_to_file, n)
    run_test(test_align_file_to_file_with_params, n * 2)
    run_test(test_align_inmem, n)
    run_test(test_align_inmem_accurate, n)
    run_test(test_align_from_file, n)
    run_test(test_align_from_file_accurate, n)
    run_test(test_align_dna_types, n)
    run_test(test_compare, n)
    run_test(test_compare_detailed, n)
    run_test(test_compare_detailed_all_cols, n)
    run_test(test_full_benchmark_loop, n)
    run_test(test_full_benchmark_loop_accurate, n)
    run_test(test_generate_sequences, n)
    run_test(test_multiple_datasets, n * 2)

    print(f"\n=== Results: {n_passed} passed, {n_failed} failed ===", flush=True)
    sys.exit(1 if n_failed > 0 else 0)


if __name__ == "__main__":
    main()
