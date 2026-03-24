"""Python-level memory stress test for kalign.

Exercises all Python API paths repeatedly to surface memory bugs
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
    """align_file_to_file with various parameter combinations."""
    import kalign
    params_list = [
        {"vsm_amax": 0.0},
        {"vsm_amax": 2.0},
        {"vsm_amax": 2.0, "refine": "confident"},
        {"vsm_amax": 2.0, "realign": 1},
        {"consistency": 3},
        {"consistency": 5, "vsm_amax": 2.0},
        {"seq_weights": 1.0, "vsm_amax": 2.0},
    ]
    for i in range(n):
        params = params_list[i % len(params_list)]
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=True) as tmp:
            kalign.align_file_to_file(str(UNALIGNED), tmp.name, **params)
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


def test_align_inmem_ensemble(n):
    """In-memory alignment with ensemble."""
    import kalign
    seqs = [
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQD",
        "MKQAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV",
        "MKTVYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEV",
    ]
    for i in range(n):
        result = kalign.align(seqs, ensemble=3, mode="fast")
        # ensemble returns (sequences, confidence)
        assert isinstance(result, tuple)
        gc.collect()


def test_align_from_file(n):
    """align_from_file loop."""
    import kalign
    for i in range(n):
        result = kalign.align_from_file(str(UNALIGNED), mode="fast")
        assert len(result.names) > 0
        gc.collect()


def test_align_from_file_ensemble(n):
    """align_from_file with ensemble."""
    import kalign
    for i in range(n):
        result = kalign.align_from_file(str(UNALIGNED), ensemble=3, mode="fast")
        assert result.column_confidence is not None
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


def test_full_benchmark_loop_precise(n):
    """Full loop with ensemble+realign (the precise/expensive path)."""
    import kalign
    for i in range(n):
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=True) as tmp:
            kalign.align_file_to_file(str(UNALIGNED), tmp.name,
                                       ensemble=3, realign=1,
                                       vsm_amax=2.0, refine="confident",
                                       mode="fast")
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
    run_test(test_align_inmem_ensemble, n)
    run_test(test_align_from_file, n)
    run_test(test_align_from_file_ensemble, n)
    run_test(test_compare, n)
    run_test(test_compare_detailed, n)
    run_test(test_compare_detailed_all_cols, n)
    run_test(test_full_benchmark_loop, n)
    run_test(test_full_benchmark_loop_precise, n)
    run_test(test_generate_sequences, n)
    run_test(test_multiple_datasets, n * 2)

    print(f"\n=== Results: {n_passed} passed, {n_failed} failed ===", flush=True)
    sys.exit(1 if n_failed > 0 else 0)


if __name__ == "__main__":
    main()
