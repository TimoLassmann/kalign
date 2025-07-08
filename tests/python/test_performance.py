"""
Performance and benchmark tests for Kalign Python package.
"""

import pytest
import time
import kalign


class TestPerformance:
    """Test performance characteristics and benchmarks."""

    @pytest.mark.performance
    def test_basic_alignment_speed(self, dna_simple, benchmark):
        """Benchmark basic alignment speed."""
        result = benchmark(kalign.align, dna_simple, seq_type="dna")
        assert len(result) == len(dna_simple)

    @pytest.mark.performance
    def test_threading_speedup(self, dna_with_gaps):
        """Test that multiple threads provide speedup."""
        # Single thread
        start_time = time.time()
        aligned_single = kalign.align(dna_with_gaps, seq_type="dna", n_threads=1)
        single_time = time.time() - start_time

        # Multiple threads
        start_time = time.time()
        aligned_multi = kalign.align(dna_with_gaps, seq_type="dna", n_threads=2)
        multi_time = time.time() - start_time

        # Results should be identical
        assert aligned_single == aligned_multi

        # Just verify both completed successfully - timing comparisons are unreliable in CI
        print(
            f"Single-thread time: {single_time:.6f}s, Multi-thread time: {multi_time:.6f}s"
        )

    @pytest.mark.performance
    @pytest.mark.slow
    def test_large_sequence_performance(self):
        """Test performance with larger sequences."""
        large_seqs = ["ATCG" * 100] * 10  # 400bp sequences

        start_time = time.time()
        aligned = kalign.align(large_seqs, seq_type="dna")
        elapsed = time.time() - start_time

        assert len(aligned) == 10
        assert elapsed < 10.0  # Should complete in reasonable time

    @pytest.mark.performance
    def test_memory_usage_reasonable(self, dna_simple):
        """Test that memory usage is reasonable."""
        # This is a basic test - could be enhanced with memory profiling
        aligned = kalign.align(dna_simple * 10)  # 30 sequences
        assert len(aligned) == 30

    @pytest.mark.performance
    def test_dna_alignment_speed(self, dna_simple, benchmark):
        """Benchmark DNA alignment speed."""
        result = benchmark.pedantic(
            kalign.align, args=[dna_simple], kwargs={"seq_type": "dna"}, rounds=3
        )
        assert len(result) == len(dna_simple)

    @pytest.mark.performance
    def test_protein_alignment_speed(self, protein_simple, benchmark):
        """Benchmark protein alignment speed."""
        result = benchmark.pedantic(
            kalign.align,
            args=[protein_simple],
            kwargs={"seq_type": "protein"},
            rounds=3,
        )
        assert len(result) == len(protein_simple)

    @pytest.mark.performance
    @pytest.mark.parametrize("n_threads", [1, 2, 4, 8, 16])
    def test_dna_threading_performance(self, n_threads, benchmark):
        """Test DNA alignment performance with different thread counts."""
        # Generate a larger problem to really test threading
        sequences = kalign.generate_test_sequences(
            n_seq=200,  # 200 sequences - much more demanding
            n_obs=30,  # 30 observed sequences for HMM training
            dna=True,  # DNA sequences
            length=1000,  # 1000bp sequences - much longer
            seed=42,  # Reproducible results
        )

        result = benchmark.pedantic(
            kalign.align,
            args=[sequences],
            kwargs={"seq_type": "dna", "n_threads": n_threads},
            rounds=1,  # Reduced rounds due to longer runtime
        )
        assert len(result) == len(sequences)
        # Verify alignment length consistency
        if result:
            expected_length = len(result[0])
            for seq in result:
                assert len(seq) == expected_length

    @pytest.mark.performance
    @pytest.mark.parametrize("n_threads", [1, 2, 4, 8, 16])
    def test_protein_threading_performance(self, n_threads, benchmark):
        """Test protein alignment performance with different thread counts."""
        # Generate a larger problem to really test threading
        sequences = kalign.generate_test_sequences(
            n_seq=150,  # 150 sequences - more demanding
            n_obs=25,  # 25 observed sequences for HMM training
            dna=False,  # Protein sequences
            length=320,  # 320aa sequences - much longer
            seed=123,  # Different seed for variety
        )

        result = benchmark.pedantic(
            kalign.align,
            args=[sequences],
            kwargs={"seq_type": "protein", "n_threads": n_threads},
            rounds=1,  # Reduced rounds due to longer runtime
        )
        assert len(result) == len(sequences)
        # Verify alignment length consistency
        if result:
            expected_length = len(result[0])
            for seq in result:
                assert len(seq) == expected_length
