"""
Shared pytest configuration and fixtures for Kalign Python tests.
"""

import pytest
import tempfile
import os
from pathlib import Path
from typing import List, Dict, Any

# Rich table support for benchmarks
try:
    from rich.console import Console
    from rich.table import Table
    from rich.text import Text

    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False


# Test data constants
DNA_SEQUENCES_SIMPLE = ["ATCGATCG", "ATCGTCG", "ATCGATCG"]

DNA_SEQUENCES_WITH_GAPS = [
    "ATCGATCGATCG",
    "ATCGTCGATCG",
    "ATCGATCATCG",
    "ATCGATCGAGATCG",
]

RNA_SEQUENCES_SIMPLE = ["AUCGAUCG", "AUCGUCG", "AUCGAUCG"]

PROTEIN_SEQUENCES_SIMPLE = ["MKTAYIAK", "MKTAYK", "MKTAYIAK"]

PROTEIN_SEQUENCES_COMPLEX = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL",
]

INVALID_DNA_SEQUENCES = [
    "ATCG123",  # Contains numbers
    "ATCGXYZ",  # Invalid characters
    "ATCG@#$",  # Special characters
]

MIXED_CASE_SEQUENCES = ["AtCgAtCg", "atcgatcg", "ATCGATCG"]


@pytest.fixture
def dna_simple():
    """Simple DNA sequences for basic testing."""
    return DNA_SEQUENCES_SIMPLE.copy()


@pytest.fixture
def dna_with_gaps():
    """DNA sequences that will likely produce gaps."""
    return DNA_SEQUENCES_WITH_GAPS.copy()


@pytest.fixture
def rna_simple():
    """Simple RNA sequences for basic testing."""
    return RNA_SEQUENCES_SIMPLE.copy()


@pytest.fixture
def protein_simple():
    """Simple protein sequences for basic testing."""
    return PROTEIN_SEQUENCES_SIMPLE.copy()


@pytest.fixture
def protein_complex():
    """Complex protein sequences for advanced testing."""
    return PROTEIN_SEQUENCES_COMPLEX.copy()


@pytest.fixture
def invalid_dna():
    """Invalid DNA sequences for error testing."""
    return INVALID_DNA_SEQUENCES.copy()


@pytest.fixture
def mixed_case():
    """Mixed case sequences for case handling testing."""
    return MIXED_CASE_SEQUENCES.copy()


@pytest.fixture
def temp_dir():
    """Create a temporary directory for file operations."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_fasta_file(temp_dir):
    """Create a sample FASTA file for testing."""
    fasta_content = """>seq1
ATCGATCGATCG
>seq2
ATCGTCGATCG
>seq3
ATCGATCATCG
"""
    fasta_file = temp_dir / "sample.fasta"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)


@pytest.fixture
def sample_protein_fasta_file(temp_dir):
    """Create a sample protein FASTA file for testing."""
    fasta_content = """>protein1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVY
>protein2
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLH
>protein3
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRL
"""
    fasta_file = temp_dir / "proteins.fasta"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)


@pytest.fixture
def invalid_fasta_file(temp_dir):
    """Create an invalid FASTA file for error testing."""
    invalid_content = """This is not a valid FASTA file
It has no proper headers
And random content
"""
    invalid_file = temp_dir / "invalid.fasta"
    invalid_file.write_text(invalid_content)
    return str(invalid_file)


@pytest.fixture
def nonexistent_file():
    """Return path to a file that doesn't exist."""
    return "/path/that/does/not/exist.fasta"


@pytest.fixture
def test_data_dir():
    """Get the test data directory path."""
    return Path(__file__).parent / "test_data"


# Parametrize fixtures for common test scenarios
@pytest.fixture(
    params=[
        ("dna", DNA_SEQUENCES_SIMPLE),
        ("rna", RNA_SEQUENCES_SIMPLE),
        ("protein", PROTEIN_SEQUENCES_SIMPLE),
    ]
)
def sequence_type_data(request):
    """Parametrized fixture providing different sequence types."""
    seq_type, sequences = request.param
    return seq_type, sequences.copy()


@pytest.fixture(params=[1, 2, 4])
def thread_count(request):
    """Parametrized fixture for testing different thread counts."""
    return request.param


@pytest.fixture(
    params=[
        ("auto", None),
        ("dna", "dna"),
        ("rna", "rna"),
        ("protein", "protein"),
        ("divergent", "divergent"),
        ("internal", "internal"),
    ]
)
def sequence_type_spec(request):
    """Parametrized fixture for testing sequence type specifications."""
    return request.param


# Pytest configuration
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line("markers", "slow: mark test as slow running")
    config.addinivalue_line("markers", "integration: mark test as integration test")
    config.addinivalue_line("markers", "performance: mark test as performance test")
    config.addinivalue_line("markers", "file_io: mark test as requiring file I/O")


def pytest_collection_modifyitems(config, items):
    """Add markers to tests based on their names."""
    for item in items:
        # Mark slow tests
        if "slow" in item.nodeid or "large" in item.nodeid:
            item.add_marker(pytest.mark.slow)

        # Mark integration tests
        if "integration" in item.nodeid:
            item.add_marker(pytest.mark.integration)

        # Mark performance tests
        if "performance" in item.nodeid or "benchmark" in item.nodeid:
            item.add_marker(pytest.mark.performance)

        # Mark file I/O tests
        if "file" in item.nodeid:
            item.add_marker(pytest.mark.file_io)


# Helper functions for tests
def assert_valid_alignment(sequences: List[str], aligned: List[str]) -> None:
    """Assert that an alignment result is valid."""
    # Same number of sequences
    assert len(aligned) == len(sequences), "Number of sequences changed"

    # All aligned sequences have same length
    if aligned:
        expected_length = len(aligned[0])
        for i, seq in enumerate(aligned):
            assert len(seq) == expected_length, f"Sequence {i} has different length"

    # No empty sequences (unless input was empty)
    for seq in aligned:
        if any(
            original for original in sequences if original
        ):  # If any input was non-empty
            assert seq, "Aligned sequence is empty"


def assert_alignment_preserves_characters(
    original: List[str], aligned: List[str]
) -> None:
    """Assert that alignment preserves the original characters (just adds gaps)."""
    for orig, align in zip(original, aligned):
        # Remove gaps from aligned sequence
        degapped = align.replace("-", "")
        assert degapped == orig, f"Original sequence {orig} not preserved in {align}"


def calculate_identity(seq1: str, seq2: str) -> float:
    """Calculate sequence identity percentage between two sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have same length")

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != "-" and b != "-")
    valid_positions = sum(1 for a, b in zip(seq1, seq2) if a != "-" and b != "-")

    return (matches / valid_positions * 100) if valid_positions > 0 else 0.0


# Expected results for regression testing
EXPECTED_ALIGNMENTS = {
    "dna_simple": {
        "sequences": DNA_SEQUENCES_SIMPLE,
        "expected_length": 8,  # This might need adjustment after testing
        "min_identity": 75.0,  # Minimum expected identity percentage
    },
    "protein_simple": {
        "sequences": PROTEIN_SEQUENCES_SIMPLE,
        "expected_length": 8,  # This might need adjustment after testing
        "min_identity": 85.0,  # Minimum expected identity percentage
    },
}


@pytest.fixture
def expected_results():
    """Expected alignment results for regression testing."""
    return EXPECTED_ALIGNMENTS.copy()


# Note: pytest_benchmark_update_json hook removed to avoid plugin validation errors
# when pytest-benchmark is not installed. If you need benchmark display functionality,
# install pytest-benchmark: pip install pytest-benchmark
