"""
Kalign - Fast multiple sequence alignment

This package provides Python bindings for the Kalign multiple sequence alignment program.
Kalign is a fast and accurate multiple sequence alignment tool for biological sequences.
"""

import os
import threading
from importlib import import_module
from typing import Any, List, Literal, Optional, Union

from . import _core, io, utils

__version__ = "3.4.5"
__author__ = "Timo Lassmann"
__email__ = "timo.lassmann@telethonkids.org.au"

# Re-export constants for convenience
DNA = _core.DNA
DNA_INTERNAL = _core.DNA_INTERNAL
RNA = _core.RNA
PROTEIN = _core.PROTEIN
PROTEIN_DIVERGENT = _core.PROTEIN_DIVERGENT
AUTO = _core.AUTO

# Global thread control
_thread_local = threading.local()
_default_threads = 1


def align(
    sequences: List[str],
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: Optional[int] = None,
    fmt: Literal["plain", "biopython", "skbio"] = "plain",
    ids: Optional[List[str]] = None,
) -> Union[List[str], Any]:
    """
    Multiple sequence alignment via Kalign.

    Parameters
    ----------
    sequences : list of str
        List of sequences to align. Sequences should be provided as strings
        containing the sequence characters (e.g., 'ATCG' for DNA, 'ACGU' for RNA,
        or amino acid codes for proteins).
    seq_type : str or int, optional
        Sequence type specification. Can be:
        - "auto" or AUTO: Auto-detect sequence type (default)
        - "dna" or DNA: DNA sequences
        - "rna" or RNA: RNA sequences
        - "protein" or PROTEIN: Protein sequences
        - "divergent" or PROTEIN_DIVERGENT: Divergent protein sequences
        - "internal" or DNA_INTERNAL: DNA with internal gap preference
    gap_open : float, optional
        Gap opening penalty. If None, uses Kalign defaults based on sequence type.
    gap_extend : float, optional
        Gap extension penalty. If None, uses Kalign defaults based on sequence type.
    terminal_gap_extend : float, optional
        Terminal gap extension penalty. If None, uses Kalign defaults based on sequence type.
    n_threads : int, optional
        Number of threads to use for alignment. If None, uses global default.
    fmt : {'plain', 'biopython', 'skbio'}, default 'plain'
        Choose return-object flavour:
        - 'plain': list of aligned sequences (fastest)
        - 'biopython': Bio.Align.MultipleSeqAlignment object
        - 'skbio': skbio.TabularMSA object
    ids : list of str, optional
        Sequence IDs (used only for Biopython / scikit-bio objects).
        If None, generates 'seq0', 'seq1', etc.

    Returns
    -------
    list of str | Bio.Align.MultipleSeqAlignment | skbio.TabularMSA
        Aligned sequences. Return type depends on `fmt` parameter.

    Raises
    ------
    ValueError
        If input sequences are empty or invalid
    RuntimeError
        If alignment fails
    ImportError
        If Biopython or scikit-bio are requested but not installed

    Examples
    --------
    >>> import kalign
    >>> sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]

    # 1) Plain list (default)
    >>> aligned = kalign.align(sequences)
    >>> print(aligned)
    ['ATCGATCGATCG', 'ATCG-TCGATCG', 'ATCGATC-ATCG']

    # 2) Biopython object
    >>> aln_bp = kalign.align(sequences, fmt="biopython", ids=["s1","s2","s3"])
    >>> print(type(aln_bp))
    <class 'Bio.Align.MultipleSeqAlignment'>

    # 3) scikit-bio object
    >>> aln_sk = kalign.align(sequences, fmt="skbio")
    >>> print(type(aln_sk))
    <class 'skbio.alignment._tabular_msa.TabularMSA'>
    """

    # Input validation - replicate C CLI robustness
    if not sequences:
        raise ValueError("No sequences were found in the input")

    if len(sequences) < 2:
        if len(sequences) == 0:
            raise ValueError("No sequences were found in the input")
        elif len(sequences) == 1:
            raise ValueError(
                "Only 1 sequence was found in the input - at least 2 sequences are required for alignment"
            )

    if not all(isinstance(seq, str) for seq in sequences):
        raise ValueError("All sequences must be strings")

    # Check for empty or whitespace-only sequences
    empty_sequences = []
    for i, seq in enumerate(sequences):
        if not seq or not seq.strip():
            empty_sequences.append(i)

    if empty_sequences:
        if len(empty_sequences) == 1:
            raise ValueError(
                f"Sequence at index {empty_sequences[0]} is empty or contains only whitespace"
            )
        else:
            raise ValueError(
                f"Sequences at indices {empty_sequences} are empty or contain only whitespace"
            )

    # Check for valid sequence characters (basic validation)
    for i, seq in enumerate(sequences):
        # Remove common whitespace and check if anything remains
        cleaned_seq = "".join(seq.split())
        if len(cleaned_seq) == 0:
            raise ValueError(
                f"Sequence at index {i} contains only whitespace characters"
            )

        # Check for obviously invalid characters (control characters, etc)
        if any(ord(char) < 32 for char in cleaned_seq if char not in "\t\n\r"):
            raise ValueError(
                f"Sequence at index {i} contains invalid control characters"
            )

        # Check for digits and other problematic characters that cause platform-specific segfaults
        invalid_chars = set(char for char in cleaned_seq if char.isdigit())
        if invalid_chars:
            raise ValueError(
                f"Sequence at index {i} contains invalid characters: {sorted(invalid_chars)}. "
                f"Sequences should only contain valid biological sequence characters."
            )

    # Warn about very short sequences (like C CLI warnings)
    very_short_sequences = [
        i for i, seq in enumerate(sequences) if len(seq.strip()) < 3
    ]
    if very_short_sequences and len(very_short_sequences) > len(sequences) * 0.5:
        import warnings

        warnings.warn(
            f"Many sequences are very short (< 3 characters). This may affect alignment quality.",
            UserWarning,
            stacklevel=2,
        )

    # Convert string sequence types to integers
    seq_type_map = {
        "auto": AUTO,
        "dna": DNA,
        "rna": RNA,
        "protein": PROTEIN,
        "divergent": PROTEIN_DIVERGENT,
        "internal": DNA_INTERNAL,
    }

    if isinstance(seq_type, str):
        seq_type_lower = seq_type.lower()
        if seq_type_lower not in seq_type_map:
            raise ValueError(
                f"Invalid seq_type: {seq_type}. Must be one of: {list(seq_type_map.keys())}"
            )
        seq_type_int = seq_type_map[seq_type_lower]
    else:
        seq_type_int = seq_type

    # Parameter validation and defaults
    if gap_open is None:
        gap_open = -1.0
    elif not isinstance(gap_open, (int, float)):
        raise ValueError("gap_open must be a number")
    elif gap_open > 0:
        raise ValueError("gap_open must be negative or zero (penalty values)")

    if gap_extend is None:
        gap_extend = -1.0
    elif not isinstance(gap_extend, (int, float)):
        raise ValueError("gap_extend must be a number")
    elif gap_extend > 0:
        raise ValueError("gap_extend must be negative or zero (penalty values)")

    if terminal_gap_extend is None:
        terminal_gap_extend = -1.0
    elif not isinstance(terminal_gap_extend, (int, float)):
        raise ValueError("terminal_gap_extend must be a number")
    elif terminal_gap_extend > 0:
        raise ValueError(
            "terminal_gap_extend must be negative or zero (penalty values)"
        )

    # Handle thread count
    if n_threads is None:
        n_threads = get_num_threads()
    elif not isinstance(n_threads, int):
        raise ValueError("n_threads must be an integer")
    elif n_threads < 1:
        raise ValueError("n_threads must be at least 1")
    elif n_threads > 1024:  # reasonable upper limit
        import warnings

        warnings.warn(
            f"Very high thread count ({n_threads}) may not improve performance and could cause issues.",
            UserWarning,
            stacklevel=2,
        )

    # Call the C++ binding for core alignment
    try:
        aligned_seqs = _core.align(
            sequences,
            seq_type_int,
            gap_open,
            gap_extend,
            terminal_gap_extend,
            n_threads,
        )
    except Exception as e:
        raise RuntimeError(f"Alignment failed: {str(e)}")

    # Validate IDs if provided (applies to all formats)
    if ids is not None and len(ids) != len(aligned_seqs):
        raise ValueError(
            f"Number of IDs ({len(ids)}) must match number of sequences ({len(aligned_seqs)})"
        )

    # Handle different return formats
    if fmt == "plain":
        return aligned_seqs

    # Generate IDs if not provided (only for ecosystem formats)
    if ids is None:
        ids = [f"seq{i}" for i in range(len(aligned_seqs))]

    if fmt == "biopython":
        try:
            MultipleSeqAlignment = import_module("Bio.Align").MultipleSeqAlignment
            SeqRecord = import_module("Bio.SeqRecord").SeqRecord
            Seq = import_module("Bio.Seq").Seq
        except ModuleNotFoundError as e:
            raise ImportError(
                "Biopython not installed. Run: pip install kalign[biopython]"
            ) from e
        return MultipleSeqAlignment(
            [SeqRecord(Seq(s), id=i) for s, i in zip(aligned_seqs, ids)]
        )

    if fmt == "skbio":
        try:
            TabularMSA = import_module("skbio").TabularMSA
            SkbioDNA = import_module("skbio.sequence").DNA
        except ModuleNotFoundError as e:
            raise ImportError(
                "scikit-bio not installed. Run: pip install kalign[skbio]"
            ) from e
        return TabularMSA(
            [SkbioDNA(s, metadata={"id": i}) for s, i in zip(aligned_seqs, ids)]
        )

    raise ValueError(f"Unknown fmt='{fmt}' (expected 'plain', 'biopython', 'skbio')")


def set_num_threads(n: int) -> None:
    """
    Set the default number of threads for alignment operations.

    This affects all future calls to align() that don't explicitly specify n_threads.
    The setting is thread-local, so different threads can have different defaults.

    Parameters
    ----------
    n : int
        Number of threads to use. Must be at least 1.

    Raises
    ------
    ValueError
        If n is less than 1

    Examples
    --------
    >>> import kalign
    >>> kalign.set_num_threads(4)
    >>> aligned = kalign.align(sequences)  # Uses 4 threads
    """
    global _default_threads
    if n < 1:
        raise ValueError("Number of threads must be at least 1")

    # Use thread-local storage for thread safety
    _thread_local.num_threads = n
    # Also update global default for new threads
    _default_threads = n


def get_num_threads() -> int:
    """
    Get the current default number of threads for alignment operations.

    Returns
    -------
    int
        Current default number of threads

    Examples
    --------
    >>> import kalign
    >>> kalign.get_num_threads()
    1
    >>> kalign.set_num_threads(8)
    >>> kalign.get_num_threads()
    8
    """
    return getattr(_thread_local, "num_threads", _default_threads)


def align_from_file(
    input_file: str,
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: int = 1,
) -> List[str]:
    """
    Align sequences from a file using Kalign.

    Parameters
    ----------
    input_file : str
        Path to input file containing sequences. Supported formats:
        FASTA, MSF, Clustal, aligned FASTA.
    seq_type : str or int, optional
        Sequence type specification (same as align function)
    gap_open : float, optional
        Gap opening penalty
    gap_extend : float, optional
        Gap extension penalty
    terminal_gap_extend : float, optional
        Terminal gap extension penalty
    n_threads : int, optional
        Number of threads to use for alignment

    Returns
    -------
    list of str
        List of aligned sequences

    Raises
    ------
    FileNotFoundError
        If input file doesn't exist
    RuntimeError
        If alignment fails
    """

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Convert string sequence types to integers
    seq_type_map = {
        "auto": AUTO,
        "dna": DNA,
        "rna": RNA,
        "protein": PROTEIN,
        "divergent": PROTEIN_DIVERGENT,
        "internal": DNA_INTERNAL,
    }

    if isinstance(seq_type, str):
        seq_type_lower = seq_type.lower()
        if seq_type_lower not in seq_type_map:
            raise ValueError(
                f"Invalid seq_type: {seq_type}. Must be one of: {list(seq_type_map.keys())}"
            )
        seq_type_int = seq_type_map[seq_type_lower]
    else:
        seq_type_int = seq_type

    # Use defaults if not specified
    if gap_open is None:
        gap_open = -1.0
    if gap_extend is None:
        gap_extend = -1.0
    if terminal_gap_extend is None:
        terminal_gap_extend = -1.0

    # Validate thread count
    if n_threads < 1:
        raise ValueError("n_threads must be at least 1")

    # Call the C++ binding
    try:
        aligned_seqs = _core.align_from_file(
            input_file,
            seq_type_int,
            gap_open,
            gap_extend,
            terminal_gap_extend,
            n_threads,
        )
        return aligned_seqs
    except Exception as e:
        raise RuntimeError(f"Alignment failed: {str(e)}")


def write_alignment(
    sequences: List[str],
    output_file: str,
    format: str = "fasta",
    ids: Optional[List[str]] = None,
) -> None:
    """
    Write aligned sequences to a file.

    Parameters
    ----------
    sequences : list of str
        List of aligned sequences
    output_file : str
        Path to output file
    format : str, optional
        Output format: "fasta", "clustal", "stockholm", "phylip" (default: "fasta")
    ids : list of str, optional
        Sequence IDs. If None, generates seq0, seq1, etc.

    Raises
    ------
    ValueError
        If invalid format or empty sequence list
    ImportError
        If Biopython is not installed for non-FASTA formats

    Examples
    --------
    >>> aligned = kalign.align(sequences)
    >>> kalign.write_alignment(aligned, "output.fasta")
    >>> kalign.write_alignment(aligned, "output.aln", format="clustal", ids=["seq1", "seq2"])
    """

    if not sequences:
        raise ValueError("Empty sequence list provided")

    format_lower = format.lower()

    # Map format aliases
    format_map = {
        "fasta": "fasta",
        "fa": "fasta",
        "clustal": "clustal",
        "aln": "clustal",
        "stockholm": "stockholm",
        "sto": "stockholm",
        "phylip": "phylip",
        "phy": "phylip",
    }

    if format_lower not in format_map:
        raise ValueError(
            f"Invalid format: {format}. Must be one of: fasta, clustal, stockholm, phylip"
        )

    mapped_format = format_map[format_lower]

    # Use appropriate writer from io module (lazy import to avoid circular imports)
    from . import io

    if mapped_format == "fasta":
        io.write_fasta(sequences, output_file, ids=ids)
    elif mapped_format == "clustal":
        io.write_clustal(sequences, output_file, ids=ids)
    elif mapped_format == "stockholm":
        io.write_stockholm(sequences, output_file, ids=ids)
    elif mapped_format == "phylip":
        io.write_phylip(sequences, output_file, ids=ids)


def generate_test_sequences(
    n_seq: int, n_obs: int, dna: bool, length: int, seed: int = 42
) -> List[str]:
    """
    Generate test sequences using DSSim HMM-based simulator.

    This function uses the DSSim (Dynamic Sequence Simulator) from the Kalign
    test suite to generate realistic evolutionary sequence data for testing
    and benchmarking purposes.

    Parameters
    ----------
    n_seq : int
        Number of sequences to generate
    n_obs : int
        Number of observed sequences for training the HMM (typically 20-50)
    dna : bool
        True to generate DNA sequences, False to generate protein sequences
    length : int
        Target sequence length
    seed : int, optional
        Random seed for reproducible results (default: 42)

    Returns
    -------
    list of str
        List of generated sequences

    Raises
    ------
    RuntimeError
        If sequence generation fails

    Examples
    --------
    >>> import kalign
    >>> # Generate 100 DNA sequences of length 200
    >>> dna_seqs = kalign.generate_test_sequences(100, 20, True, 200)
    >>> len(dna_seqs)
    100
    >>> len(dna_seqs[0])
    200

    >>> # Generate 50 protein sequences of length 150
    >>> protein_seqs = kalign.generate_test_sequences(50, 30, False, 150)
    >>> len(protein_seqs)
    50
    """

    if n_seq < 1:
        raise ValueError("n_seq must be at least 1")
    if n_obs < 1:
        raise ValueError("n_obs must be at least 1")
    if length < 1:
        raise ValueError("length must be at least 1")

    try:
        sequences = _core.generate_test_sequences(n_seq, n_obs, dna, length, seed)
        return sequences
    except Exception as e:
        raise RuntimeError(f"Test sequence generation failed: {str(e)}")


# Convenience aliases
kalign = align  # For backward compatibility or alternative naming


__all__ = [
    "align",
    "align_from_file",
    "write_alignment",
    "generate_test_sequences",
    "set_num_threads",
    "get_num_threads",
    "kalign",
    "DNA",
    "DNA_INTERNAL",
    "RNA",
    "PROTEIN",
    "PROTEIN_DIVERGENT",
    "AUTO",
    "__version__",
    "__author__",
    "__email__",
    "io",
    "utils",
]
