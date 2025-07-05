"""
Kalign - Fast multiple sequence alignment

This package provides Python bindings for the Kalign multiple sequence alignment program.
Kalign is a fast and accurate multiple sequence alignment tool for biological sequences.
"""

from typing import List, Optional, Union
import os
from . import _core

__version__ = "3.4.1"
__author__ = "Timo Lassmann"
__email__ = "timo.lassmann@telethonkids.org.au"

# Re-export constants for convenience
DNA = _core.DNA
DNA_INTERNAL = _core.DNA_INTERNAL
RNA = _core.RNA
PROTEIN = _core.PROTEIN
PROTEIN_DIVERGENT = _core.PROTEIN_DIVERGENT
AUTO = _core.AUTO


def align(
    sequences: List[str],
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: int = 1
) -> List[str]:
    """
    Align a list of sequences using Kalign.
    
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
        Number of threads to use for alignment (default: 1)
        
    Returns
    -------
    list of str
        List of aligned sequences in the same order as input sequences.
        All sequences will have the same length after alignment.
        
    Raises
    ------
    ValueError
        If input sequences are empty or invalid
    RuntimeError
        If alignment fails
        
    Examples
    --------
    >>> import kalign
    >>> sequences = [
    ...     "ATCGATCGATCG",
    ...     "ATCGTCGATCG", 
    ...     "ATCGATCATCG"
    ... ]
    >>> aligned = kalign.align(sequences, seq_type="dna")
    >>> for seq in aligned:
    ...     print(seq)
    ATCGATCGATCG
    ATCG-TCGATCG
    ATCGATC-ATCG
    
    >>> # Protein alignment with custom parameters
    >>> protein_seqs = [
    ...     "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL",
    ...     "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL"
    ... ]
    >>> aligned = kalign.align(protein_seqs, seq_type="protein", gap_open=-10.0, gap_extend=-1.0)
    """
    
    if not sequences:
        raise ValueError("Empty sequence list provided")
    
    if not all(isinstance(seq, str) for seq in sequences):
        raise ValueError("All sequences must be strings")
    
    if not all(seq.strip() for seq in sequences):
        raise ValueError("All sequences must be non-empty")
    
    # Convert string sequence types to integers
    seq_type_map = {
        "auto": AUTO,
        "dna": DNA,
        "rna": RNA,
        "protein": PROTEIN,
        "divergent": PROTEIN_DIVERGENT,
        "internal": DNA_INTERNAL
    }
    
    if isinstance(seq_type, str):
        seq_type_lower = seq_type.lower()
        if seq_type_lower not in seq_type_map:
            raise ValueError(f"Invalid seq_type: {seq_type}. Must be one of: {list(seq_type_map.keys())}")
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
        aligned_seqs = _core.align(
            sequences,
            seq_type_int,
            gap_open,
            gap_extend,
            terminal_gap_extend,
            n_threads
        )
        return aligned_seqs
    except Exception as e:
        raise RuntimeError(f"Alignment failed: {str(e)}")


def align_from_file(
    input_file: str,
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: int = 1
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
        "internal": DNA_INTERNAL
    }
    
    if isinstance(seq_type, str):
        seq_type_lower = seq_type.lower()
        if seq_type_lower not in seq_type_map:
            raise ValueError(f"Invalid seq_type: {seq_type}. Must be one of: {list(seq_type_map.keys())}")
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
            n_threads
        )
        return aligned_seqs
    except Exception as e:
        raise RuntimeError(f"Alignment failed: {str(e)}")


def write_alignment(
    sequences: List[str],
    output_file: str,
    format: str = "fasta"
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
        Output format: "fasta", "msf", or "clustal" (default: "fasta")
        
    Raises
    ------
    RuntimeError
        If writing fails
    """
    
    if not sequences:
        raise ValueError("Empty sequence list provided")
    
    valid_formats = ["fasta", "msf", "clustal"]
    if format.lower() not in valid_formats:
        raise ValueError(f"Invalid format: {format}. Must be one of: {valid_formats}")
    
    try:
        _core.write_alignment(sequences, output_file, format.lower())
    except Exception as e:
        raise RuntimeError(f"Writing alignment failed: {str(e)}")


def generate_test_sequences(
    n_seq: int,
    n_obs: int,
    dna: bool,
    length: int,
    seed: int = 42
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
    "kalign",
    "DNA",
    "DNA_INTERNAL", 
    "RNA",
    "PROTEIN",
    "PROTEIN_DIVERGENT",
    "AUTO",
    "__version__",
    "__author__",
    "__email__"
]