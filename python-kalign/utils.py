"""
Utility functions for working with multiple sequence alignments.

This module provides helpful utilities for analyzing and manipulating
alignments produced by Kalign.
"""

from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np


def to_array(alignment: List[str]) -> np.ndarray:
    """
    Convert alignment to NumPy character array.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences

    Returns
    -------
    numpy.ndarray
        2D character array with shape (n_sequences, alignment_length)

    Examples
    --------
    >>> import kalign
    >>> aligned = kalign.align(["ATCG", "ATCG", "ATGG"])
    >>> arr = kalign.utils.to_array(aligned)
    >>> print(arr.shape)
    (3, 4)
    """
    if not alignment:
        raise ValueError("Empty alignment provided")

    # Verify all sequences have same length
    lengths = [len(seq) for seq in alignment]
    if len(set(lengths)) > 1:
        raise ValueError("All sequences in alignment must have the same length")

    # Convert to numpy array
    return np.array([list(seq) for seq in alignment], dtype="U1")


def alignment_stats(alignment: List[str]) -> Dict[str, float]:
    """
    Calculate basic statistics for a multiple sequence alignment.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences

    Returns
    -------
    dict
        Dictionary containing alignment statistics:
        - length: alignment length
        - n_sequences: number of sequences
        - gap_fraction: fraction of positions that are gaps
        - conservation: fraction of positions that are fully conserved
        - identity: average pairwise identity

    Examples
    --------
    >>> import kalign
    >>> aligned = kalign.align(["ATCG", "ATCG", "ATGG"])
    >>> stats = kalign.utils.alignment_stats(aligned)
    >>> print(f"Conservation: {stats['conservation']:.2f}")
    Conservation: 0.75
    """
    if not alignment:
        raise ValueError("Empty alignment provided")

    # Convert to array for easier analysis
    arr = to_array(alignment)
    n_sequences, length = arr.shape

    # Calculate gap fraction
    n_gaps = np.sum(arr == "-")
    gap_fraction = n_gaps / (n_sequences * length)

    # Calculate conservation (fully conserved positions)
    conserved_positions = 0
    for col in range(length):
        column = arr[:, col]
        # Remove gaps for conservation calculation
        non_gap = column[column != "-"]
        if len(non_gap) > 0 and len(set(non_gap)) == 1:
            conserved_positions += 1

    conservation = conserved_positions / length

    # Calculate average pairwise identity
    total_comparisons = 0
    total_matches = 0

    for i in range(n_sequences):
        for j in range(i + 1, n_sequences):
            seq1, seq2 = arr[i], arr[j]
            # Only compare non-gap positions
            valid_positions = (seq1 != "-") & (seq2 != "-")
            if np.sum(valid_positions) > 0:
                matches = np.sum(seq1[valid_positions] == seq2[valid_positions])
                total_matches += matches
                total_comparisons += np.sum(valid_positions)

    identity = total_matches / total_comparisons if total_comparisons > 0 else 0.0

    return {
        "length": length,
        "n_sequences": n_sequences,
        "gap_fraction": gap_fraction,
        "conservation": conservation,
        "identity": identity,
    }


def consensus_sequence(alignment: List[str], threshold: float = 0.5) -> str:
    """
    Generate consensus sequence from alignment.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences
    threshold : float, optional
        Minimum fraction of sequences that must have the same character
        for it to be included in consensus (default: 0.5)

    Returns
    -------
    str
        Consensus sequence. Positions where no character meets the threshold
        are represented as 'N' (nucleotides) or 'X' (amino acids).

    Examples
    --------
    >>> import kalign
    >>> aligned = ["ATCG", "ATCG", "ATGG"]
    >>> consensus = kalign.utils.consensus_sequence(aligned)
    >>> print(consensus)
    ATCG
    """
    if not alignment:
        raise ValueError("Empty alignment provided")

    if not 0 <= threshold <= 1:
        raise ValueError("Threshold must be between 0 and 1")

    arr = to_array(alignment)
    n_sequences, length = arr.shape
    consensus = []

    # Detect sequence type (crude heuristic)
    all_chars = set("".join(alignment).upper().replace("-", ""))
    is_nucleotide = all_chars.issubset(set("ATCGUN"))
    ambiguous_char = "N" if is_nucleotide else "X"

    for col in range(length):
        column = arr[:, col]
        # Remove gaps for consensus calculation
        non_gap = column[column != "-"]

        if len(non_gap) == 0:
            consensus.append("-")
            continue

        # Count character frequencies
        char_counts = Counter(non_gap)
        most_common_char, count = char_counts.most_common(1)[0]

        # Check if most common character meets threshold
        if count / len(non_gap) >= threshold:
            consensus.append(most_common_char)
        else:
            consensus.append(ambiguous_char)

    return "".join(consensus)


def remove_gap_columns(alignment: List[str], threshold: float = 1.0) -> List[str]:
    """
    Remove columns with gaps from alignment.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences
    threshold : float, optional
        Remove columns where fraction of gaps >= threshold (default: 1.0)
        Default removes only all-gap columns.

    Returns
    -------
    list of str
        Alignment with gap columns removed

    Examples
    --------
    >>> aligned = ["AT-G", "AT-G", "AT-G"]
    >>> nogaps = kalign.utils.remove_gap_columns(aligned)
    >>> print(nogaps)
    ['ATG', 'ATG', 'ATG']
    """
    if not alignment:
        raise ValueError("Empty alignment provided")

    if not 0 <= threshold <= 1:
        raise ValueError("Threshold must be between 0 and 1")

    arr = to_array(alignment)
    n_sequences, length = arr.shape

    # Identify columns to keep
    keep_columns = []
    for col in range(length):
        column = arr[:, col]
        gap_fraction = np.sum(column == "-") / n_sequences
        if gap_fraction < threshold:
            keep_columns.append(col)

    # Extract kept columns
    if not keep_columns:
        return [""] * n_sequences  # All columns removed

    filtered_arr = arr[:, keep_columns]
    return ["".join(seq) for seq in filtered_arr]


def pairwise_identity_matrix(alignment: List[str]) -> np.ndarray:
    """
    Calculate pairwise identity matrix for alignment.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences

    Returns
    -------
    numpy.ndarray
        Symmetric matrix of pairwise identities (shape: n_sequences x n_sequences)

    Examples
    --------
    >>> aligned = kalign.align(["ATCG", "ATCG", "ATGG"])
    >>> matrix = kalign.utils.pairwise_identity_matrix(aligned)
    >>> print(matrix[0, 1])  # Identity between first two sequences
    1.0
    """
    if not alignment:
        raise ValueError("Empty alignment provided")

    arr = to_array(alignment)
    n_sequences = arr.shape[0]
    identity_matrix = np.zeros((n_sequences, n_sequences))

    for i in range(n_sequences):
        identity_matrix[i, i] = 1.0  # Self-identity is always 1

        for j in range(i + 1, n_sequences):
            seq1, seq2 = arr[i], arr[j]
            # Only compare non-gap positions in both sequences
            valid_positions = (seq1 != "-") & (seq2 != "-")

            if np.sum(valid_positions) > 0:
                matches = np.sum(seq1[valid_positions] == seq2[valid_positions])
                identity = matches / np.sum(valid_positions)
            else:
                identity = 0.0

            identity_matrix[i, j] = identity
            identity_matrix[j, i] = identity  # Symmetric

    return identity_matrix


def trim_alignment(
    alignment: List[str], start: Optional[int] = None, end: Optional[int] = None
) -> List[str]:
    """
    Trim alignment to specified positions.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences
    start : int, optional
        Start position (0-based, inclusive). If None, starts from beginning.
    end : int, optional
        End position (0-based, exclusive). If None, goes to end.

    Returns
    -------
    list of str
        Trimmed alignment

    Examples
    --------
    >>> aligned = ["ATCGATCG", "ATCGATCG"]
    >>> trimmed = kalign.utils.trim_alignment(aligned, start=2, end=6)
    >>> print(trimmed)
    ['CGAT', 'CGAT']
    """
    if not alignment:
        raise ValueError("Empty alignment provided")

    length = len(alignment[0])

    if start is None:
        start = 0
    if end is None:
        end = length

    if start < 0:
        start = max(0, length + start)
    if end < 0:
        end = max(0, length + end)

    if start >= end:
        raise ValueError("Start position must be less than end position")

    return [seq[start:end] for seq in alignment]
