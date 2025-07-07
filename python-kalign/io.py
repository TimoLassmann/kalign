"""
I/O helper functions for reading and writing sequence alignments.

This module provides convenient functions for reading sequences from files
and writing alignments in various formats, with optional Biopython integration.
"""

import os
from pathlib import Path
from typing import List, Optional, TextIO, Tuple, Union


def read_fasta(path: Union[str, Path]) -> List[str]:
    """
    Read FASTA file and return sequences as list of strings.

    Parameters
    ----------
    path : str or Path
        Path to FASTA file

    Returns
    -------
    list of str
        List of sequences (without headers)

    Raises
    ------
    ImportError
        If Biopython is not installed
    FileNotFoundError
        If file doesn't exist

    Examples
    --------
    >>> sequences = kalign.io.read_fasta("sequences.fasta")
    >>> aligned = kalign.align(sequences)
    """
    try:
        from Bio import SeqIO
    except ImportError as e:
        raise ImportError(
            "Biopython required for FASTA I/O. Run: pip install kalign[io]"
        ) from e

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    with open(path, "r") as handle:
        return [str(record.seq) for record in SeqIO.parse(handle, "fasta")]


def read_sequences(
    path: Union[str, Path], format: str = "auto"
) -> Tuple[List[str], List[str]]:
    """
    Read sequences from file and return sequences with their IDs.

    Parameters
    ----------
    path : str or Path
        Path to sequence file
    format : str, optional
        File format. Options: "auto", "fasta", "genbank", "embl", "swiss".
        Default "auto" attempts to detect format from file extension.

    Returns
    -------
    sequences : list of str
        List of sequences
    ids : list of str
        List of sequence IDs

    Raises
    ------
    ImportError
        If Biopython is not installed
    FileNotFoundError
        If file doesn't exist
    ValueError
        If format is unsupported or auto-detection fails

    Examples
    --------
    >>> sequences, ids = kalign.io.read_sequences("data.fasta")
    >>> aligned = kalign.align(sequences, fmt="biopython", ids=ids)
    """
    try:
        from Bio import SeqIO
    except ImportError as e:
        raise ImportError(
            "Biopython required for sequence I/O. Run: pip install kalign[io]"
        ) from e

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # Auto-detect format from extension
    if format == "auto":
        ext = path.suffix.lower()
        format_map = {
            ".fasta": "fasta",
            ".fa": "fasta",
            ".fas": "fasta",
            ".fna": "fasta",
            ".ffn": "fasta",
            ".faa": "fasta",
            ".gb": "genbank",
            ".gbk": "genbank",
            ".embl": "embl",
            ".swiss": "swiss-prot",
        }
        format = format_map.get(ext, "fasta")  # Default to FASTA

    try:
        with open(path, "r") as handle:
            records = list(SeqIO.parse(handle, format))

        sequences = [str(record.seq) for record in records]
        ids = [record.id for record in records]

        return sequences, ids

    except Exception as e:
        raise ValueError(f"Failed to parse file as {format}: {e}")


def write_fasta(
    alignment: List[str],
    path: Union[str, Path, TextIO],
    ids: Optional[List[str]] = None,
    line_length: int = 80,
) -> None:
    """
    Write aligned sequences to FASTA format.

    This function doesn't require Biopython and is used as the default
    writer for FASTA format.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences
    path : str, Path, or file-like object
        Output path or file handle
    ids : list of str, optional
        Sequence IDs. If None, generates seq0, seq1, etc.
    line_length : int, optional
        Maximum line length for sequence lines (default: 80)

    Examples
    --------
    >>> aligned = kalign.align(sequences)
    >>> kalign.io.write_fasta(aligned, "output.fasta", ids=["seq1", "seq2"])
    """
    if not alignment:
        raise ValueError("Empty alignment provided")

    if ids is None:
        ids = [f"seq{i}" for i in range(len(alignment))]
    elif len(ids) != len(alignment):
        raise ValueError(
            f"Number of IDs ({len(ids)}) must match alignment length ({len(alignment)})"
        )

    def write_to_handle(handle):
        for seq_id, seq in zip(ids, alignment):
            handle.write(f">{seq_id}\n")
            # Write sequence with line wrapping
            for i in range(0, len(seq), line_length):
                handle.write(seq[i : i + line_length] + "\n")

    if hasattr(path, "write"):
        # File-like object
        write_to_handle(path)
    else:
        # File path
        with open(path, "w") as handle:
            write_to_handle(handle)


def write_clustal(
    alignment: List[str],
    path: Union[str, Path, TextIO],
    ids: Optional[List[str]] = None,
) -> None:
    """
    Write aligned sequences to Clustal format.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences
    path : str, Path, or file-like object
        Output path or file handle
    ids : list of str, optional
        Sequence IDs. If None, generates seq0, seq1, etc.

    Examples
    --------
    >>> aligned = kalign.align(sequences)
    >>> kalign.io.write_clustal(aligned, "output.aln")
    """
    try:
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError as e:
        raise ImportError(
            "Biopython required for Clustal I/O. Run: pip install kalign[io]"
        ) from e

    if not alignment:
        raise ValueError("Empty alignment provided")

    if ids is None:
        ids = [f"seq{i}" for i in range(len(alignment))]
    elif len(ids) != len(alignment):
        raise ValueError(
            f"Number of IDs ({len(ids)}) must match alignment length ({len(alignment)})"
        )

    # Create Biopython alignment object
    records = [SeqRecord(Seq(seq), id=seq_id) for seq, seq_id in zip(alignment, ids)]
    msa = MultipleSeqAlignment(records)

    if hasattr(path, "write"):
        # File-like object
        AlignIO.write(msa, path, "clustal")
    else:
        # File path
        with open(path, "w") as handle:
            AlignIO.write(msa, handle, "clustal")


def write_stockholm(
    alignment: List[str],
    path: Union[str, Path, TextIO],
    ids: Optional[List[str]] = None,
) -> None:
    """
    Write aligned sequences to Stockholm format.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences
    path : str, Path, or file-like object
        Output path or file handle
    ids : list of str, optional
        Sequence IDs. If None, generates seq0, seq1, etc.

    Examples
    --------
    >>> aligned = kalign.align(sequences)
    >>> kalign.io.write_stockholm(aligned, "output.sto")
    """
    try:
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError as e:
        raise ImportError(
            "Biopython required for Stockholm I/O. Run: pip install kalign[io]"
        ) from e

    if not alignment:
        raise ValueError("Empty alignment provided")

    if ids is None:
        ids = [f"seq{i}" for i in range(len(alignment))]
    elif len(ids) != len(alignment):
        raise ValueError(
            f"Number of IDs ({len(ids)}) must match alignment length ({len(alignment)})"
        )

    # Create Biopython alignment object
    records = [SeqRecord(Seq(seq), id=seq_id) for seq, seq_id in zip(alignment, ids)]
    msa = MultipleSeqAlignment(records)

    if hasattr(path, "write"):
        # File-like object
        AlignIO.write(msa, path, "stockholm")
    else:
        # File path
        with open(path, "w") as handle:
            AlignIO.write(msa, handle, "stockholm")


def write_phylip(
    alignment: List[str],
    path: Union[str, Path, TextIO],
    ids: Optional[List[str]] = None,
    interleaved: bool = False,
) -> None:
    """
    Write aligned sequences to PHYLIP format.

    Parameters
    ----------
    alignment : list of str
        List of aligned sequences
    path : str, Path, or file-like object
        Output path or file handle
    ids : list of str, optional
        Sequence IDs. If None, generates seq0, seq1, etc.
    interleaved : bool, optional
        Whether to use interleaved format (default: False, sequential)

    Examples
    --------
    >>> aligned = kalign.align(sequences)
    >>> kalign.io.write_phylip(aligned, "output.phy")
    """
    try:
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError as e:
        raise ImportError(
            "Biopython required for PHYLIP I/O. Run: pip install kalign[io]"
        ) from e

    if not alignment:
        raise ValueError("Empty alignment provided")

    if ids is None:
        ids = [f"seq{i}" for i in range(len(alignment))]
    elif len(ids) != len(alignment):
        raise ValueError(
            f"Number of IDs ({len(ids)}) must match alignment length ({len(alignment)})"
        )

    # Create Biopython alignment object
    records = [SeqRecord(Seq(seq), id=seq_id) for seq, seq_id in zip(alignment, ids)]
    msa = MultipleSeqAlignment(records)

    format_name = "phylip" if interleaved else "phylip-sequential"

    if hasattr(path, "write"):
        # File-like object
        AlignIO.write(msa, path, format_name)
    else:
        # File path
        with open(path, "w") as handle:
            AlignIO.write(msa, handle, format_name)
