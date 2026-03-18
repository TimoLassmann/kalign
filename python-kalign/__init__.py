"""
Kalign - Fast multiple sequence alignment

This package provides Python bindings for the Kalign multiple sequence alignment program.
Kalign is a fast and accurate multiple sequence alignment tool for biological sequences.
"""

import os
import threading
import warnings
from importlib import import_module
from typing import Any, List, Literal, Optional, Union


class AlignedSequences:
    """Result of aligning sequences from a file, preserving sequence names.

    Supports unpacking as ``names, sequences = result`` for backward
    compatibility, while also exposing optional confidence attributes.
    """

    __slots__ = ("names", "sequences", "column_confidence", "residue_confidence")

    def __init__(
        self,
        names: List[str],
        sequences: List[str],
        column_confidence: Optional[List[float]] = None,
        residue_confidence: Optional[List[List[float]]] = None,
    ):
        self.names = names
        self.sequences = sequences
        self.column_confidence = column_confidence
        self.residue_confidence = residue_confidence

    # Backward-compatible 2-tuple unpacking: names, sequences = result
    def __iter__(self):
        return iter((self.names, self.sequences))

    def __len__(self):
        return 2

    def __getitem__(self, index):
        return (self.names, self.sequences)[index]

    def __repr__(self):
        return (
            f"AlignedSequences(names={self.names!r}, sequences={self.sequences!r}, "
            f"column_confidence={'[...]' if self.column_confidence else None}, "
            f"residue_confidence={'[...]' if self.residue_confidence else None})"
        )


from importlib.metadata import version

from . import _core, io, utils

__version__ = version("kalign-python")
__author__ = "Timo Lassmann"
__email__ = "timolassmann@icloud.com"

# Re-export constants for convenience
DNA = _core.DNA
DNA_INTERNAL = _core.DNA_INTERNAL
RNA = _core.RNA
PROTEIN = _core.PROTEIN
PROTEIN_PFASUM43 = _core.PROTEIN_PFASUM43
PROTEIN_PFASUM60 = _core.PROTEIN_PFASUM60
PROTEIN_PFASUM_AUTO = _core.PROTEIN_PFASUM_AUTO
PROTEIN_DIVERGENT = _core.PROTEIN_DIVERGENT
PROTEIN_CORBLOSUM66 = _core.PROTEIN_CORBLOSUM66
AUTO = _core.AUTO

# Refinement mode constants (used by ensemble_custom_file_to_file optimizer path)
REFINE_NONE = _core.REFINE_NONE
REFINE_ALL = _core.REFINE_ALL
REFINE_CONFIDENT = _core.REFINE_CONFIDENT
REFINE_INLINE = _core.REFINE_INLINE

# Mode constants
MODE_DEFAULT = "default"
MODE_FAST = "fast"
MODE_ACCURATE = "accurate"
MODE_PRECISE = "precise"  # deprecated alias for "accurate"

# Valid preset modes (resolved by C library)
_PRESET_MODES = {"fast", "default", "recall", "accurate"}

# Global thread control
_thread_local = threading.local()
_default_threads = 1

# Sequence type string->int mapping (shared by all entry points)
_SEQ_TYPE_MAP = {
    "auto": AUTO,
    "dna": DNA,
    "rna": RNA,
    "protein": PROTEIN,
    "pfasum43": PROTEIN_PFASUM43,
    "pfasum60": PROTEIN_PFASUM60,
    "pfasum": PROTEIN_PFASUM_AUTO,
    "divergent": PROTEIN_DIVERGENT,
    "internal": DNA_INTERNAL,
}


def _resolve_seq_type(seq_type):
    """Convert string or int seq_type to integer constant."""
    if isinstance(seq_type, int):
        return seq_type
    lower = seq_type.lower()
    if lower not in _SEQ_TYPE_MAP:
        raise ValueError(
            f"Invalid seq_type: {seq_type}. Must be one of: {list(_SEQ_TYPE_MAP.keys())}"
        )
    return _SEQ_TYPE_MAP[lower]


def _resolve_mode_name(mode):
    """Normalize mode name, handling 'precise' -> 'accurate' alias."""
    if mode is None:
        return "default"
    lower = mode.lower()
    if lower == "precise":
        warnings.warn(
            'mode="precise" is deprecated, use mode="accurate" instead.',
            DeprecationWarning,
            stacklevel=3,
        )
        return "accurate"
    if lower not in _PRESET_MODES:
        raise ValueError(
            f"Invalid mode: {mode!r}. Must be one of: 'default', 'fast', 'accurate'"
        )
    return lower


def _conf_to_pp(conf: float) -> str:
    """Convert a confidence value [0..1] to HMMER-style PP character."""
    if conf >= 0.95:
        return "*"
    return str(int(conf * 10))


def _confidence_to_pp_string(seq: str, confidences: list) -> str:
    """Convert per-residue confidence array to PP string."""
    pp = []
    for ch, conf in zip(seq, confidences):
        if ch == "-" or ch == ".":
            pp.append(".")
        else:
            pp.append(_conf_to_pp(conf))
    return "".join(pp)


def _parse_refine_mode(refine):
    """Convert string or int refine mode to integer constant."""
    if isinstance(refine, int):
        return refine
    refine_map = {
        "none": REFINE_NONE,
        "all": REFINE_ALL,
        "confident": REFINE_CONFIDENT,
        "inline": REFINE_INLINE,
    }
    refine_lower = refine.lower()
    if refine_lower not in refine_map:
        raise ValueError(
            f"Invalid refine mode: {refine}. Must be one of: {list(refine_map.keys())}"
        )
    return refine_map[refine_lower]


def _infer_skbio_type(sequences, skbio_seq):
    """Infer the appropriate skbio sequence class from raw sequence content."""
    chars = set()
    for seq in sequences:
        chars.update(seq.upper().replace("-", "").replace(".", ""))
    dna_chars = set("ACGTNRYSWKMBDHV")
    rna_chars = set("ACGUNRYSWKMBDHV")
    if "U" in chars and "T" not in chars and chars <= rna_chars:
        return skbio_seq.RNA
    if chars <= dna_chars:
        return skbio_seq.DNA
    return skbio_seq.Protein


def _validate_sequences(sequences):
    """Validate input sequences, raising on empty/invalid input."""
    if not sequences:
        raise ValueError("No sequences were found in the input")

    if len(sequences) == 1:
        raise ValueError(
            "Only 1 sequence was found in the input - at least 2 sequences are required for alignment"
        )

    if not all(isinstance(seq, str) for seq in sequences):
        raise ValueError("All sequences must be strings")

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

    for i, seq in enumerate(sequences):
        cleaned_seq = "".join(seq.split())
        if len(cleaned_seq) == 0:
            raise ValueError(
                f"Sequence at index {i} contains only whitespace characters"
            )
        if any(ord(char) < 32 for char in cleaned_seq if char not in "\t\n\r"):
            raise ValueError(
                f"Sequence at index {i} contains invalid control characters"
            )
        invalid_chars = set(char for char in cleaned_seq if char.isdigit())
        if invalid_chars:
            raise ValueError(
                f"Sequence at index {i} contains invalid characters: {sorted(invalid_chars)}. "
                f"Sequences should only contain valid biological sequence characters."
            )

    very_short_sequences = [
        i for i, seq in enumerate(sequences) if len(seq.strip()) < 3
    ]
    if very_short_sequences and len(very_short_sequences) > len(sequences) * 0.5:
        warnings.warn(
            "Many sequences are very short (< 3 characters). This may affect alignment quality.",
            UserWarning,
            stacklevel=3,
        )


def align(
    sequences: List[str],
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: Optional[int] = None,
    mode: Optional[str] = None,
    fmt: Literal["plain", "biopython", "skbio"] = "plain",
    ids: Optional[List[str]] = None,
) -> Union[List[str], Any]:
    """
    Multiple sequence alignment via Kalign.

    Parameters
    ----------
    sequences : list of str
        List of sequences to align.
    seq_type : str or int, optional
        Sequence type: "auto", "dna", "rna", "protein" (default: "auto")
    gap_open : float, optional
        Gap opening penalty. When set, mode is ignored and "fast" preset
        is used with the specified penalty.
    gap_extend : float, optional
        Gap extension penalty.
    terminal_gap_extend : float, optional
        Terminal gap extension penalty.
    n_threads : int, optional
        Number of threads (default: global setting).
    mode : str, optional
        Preset mode: "fast", "default", "accurate" (default: "default").
    fmt : {'plain', 'biopython', 'skbio'}, default 'plain'
        Return format.
    ids : list of str, optional
        Sequence IDs (for Biopython/scikit-bio formats).

    Returns
    -------
    list of str | Bio.Align.MultipleSeqAlignment | skbio.TabularMSA
    """
    _validate_sequences(sequences)

    seq_type_int = _resolve_seq_type(seq_type)

    # Parameter validation for gap penalties
    if gap_open is not None:
        if not isinstance(gap_open, (int, float)):
            raise ValueError("gap_open must be a number")
        if gap_open < 0:
            raise ValueError("gap_open must be a positive number (penalty value)")
    if gap_extend is not None:
        if not isinstance(gap_extend, (int, float)):
            raise ValueError("gap_extend must be a number")
        if gap_extend < 0:
            raise ValueError("gap_extend must be a positive number (penalty value)")
    if terminal_gap_extend is not None:
        if not isinstance(terminal_gap_extend, (int, float)):
            raise ValueError("terminal_gap_extend must be a number")
        if terminal_gap_extend < 0:
            raise ValueError(
                "terminal_gap_extend must be a positive number (penalty value)"
            )

    # Handle thread count
    if n_threads is None:
        n_threads = get_num_threads()
    elif not isinstance(n_threads, int):
        raise ValueError("n_threads must be an integer")
    elif n_threads < 1:
        raise ValueError("n_threads must be at least 1")
    elif n_threads > 1024:
        warnings.warn(
            f"Very high thread count ({n_threads}) may not improve performance.",
            UserWarning,
            stacklevel=2,
        )

    # Gap penalty override rule: if any gap penalty is set, use "fast" preset
    has_gap_override = (
        gap_open is not None or gap_extend is not None
        or terminal_gap_extend is not None
    )
    if has_gap_override:
        effective_mode = "fast"
    else:
        effective_mode = _resolve_mode_name(mode)

    confidence_data = None
    result = _core.align_mode(
        sequences, effective_mode, seq_type_int,
        gap_open if gap_open is not None else -1.0,
        gap_extend if gap_extend is not None else -1.0,
        terminal_gap_extend if terminal_gap_extend is not None else -1.0,
        n_threads,
    )
    if isinstance(result, tuple):
        aligned_seqs = result[0]
        confidence_data = result[1]
    else:
        aligned_seqs = result

    return _format_result(
        aligned_seqs, sequences, seq_type_int, confidence_data, fmt, ids
    )


def _format_result(aligned_seqs, sequences, seq_type_int, confidence_data, fmt, ids):
    """Format alignment result based on requested format."""
    if ids is not None and len(ids) != len(aligned_seqs):
        raise ValueError(
            f"Number of IDs ({len(ids)}) must match number of sequences ({len(aligned_seqs)})"
        )

    if fmt == "plain":
        return aligned_seqs

    if ids is None:
        ids = [f"seq{i}" for i in range(len(aligned_seqs))]

    if fmt == "biopython":
        try:
            MultipleSeqAlignment = import_module("Bio.Align").MultipleSeqAlignment
            SeqRecord = import_module("Bio.SeqRecord").SeqRecord
            Seq = import_module("Bio.Seq").Seq
        except ModuleNotFoundError as e:
            raise ImportError(
                "Biopython not installed. Run: pip install kalign-python[biopython]"
            ) from e

        records = []
        for idx, (s, i) in enumerate(zip(aligned_seqs, ids)):
            rec = SeqRecord(Seq(s), id=i)
            if confidence_data is not None:
                res_conf = confidence_data["residue_confidence"]
                if idx < len(res_conf) and len(res_conf[idx]) == len(s):
                    pp_str = _confidence_to_pp_string(s, res_conf[idx])
                    rec.letter_annotations["posterior_probability"] = pp_str
            records.append(rec)
        return MultipleSeqAlignment(records)

    if fmt == "skbio":
        try:
            skbio_mod = import_module("skbio")
            skbio_seq = import_module("skbio.sequence")
            TabularMSA = skbio_mod.TabularMSA
        except ModuleNotFoundError as e:
            raise ImportError(
                "scikit-bio not installed. Run: pip install kalign-python[skbio]"
            ) from e

        skbio_type_map = {
            DNA: skbio_seq.DNA,
            DNA_INTERNAL: skbio_seq.DNA,
            RNA: skbio_seq.RNA,
            PROTEIN: skbio_seq.Protein,
            PROTEIN_DIVERGENT: skbio_seq.Protein,
        }
        if seq_type_int in skbio_type_map:
            SeqClass = skbio_type_map[seq_type_int]
        else:
            SeqClass = _infer_skbio_type(sequences, skbio_seq)

        return TabularMSA(
            [SeqClass(s, metadata={"id": i}) for s, i in zip(aligned_seqs, ids)]
        )

    raise ValueError(f"Unknown fmt='{fmt}' (expected 'plain', 'biopython', 'skbio')")


def set_num_threads(n: int) -> None:
    """Set the default number of threads for alignment operations."""
    global _default_threads
    if n < 1:
        raise ValueError("Number of threads must be at least 1")
    _thread_local.num_threads = n
    _default_threads = n


def get_num_threads() -> int:
    """Get the current default number of threads."""
    return getattr(_thread_local, "num_threads", _default_threads)


def align_from_file(
    input_file: str,
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: Optional[int] = None,
    mode: Optional[str] = None,
) -> AlignedSequences:
    """
    Align sequences from a file using Kalign.

    Parameters
    ----------
    input_file : str
        Path to input file (FASTA, MSF, Clustal).
    seq_type : str or int, optional
        Sequence type (default: "auto").
    gap_open, gap_extend, terminal_gap_extend : float, optional
        Gap penalties. When set, mode is ignored and "fast" preset is used.
    n_threads : int, optional
        Number of threads.
    mode : str, optional
        Preset mode: "fast", "default", "accurate" (default: "default").

    Returns
    -------
    AlignedSequences
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    seq_type_int = _resolve_seq_type(seq_type)

    if n_threads is None:
        n_threads = get_num_threads()
    if n_threads < 1:
        raise ValueError("n_threads must be at least 1")

    has_gap_override = (
        gap_open is not None or gap_extend is not None
        or terminal_gap_extend is not None
    )
    if has_gap_override:
        effective_mode = "fast"
    else:
        effective_mode = _resolve_mode_name(mode)

    result = _core.align_from_file_mode(
        input_file, effective_mode, seq_type_int,
        gap_open if gap_open is not None else -1.0,
        gap_extend if gap_extend is not None else -1.0,
        terminal_gap_extend if terminal_gap_extend is not None else -1.0,
        n_threads,
    )
    if len(result) == 3:
        names, sequences, conf = result
        col_conf = list(conf["column_confidence"])
        res_conf = [list(row) for row in conf["residue_confidence"]]
        return AlignedSequences(
            names=names, sequences=sequences,
            column_confidence=col_conf, residue_confidence=res_conf,
        )
    else:
        names, sequences = result
        return AlignedSequences(names=names, sequences=sequences)


def write_alignment(
    sequences: List[str],
    output_file: str,
    format: str = "fasta",
    ids: Optional[List[str]] = None,
    column_confidence: Optional[List[float]] = None,
    residue_confidence: Optional[List[List[float]]] = None,
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
    """
    if not sequences:
        raise ValueError("Empty sequence list provided")

    format_lower = format.lower()
    format_map = {
        "fasta": "fasta", "fa": "fasta",
        "clustal": "clustal", "aln": "clustal",
        "stockholm": "stockholm", "sto": "stockholm",
        "phylip": "phylip", "phy": "phylip",
    }
    if format_lower not in format_map:
        raise ValueError(
            f"Invalid format: {format}. Must be one of: fasta, clustal, stockholm, phylip"
        )
    mapped_format = format_map[format_lower]

    from . import io as _io
    if mapped_format == "fasta":
        _io.write_fasta(sequences, output_file, ids=ids)
    elif mapped_format == "clustal":
        _io.write_clustal(sequences, output_file, ids=ids)
    elif mapped_format == "stockholm":
        _io.write_stockholm(
            sequences, output_file, ids=ids,
            column_confidence=column_confidence,
            residue_confidence=residue_confidence,
        )
    elif mapped_format == "phylip":
        _io.write_phylip(sequences, output_file, ids=ids)


def generate_test_sequences(
    n_seq: int, n_obs: int, dna: bool, length: int, seed: int = 42
) -> List[str]:
    """Generate test sequences using DSSim HMM-based simulator."""
    if n_seq < 1:
        raise ValueError("n_seq must be at least 1")
    if n_obs < 1:
        raise ValueError("n_obs must be at least 1")
    if length < 1:
        raise ValueError("length must be at least 1")

    sequences = _core.generate_test_sequences(n_seq, n_obs, dna, length, seed)
    return sequences


def compare(reference_file: str, test_file: str) -> float:
    """
    Compare two multiple sequence alignments and return SP score.

    Parameters
    ----------
    reference_file : str
        Path to reference alignment file
    test_file : str
        Path to test alignment file

    Returns
    -------
    float
        SP score (0-100)
    """
    if not os.path.exists(reference_file):
        raise FileNotFoundError(f"Reference file not found: {reference_file}")
    if not os.path.exists(test_file):
        raise FileNotFoundError(f"Test file not found: {test_file}")
    return _core.compare(reference_file, test_file)


def compare_detailed(
    reference_file: str,
    test_file: str,
    max_gap_frac: float = 0.2,
    column_mask: Optional[List[int]] = None,
) -> dict:
    """
    Compare two MSAs returning detailed POAR scores (recall, precision, F1, TC).

    Parameters
    ----------
    reference_file : str
        Path to reference alignment file
    test_file : str
        Path to test alignment file
    max_gap_frac : float, optional
        Max gap fraction for scored columns (default: 0.2).
    column_mask : list of int, optional
        Explicit binary mask for column scoring.

    Returns
    -------
    dict
        Keys: recall, precision, f1, tc, ref_pairs, test_pairs, common_pairs
    """
    if not os.path.exists(reference_file):
        raise FileNotFoundError(f"Reference file not found: {reference_file}")
    if not os.path.exists(test_file):
        raise FileNotFoundError(f"Test file not found: {test_file}")

    if column_mask is not None:
        return _core.compare_detailed_with_mask(reference_file, test_file, column_mask)
    return _core.compare_detailed(reference_file, test_file, max_gap_frac)


def align_file_to_file(
    input_file: str,
    output_file: str,
    format: str = "fasta",
    seq_type: Union[str, int] = "auto",
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    terminal_gap_extend: Optional[float] = None,
    n_threads: Optional[int] = None,
    mode: Optional[str] = None,
) -> None:
    """
    Align sequences from input file and write result to output file.

    Parameters
    ----------
    input_file : str
        Path to input file
    output_file : str
        Path to output alignment file
    format : str, optional
        Output format: "fasta", "msf", "clu" (default: "fasta")
    seq_type : str or int, optional
        Sequence type (default: "auto")
    gap_open, gap_extend, terminal_gap_extend : float, optional
        Gap penalties. When set, mode is ignored and "fast" preset is used.
    n_threads : int, optional
        Number of threads.
    mode : str, optional
        Preset mode: "fast", "default", "accurate" (default: "default").
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    seq_type_int = _resolve_seq_type(seq_type)

    if n_threads is None:
        n_threads = get_num_threads()

    has_gap_override = (
        gap_open is not None or gap_extend is not None
        or terminal_gap_extend is not None
    )
    if has_gap_override:
        effective_mode = "fast"
    else:
        effective_mode = _resolve_mode_name(mode)

    _core.align_file_to_file_mode(
        input_file, output_file, effective_mode, format, n_threads,
        seq_type_int,
        gap_open if gap_open is not None else -1.0,
        gap_extend if gap_extend is not None else -1.0,
        terminal_gap_extend if terminal_gap_extend is not None else -1.0,
    )


# Convenience aliases
kalign = align


__all__ = [
    "align",
    "align_from_file",
    "align_file_to_file",
    "compare",
    "compare_detailed",
    "write_alignment",
    "generate_test_sequences",
    "set_num_threads",
    "get_num_threads",
    "kalign",
    "AlignedSequences",
    "DNA",
    "DNA_INTERNAL",
    "RNA",
    "PROTEIN",
    "PROTEIN_PFASUM43",
    "PROTEIN_PFASUM60",
    "PROTEIN_PFASUM_AUTO",
    "PROTEIN_DIVERGENT",
    "PROTEIN_CORBLOSUM66",
    "AUTO",
    "REFINE_NONE",
    "REFINE_ALL",
    "REFINE_CONFIDENT",
    "REFINE_INLINE",
    "MODE_DEFAULT",
    "MODE_FAST",
    "MODE_ACCURATE",
    "MODE_PRECISE",
    "__version__",
    "__author__",
    "__email__",
    "io",
    "utils",
]
