"""
Command-line interface for the Python Kalign package.

This CLI uses the compiled Python extension module shipped with the package
instead of shelling out to a separately-installed `kalign` binary.

Installed as ``kalign-py`` so it does not shadow the C binary on PATH.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as dist_version
from pathlib import Path
from typing import Optional


def _resolve_version() -> str:
    for dist_name in ("kalign-python",):
        try:
            return dist_version(dist_name)
        except PackageNotFoundError:
            continue

    try:
        from . import __version__

        return __version__
    except Exception:
        return "unknown"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="kalign-py",
        description="Multiple sequence alignment via the Kalign Python package",
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input sequence file path, or '-' to read from stdin.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output file path, or '-' to write to stdout (default).",
    )
    parser.add_argument(
        "--format",
        default="fasta",
        help=(
            "Output format: fasta (fa), clustal (aln), stockholm (sto), "
            "phylip (phy). Case-insensitive (default: fasta). "
            "Confidence annotations are written only for stockholm."
        ),
    )
    parser.add_argument(
        "--type",
        default="auto",
        dest="seq_type",
        choices=["auto", "protein", "dna", "rna"],
        help="Sequence type: protein, dna or rna overrides detection (default: auto).",
    )
    parser.add_argument(
        "--mode",
        default="fast",
        choices=["fast", "default", "recall", "accurate"],
        help=(
            "fast: single run (default; alias: default). "
            "accurate: fitted 5-member ensemble with consensus and per-column "
            "confidence, at most 5,000 sequences (alias: recall). "
            "--gpo, --gpe and --tgpe can only be used with fast."
        ),
    )
    parser.add_argument(
        "--gpo",
        type=float,
        default=None,
        help="Gap open penalty. Fast mode only.",
    )
    parser.add_argument(
        "--gpe",
        type=float,
        default=None,
        help="Gap extension penalty. Fast mode only.",
    )
    parser.add_argument(
        "--tgpe",
        type=float,
        default=None,
        help="Terminal gap extension penalty. Fast mode only.",
    )
    parser.add_argument(
        "-n",
        "--nthreads",
        type=int,
        default=1,
        help=(
            "Number of threads. The default is 1 (the kalign command-line "
            "program uses more by default). Use -n to request more; the "
            "count is capped at the number of cores."
        ),
    )

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {_resolve_version()}",
        help="Print version and exit",
    )
    return parser


def _write_stdout(
    names: list[str],
    sequences: list[str],
    fmt: str,
    column_confidence: Optional[list[float]] = None,
    residue_confidence: Optional[list[list[float]]] = None,
) -> None:
    import kalign

    fmt_lower = fmt.lower()
    if fmt_lower in {"fasta", "fa"}:
        kalign.io.write_fasta(sequences, sys.stdout, ids=names)
        return

    if fmt_lower in {"clustal", "aln"}:
        kalign.io.write_clustal(sequences, sys.stdout, ids=names)
        return

    if fmt_lower in {"stockholm", "sto"}:
        kalign.io.write_stockholm(
            sequences,
            sys.stdout,
            ids=names,
            column_confidence=column_confidence,
            residue_confidence=residue_confidence,
        )
        return

    if fmt_lower in {"phylip", "phy"}:
        kalign.io.write_phylip(sequences, sys.stdout, ids=names)
        return

    raise ValueError(
        f"Invalid format: {fmt}. Must be one of: fasta, clustal, stockholm, phylip"
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    import kalign

    input_path = args.input
    tmp_path: Optional[Path] = None
    if input_path == "-":
        with tempfile.NamedTemporaryFile(
            prefix="kalign-", suffix=".fa", delete=False
        ) as tmp:
            tmp.write(sys.stdin.buffer.read())
            tmp_path = Path(tmp.name)
        input_path = str(tmp_path)

    try:
        result = kalign.align_from_file(
            input_path,
            seq_type=args.seq_type,
            gap_open=args.gpo,
            gap_extend=args.gpe,
            terminal_gap_extend=args.tgpe,
            n_threads=args.nthreads,
            mode=args.mode,
        )

        write_kwargs = dict(
            format=args.format,
            ids=result.names,
        )
        if result.column_confidence is not None:
            write_kwargs["column_confidence"] = result.column_confidence
        if result.residue_confidence is not None:
            write_kwargs["residue_confidence"] = result.residue_confidence

        if args.output == "-":
            _write_stdout(
                result.names,
                result.sequences,
                args.format,
                column_confidence=result.column_confidence,
                residue_confidence=result.residue_confidence,
            )
        else:
            kalign.write_alignment(
                result.sequences,
                args.output,
                **write_kwargs,
            )

        return 0
    except KeyboardInterrupt:
        return 130
    except Exception as exc:
        print(f"kalign-py: error: {exc}", file=sys.stderr)
        return 2
    finally:
        if tmp_path is not None:
            try:
                tmp_path.unlink()
            except OSError:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
