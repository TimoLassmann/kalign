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
    # kalign-test is the name of the test distribution, which may be installed in test environments. If it's present, use its version; otherwise, fall back to the main kalign distribution. If neither is found, try to import __version__ from the package, and if that fails, return "unknown".
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
        help="Output format: fasta, clustal, stockholm, phylip (default: fasta).",
    )
    parser.add_argument(
        "--type",
        default="auto",
        dest="seq_type",
        help="Sequence type: auto, dna, rna, internal, protein, divergent (default: auto).",
    )
    parser.add_argument(
        "--gpo",
        type=float,
        default=None,
        help="Gap open penalty (default: Kalign internal defaults).",
    )
    parser.add_argument(
        "--gpe",
        type=float,
        default=None,
        help="Gap extension penalty (default: Kalign internal defaults).",
    )
    parser.add_argument(
        "--tgpe",
        type=float,
        default=None,
        help="Terminal gap extension penalty (default: Kalign internal defaults).",
    )
    parser.add_argument(
        "-n",
        "--nthreads",
        type=int,
        default=1,
        help="Number of threads to use (default: 1).",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {_resolve_version()}",
        help="Print version and exit",
    )
    return parser


def _write_stdout(names: list[str], sequences: list[str], fmt: str) -> None:
    import kalign

    fmt_lower = fmt.lower()
    if fmt_lower in {"fasta", "fa"}:
        kalign.io.write_fasta(sequences, sys.stdout, ids=names)
        return

    if fmt_lower in {"clustal", "aln"}:
        kalign.io.write_clustal(sequences, sys.stdout, ids=names)
        return

    if fmt_lower in {"stockholm", "sto"}:
        kalign.io.write_stockholm(sequences, sys.stdout, ids=names)
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
        )

        if args.output == "-":
            _write_stdout(result.names, result.sequences, args.format)
        else:
            kalign.write_alignment(
                result.sequences,
                args.output,
                format=args.format,
                ids=result.names,
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
