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
        "--refine",
        default="confident",
        help="Refinement mode: none, all, confident (default: confident).",
    )
    parser.add_argument(
        "--adaptive-budget",
        action="store_true",
        default=False,
        help="Scale refinement trial count by uncertainty.",
    )
    ens = parser.add_argument_group(
        "ensemble options",
        "These options only take effect when --ensemble is used.",
    )
    ens.add_argument(
        "--ensemble",
        type=int,
        default=0,
        help="Number of ensemble runs (default: 0 = off). Try 3-5 for better accuracy.",
    )
    ens.add_argument(
        "--ensemble-seed",
        type=int,
        default=42,
        help="RNG seed for ensemble (default: 42).",
    )
    ens.add_argument(
        "--min-support",
        type=int,
        default=0,
        help="Explicit consensus threshold (default: 0 = auto).",
    )
    ens.add_argument(
        "--save-poar",
        default=None,
        help="Save POAR consensus table to file for later re-thresholding.",
    )
    ens.add_argument(
        "--load-poar",
        default=None,
        help="Load POAR consensus table from file (skip alignment, just re-threshold).",
    )

    adv = parser.add_argument_group("advanced options")
    adv.add_argument(
        "--dist-scale",
        type=float,
        default=0.0,
        help="Distance scaling parameter (default: 0.0).",
    )
    adv.add_argument(
        "--vsm-amax",
        type=float,
        default=None,
        help="Variable Scoring Matrix a_max (default: Kalign internal defaults).",
    )
    adv.add_argument(
        "--realign",
        type=int,
        default=0,
        help="Realignment iterations (default: 0 = off).",
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
        kwargs = dict(
            seq_type=args.seq_type,
            gap_open=args.gpo,
            gap_extend=args.gpe,
            terminal_gap_extend=args.tgpe,
            n_threads=args.nthreads,
            refine=args.refine,
            adaptive_budget=args.adaptive_budget,
            ensemble=args.ensemble,
            ensemble_seed=args.ensemble_seed,
            min_support=args.min_support,
            realign=args.realign,
            dist_scale=args.dist_scale,
        )
        if args.vsm_amax is not None:
            kwargs["vsm_amax"] = args.vsm_amax
        if args.save_poar is not None:
            kwargs["save_poar"] = args.save_poar
        if args.load_poar is not None:
            kwargs["load_poar"] = args.load_poar

        result = kalign.align_from_file(input_path, **kwargs)

        # Pass confidence data to Stockholm writer when available
        write_kwargs = dict(
            format=args.format,
            ids=result.names,
        )
        if result.column_confidence is not None:
            write_kwargs["column_confidence"] = result.column_confidence
        if result.residue_confidence is not None:
            write_kwargs["residue_confidence"] = result.residue_confidence

        if args.output == "-":
            _write_stdout(result.names, result.sequences, args.format)
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
