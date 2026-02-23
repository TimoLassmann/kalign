"""INDELible sequence simulation and tree generation via dendropy.

Provides utilities for generating simulated multiple sequence alignments
with known true alignments, for benchmarking alignment accuracy.
"""

from __future__ import annotations

import itertools
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, List, Optional


@dataclass
class SimulatedDataset:
    """Result of a single INDELible simulation run."""

    true_alignment: Path  # FASTA with gaps (true alignment)
    unaligned: Path  # FASTA without gaps (input for aligners)
    true_tree: Path  # Newick file
    site_classes: List[int]  # Per-site class (0=neutral, 1=positive selection)
    params: dict = field(default_factory=dict)  # All simulation parameters


# ---------------------------------------------------------------------------
# Tree generation
# ---------------------------------------------------------------------------


def random_birth_death_tree(
    n_taxa: int,
    target_depth: float,
    seed: int = 42,
) -> str:
    """Generate a random birth-death tree with dendropy, scale to target depth.

    Uses ``dendropy.simulate.treesim.birth_death_tree()`` with
    birth_rate=1.0, death_rate=0.5.  Then scales all branch lengths so the
    maximum root-to-tip distance equals *target_depth*.

    Returns a Newick string.
    """
    import random

    from dendropy.simulate import treesim

    tree = treesim.birth_death_tree(
        birth_rate=1.0,
        death_rate=0.5,
        num_extant_tips=n_taxa,
        rng=random.Random(seed),
    )

    # Compute the maximum root-to-tip distance.
    max_dist = 0.0
    for leaf in tree.leaf_node_iter():
        dist = leaf.distance_from_root()
        if dist > max_dist:
            max_dist = dist

    if max_dist > 0.0:
        scale_factor = target_depth / max_dist
        for edge in tree.preorder_edge_iter():
            if edge.length is not None:
                edge.length *= scale_factor

    nwk = tree.as_string(schema="newick").strip()
    # Strip dendropy annotations like [&R] that INDELible doesn't understand,
    # and remove the root branch length (e.g. "...):0.5;" → "...);")
    # which INDELible rejects.
    import re

    nwk = re.sub(r"\[&[^\]]*\]", "", nwk).strip()
    nwk = re.sub(r"\):[0-9eE.+-]+;$", ");", nwk)
    # INDELible cannot parse scientific notation (e.g. 1.23e-10) in branch
    # lengths.  Reformat all `:NUMBER` tokens to fixed-point decimal.
    def _fix_brlen(m: re.Match) -> str:
        return f":{float(m.group(1)):.10f}"
    nwk = re.sub(r":([0-9eE.+-]+)", _fix_brlen, nwk)
    return nwk


# ---------------------------------------------------------------------------
# INDELible control-file writers
# ---------------------------------------------------------------------------


def _write_codon_control(
    output_dir: Path,
    tree: str,
    n_codons: int = 400,
    kappa: float = 2.0,
    indel_rate: float = 0.05,
    psel_fraction: float = 0.0,
    omega_positive: float = 3.0,
    seed: int = 42,
) -> Path:
    """Write INDELible control file for codon simulation.

    If *psel_fraction* > 0, uses M1a (two dN/dS categories: purifying +
    positive selection).  Otherwise uses M0 (single omega, purifying only).

    INDELible ``[submodel]`` parameter counts determine the codon model:
    - 2 params → M0: ``kappa omega``
    - 4 params → M1a: ``kappa p0 omega0 omega1``

    Returns the path to the written ``control.txt``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build the model block.
    # Max indel length (in codons).  Without a bound, INDELible's power-law
    # sampler can hang when indel_rate is high.
    max_indel_codons = 10

    if psel_fraction > 0.0:
        # M1a: two categories — purifying + positive selection
        p0 = 1.0 - psel_fraction  # proportion in purifying class
        omega_purifying = 0.2
        model_block = (
            f"[MODEL] simmodel\n"
            f"  [submodel] {kappa} {p0} {omega_purifying} {omega_positive}\n"
            f"  [indelmodel] POW 1.7 {max_indel_codons}\n"
            f"  [indelrate] {indel_rate}\n"
        )
    else:
        # M0: single omega (purifying only)
        omega_purifying = 0.2
        model_block = (
            f"[MODEL] simmodel\n"
            f"  [submodel] {kappa} {omega_purifying}\n"
            f"  [indelmodel] POW 1.7 {max_indel_codons}\n"
            f"  [indelrate] {indel_rate}\n"
        )

    control = (
        f"[TYPE] CODON 1\n"
        f"\n"
        f"[SETTINGS]\n"
        f"  [randomseed] {seed}\n"
        f"  [printrates] TRUE\n"
        f"\n"
        f"{model_block}\n"
        f"[TREE] t1 {tree}\n"
        f"\n"
        f"[PARTITIONS] p1\n"
        f"  [t1 simmodel {n_codons}]\n"
        f"\n"
        f"[EVOLVE] p1 1 output\n"
    )

    control_path = output_dir / "control.txt"
    control_path.write_text(control)
    return control_path


def _write_protein_control(
    output_dir: Path,
    tree: str,
    seq_length: int = 300,
    indel_rate: float = 0.05,
    indel_length_mean: float = 2.0,
    seed: int = 42,
) -> Path:
    """Write INDELible control file for protein simulation (WAG+G4).

    Returns the path to the written ``control.txt``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Max indel length is set from the mean; INDELible uses a power-law
    # distribution with an upper bound.
    max_indel_length = int(round(indel_length_mean * 5))
    if max_indel_length < 1:
        max_indel_length = 1

    control = (
        f"[TYPE] AMINOACID 1\n"
        f"\n"
        f"[SETTINGS]\n"
        f"  [randomseed] {seed}\n"
        f"\n"
        f"[MODEL] simmodel\n"
        f"  [submodel] WAG\n"
        f"  [rates] 0 0.5 4\n"
        f"  [indelmodel] POW 1.7 {max_indel_length}\n"
        f"  [indelrate] {indel_rate}\n"
        f"\n"
        f"[TREE] t1 {tree}\n"
        f"\n"
        f"[PARTITIONS] p1\n"
        f"  [t1 simmodel {seq_length}]\n"
        f"\n"
        f"[EVOLVE] p1 1 output\n"
    )

    control_path = output_dir / "control.txt"
    control_path.write_text(control)
    return control_path


# ---------------------------------------------------------------------------
# INDELible execution and output parsing
# ---------------------------------------------------------------------------


def _run_indelible(control_dir: Path) -> None:
    """Run INDELible in the given directory.

    Raises ``RuntimeError`` if the executable is not found or if it exits
    with a non-zero return code.
    """
    indelible_bin = shutil.which("indelible") or shutil.which("INDELible")
    if indelible_bin is None:
        raise RuntimeError(
            "INDELible is not installed or not on PATH. "
            "Please install INDELible (http://abacus.gene.ucl.ac.uk/software/indelible/) "
            "and ensure the 'indelible' binary is accessible."
        )

    try:
        result = subprocess.run(
            [indelible_bin, "control.txt"],
            cwd=str(control_dir),
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"INDELible timed out after 120s in {control_dir}"
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"INDELible failed (exit code {result.returncode}).\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def _parse_indelible_output(
    output_dir: Path,
) -> tuple[list[str], list[str], list[str], list[str]]:
    """Parse INDELible output files.

    Returns ``(true_names, true_sequences_with_gaps, unaligned_names,
    unaligned_sequences)``.

    INDELible produces (original v1.03 naming / matsengrp fork naming):
    - ``output_TRUE_1.phy`` or ``output_TRUE.phy`` — true alignment (relaxed PHYLIP)
    - ``output_1.fas`` or ``output.fas`` — unaligned sequences (FASTA)
    """
    output_dir = Path(output_dir)

    # -- Parse the true alignment (relaxed PHYLIP) --------------------------
    # Handle both original INDELible naming (_1 suffix) and matsengrp fork (no suffix)
    true_phy = output_dir / "output_TRUE_1.phy"
    if not true_phy.exists():
        true_phy = output_dir / "output_TRUE.phy"
    if not true_phy.exists():
        raise FileNotFoundError(
            f"Expected INDELible true-alignment file not found: {output_dir}/output_TRUE{{_1,}}.phy"
        )

    true_names: list[str] = []
    true_seqs: list[str] = []
    with open(true_phy) as fh:
        header = fh.readline().strip().split()
        # First line: "n_seqs seq_length"
        _n_seqs = int(header[0])
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                true_names.append(parts[0])
                true_seqs.append(parts[1].replace(" ", ""))

    # -- Parse the unaligned FASTA ------------------------------------------
    unaln_fas = output_dir / "output_1.fas"
    if not unaln_fas.exists():
        unaln_fas = output_dir / "output.fas"
    if not unaln_fas.exists():
        raise FileNotFoundError(
            f"Expected INDELible unaligned FASTA file not found: {output_dir}/output{{_1,}}.fas"
        )

    unaln_names: list[str] = []
    unaln_seqs: list[str] = []
    current_name: Optional[str] = None
    current_seq_parts: list[str] = []
    with open(unaln_fas) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_name is not None:
                    unaln_names.append(current_name)
                    unaln_seqs.append("".join(current_seq_parts))
                current_name = line[1:].strip()
                current_seq_parts = []
            else:
                current_seq_parts.append(line.strip())
        # Last sequence
        if current_name is not None:
            unaln_names.append(current_name)
            unaln_seqs.append("".join(current_seq_parts))

    return true_names, true_seqs, unaln_names, unaln_seqs


def _parse_site_classes(output_dir: Path, n_codons: int) -> List[int]:
    """Parse site class assignments from INDELible rates output.

    Returns a list of ints: 0 = neutral, 1 = positive selection.
    If no rates file is found, returns a list of zeros.

    INDELible writes per-site rate information to ``output_1_RATES.txt``
    when ``[printrates] TRUE`` is set.  The file has columns:
    Site, Class, Rate (tab- or space-separated, with a header line).
    For M8, class index equal to the last class is the positive-selection
    class.
    """
    rates_file = output_dir / "output_1_RATES.txt"
    if not rates_file.exists():
        rates_file = output_dir / "output_RATES.txt"
    if not rates_file.exists():
        return [0] * n_codons

    classes: list[int] = []
    max_class = 0
    with open(rates_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("Site") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    cls = int(parts[1])
                except ValueError:
                    continue
                classes.append(cls)
                if cls > max_class:
                    max_class = cls

    # In M8, the highest class index is the positive-selection class.
    # Map: highest class -> 1 (positive selection), everything else -> 0.
    if max_class == 0:
        # All neutral (M7 or single class)
        return [0] * len(classes) if classes else [0] * n_codons

    site_classes = [1 if c == max_class else 0 for c in classes]

    # Pad or truncate to n_codons (rates file should match, but be safe).
    if len(site_classes) < n_codons:
        site_classes.extend([0] * (n_codons - len(site_classes)))
    elif len(site_classes) > n_codons:
        site_classes = site_classes[:n_codons]

    return site_classes


# ---------------------------------------------------------------------------
# Helper: strip gaps
# ---------------------------------------------------------------------------


def strip_gaps(sequences: list[str]) -> list[str]:
    """Remove all gap characters from sequences."""
    return [s.replace("-", "").replace(".", "") for s in sequences]


# ---------------------------------------------------------------------------
# FASTA I/O helpers
# ---------------------------------------------------------------------------


def _write_fasta(path: Path, names: list[str], sequences: list[str]) -> None:
    """Write sequences in FASTA format."""
    with open(path, "w") as fh:
        for name, seq in zip(names, sequences):
            fh.write(f">{name}\n{seq}\n")


# ---------------------------------------------------------------------------
# Main simulation function
# ---------------------------------------------------------------------------


def generate_indelible_dataset(
    tree: str,
    model: str,
    seq_length: int,
    indel_rate: float,
    indel_length_mean: float = 2.0,
    output_dir: Optional[Path] = None,
    seed: int = 42,
    psel_fraction: float = 0.0,
    omega_positive: float = 3.0,
    kappa: float = 2.0,
) -> SimulatedDataset:
    """Generate a simulated dataset using INDELible.

    Parameters
    ----------
    tree : str
        Newick tree string.
    model : str
        Substitution model: ``"WAG"``, ``"M7"``, or ``"M8"``.
    seq_length : int
        Number of sites (amino acids) or codons depending on model.
    indel_rate : float
        Insertion/deletion rate.
    indel_length_mean : float
        Mean indel length (protein models only).
    output_dir : Path or None
        Directory for output files.  If *None*, a temporary directory is used.
    seed : int
        Random seed for reproducibility.
    psel_fraction : float
        Fraction of sites under positive selection (codon models only).
    omega_positive : float
        Omega for the positive-selection class (codon models only).
    kappa : float
        Transition/transversion ratio (codon models only).

    Returns
    -------
    SimulatedDataset
        Paths to output files and simulation metadata.
    """
    use_temp = output_dir is None
    if use_temp:
        tmpdir = tempfile.mkdtemp(prefix="indelible_")
        output_dir = Path(tmpdir)
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Store all params for provenance.
    params = dict(
        model=model,
        seq_length=seq_length,
        indel_rate=indel_rate,
        indel_length_mean=indel_length_mean,
        seed=seed,
        psel_fraction=psel_fraction,
        omega_positive=omega_positive,
        kappa=kappa,
    )

    tree_path = output_dir / "tree.nwk"
    true_aln_path = output_dir / "true_alignment.fasta"
    unaligned_path = output_dir / "unaligned.fasta"

    # Skip INDELible if output files already exist (deterministic given params)
    if (
        not use_temp
        and tree_path.exists()
        and true_aln_path.exists()
        and unaligned_path.exists()
        and true_aln_path.stat().st_size > 0
    ):
        model_upper = model.upper()
        if model_upper in ("M7", "M8"):
            site_classes = _parse_site_classes(output_dir, seq_length)
        else:
            site_classes = [0] * seq_length
        return SimulatedDataset(
            true_alignment=true_aln_path,
            unaligned=unaligned_path,
            true_tree=tree_path,
            site_classes=site_classes,
            params=params,
        )

    # Write the tree to a file for reference.
    tree_path.write_text(tree + "\n")

    # Write control file.
    model_upper = model.upper()
    if model_upper in ("M7", "M8"):
        effective_psel = psel_fraction if model_upper == "M8" else 0.0
        _write_codon_control(
            output_dir=output_dir,
            tree=tree,
            n_codons=seq_length,
            kappa=kappa,
            indel_rate=indel_rate,
            psel_fraction=effective_psel,
            omega_positive=omega_positive,
            seed=seed,
        )
        n_codons = seq_length
    elif model_upper == "WAG":
        _write_protein_control(
            output_dir=output_dir,
            tree=tree,
            seq_length=seq_length,
            indel_rate=indel_rate,
            indel_length_mean=indel_length_mean,
            seed=seed,
        )
        n_codons = 0
    else:
        raise ValueError(f"Unknown model: {model!r}. Expected 'WAG', 'M7', or 'M8'.")

    # Run INDELible.
    _run_indelible(output_dir)

    # Parse output.
    true_names, true_seqs, unaln_names, unaln_seqs = _parse_indelible_output(
        output_dir
    )

    # Parse site classes (codon models only).
    if model_upper in ("M7", "M8"):
        site_classes = _parse_site_classes(output_dir, n_codons)
    else:
        site_classes = [0] * seq_length

    # Write FASTA outputs.
    true_aln_path = output_dir / "true_alignment.fasta"
    unaligned_path = output_dir / "unaligned.fasta"

    _write_fasta(true_aln_path, true_names, true_seqs)
    _write_fasta(unaligned_path, unaln_names, unaln_seqs)

    return SimulatedDataset(
        true_alignment=true_aln_path,
        unaligned=unaligned_path,
        true_tree=tree_path,
        site_classes=site_classes,
        params=params,
    )


# ---------------------------------------------------------------------------
# Batch simulation: parameter grids
# ---------------------------------------------------------------------------


@dataclass
class SimulationGrid:
    """Full factorial simulation parameter grid."""

    n_taxa: List[int]
    tree_depths: List[float]
    indel_rates: List[float]
    replicates: int = 10
    # Codon-specific
    psel_fractions: List[float] = field(default_factory=lambda: [0.0])
    omega_positive: float = 3.0
    kappa: float = 2.0
    n_codons: int = 400
    # Protein-specific
    seq_length: int = 300
    indel_length_means: List[float] = field(default_factory=lambda: [2.0])


# -- Full grids (for final publication runs with --full) --------------------

CODON_GRID_FULL = SimulationGrid(
    n_taxa=[32, 64, 128],
    tree_depths=[0.5, 1.0, 2.0],
    indel_rates=[0.01, 0.05, 0.10],
    psel_fractions=[0.10, 0.20, 0.30],
    replicates=10,
    n_codons=400,
)

PROTEIN_GRID_FULL = SimulationGrid(
    n_taxa=[16, 32, 64],
    tree_depths=[0.5, 1.0, 2.0, 4.0],
    indel_rates=[0.02, 0.05, 0.10, 0.20],
    indel_length_means=[2.0, 5.0],
    replicates=20,
)

# -- Slim grids (default — faster iteration, still publication quality) -----
# Protein: 2×4×4×2×3 = 192 sims
# Codon:   3×3×3×3×3 = 243 sims
#
# Codon grid rationale: FUBAR is a Bayesian site-level test that needs
# sufficient taxa for statistical power (≥32) and enough sites under
# selection (psel ≥ 0.10) to produce meaningful F1 scores.  Previous
# grid (n_taxa=[16,32,64], psel=[0.05,0.10,0.20]) yielded F1 < 0.09
# for ALL methods including the true alignment, indicating the test was
# underpowered rather than aligners being poor.

PROTEIN_GRID = SimulationGrid(
    n_taxa=[16, 32],
    tree_depths=[0.5, 1.0, 2.0, 4.0],
    indel_rates=[0.02, 0.05, 0.10, 0.20],
    indel_length_means=[2.0, 5.0],
    replicates=3,
)

CODON_GRID = SimulationGrid(
    n_taxa=[32, 64, 128],
    tree_depths=[0.5, 1.0, 2.0],
    indel_rates=[0.01, 0.05, 0.10],
    psel_fractions=[0.10, 0.20, 0.30],
    replicates=3,
    n_codons=400,
)


def iter_simulation_params(grid: SimulationGrid, model: str) -> Iterator[dict]:
    """Yield parameter dicts for each simulation in the grid.

    Each dict has keys matching :func:`generate_indelible_dataset` kwargs
    plus a ``sim_id`` string that uniquely identifies the parameter
    combination and replicate number.
    """
    model_upper = model.upper()

    if model_upper == "WAG":
        # Protein grid: iterate over n_taxa x tree_depths x indel_rates
        #   x indel_length_means x replicates
        combos = itertools.product(
            grid.n_taxa,
            grid.tree_depths,
            grid.indel_rates,
            grid.indel_length_means,
            range(grid.replicates),
        )
        for n_taxa, depth, indel_rate, indel_len, rep in combos:
            sim_id = (
                f"{model_upper}_t{n_taxa}_d{depth}_ir{indel_rate}"
                f"_il{indel_len}_r{rep}"
            )
            yield dict(
                sim_id=sim_id,
                model=model_upper,
                n_taxa=n_taxa,
                target_depth=depth,
                seq_length=grid.seq_length,
                indel_rate=indel_rate,
                indel_length_mean=indel_len,
                seed=42 + rep,
                psel_fraction=0.0,
                omega_positive=grid.omega_positive,
                kappa=grid.kappa,
            )

    elif model_upper in ("M7", "M8"):
        # Codon grid: iterate over n_taxa x tree_depths x indel_rates
        #   x psel_fractions x replicates
        combos = itertools.product(
            grid.n_taxa,
            grid.tree_depths,
            grid.indel_rates,
            grid.psel_fractions,
            range(grid.replicates),
        )
        for n_taxa, depth, indel_rate, psel_frac, rep in combos:
            sim_id = (
                f"{model_upper}_t{n_taxa}_d{depth}_ir{indel_rate}"
                f"_ps{psel_frac}_r{rep}"
            )
            yield dict(
                sim_id=sim_id,
                model=model_upper,
                n_taxa=n_taxa,
                target_depth=depth,
                seq_length=grid.n_codons,
                indel_rate=indel_rate,
                indel_length_mean=2.0,
                seed=42 + rep,
                psel_fraction=psel_frac,
                omega_positive=grid.omega_positive,
                kappa=grid.kappa,
            )
    else:
        raise ValueError(f"Unknown model: {model!r}. Expected 'WAG', 'M7', or 'M8'.")
