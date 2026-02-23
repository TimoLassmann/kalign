"""Integration tests for downstream benchmark package.

Verifies that all modules import cleanly, dataclasses have correct fields,
pipeline functions exist with correct signatures, the CLI works, and
cross-module interactions are correct.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pytest


# ======================================================================
# Module import tests
# ======================================================================


class TestModuleImports:
    """Verify every module in the downstream package imports without error."""

    def test_import_package(self):
        from benchmarks.downstream import __init__  # noqa: F401

    def test_import_provenance(self):
        from benchmarks.downstream import provenance  # noqa: F401

    def test_import_utils(self):
        from benchmarks.downstream import utils  # noqa: F401

    def test_import_simulation(self):
        from benchmarks.downstream import simulation  # noqa: F401

    def test_import_calibration(self):
        from benchmarks.downstream import calibration  # noqa: F401

    def test_import_positive_selection(self):
        from benchmarks.downstream import positive_selection  # noqa: F401

    def test_import_phylo_accuracy(self):
        from benchmarks.downstream import phylo_accuracy  # noqa: F401

    def test_import_hmmer_detection(self):
        from benchmarks.downstream import hmmer_detection  # noqa: F401

    def test_import_figures(self):
        from benchmarks.downstream import figures  # noqa: F401

    def test_import_main(self):
        from benchmarks.downstream import __main__  # noqa: F401


# ======================================================================
# Dataclass field tests
# ======================================================================


class TestDataclassFields:
    """Verify all dataclasses have the expected fields."""

    @staticmethod
    def _field_names(cls):
        import dataclasses

        return {f.name for f in dataclasses.fields(cls)}

    def test_provenance_fields(self):
        from benchmarks.downstream.provenance import Provenance

        fields = self._field_names(Provenance)
        for f in [
            "timestamp", "kalign_version", "kalign_commit",
            "container_image", "hostname", "cpu_model", "cpu_cores",
            "ram_gb", "os_version", "python_version",
            "tool_versions", "parameters",
        ]:
            assert f in fields, f"Provenance missing field: {f}"

    def test_case_result_fields(self):
        from benchmarks.downstream.calibration import CaseResult

        fields = self._field_names(CaseResult)
        for f in [
            "sim_id", "method", "predicted_confidence", "actual_correct",
            "sp_score", "tc_score", "wall_time", "peak_memory_mb",
        ]:
            assert f in fields, f"CaseResult missing field: {f}"

    def test_selection_case_result_fields(self):
        from benchmarks.downstream.positive_selection import SelectionCaseResult

        fields = self._field_names(SelectionCaseResult)
        for f in [
            "sim_id", "method", "true_positives", "false_positives",
            "false_negatives", "true_negatives", "precision", "recall",
            "f1", "sp_score",
        ]:
            assert f in fields, f"SelectionCaseResult missing field: {f}"

    def test_phylo_case_result_fields(self):
        from benchmarks.downstream.phylo_accuracy import PhyloCaseResult

        fields = self._field_names(PhyloCaseResult)
        for f in [
            "sim_id", "method", "nrf", "branch_score_dist", "sp_score",
            "wall_time_align", "wall_time_iqtree",
        ]:
            assert f in fields, f"PhyloCaseResult missing field: {f}"

    def test_hmmer_case_result_fields(self):
        from benchmarks.downstream.hmmer_detection import HmmerCaseResult

        fields = self._field_names(HmmerCaseResult)
        for f in [
            "family_id", "method", "true_positives", "false_positives",
            "false_negatives", "sensitivity",
        ]:
            assert f in fields, f"HmmerCaseResult missing field: {f}"

    def test_simulated_dataset_fields(self):
        from benchmarks.downstream.simulation import SimulatedDataset

        fields = self._field_names(SimulatedDataset)
        for f in ["true_alignment", "unaligned", "true_tree", "site_classes", "params"]:
            assert f in fields, f"SimulatedDataset missing field: {f}"


# ======================================================================
# Public API tests
# ======================================================================


class TestPublicAPI:
    """Verify each module exposes the expected public functions."""

    def test_provenance_api(self):
        from benchmarks.downstream import provenance

        assert callable(provenance.collect_provenance)
        assert callable(provenance.collect_tool_versions)
        assert callable(provenance.result_path)
        assert callable(provenance.update_latest_symlink)
        assert callable(provenance.load_latest_results)

    def test_utils_api(self):
        from benchmarks.downstream import utils

        assert callable(utils.parse_fasta)
        assert callable(utils.write_fasta)
        assert callable(utils.mask_alignment_by_confidence)
        assert callable(utils.write_site_weights)
        assert callable(utils.alignment_accuracy)
        assert callable(utils.compare_trees)
        assert callable(utils.bootstrap_ci)
        assert callable(utils.wilcoxon_paired)
        assert callable(utils.holm_bonferroni)
        assert callable(utils.run_method)
        assert isinstance(utils.METHODS, dict)
        assert isinstance(utils.METHOD_COLORS, dict)

    def test_simulation_api(self):
        from benchmarks.downstream import simulation

        assert callable(simulation.random_birth_death_tree)
        assert callable(simulation.generate_indelible_dataset)
        assert callable(simulation.iter_simulation_params)
        assert callable(simulation.strip_gaps)

    def test_calibration_api(self):
        from benchmarks.downstream import calibration

        assert callable(calibration.brier_score)
        assert callable(calibration.expected_calibration_error)
        assert callable(calibration.calibration_curve)
        assert callable(calibration.run_calibration_case)
        assert callable(calibration.run_pipeline)
        assert callable(calibration.load_results)

    def test_positive_selection_api(self):
        from benchmarks.downstream import positive_selection

        assert callable(positive_selection.run_selection_case)
        assert callable(positive_selection.run_pipeline)
        assert callable(positive_selection.load_results)

    def test_phylo_accuracy_api(self):
        from benchmarks.downstream import phylo_accuracy

        assert callable(phylo_accuracy.run_phylo_case)
        assert callable(phylo_accuracy.run_pipeline)
        assert callable(phylo_accuracy.load_results)

    def test_hmmer_detection_api(self):
        from benchmarks.downstream import hmmer_detection

        assert callable(hmmer_detection.run_hmmer_case)
        assert callable(hmmer_detection.run_pipeline)
        assert callable(hmmer_detection.load_results)
        assert isinstance(hmmer_detection.PFAM_FAMILIES, list)
        assert len(hmmer_detection.PFAM_FAMILIES) == 50

    def test_figures_api(self):
        from benchmarks.downstream import figures

        assert callable(figures.figure_calibration)
        assert callable(figures.figure_positive_selection)
        assert callable(figures.figure_phylo_accuracy)
        assert callable(figures.figure_hmmer_detection)
        assert callable(figures.figure_speed_comparison)
        assert callable(figures.figure_summary_heatmap)
        assert callable(figures.generate_all_figures)


# ======================================================================
# Calibration metric tests
# ======================================================================


class TestCalibrationMetrics:
    """Test calibration scoring functions with known inputs."""

    def test_brier_score_perfect(self):
        from benchmarks.downstream.calibration import brier_score

        # Perfect predictions: confidence matches outcome exactly
        assert brier_score([1.0, 0.0, 1.0], [1, 0, 1]) == pytest.approx(0.0)

    def test_brier_score_worst(self):
        from benchmarks.downstream.calibration import brier_score

        # Worst predictions: fully confident but wrong
        assert brier_score([1.0, 1.0], [0, 0]) == pytest.approx(1.0)

    def test_brier_score_uniform(self):
        from benchmarks.downstream.calibration import brier_score

        # Uniform 0.5 confidence, half correct
        score = brier_score([0.5, 0.5, 0.5, 0.5], [1, 0, 1, 0])
        assert score == pytest.approx(0.25)

    def test_brier_score_empty(self):
        from benchmarks.downstream.calibration import brier_score

        assert brier_score([], []) == pytest.approx(0.0)

    def test_ece_perfect(self):
        from benchmarks.downstream.calibration import expected_calibration_error

        # Perfect calibration: 100% confidence, all correct
        assert expected_calibration_error([1.0, 1.0], [1, 1]) == pytest.approx(
            0.0, abs=0.01
        )

    def test_calibration_curve_shape(self):
        from benchmarks.downstream.calibration import calibration_curve

        pred = [0.1, 0.2, 0.3, 0.8, 0.9, 0.95]
        actual = [0, 0, 1, 1, 1, 1]
        centers, fracs, counts = calibration_curve(pred, actual, n_bins=5)
        assert len(centers) == 5
        assert len(fracs) == 5
        assert len(counts) == 5


# ======================================================================
# Utils integration tests
# ======================================================================


class TestUtilsIntegration:
    """Test utility functions with realistic data."""

    def test_fasta_roundtrip_with_long_sequences(self):
        from benchmarks.downstream.utils import parse_fasta, write_fasta

        names = [f"seq{i}" for i in range(10)]
        seqs = ["ACDEFGHIKLMNPQRSTVWY" * 50 for _ in range(10)]  # 1000 chars each

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.fasta"
            write_fasta(names, seqs, path)
            loaded_names, loaded_seqs = parse_fasta(path)
            assert loaded_names == names
            assert loaded_seqs == seqs

    def test_mask_then_write(self):
        from benchmarks.downstream.utils import (
            mask_alignment_by_confidence,
            parse_fasta,
            write_fasta,
        )

        names = ["s1", "s2"]
        seqs = ["ABCDEFGH", "IJKLMNOP"]
        conf = [0.9, 0.1, 0.8, 0.2, 0.7, 0.3, 0.6, 0.4]

        masked, n_kept = mask_alignment_by_confidence(seqs, conf, 0.5)
        assert n_kept == 4

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "masked.fasta"
            write_fasta(names, masked, path)
            loaded_names, loaded_seqs = parse_fasta(path)
            assert loaded_names == names
            assert all(len(s) == 4 for s in loaded_seqs)

    def test_methods_registry_consistency(self):
        """Every method in METHODS should have a matching color."""
        from benchmarks.downstream.utils import METHOD_COLORS, METHODS

        for name in METHODS:
            assert name in METHOD_COLORS, (
                f"Method {name!r} in METHODS but not in METHOD_COLORS"
            )


# ======================================================================
# Simulation grid tests
# ======================================================================


class TestSimulationGrids:
    """Test simulation parameter generation."""

    def test_codon_grid_count(self):
        from benchmarks.downstream.simulation import CODON_GRID, iter_simulation_params

        params = list(iter_simulation_params(CODON_GRID, "M8"))
        # 3 taxa x 3 depths x 3 indel_rates x 3 psel_fracs x 3 reps = 243
        assert len(params) == 243

    def test_protein_grid_count(self):
        from benchmarks.downstream.simulation import (
            PROTEIN_GRID,
            iter_simulation_params,
        )

        params = list(iter_simulation_params(PROTEIN_GRID, "WAG"))
        # 2 taxa x 4 depths x 4 indel_rates x 2 indel_lengths x 3 reps = 192
        assert len(params) == 192

    def test_sim_ids_unique(self):
        from benchmarks.downstream.simulation import CODON_GRID, iter_simulation_params

        params = list(iter_simulation_params(CODON_GRID, "M8"))
        sim_ids = [p["sim_id"] for p in params]
        assert len(sim_ids) == len(set(sim_ids)), "sim_ids are not unique"


# ======================================================================
# CLI tests
# ======================================================================


class TestCLI:
    """Test the __main__.py CLI."""

    def test_help_flag(self):
        result = subprocess.run(
            ["uv", "run", "python", "-m", "benchmarks.downstream", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "downstream" in result.stdout.lower()
        assert "--quick" in result.stdout
        assert "--all" in result.stdout
        assert "--figures" in result.stdout

    def test_no_args_exits_nonzero(self):
        result = subprocess.run(
            ["uv", "run", "python", "-m", "benchmarks.downstream"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode != 0  # No action specified


# ======================================================================
# Provenance integration
# ======================================================================


class TestProvenanceIntegration:
    """Test provenance works end-to-end with result file management."""

    def test_full_provenance_cycle(self):
        from dataclasses import asdict

        from benchmarks.downstream.provenance import (
            collect_provenance,
            load_latest_results,
            result_path,
            update_latest_symlink,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Collect provenance
            prov = collect_provenance({"pipeline": "test", "quick": True})
            prov_dict = asdict(prov)

            # Generate result path
            path = result_path(Path(tmpdir), "test_pipeline")

            # Write result
            result_data = {
                "provenance": prov_dict,
                "cases": [{"method": "kalign", "sp_score": 0.85}],
                "summary": {"kalign": {"mean_sp": 0.85}},
            }
            path.write_text(json.dumps(result_data, indent=2))

            # Create symlink
            update_latest_symlink(path)

            # Load via symlink
            loaded = load_latest_results(path.parent)
            assert loaded["provenance"]["cpu_cores"] > 0
            assert loaded["cases"][0]["sp_score"] == 0.85

    def test_kalign_version_detected(self):
        from benchmarks.downstream.provenance import collect_provenance

        prov = collect_provenance({})
        # kalign should be installed in the dev environment
        assert prov.kalign_version != "not installed"
