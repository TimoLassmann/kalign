"""Tests for downstream benchmark foundation modules.

Verifies provenance, utils, and simulation modules can be imported
and their core functions work correctly without external tools.
"""

import json
import os
import tempfile
from pathlib import Path

import pytest

pytest.importorskip("benchmarks", reason="benchmarks package not on sys.path (run from repo root)")


# ======================================================================
# provenance.py
# ======================================================================


class TestProvenance:
    """Tests for benchmarks.downstream.provenance."""

    def test_import(self):
        from benchmarks.downstream import provenance

        assert hasattr(provenance, "Provenance")
        assert hasattr(provenance, "collect_provenance")
        assert hasattr(provenance, "collect_tool_versions")
        assert hasattr(provenance, "result_path")
        assert hasattr(provenance, "update_latest_symlink")
        assert hasattr(provenance, "load_latest_results")

    def test_provenance_dataclass_fields(self):
        from benchmarks.downstream.provenance import Provenance

        import dataclasses

        fields = {f.name for f in dataclasses.fields(Provenance)}
        expected = {
            "timestamp",
            "kalign_version",
            "kalign_commit",
            "container_image",
            "hostname",
            "cpu_model",
            "cpu_cores",
            "ram_gb",
            "os_version",
            "python_version",
            "tool_versions",
            "parameters",
        }
        assert expected == fields

    def test_collect_tool_versions_returns_dict(self):
        from benchmarks.downstream.provenance import collect_tool_versions

        versions = collect_tool_versions()
        assert isinstance(versions, dict)
        # Should have entries for all known tools
        for tool in ["kalign", "mafft", "muscle", "clustalo", "hmmer", "iqtree"]:
            assert tool in versions
            assert isinstance(versions[tool], str)

    def test_collect_provenance_structure(self):
        from benchmarks.downstream.provenance import Provenance, collect_provenance

        prov = collect_provenance({"test_param": "value"})
        assert isinstance(prov, Provenance)
        assert prov.cpu_cores > 0
        assert prov.ram_gb > 0
        assert prov.python_version  # non-empty
        assert prov.hostname  # non-empty
        assert prov.parameters == {"test_param": "value"}

    def test_provenance_json_roundtrip(self):
        from dataclasses import asdict

        from benchmarks.downstream.provenance import collect_provenance

        prov = collect_provenance({"key": 42})
        data = asdict(prov)
        # Must be JSON-serialisable
        text = json.dumps(data)
        loaded = json.loads(text)
        assert loaded["parameters"] == {"key": 42}
        assert loaded["cpu_cores"] > 0

    def test_result_path_generation(self):
        from benchmarks.downstream.provenance import result_path

        with tempfile.TemporaryDirectory() as tmpdir:
            path = result_path(Path(tmpdir), "calibration")
            assert path.parent.name == "calibration"
            assert path.name.startswith("run_")
            assert path.name.endswith(".json")
            assert path.parent.exists()

    def test_symlink_management(self):
        from benchmarks.downstream.provenance import (
            load_latest_results,
            update_latest_symlink,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a fake result file
            result_file = Path(tmpdir) / "run_test.json"
            result_file.write_text(json.dumps({"score": 42}))

            # Create symlink
            update_latest_symlink(result_file)

            link = Path(tmpdir) / "latest.json"
            assert link.is_symlink()

            # Load through symlink
            data = load_latest_results(Path(tmpdir))
            assert data == {"score": 42}

    def test_symlink_update_replaces(self):
        from benchmarks.downstream.provenance import update_latest_symlink

        with tempfile.TemporaryDirectory() as tmpdir:
            # First result
            f1 = Path(tmpdir) / "run_1.json"
            f1.write_text(json.dumps({"v": 1}))
            update_latest_symlink(f1)

            # Second result replaces
            f2 = Path(tmpdir) / "run_2.json"
            f2.write_text(json.dumps({"v": 2}))
            update_latest_symlink(f2)

            link = Path(tmpdir) / "latest.json"
            data = json.loads(link.read_text())
            assert data == {"v": 2}


# ======================================================================
# utils.py
# ======================================================================


class TestUtils:
    """Tests for benchmarks.downstream.utils."""

    def test_import(self):
        from benchmarks.downstream import utils

        assert hasattr(utils, "AlignResult")
        assert hasattr(utils, "parse_fasta")
        assert hasattr(utils, "write_fasta")
        assert hasattr(utils, "mask_alignment_by_confidence")
        assert hasattr(utils, "write_site_weights")
        assert hasattr(utils, "holm_bonferroni")
        assert hasattr(utils, "METHOD_COLORS")
        assert hasattr(utils, "METHODS")
        assert hasattr(utils, "run_method")

    def test_align_result_fields(self):
        from benchmarks.downstream.utils import AlignResult

        r = AlignResult(
            sequences=["ACGT", "AC-T"],
            names=["s1", "s2"],
            column_confidence=[0.9, 0.8, 0.7, 0.6],
            residue_confidence=[[0.9, 0.8, 0.7, 0.6], [0.8, 0.7, 0.0, 0.5]],
            wall_time=1.5,
            peak_memory_mb=100.0,
        )
        assert len(r.sequences) == 2
        assert r.wall_time == 1.5

    def test_fasta_roundtrip(self):
        from benchmarks.downstream.utils import parse_fasta, write_fasta

        names = ["seq1", "seq2", "seq3"]
        seqs = ["ACGTACGT" * 20, "TGCATGCA" * 20, "AAAA"]

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "test.fasta"
            write_fasta(names, seqs, fasta_path)

            loaded_names, loaded_seqs = parse_fasta(fasta_path)
            assert loaded_names == names
            assert loaded_seqs == seqs

    def test_mask_alignment_basic(self):
        from benchmarks.downstream.utils import mask_alignment_by_confidence

        seqs = ["ABCD", "EFGH"]
        conf = [0.9, 0.3, 0.8, 0.2]
        masked, n_kept = mask_alignment_by_confidence(seqs, conf, 0.5)
        assert n_kept == 2
        assert masked == ["AC", "EG"]

    def test_mask_alignment_empty(self):
        from benchmarks.downstream.utils import mask_alignment_by_confidence

        masked, n_kept = mask_alignment_by_confidence([], [], 0.5)
        assert n_kept == 0
        assert masked == []

    def test_mask_all_filtered(self):
        from benchmarks.downstream.utils import mask_alignment_by_confidence

        seqs = ["AB", "CD"]
        conf = [0.1, 0.2]
        masked, n_kept = mask_alignment_by_confidence(seqs, conf, 0.5)
        assert n_kept == 0
        assert masked == ["", ""]

    def test_write_site_weights(self):
        from benchmarks.downstream.utils import write_site_weights

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "weights.txt"
            conf = [0.95, 0.3, 1.5, -0.1]  # includes out-of-range
            write_site_weights(conf, path)
            lines = path.read_text().strip().split("\n")
            assert len(lines) == 4
            vals = [float(x) for x in lines]
            assert vals[0] == pytest.approx(0.95, abs=1e-5)
            assert vals[1] == pytest.approx(0.3, abs=1e-5)
            assert vals[2] == pytest.approx(1.0, abs=1e-5)  # clamped
            assert vals[3] == pytest.approx(0.0, abs=1e-5)  # clamped

    def test_holm_bonferroni_basic(self):
        from benchmarks.downstream.utils import holm_bonferroni

        # Simple case: 3 p-values
        adjusted = holm_bonferroni([0.01, 0.04, 0.03])
        assert len(adjusted) == 3
        # All adjusted values should be >= original
        for adj, orig in zip(adjusted, [0.01, 0.04, 0.03]):
            assert adj >= orig
        # All adjusted values should be <= 1.0
        for adj in adjusted:
            assert adj <= 1.0

    def test_holm_bonferroni_empty(self):
        from benchmarks.downstream.utils import holm_bonferroni

        assert holm_bonferroni([]) == []

    def test_holm_bonferroni_single(self):
        from benchmarks.downstream.utils import holm_bonferroni

        result = holm_bonferroni([0.05])
        assert len(result) == 1
        assert result[0] == pytest.approx(0.05)

    def test_method_colors_has_entries(self):
        from benchmarks.downstream.utils import METHOD_COLORS

        assert "kalign" in METHOD_COLORS
        assert "mafft" in METHOD_COLORS
        assert "guidance2_mafft" in METHOD_COLORS

    def test_methods_registry_structure(self):
        from benchmarks.downstream.utils import METHODS

        assert "kalign" in METHODS
        assert "kalign_ens3" in METHODS
        assert "mafft" in METHODS
        assert "muscle" in METHODS
        assert "clustalo" in METHODS
        assert "guidance2_mafft" in METHODS
        # All entries have "fn"
        for name, cfg in METHODS.items():
            assert "fn" in cfg, f"METHODS[{name!r}] missing 'fn'"
            assert callable(cfg["fn"]), f"METHODS[{name!r}]['fn'] not callable"

    def test_run_method_unknown_raises(self):
        from benchmarks.downstream.utils import run_method

        with pytest.raises(ValueError, match="Unknown method"):
            run_method("nonexistent_aligner", Path("/tmp/x.fa"), Path("/tmp"))


# ======================================================================
# simulation.py
# ======================================================================


class TestSimulation:
    """Tests for benchmarks.downstream.simulation."""

    def test_import(self):
        from benchmarks.downstream import simulation

        assert hasattr(simulation, "SimulatedDataset")
        assert hasattr(simulation, "random_birth_death_tree")
        assert hasattr(simulation, "generate_indelible_dataset")
        assert hasattr(simulation, "CODON_GRID")
        assert hasattr(simulation, "PROTEIN_GRID")
        assert hasattr(simulation, "iter_simulation_params")
        assert hasattr(simulation, "strip_gaps")

    def test_simulated_dataset_fields(self):
        import dataclasses

        from benchmarks.downstream.simulation import SimulatedDataset

        fields = {f.name for f in dataclasses.fields(SimulatedDataset)}
        expected = {"true_alignment", "unaligned", "true_tree", "site_classes", "params"}
        assert expected == fields

    def test_strip_gaps(self):
        from benchmarks.downstream.simulation import strip_gaps

        assert strip_gaps(["A-C.G-T", "--ACGT", "ACGT"]) == [
            "ACGT",
            "ACGT",
            "ACGT",
        ]

    def test_strip_gaps_empty(self):
        from benchmarks.downstream.simulation import strip_gaps

        assert strip_gaps([]) == []
        assert strip_gaps(["---"]) == [""]

    def test_codon_grid_structure(self):
        from benchmarks.downstream.simulation import CODON_GRID

        assert len(CODON_GRID.n_taxa) == 3
        assert len(CODON_GRID.tree_depths) == 3
        assert len(CODON_GRID.indel_rates) == 3
        assert CODON_GRID.replicates == 10

    def test_protein_grid_structure(self):
        from benchmarks.downstream.simulation import PROTEIN_GRID

        assert len(PROTEIN_GRID.n_taxa) == 3
        assert len(PROTEIN_GRID.tree_depths) == 4
        assert len(PROTEIN_GRID.indel_rates) == 3
        assert PROTEIN_GRID.replicates == 20

    def test_iter_simulation_params_wag(self):
        from benchmarks.downstream.simulation import SimulationGrid, iter_simulation_params

        small_grid = SimulationGrid(
            n_taxa=[8],
            tree_depths=[0.5],
            indel_rates=[0.05],
            indel_length_means=[2.0],
            replicates=3,
        )
        params = list(iter_simulation_params(small_grid, "WAG"))
        assert len(params) == 3  # 1x1x1x1x3
        for p in params:
            assert "sim_id" in p
            assert p["model"] == "WAG"
            assert p["n_taxa"] == 8
            assert "WAG" in p["sim_id"]

    def test_iter_simulation_params_m8(self):
        from benchmarks.downstream.simulation import SimulationGrid, iter_simulation_params

        small_grid = SimulationGrid(
            n_taxa=[16],
            tree_depths=[0.7],
            indel_rates=[0.05],
            psel_fractions=[0.0, 0.10],
            replicates=2,
        )
        params = list(iter_simulation_params(small_grid, "M8"))
        assert len(params) == 4  # 1x1x1x2x2
        for p in params:
            assert p["model"] == "M8"

    def test_iter_simulation_params_unknown_model(self):
        from benchmarks.downstream.simulation import SimulationGrid, iter_simulation_params

        small_grid = SimulationGrid(n_taxa=[8], tree_depths=[0.5], indel_rates=[0.05])
        with pytest.raises(ValueError, match="Unknown model"):
            list(iter_simulation_params(small_grid, "JTT"))

    def test_indelible_control_file_codon(self):
        """Test that codon control file is written correctly."""
        from benchmarks.downstream.simulation import _write_codon_control

        with tempfile.TemporaryDirectory() as tmpdir:
            tree = "((A:0.1,B:0.1):0.2,C:0.3);"
            path = _write_codon_control(
                Path(tmpdir), tree, n_codons=100, seed=123
            )
            assert path.exists()
            text = path.read_text()
            assert "[TYPE] CODON 1" in text
            assert "[randomseed] 123" in text
            assert "simmodel" in text
            assert "[EVOLVE]" in text

    def test_indelible_control_file_protein(self):
        """Test that protein control file is written correctly."""
        from benchmarks.downstream.simulation import _write_protein_control

        with tempfile.TemporaryDirectory() as tmpdir:
            tree = "((A:0.1,B:0.1):0.2,C:0.3);"
            path = _write_protein_control(
                Path(tmpdir), tree, seq_length=200, seed=456
            )
            assert path.exists()
            text = path.read_text()
            assert "[TYPE] AMINOACID 1" in text
            assert "[randomseed] 456" in text
            assert "WAG" in text

    def test_parse_indelible_output_phylip(self):
        """Test PHYLIP parsing with synthetic data."""
        from benchmarks.downstream.simulation import _parse_indelible_output

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write synthetic PHYLIP true alignment
            phy_path = tmpdir / "output_TRUE_1.phy"
            phy_path.write_text(
                " 3 10\n"
                "T1        ACGTACGTAC\n"
                "T2        ACGT--GTAC\n"
                "T3        AC--ACGTAC\n"
            )

            # Write synthetic FASTA unaligned
            fas_path = tmpdir / "output_1.fas"
            fas_path.write_text(
                ">T1\nACGTACGTAC\n>T2\nACGTGTAC\n>T3\nACACGTAC\n"
            )

            true_names, true_seqs, unaln_names, unaln_seqs = (
                _parse_indelible_output(tmpdir)
            )
            assert true_names == ["T1", "T2", "T3"]
            assert true_seqs[0] == "ACGTACGTAC"
            assert true_seqs[1] == "ACGT--GTAC"
            assert unaln_names == ["T1", "T2", "T3"]
            assert unaln_seqs[1] == "ACGTGTAC"

    def test_parse_site_classes_empty(self):
        """No rates file → all zeros."""
        from benchmarks.downstream.simulation import _parse_site_classes

        with tempfile.TemporaryDirectory() as tmpdir:
            classes = _parse_site_classes(Path(tmpdir), 10)
            assert classes == [0] * 10

    def test_parse_site_classes_with_data(self):
        """Synthetic rates file parsing."""
        from benchmarks.downstream.simulation import _parse_site_classes

        with tempfile.TemporaryDirectory() as tmpdir:
            rates = Path(tmpdir) / "output_1_RATES.txt"
            rates.write_text(
                "Site\tClass\tRate\n"
                "1\t0\t0.5\n"
                "2\t0\t0.3\n"
                "3\t1\t2.5\n"
                "4\t0\t0.1\n"
                "5\t1\t3.0\n"
            )
            classes = _parse_site_classes(Path(tmpdir), 5)
            assert classes == [0, 0, 1, 0, 1]


# ======================================================================
# Cross-module integration
# ======================================================================


class TestCrossModule:
    """Tests that modules work together correctly."""

    def test_provenance_in_result_json(self):
        """provenance + utils.write_fasta can produce a complete result."""
        from dataclasses import asdict

        from benchmarks.downstream.provenance import collect_provenance, result_path
        from benchmarks.downstream.utils import write_fasta

        with tempfile.TemporaryDirectory() as tmpdir:
            # Generate result path
            path = result_path(Path(tmpdir), "test_pipeline")
            assert path.parent.exists()

            # Collect provenance
            prov = collect_provenance({"method": "kalign_ens3"})
            prov_dict = asdict(prov)

            # Write a fake result
            result = {
                "provenance": prov_dict,
                "cases": [{"sp_score": 0.85, "tc_score": 0.72}],
            }
            path.write_text(json.dumps(result, indent=2))
            assert path.exists()

            # Verify it can be read back
            loaded = json.loads(path.read_text())
            assert "provenance" in loaded
            assert loaded["provenance"]["cpu_cores"] > 0

    def test_utils_fasta_compatible_with_simulation_strip_gaps(self):
        """write_fasta → parse_fasta → strip_gaps chain works."""
        from benchmarks.downstream.simulation import strip_gaps
        from benchmarks.downstream.utils import parse_fasta, write_fasta

        names = ["s1", "s2"]
        seqs = ["AC-GT-A", "A--GTCA"]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "aligned.fasta"
            write_fasta(names, seqs, path)
            loaded_names, loaded_seqs = parse_fasta(path)
            assert loaded_names == names
            assert loaded_seqs == seqs

            # Strip gaps
            ungapped = strip_gaps(loaded_seqs)
            assert ungapped == ["ACGTA", "AGTCA"]
