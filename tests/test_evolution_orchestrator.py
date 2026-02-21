"""
Unit tests for evolution_orchestrator: fitness, EvolutionGym, mutations, HEPN indices.
No GPU/Protenix/PXDesign. Mocks subprocess for run_protenix_inference.
"""
import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

# Import orchestrator components (avoid loading full loop)
import sys
scripts_path = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(scripts_path))

# Import only what we need - avoid running main
import importlib.util
spec = importlib.util.spec_from_file_location(
    "evolution_orchestrator",
    scripts_path / "evolution_orchestrator.py",
)
orch = importlib.util.module_from_spec(spec)
# Prevent side effects during import (e.g. reading METADATA_FILE)
with patch.dict("os.environ", {}):
    spec.loader.exec_module(orch)

compute_fitness = orch.compute_fitness
EvolutionGym = orch.EvolutionGym
get_catalytic_histidine_indices = orch.get_catalytic_histidine_indices
extract_mutations = orch.extract_mutations
build_metadata_override_for_evolved = orch.build_metadata_override_for_evolved


class TestComputeFitness:
    """Test compute_fitness."""

    def test_high_off_low_on_is_good(self):
        # OFF=30A, ON=10A, shift=20. Good ternary.
        f = compute_fitness(30, 10, 0.9, 0.85, is_full_ternary=True)
        assert f > 0

    def test_low_off_is_bad(self):
        f = compute_fitness(10, 8, 0.9, 0.85, is_full_ternary=True)
        assert f < compute_fitness(30, 8, 0.9, 0.85, is_full_ternary=True)

    def test_specificity_penalty(self):
        """Activity on off-targets should reduce fitness."""
        base = compute_fitness(30, 10, 0.9, 0.85, is_full_ternary=True, offtarget_by_mismatch=None)
        with_penalty = compute_fitness(
            30, 10, 0.9, 0.85, is_full_ternary=True,
            offtarget_by_mismatch={1: 20.0, 2: 15.0, 3: 10.0}  # 3mm active at 10A is bad
        )
        assert with_penalty < base

    def test_ternary_multiplier(self):
        f_mini = compute_fitness(30, 10, 0.9, 0.85, is_full_ternary=False)
        f_full = compute_fitness(30, 10, 0.9, 0.85, is_full_ternary=True)
        assert f_full == 2.0 * f_mini


class TestEvolutionGym:
    """Test EvolutionGym."""

    def test_register_and_bias(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        gym = EvolutionGym()
        gym.register_evaluation("v1", ["10_P", "20_A"], 30, 10, 0.9, is_full_ternary=True)
        gym.register_evaluation("v2", ["10_P", "30_G"], 28, 11, 0.85, is_full_ternary=True)
        assert len(gym.generation_history) == 2
        assert "10_P" in gym.mutation_weights
        bias_file = gym.generate_mpnn_bias_matrix(1)
        assert Path(bias_file).exists()
        data = json.loads(Path(bias_file).read_text())
        assert "10" in data

    def test_bias_clipping(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        gym = EvolutionGym()
        # Extremely good fitness to test clipping
        for i in range(20):
            gym.register_evaluation(f"v{i}", ["100_X"], 50, 5, 0.99, is_full_ternary=True)
        bias_file = gym.generate_mpnn_bias_matrix(1)
        data = json.loads(Path(bias_file).read_text())
        assert data["100"]["X"] <= 5.0 and data["100"]["X"] >= -5.0


class TestGetCatalyticHistidineIndices:
    """Test get_catalytic_histidine_indices."""

    def test_returns_indices_for_two_motifs(self, variant_fasta):
        h1, h2 = get_catalytic_histidine_indices(variant_fasta)
        assert h1 is not None and h2 is not None
        assert h1 < h2

    def test_returns_none_for_one_motif(self, tmpdir):
        p = tmpdir / "one_hepn.fasta"
        p.write_text(">x\nMAAAAARAILXHGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
        h1, h2 = get_catalytic_histidine_indices(str(p))
        assert h1 is None and h2 is None

    def test_returns_none_for_no_motif(self, tmpdir):
        p = tmpdir / "no_hepn.fasta"
        p.write_text(">x\nMGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
        h1, h2 = get_catalytic_histidine_indices(str(p))
        assert h1 is None and h2 is None


class TestExtractMutations:
    """Test extract_mutations."""

    def test_single_substitution(self, variant_fasta, baseline_json, tmpdir):
        # variant has P at position 12 (1-based) vs baseline G
        base_dir = Path(baseline_json).parent
        with patch.object(orch, "BASE_JSON_DIR", str(base_dir)):
            muts = extract_mutations("baseline_id", variant_fasta)
        assert len(muts) >= 1
        assert any("_P" in m or "12_" in m for m in muts)

    def test_uses_baseline_fasta_when_provided(self, variant_fasta, tmpdir):
        baseline_fasta = tmpdir / "baseline.fasta"
        baseline_fasta.write_text(">baseline\nMAAAAARAILXHGGGGGGGGGGGGGGGGGGGRVVVXHGGGGGGGG\n")
        variant = tmpdir / "v.fasta"
        variant.write_text(">v\nMAAAAARAILXHPGGGGGGGGGGGGGGGGGRVVVXHGGGGGGGG\n")  # 12 G->P
        muts = extract_mutations("x", str(variant), baseline_fasta_path=str(baseline_fasta))
        assert any("12_P" in m for m in muts)

    def test_detects_deletion(self, tmpdir):
        base = tmpdir / "base.fasta"
        base.write_text(">b\nMABCDEF\n")
        var = tmpdir / "var.fasta"
        var.write_text(">v\nMABDEF\n")  # C deleted at 4
        muts = extract_mutations("b", str(var), baseline_fasta_path=str(base))
        assert any("del" in m for m in muts)


class TestBuildMetadataOverride:
    """Test build_metadata_override_for_evolved."""

    def test_returns_override_dict(self, tmpdir, sample_metadata_json):
        baseline_fasta = tmpdir / "evolved.fasta"
        baseline_fasta.write_text(">evolved\nMAAAAARAILXHGGGGGGGGGGGGGGGGGGGRVVVXHGGGGGGGG\n")
        meta = json.loads(Path(sample_metadata_json).read_text())
        override = build_metadata_override_for_evolved(
            "evolved_id", str(baseline_fasta), "test_cas13", meta
        )
        assert override is not None
        assert "evolved_id" in override
        assert "domains" in override["evolved_id"]
        assert "HEPN1" in override["evolved_id"]["domains"]

    def test_returns_none_for_bad_seq(self, tmpdir, sample_metadata_json):
        bad_fasta = tmpdir / "bad.fasta"
        bad_fasta.write_text(">bad\nMGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
        meta = json.loads(Path(sample_metadata_json).read_text())
        override = build_metadata_override_for_evolved("bad", str(bad_fasta), "test_cas13", meta)
        assert override is None
