"""
Unit tests for evolution_orchestrator: fitness, EvolutionGym, mutations, HEPN indices,
worker isolation, and GPU lock behavior.
No GPU/Protenix/PXDesign. Mocks subprocess for run_protenix_inference.
"""
import json
import os
import threading
import time
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
_DummyLock = orch._DummyLock
_update_population = orch._update_population
_tournament_select = orch._tournament_select


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
    """Test EvolutionGym with relative fitness normalization."""

    def test_register_flush_and_bias(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        gym = EvolutionGym()
        gym.register_evaluation("v1", ["10_P", "20_A"], 30, 10, 0.9, is_full_ternary=True)
        gym.register_evaluation("v2", ["10_P", "30_G"], 28, 11, 0.85, is_full_ternary=True)
        gym.flush_generation()
        assert len(gym.generation_history) == 2
        assert "10_P" in gym.mutation_weights
        bias_file = gym.generate_mpnn_bias_matrix(1)
        assert Path(bias_file).exists()
        data = json.loads(Path(bias_file).read_text())
        assert "10" in data

    def test_relative_fitness_produces_positive_weights(self, tmpdir, monkeypatch):
        """Even with deeply negative absolute fitness, the best variant should
        get positive mutation weights via relative normalization."""
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        gym = EvolutionGym()
        gym.register_evaluation("good", ["50_L"], 40, 30, 0.13, is_full_ternary=False)
        gym.register_evaluation("bad", ["50_A"], 40, 50, 0.11, is_full_ternary=False)
        gym.flush_generation()
        assert gym.mutation_weights["50_L"] > 0, "Best variant's mutation should have positive weight"
        assert gym.mutation_weights["50_A"] < 0, "Worst variant's mutation should have negative weight"

    def test_baseline_reference_shifts_center(self, tmpdir, monkeypatch):
        """When baseline fitness is set, weights are relative to that reference."""
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        gym = EvolutionGym()
        gym.set_baseline_fitness(-50.0)
        gym.register_evaluation("better", ["10_W"], 30, 10, 0.9, is_full_ternary=True)
        gym.flush_generation()
        assert gym.mutation_weights["10_W"] > 0, "Variant better than baseline should be positive"

    def test_bias_clipping(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        gym = EvolutionGym()
        for i in range(20):
            gym.register_evaluation(f"v{i}", ["100_X"], 50, 5, 0.99, is_full_ternary=True)
        gym.flush_generation()
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
        # variant has P at position 13 (1-based) vs baseline G
        base_dir = Path(baseline_json).parent
        with patch.object(orch, "BASE_JSON_DIR", str(base_dir)):
            muts = extract_mutations("baseline_id", variant_fasta)
        assert len(muts) >= 1
        assert any("_P" in m or "13_" in m for m in muts)

    def test_uses_baseline_fasta_when_provided(self, variant_fasta, tmpdir):
        pad = "G" * 200
        base_seq = f"M{'A' * 109}RAILXH{pad}{'A' * 94}RVVVXH{'G' * 100}"
        var_seq = base_seq[:112] + "P" + base_seq[113:]
        baseline_fasta = tmpdir / "baseline.fasta"
        baseline_fasta.write_text(f">baseline\n{base_seq}\n")
        variant = tmpdir / "v.fasta"
        variant.write_text(f">v\n{var_seq}\n")
        muts = extract_mutations("x", str(variant), baseline_fasta_path=str(baseline_fasta))
        assert any("_P" in m for m in muts)

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
        pad = "G" * 200
        evolved_seq = f"M{'A' * 109}RAILXH{pad}{'A' * 94}RVVVXH{'G' * 100}"
        baseline_fasta = tmpdir / "evolved.fasta"
        baseline_fasta.write_text(f">evolved\n{evolved_seq}\n")
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
        bad_fasta.write_text(f">bad\nM{'G' * 499}\n")
        meta = json.loads(Path(sample_metadata_json).read_text())
        override = build_metadata_override_for_evolved("bad", str(bad_fasta), "test_cas13", meta)
        assert override is None


class TestWorkerIsolation:
    """Test that EvolutionGym with different worker_ids creates isolated directories."""

    def test_separate_gym_dirs(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        monkeypatch.setattr(orch, "NUM_WORKERS", 3)
        gyms = [EvolutionGym(worker_id=i) for i in range(3)]
        dirs = [g.gym_dir for g in gyms]
        assert len(set(dirs)) == 3, "Each worker should get a unique gym directory"
        for d in dirs:
            assert os.path.isdir(d)

    def test_bias_files_isolated(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        monkeypatch.setattr(orch, "NUM_WORKERS", 2)
        g0 = EvolutionGym(worker_id=0)
        g1 = EvolutionGym(worker_id=1)
        g0.register_evaluation("v0", ["10_P"], 30, 10, 0.9, is_full_ternary=True)
        g0.flush_generation()
        g1.register_evaluation("v1", ["20_A"], 25, 12, 0.8, is_full_ternary=True)
        g1.flush_generation()
        bf0 = g0.generate_mpnn_bias_matrix(1)
        bf1 = g1.generate_mpnn_bias_matrix(1)
        assert bf0 != bf1, "Bias files should be in different directories"
        assert os.path.isfile(bf0)
        assert os.path.isfile(bf1)

    def test_single_worker_uses_base_dir(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        monkeypatch.setattr(orch, "NUM_WORKERS", 1)
        g = EvolutionGym(worker_id=0)
        assert g.gym_dir == str(tmpdir / "gym"), "Single-worker mode should use base GYM_DIR"


class TestGpuLockBehavior:
    """Test that the GPU lock serializes access correctly."""

    def test_dummy_lock_is_reentrant(self):
        lock = _DummyLock()
        with lock:
            with lock:
                pass  # Should not deadlock

    def test_dummy_lock_context_manager(self):
        lock = _DummyLock()
        with lock:
            result = 42
        assert result == 42

    def test_lock_serializes_access(self):
        """Verify a real threading lock prevents concurrent access."""
        lock = threading.Lock()
        concurrent_count = []
        active = [0]

        def worker():
            with lock:
                active[0] += 1
                concurrent_count.append(active[0])
                time.sleep(0.05)
                active[0] -= 1

        threads = [threading.Thread(target=worker) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert max(concurrent_count) == 1, "Lock should serialize to at most 1 concurrent entry"


class TestAggregateResults:
    """Test the aggregation of per-worker RL datasets."""

    def test_merge_worker_datasets(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        monkeypatch.setattr(orch, "RL_TRAINING_DATASET", str(tmpdir / "gym" / "rl_training_dataset.jsonl"))
        w0 = tmpdir / "gym" / "worker_0"
        w1 = tmpdir / "gym" / "worker_1"
        w0.mkdir(parents=True)
        w1.mkdir(parents=True)
        (w0 / "rl_training_dataset.jsonl").write_text(
            json.dumps({"variant_id": "v0", "fitness": 1.0}) + "\n"
            + json.dumps({"variant_id": "v1", "fitness": 2.0}) + "\n"
        )
        (w1 / "rl_training_dataset.jsonl").write_text(
            json.dumps({"variant_id": "v2", "fitness": 3.0}) + "\n"
        )
        orch._aggregate_worker_results()
        merged = (tmpdir / "gym" / "rl_training_dataset.jsonl").read_text()
        lines = [l for l in merged.strip().split("\n") if l]
        assert len(lines) == 3
        ids = [json.loads(l)["variant_id"] for l in lines]
        assert set(ids) == {"v0", "v1", "v2"}

    def test_noop_when_no_workers(self, tmpdir, monkeypatch):
        monkeypatch.setattr(orch, "GYM_DIR", str(tmpdir / "gym"))
        monkeypatch.setattr(orch, "RL_TRAINING_DATASET", str(tmpdir / "gym" / "rl_training_dataset.jsonl"))
        (tmpdir / "gym").mkdir(parents=True, exist_ok=True)
        orch._aggregate_worker_results()
        assert not (tmpdir / "gym" / "rl_training_dataset.jsonl").exists()


def _make_individual(name, fitness, **kwargs):
    """Helper to create a population individual dict for tests."""
    d = {"name": name, "fasta": f"/tmp/{name}.fasta", "fitness": fitness,
         "off": 30.0, "on": 10.0, "iptm": 0.8, "af2_ig": 0.7,
         "hf_pdb": None, "crrna_lid": "test", "offtarget": None}
    d.update(kwargs)
    return d


class TestPopulationManagement:
    """Test _update_population merges and truncates correctly."""

    def test_keeps_top_k(self):
        pop = [_make_individual("a", 10.0), _make_individual("b", 5.0)]
        new = [_make_individual("c", 15.0), _make_individual("d", 1.0)]
        result = _update_population(pop, new, max_size=3)
        assert len(result) == 3
        assert result[0]["name"] == "c"
        assert result[1]["name"] == "a"
        assert result[2]["name"] == "b"

    def test_empty_population_fills_from_results(self):
        result = _update_population([], [_make_individual("x", 7.0)], max_size=3)
        assert len(result) == 1
        assert result[0]["name"] == "x"

    def test_existing_population_not_displaced_by_worse(self):
        pop = [_make_individual("a", 20.0), _make_individual("b", 15.0)]
        new = [_make_individual("c", 1.0)]
        result = _update_population(pop, new, max_size=2)
        names = [r["name"] for r in result]
        assert "c" not in names

    def test_sorted_descending(self):
        pop = [_make_individual("low", 1.0)]
        new = [_make_individual("high", 100.0), _make_individual("mid", 50.0)]
        result = _update_population(pop, new, max_size=5)
        fitnesses = [r["fitness"] for r in result]
        assert fitnesses == sorted(fitnesses, reverse=True)


class TestTournamentSelection:
    """Test _tournament_select picks from population with fitness pressure."""

    def test_single_member(self):
        pop = [_make_individual("only", 5.0)]
        winner = _tournament_select(pop, k=2)
        assert winner["name"] == "only"

    def test_always_picks_best_when_k_equals_population(self):
        pop = [_make_individual("best", 100.0), _make_individual("worst", 1.0)]
        for _ in range(20):
            winner = _tournament_select(pop, k=2)
            assert winner["name"] == "best"

    def test_returns_valid_member(self):
        pop = [_make_individual(f"v{i}", float(i)) for i in range(5)]
        names = {p["name"] for p in pop}
        for _ in range(30):
            winner = _tournament_select(pop, k=2)
            assert winner["name"] in names
