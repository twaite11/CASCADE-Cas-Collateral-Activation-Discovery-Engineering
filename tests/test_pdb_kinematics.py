"""
Unit tests for pdb_kinematics: HEPN distance calculation and Protenix score extraction.
No GPU required.
"""
import json
import pytest
from pathlib import Path

from utils.pdb_kinematics import calculate_hepn_shift, extract_protenix_scores


class TestExtractProtenixScores:
    """Test extract_protenix_scores."""

    def test_extracts_all_fields(self, protenix_summary_json):
        scores = extract_protenix_scores(protenix_summary_json)
        assert scores["iptm"] == 0.92
        assert scores["ptm"] == 0.88
        assert scores["ranking_score"] == 0.9
        assert scores["af2_ig"] == 0.85

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError, match="not found"):
            extract_protenix_scores("/nonexistent/path/summary.json")

    def test_handles_missing_keys(self, tmpdir):
        p = tmpdir / "minimal.json"
        p.write_text("{}")
        scores = extract_protenix_scores(str(p))
        assert scores["iptm"] == 0.0
        assert scores["af2_ig"] == 0.0

    def test_af2_ig_score_fallback(self, tmpdir):
        """af2_ig can be stored as af2_ig_score."""
        p = tmpdir / "alt_key.json"
        p.write_text('{"iptm": 0.5, "af2_ig_score": 0.77}')
        scores = extract_protenix_scores(str(p))
        assert scores["af2_ig"] == 0.77


class TestCalculateHepnShift:
    """Test calculate_hepn_shift."""

    def test_returns_30_angstroms(self, minimal_pdb):
        dist = calculate_hepn_shift(minimal_pdb, 10, 50)
        assert abs(dist - 30.0) < 0.01

    def test_returns_15_angstroms(self, minimal_pdb_15a):
        dist = calculate_hepn_shift(minimal_pdb_15a, 10, 50)
        assert abs(dist - 15.0) < 0.01

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError, match="not found"):
            calculate_hepn_shift("/nonexistent/model.pdb", 10, 50)

    def test_invalid_residue_raises(self, minimal_pdb):
        with pytest.raises(ValueError, match="Could not find"):
            calculate_hepn_shift(minimal_pdb, 999, 50)

    def test_custom_chain_id(self, minimal_pdb):
        """Default chain A works; other chains would need different PDB."""
        dist = calculate_hepn_shift(minimal_pdb, 10, 50, protein_chain_id="A")
        assert abs(dist - 30.0) < 0.01
