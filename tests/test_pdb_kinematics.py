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

    def test_af2_ig_derived_from_chain_pair_iptm(self, tmpdir):
        """When af2_ig key is missing, derive from chain_pair_iptm off-diagonal."""
        data = {
            "iptm": 0.85,
            "ptm": 0.80,
            "ranking_score": 0.82,
            "chain_pair_iptm": [
                [1.0, 0.7, 0.6],
                [0.7, 1.0, 0.5],
                [0.6, 0.5, 1.0],
            ],
        }
        p = tmpdir / "chain_pair.json"
        p.write_text(json.dumps(data))
        scores = extract_protenix_scores(str(p))
        expected = (0.7 + 0.6 + 0.7 + 0.5 + 0.6 + 0.5) / 6.0
        assert abs(scores["af2_ig"] - expected) < 0.001

    def test_af2_ig_iptm_ptm_proxy_fallback(self, tmpdir):
        """When no chain_pair_iptm, fall back to 0.8*iptm + 0.2*ptm."""
        data = {"iptm": 0.90, "ptm": 0.80}
        p = tmpdir / "proxy.json"
        p.write_text(json.dumps(data))
        scores = extract_protenix_scores(str(p))
        expected = 0.8 * 0.90 + 0.2 * 0.80
        assert abs(scores["af2_ig"] - expected) < 0.001


class TestCalculateHepnShift:
    """Test calculate_hepn_shift — uses NE2 side-chain atoms when available."""

    def test_uses_ne2_returns_26_angstroms(self, minimal_pdb):
        """With NE2 atoms present, should measure NE2-NE2 = 26 A (not CA-CA = 30 A)."""
        dist = calculate_hepn_shift(minimal_pdb, 10, 50)
        assert abs(dist - 26.0) < 0.01

    def test_uses_ne2_returns_11_angstroms(self, minimal_pdb_15a):
        """NE2-NE2 = 11 A (not CA-CA = 15 A)."""
        dist = calculate_hepn_shift(minimal_pdb_15a, 10, 50)
        assert abs(dist - 11.0) < 0.01

    def test_falls_back_to_ca(self, minimal_pdb_ca_only):
        """When no side-chain atoms exist, should fall back to CA-CA = 30 A."""
        dist = calculate_hepn_shift(minimal_pdb_ca_only, 10, 50)
        assert abs(dist - 30.0) < 0.01

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError, match="not found"):
            calculate_hepn_shift("/nonexistent/model.pdb", 10, 50)

    def test_invalid_residue_raises(self, minimal_pdb):
        with pytest.raises(ValueError, match="Could not find"):
            calculate_hepn_shift(minimal_pdb, 999, 50)

    def test_custom_chain_id(self, minimal_pdb):
        dist = calculate_hepn_shift(minimal_pdb, 10, 50, protein_chain_id="A")
        assert abs(dist - 26.0) < 0.01
