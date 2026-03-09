"""
Unit tests for 03_pxdesign_wrapper: generate_frozen_rec_config and X-resolution.
No PXDesign/GPU - we do NOT call run_pxdesign_generation.

All coordinates are 0-based Python-slice convention:
    linker1_len = hepn1_start - rec_end  (no -1)
    linker2_len = hepn2_start - hepn1_end  (no -1)
"""
import json
import pytest
from pathlib import Path

import sys
scripts_path = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(scripts_path))

from importlib.util import spec_from_file_location, module_from_spec
spec = spec_from_file_location("pxdesign_wrapper", scripts_path / "03_pxdesign_wrapper.py")
pxd = module_from_spec(spec)
spec.loader.exec_module(pxd)

generate_frozen_rec_config = pxd.generate_frozen_rec_config
_resolve_unknown_residues = pxd._resolve_unknown_residues
_apply_bias_to_sequence = pxd._apply_bias_to_sequence


class TestGenerateFrozenRecConfig:
    """Test generate_frozen_rec_config (returns coords dict)."""

    def test_returns_coords_dict(self, sample_metadata_json):
        coords = generate_frozen_rec_config(sample_metadata_json, "test_cas13")
        assert isinstance(coords, dict)
        assert "rec_end" in coords
        assert "binder_length" in coords
        assert "hepn1_start" in coords
        assert "hepn2_end" in coords
        assert "linker1_len" in coords
        assert "linker2_len" in coords
        assert coords["rec_end"] >= 0
        assert coords["binder_length"] >= 20  # MIN_BINDER_LENGTH

    def test_respects_hepn1_boundary(self, sample_metadata_json, tmpdir):
        """REC should end 10aa before HEPN1 start. binder_length = linkers only (0-based)."""
        meta = json.loads(Path(sample_metadata_json).read_text())
        meta["test_cas13"]["domains"]["HEPN1"] = {"start": 360, "end": 400}
        meta["test_cas13"]["domains"]["HEPN2"] = {"start": 420, "end": 480}
        meta["test_cas13"]["sequence_length"] = 500
        meta_path = tmpdir / "meta.json"
        meta_path.write_text(json.dumps(meta))
        coords = generate_frozen_rec_config(str(meta_path), "test_cas13")
        assert coords["rec_end"] == 350   # 360 - 10
        assert coords["linker1_len"] == 10  # 360 - 350  (0-based, no -1)
        assert coords["linker2_len"] == 20  # 420 - 400  (0-based, no -1)
        assert coords["binder_length"] == 30  # 10 + 20

    def test_metadata_override(self, tmpdir):
        """When metadata_override is provided, use it instead of file."""
        override = {
            "custom_variant": {
                "domains": {"HEPN1": {"start": 200, "end": 280}, "HEPN2": {"start": 300, "end": 380}},
                "sequence_length": 400,
            }
        }
        coords = generate_frozen_rec_config(
            "/nonexistent.json", "custom_variant",
            metadata_override=override
        )
        assert coords["rec_end"] == 190   # 200 - 10
        assert coords["linker1_len"] == 10  # 200 - 190
        assert coords["linker2_len"] == 20  # 300 - 280
        assert coords["binder_length"] == 30  # 10 + 20

    def test_missing_variant_raises(self, sample_metadata_json):
        with pytest.raises(ValueError, match="not found"):
            generate_frozen_rec_config(sample_metadata_json, "nonexistent_id")


class TestResolveUnknownResidues:
    """Test _resolve_unknown_residues."""

    def test_no_x_unchanged(self):
        assert _resolve_unknown_residues("ACDEFGH") == "ACDEFGH"

    def test_x_replaced_by_baseline(self):
        assert _resolve_unknown_residues("AXDEFGH", "ACDEFGH") == "ACDEFGH"

    def test_x_replaced_by_glycine_when_no_baseline(self):
        assert _resolve_unknown_residues("AXD") == "AGD"

    def test_multiple_x(self):
        result = _resolve_unknown_residues("XXXX", "ACDE")
        assert result == "ACDE"

    def test_x_at_end_past_baseline(self):
        result = _resolve_unknown_residues("ACX", "AC")
        assert result == "ACG"  # Past baseline → Glycine


class TestApplyBiasToSequence:
    """Test _apply_bias_to_sequence."""

    def test_no_bias_file(self):
        seq = "AAAAAAAAAA"
        assert _apply_bias_to_sequence(seq, None, {}) == seq

    def test_bias_applied_in_linker(self, tmpdir):
        bias_file = tmpdir / "bias.json"
        # Position 3 (1-based) = index 2 → in linker1 region
        bias_file.write_text(json.dumps({"3": {"W": 3.0}}))
        coords = {"rec_end": 0, "hepn1_start": 5, "hepn1_end": 8, "hepn2_start": 10}
        seq = "AAAAAHHHAABB"
        result = _apply_bias_to_sequence(seq, str(bias_file), coords)
        assert result[2] == "W"

    def test_bias_not_applied_in_hepn(self, tmpdir):
        bias_file = tmpdir / "bias.json"
        # Position 7 (1-based) = index 6 → in HEPN1 region [5:8]
        bias_file.write_text(json.dumps({"7": {"W": 3.0}}))
        coords = {"rec_end": 0, "hepn1_start": 5, "hepn1_end": 8, "hepn2_start": 10}
        seq = "AAAAAHHHAABB"
        result = _apply_bias_to_sequence(seq, str(bias_file), coords)
        assert result[6] == "H"  # Unchanged: inside HEPN
