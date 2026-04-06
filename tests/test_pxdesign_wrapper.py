"""
Unit tests for 03_pxdesign_wrapper: generate_frozen_rec_config, X-resolution,
sequence-level linker mutator.
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
_compact_variant_name = pxd._compact_variant_name
_generate_variants_sequence_mutator = pxd._generate_variants_sequence_mutator
_pxdesign_available = pxd._pxdesign_available


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


class TestCompactVariantNaming:
    """Ensure concise stable naming across generations."""

    def test_compact_name_shape(self):
        name = _compact_variant_name("3174363721_ORF_Score_0.928", generation_num=3, variant_index=1)
        assert name.startswith("L")
        assert "_g03_v01" in name
        assert "_variant_" not in name

    def test_compact_fallback_suffix(self):
        name = _compact_variant_name("baseline_id", generation_num=12, variant_index=0, fallback=True)
        assert name.endswith("_fb")


class TestSequenceMutator:
    """Test _generate_variants_sequence_mutator (pure-sequence linker mutator)."""

    # Realistic-ish 100-aa sequence: REC[0:40] + linker1[40:60] + HEPN1[60:80] + linker2[80:90] + HEPN2[90:100]
    COORDS = {
        "rec_end": 40,
        "hepn1_start": 60,
        "hepn1_end": 80,
        "hepn2_start": 90,
        "hepn2_end": 100,
        "linker1_len": 20,  # 60 - 40
        "linker2_len": 10,  # 90 - 80
        "binder_length": 30,
        "seq_len": 100,
    }
    SEQ = "A" * 40 + "G" * 20 + "H" * 20 + "G" * 10 + "H" * 10

    def test_returns_correct_count(self):
        variants = _generate_variants_sequence_mutator(
            self.SEQ, self.COORDS, variant_count=5, bias_json_path=None,
            generation_num=1, seed=0,
        )
        assert len(variants) == 5

    def test_preserves_domains(self):
        """REC, HEPN1, and HEPN2 must be identical to baseline in every variant."""
        variants = _generate_variants_sequence_mutator(
            self.SEQ, self.COORDS, variant_count=10, bias_json_path=None,
            generation_num=3, seed=123,
        )
        for v in variants:
            assert len(v) == len(self.SEQ), "Variant length must match baseline"
            assert v[:40] == self.SEQ[:40], "REC domain modified"
            assert v[60:80] == self.SEQ[60:80], "HEPN1 domain modified"
            assert v[90:100] == self.SEQ[90:100], "HEPN2 domain modified"

    def test_mutations_only_in_linkers(self):
        """All mutations should be in linker1 [40:60) or linker2 [80:90)."""
        variants = _generate_variants_sequence_mutator(
            self.SEQ, self.COORDS, variant_count=20, bias_json_path=None,
            generation_num=5, seed=42,
        )
        for v in variants:
            for i, (orig, new) in enumerate(zip(self.SEQ, v)):
                if orig != new:
                    assert (40 <= i < 60) or (80 <= i < 90), \
                        f"Mutation at position {i} is outside linker regions"

    def test_variants_are_diverse(self):
        """Not all variants should be identical (unless mutation rate is 0)."""
        variants = _generate_variants_sequence_mutator(
            self.SEQ, self.COORDS, variant_count=10, bias_json_path=None,
            generation_num=1, seed=42,
        )
        unique = set(variants)
        assert len(unique) > 1, "All variants are identical — no diversity"

    def test_bias_influences_mutations(self, tmpdir):
        """When bias strongly prefers W at a linker position, that position should be W."""
        bias_file = tmpdir / "bias.json"
        # Position 45 (1-based=46) is in linker1 [40:60); strong W bias
        bias_file.write_text(json.dumps({"46": {"W": 100.0}}))
        variants = _generate_variants_sequence_mutator(
            self.SEQ, self.COORDS, variant_count=50, bias_json_path=str(bias_file),
            generation_num=1, seed=0,
        )
        # At least some variants should have W at position 45
        w_count = sum(1 for v in variants if v[45] == "W")
        assert w_count > 0, "Bias for W at pos 45 had no effect"

    def test_deterministic_with_same_seed(self):
        v1 = _generate_variants_sequence_mutator(
            self.SEQ, self.COORDS, variant_count=5, bias_json_path=None,
            generation_num=1, seed=999,
        )
        v2 = _generate_variants_sequence_mutator(
            self.SEQ, self.COORDS, variant_count=5, bias_json_path=None,
            generation_num=1, seed=999,
        )
        assert v1 == v2


class TestPxdesignAvailability:
    """Test _pxdesign_available detection."""

    def test_returns_bool(self):
        result = _pxdesign_available()
        assert isinstance(result, bool)
