"""
Unit tests for 03_pxdesign_wrapper: generate_frozen_rec_config only.
No PXDesign/GPU - we do NOT call run_pxdesign_generation.


TODO: Add tests per generation loop through determining if the steering is improving the fitness score per run.
Are we steering in the right direction per bias addition after a full generative single pass through the loop and does that 
magnitude of improvement increase with each addition of bias modificiation.

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
        """REC should end 10aa before HEPN1 start. binder_length = linkers only."""
        meta = json.loads(Path(sample_metadata_json).read_text())
        meta["test_cas13"]["domains"]["HEPN1"] = {"start": 360, "end": 400}
        meta["test_cas13"]["domains"]["HEPN2"] = {"start": 420, "end": 480}
        meta["test_cas13"]["sequence_length"] = 500
        meta_path = tmpdir / "meta.json"
        meta_path.write_text(json.dumps(meta))
        coords = generate_frozen_rec_config(str(meta_path), "test_cas13")
        assert coords["rec_end"] == 350  # 360 - 10
        assert coords["linker1_len"] == 9  # 360 - 350 - 1
        assert coords["linker2_len"] == 19  # 420 - 400 - 1
        assert coords["binder_length"] == 28  # 9 + 19

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
        assert coords["rec_end"] == 190  # 200 - 10
        assert coords["linker1_len"] == 9
        assert coords["linker2_len"] == 19
        assert coords["binder_length"] == 28

    def test_missing_variant_raises(self, sample_metadata_json):
        with pytest.raises(ValueError, match="not found"):
            generate_frozen_rec_config(sample_metadata_json, "nonexistent_id")
