"""
Unit tests for 03_pxdesign_wrapper: generate_frozen_rec_config only.
No PXDesign/GPU - we do NOT call run_pxdesign_generation.
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
    """Test generate_frozen_rec_config (no subprocess)."""

    def test_creates_freeze_json(self, sample_metadata_json, tmpdir):
        path = generate_frozen_rec_config(
            sample_metadata_json, "test_cas13", str(tmpdir)
        )
        assert Path(path).exists()
        data = json.loads(Path(path).read_text())
        assert "freeze_residues" in data
        assert data["freeze_residues"].startswith("A1-A")
        assert "description" in data

    def test_respects_hepn1_boundary(self, sample_metadata_json, tmpdir):
        """REC should end 10aa before HEPN1 start."""
        meta = json.loads(Path(sample_metadata_json).read_text())
        meta["test_cas13"]["domains"]["HEPN1"]["start"] = 360
        meta_path = tmpdir / "meta.json"
        meta_path.write_text(json.dumps(meta))
        path = generate_frozen_rec_config(str(meta_path), "test_cas13", str(tmpdir))
        data = json.loads(Path(path).read_text())
        # rec_end = max(0, 360 - 10) = 350
        assert "A350" in data["freeze_residues"] or "A350" in str(data)

    def test_metadata_override(self, tmpdir):
        """When metadata_override is provided, use it instead of file."""
        override = {
            "custom_variant": {
                "domains": {"HEPN1": {"start": 200, "end": 280}},
            }
        }
        path = generate_frozen_rec_config(
            "/nonexistent.json", "custom_variant", str(tmpdir),
            metadata_override=override
        )
        data = json.loads(Path(path).read_text())
        assert "A190" in data["freeze_residues"]  # 200 - 10

    def test_missing_variant_raises(self, sample_metadata_json, tmpdir):
        with pytest.raises(ValueError, match="not found"):
            generate_frozen_rec_config(
                sample_metadata_json, "nonexistent_id", str(tmpdir)
            )
