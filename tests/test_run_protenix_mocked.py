"""
Unit test for run_protenix_inference with mocked subprocess.
Verifies the correct CLI command is built without calling Protenix.
"""
import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

import sys
scripts_path = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(scripts_path))

from utils.protenix_eval import run_protenix_inference


class TestRunProtenixInferenceMocked:
    """Test run_protenix_inference with subprocess mocked."""

    def test_builds_mini_command(self, tmpdir):
        """Verify mini model CLI args."""
        json_path = tmpdir / "test_OFF.json"
        json_path.write_text(json.dumps([{"name": "test_OFF", "sequences": []}]))
        out_dir = tmpdir / "out"
        out_dir.mkdir()

        with patch("utils.protenix_eval.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            with patch("utils.protenix_eval.glob.glob") as mock_glob:
                mock_glob.side_effect = [
                    [str(out_dir / "test_OFF" / "model.pdb")],
                    [str(out_dir / "test_OFF" / "model_summary.json")],
                ]
                run_protenix_inference(str(json_path), str(out_dir), model_tier="mini")

        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "protenix" in call_args[0]
        assert "predict" in call_args
        assert "protenix_mini_default_v0.5.0" in call_args
        assert "bf16" in call_args
        assert "enable_fusion" in str(call_args)

    def test_builds_base_command(self, tmpdir):
        """Verify base model CLI args."""
        json_path = tmpdir / "test_ON.json"
        json_path.write_text(json.dumps([{"name": "test_ON", "sequences": []}]))
        out_dir = tmpdir / "out"
        out_dir.mkdir()

        with patch("utils.protenix_eval.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            with patch("utils.protenix_eval.glob.glob") as mock_glob:
                mock_glob.side_effect = [
                    [str(out_dir / "test_ON" / "model.pdb")],
                    [str(out_dir / "test_ON" / "model_summary.json")],
                ]
                pdb_path, summary_path = run_protenix_inference(
                    str(json_path), str(out_dir), model_tier="base"
                )
        assert "model.pdb" in pdb_path
        assert "model_summary.json" in summary_path

        call_args = mock_run.call_args[0][0]
        assert "protenix_base_default_v1.0.0" in call_args
