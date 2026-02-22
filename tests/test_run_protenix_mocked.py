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
                # protenix_eval prefers CIF, then PDB (recursive then flat)
                pdb_path = str(out_dir / "test_OFF" / "model.pdb")
                summary_path = str(out_dir / "test_OFF" / "model_summary.json")
                mock_glob.side_effect = [
                    [], [],  # no cif (recursive, flat)
                    [pdb_path], [],  # pdb recursive found; summary (flat, recursive)
                    [summary_path], [],
                ]
                run_protenix_inference(str(json_path), str(out_dir), model_tier="mini")

        mock_run.assert_called()
        # May be called twice (msa + predict) or once (predict only)
        predict_calls = [c for c in mock_run.call_args_list if c[0][0] and "predict" in c[0][0]]
        assert len(predict_calls) >= 1
        call_args = predict_calls[-1][0][0]
        assert "protenix" in call_args[0]
        assert "predict" in call_args
        assert "protenix_mini_default_v0.5.0" in call_args
        assert "use_default_params" in str(call_args)

    def test_builds_base_command(self, tmpdir):
        """Verify base model CLI args."""
        json_path = tmpdir / "test_ON.json"
        json_path.write_text(json.dumps([{"name": "test_ON", "sequences": []}]))
        out_dir = tmpdir / "out"
        out_dir.mkdir()

        with patch("utils.protenix_eval.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            with patch("utils.protenix_eval.glob.glob") as mock_glob:
                pdb_path = str(out_dir / "test_ON" / "model.pdb")
                summary_path = str(out_dir / "test_ON" / "model_summary.json")
                mock_glob.side_effect = [
                    [], [], [pdb_path], [], [summary_path], [],
                ]
                struct_path, sum_path = run_protenix_inference(
                    str(json_path), str(out_dir), model_tier="base"
                )
        assert "model.pdb" in struct_path or "model.cif" in struct_path
        assert "model_summary.json" in sum_path

        call_args = mock_run.call_args[0][0]
        assert "protenix_base_default_v1.0.0" in call_args
