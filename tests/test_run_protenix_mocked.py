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
        pdb_path = str(out_dir / "test_OFF" / "model.pdb")
        summary_path = str(out_dir / "test_OFF" / "model_summary.json")

        with patch("utils.protenix_eval.subprocess.run") as mock_run, \
             patch("utils.protenix_eval._find_cached_outputs", return_value=(None, None)), \
             patch("utils.protenix_eval._find_structure_and_summary", return_value=(pdb_path, summary_path)), \
             patch("utils.protenix_eval._run_msa_step", return_value=str(json_path)):
            mock_run.return_value = MagicMock(returncode=0)
            run_protenix_inference(str(json_path), str(out_dir), model_tier="mini")

        mock_run.assert_called()
        predict_calls = [c for c in mock_run.call_args_list if c[0][0] and ("pred" in c[0][0] or "predict" in c[0][0])]
        assert len(predict_calls) >= 1
        call_args = predict_calls[-1][0][0]
        assert "protenix" in call_args[0]
        assert "pred" in call_args or "predict" in call_args
        assert "protenix_mini_default_v0.5.0" in call_args
        assert "use_default_params" in str(call_args)
        msa_idx = call_args.index("--use_msa")
        assert call_args[msa_idx + 1] == "false", "Mini tier should skip MSA"

    def test_builds_base_command(self, tmpdir):
        """Verify base model CLI args."""
        json_path = tmpdir / "test_ON.json"
        json_path.write_text(json.dumps([{"name": "test_ON", "sequences": []}]))
        out_dir = tmpdir / "out"
        out_dir.mkdir()
        pdb_path = str(out_dir / "test_ON" / "model.pdb")
        summary_path = str(out_dir / "test_ON" / "model_summary.json")

        with patch("utils.protenix_eval.subprocess.run") as mock_run, \
             patch("utils.protenix_eval._find_cached_outputs", return_value=(None, None)), \
             patch("utils.protenix_eval._find_structure_and_summary", return_value=(pdb_path, summary_path)), \
             patch("utils.protenix_eval._run_msa_step", return_value=str(json_path)):
            mock_run.return_value = MagicMock(returncode=0)
            struct_path, sum_path = run_protenix_inference(
                str(json_path), str(out_dir), model_tier="base"
            )
        assert "model.pdb" in struct_path or "model.cif" in struct_path
        assert "model_summary.json" in sum_path

        call_args = mock_run.call_args[0][0]
        assert "protenix_base_default_v1.0.0" in call_args
        msa_idx = call_args.index("--use_msa")
        assert call_args[msa_idx + 1] == "true", "Base tier should enable MSA"

    def test_cattle_prod_fallback_when_non_strict(self, tmpdir):
        """If cattle-prod fails and strict mode is off, fallback to protenix."""
        json_path = tmpdir / "test_OFF.json"
        json_path.write_text(json.dumps([{"name": "test_OFF", "sequences": []}]))
        out_dir = tmpdir / "out"
        out_dir.mkdir()
        pdb_path = str(out_dir / "test_OFF" / "model.pdb")
        summary_path = str(out_dir / "test_OFF" / "model_summary.json")

        fail = MagicMock(returncode=1, stdout="", stderr="TensorNotFound(\"relpos.linear.weight\")")
        ok = MagicMock(returncode=0, stdout="", stderr="")
        with patch("utils.protenix_eval.subprocess.run", side_effect=[fail, ok]) as mock_run, \
             patch("utils.protenix_eval._find_cached_outputs", return_value=(None, None)), \
             patch("utils.protenix_eval._find_structure_and_summary", return_value=(pdb_path, summary_path)), \
             patch("utils.protenix_eval._run_msa_step", return_value=str(json_path)), \
             patch("utils.protenix_eval._resolve_cattle_prod_checkpoint", return_value=str(tmpdir / "ckpt")), \
             patch("utils.protenix_eval.shutil.which", return_value="protenix"), \
             patch("utils.protenix_eval.EVAL_ENGINE", "cattle-prod"), \
             patch("utils.protenix_eval.EVAL_BIN", "cattle-prod"), \
             patch("utils.protenix_eval._CATTLE_PROD_STRICT", False):
            run_protenix_inference(str(json_path), str(out_dir), model_tier="mini")

        assert mock_run.call_count == 2
        first_cmd = mock_run.call_args_list[0][0][0]
        second_cmd = mock_run.call_args_list[1][0][0]
        assert first_cmd[0] == "cattle-prod"
        assert second_cmd[0] == "protenix"

    def test_cattle_prod_strict_mode_raises(self, tmpdir):
        """If strict mode is on, cattle-prod failure is raised."""
        json_path = tmpdir / "test_OFF.json"
        json_path.write_text(json.dumps([{"name": "test_OFF", "sequences": []}]))
        out_dir = tmpdir / "out"
        out_dir.mkdir()

        fail = MagicMock(returncode=1, stdout="", stderr="TensorNotFound(\"relpos.linear.weight\")")
        with patch("utils.protenix_eval.subprocess.run", return_value=fail), \
             patch("utils.protenix_eval._find_cached_outputs", return_value=(None, None)), \
             patch("utils.protenix_eval._run_msa_step", return_value=str(json_path)), \
             patch("utils.protenix_eval._resolve_cattle_prod_checkpoint", return_value=str(tmpdir / "ckpt")), \
             patch("utils.protenix_eval.EVAL_ENGINE", "cattle-prod"), \
             patch("utils.protenix_eval.EVAL_BIN", "cattle-prod"), \
             patch("utils.protenix_eval._CATTLE_PROD_STRICT", True):
            with pytest.raises(Exception):
                run_protenix_inference(str(json_path), str(out_dir), model_tier="mini")
