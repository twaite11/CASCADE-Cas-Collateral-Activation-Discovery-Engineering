"""Phase 4 determinism + pipeline integrity regression tests.

Covers:
  * B-12 seed_all / resolve_seed / --seed flag and CASCADE_SEED env
  * B-21 single CRISPR-array definition (utils.crispr_constants)
  * C-5  switch_report append lock signature
  * C-13 deterministic glob (sorted everywhere it matters)
  * C-14 MSA flag respect (env override + DB presence)
  * C-15 sys.executable (no unqualified "python")
  * C-27 GPU lock present around every Protenix call
"""
from __future__ import annotations

import importlib
import os
import random
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ===========================================================================
# B-12 seed_all
# ===========================================================================
class TestB12_SeedAll:
    def test_seed_all_is_idempotent_and_seeds_random(self):
        from utils.determinism import seed_all
        seed_all(123)
        a = [random.random() for _ in range(5)]
        seed_all(123)
        b = [random.random() for _ in range(5)]
        assert a == b, "random must be deterministic after seed_all"

    def test_seed_all_seeds_numpy(self):
        np = pytest.importorskip("numpy")
        from utils.determinism import seed_all
        seed_all(7)
        a = np.random.random(8).tolist()
        seed_all(7)
        b = np.random.random(8).tolist()
        assert a == b

    def test_resolve_seed_precedence(self, monkeypatch):
        from utils.determinism import resolve_seed, DEFAULT_SEED
        monkeypatch.delenv("CASCADE_SEED", raising=False)
        # 1) CLI wins
        assert resolve_seed(99) == 99
        # 2) env beats default
        monkeypatch.setenv("CASCADE_SEED", "777")
        assert resolve_seed(None) == 777
        # 3) default when nothing set
        monkeypatch.delenv("CASCADE_SEED", raising=False)
        assert resolve_seed(None) == DEFAULT_SEED
        # 4) bad env -> default + log
        monkeypatch.setenv("CASCADE_SEED", "not-an-int")
        assert resolve_seed(None) == DEFAULT_SEED

    def test_orchestrator_exposes_seed_cli_flag(self):
        src = (SCRIPTS / "evolution_orchestrator.py").read_text(encoding="utf-8")
        assert '"--seed"' in src
        assert "resolve_seed" in src
        assert "seed_all(" in src
        # mismatch generator must use resolved seed instead of hardcoded 42
        assert "seed=_ms_seed" in src


# ===========================================================================
# B-21 single CRISPR-array definition
# ===========================================================================
class TestB21_CrisprConstants:
    def test_canonical_values(self):
        from utils.crispr_constants import (
            CRISPR_REPEAT_MIN, CRISPR_REPEAT_MAX,
            CRISPR_SPACER_MIN, CRISPR_SPACER_MAX,
            MIN_ARRAY_UNITS,
        )
        assert CRISPR_REPEAT_MIN == 23
        assert CRISPR_REPEAT_MAX == 50
        assert CRISPR_SPACER_MIN == 15
        assert CRISPR_SPACER_MAX == 80
        assert MIN_ARRAY_UNITS == 3

    def test_fix_crrna_assignments_imports_shared(self):
        # Re-export check: fix_crrna_assignments.CRISPR_REPEAT_MIN must
        # be the same int as the canonical module's value.
        from utils import crispr_constants as c
        fix_mod = importlib.import_module("fix_crrna_assignments")
        assert fix_mod.CRISPR_REPEAT_MIN is c.CRISPR_REPEAT_MIN
        assert fix_mod.MIN_ARRAY_UNITS is c.MIN_ARRAY_UNITS

    def test_discover_crrna_uses_shared(self):
        src = (SCRIPTS / "discover_crrna.py").read_text(encoding="utf-8")
        # The legacy "len(...) >= 20" gate must be gone; 23-50 nt window used.
        assert "len(repeat_candidate) >= 20" not in src
        assert "CRISPR_REPEAT_MIN <= len(repeat_candidate)" in src
        assert "CRISPR_REPEAT_MIN <= len(repeat_seq)" in src


# ===========================================================================
# C-5 switch_report lock plumbing
# ===========================================================================
class TestC5_SwitchReportLock:
    def test_append_signature_includes_lock(self):
        src = (SCRIPTS / "evolution_orchestrator.py").read_text(encoding="utf-8")
        assert "report_lock=None" in src, "append_switch_report should accept a report_lock kwarg"
        assert "with lock_cm:" in src, "writer should serialise via the lock"

    def test_multiworker_passes_a_real_lock(self):
        src = (SCRIPTS / "evolution_orchestrator.py").read_text(encoding="utf-8")
        # manager.Lock() should be created for multi-worker runs, and the
        # report_lock=... kwarg should flow through to _run_single_lineage.
        assert "report_lock = manager.Lock()" in src
        assert "chunks[i], gpu_lock, report_lock" in src


# ===========================================================================
# C-13 deterministic glob
# ===========================================================================
class TestC13_DeterministicGlob:
    @pytest.mark.parametrize(
        "filename",
        [
            "evolution_orchestrator.py",
            "03_pxdesign_wrapper.py",
            "utils/protenix_eval.py",
        ],
    )
    def test_no_bare_glob(self, filename):
        src = (SCRIPTS / filename).read_text(encoding="utf-8")
        # Strip line comments so we don't false-match on doc strings that
        # talk about "glob.glob" in prose.
        code_lines = [
            line for line in src.splitlines()
            if not line.lstrip().startswith("#")
        ]
        for line in code_lines:
            if "= glob.glob(" in line or "= _glob.glob(" in line:
                pytest.fail(
                    f"{filename}: bare glob() result assigned without sorted(): "
                    f"{line.strip()!r}"
                )


# ===========================================================================
# C-14 MSA flag respect
# ===========================================================================
class TestC14_MsaFlag:
    def test_protenix_eval_respects_env_and_db_presence(self):
        src = (SCRIPTS / "utils" / "protenix_eval.py").read_text(encoding="utf-8")
        assert "CASCADE_DISABLE_MSA" in src
        assert "os.path.isdir(seqres_db_path)" in src
        # The old hardcoded "use_msa = model_tier != \"mini\"" must be
        # gone -- replaced with the tier + env + db gate.
        assert (
            "use_msa = (model_tier != \"mini\") and (not disable_msa) and db_available"
            in src
        )


# ===========================================================================
# C-15 sys.executable
# ===========================================================================
class TestC15_SysExecutable:
    def test_pxdesign_wrapper_uses_sys_executable(self):
        src = (SCRIPTS / "03_pxdesign_wrapper.py").read_text(encoding="utf-8")
        # The unqualified "python" subprocess invocation must be gone.
        assert '["python", parse_script' not in src
        assert "sys.executable, parse_script" in src
        assert "import sys" in src


# ===========================================================================
# C-27 GPU lock around every Protenix call (smoke-level audit)
# ===========================================================================
class TestC27_GpuLock:
    def test_every_run_protenix_call_is_inside_gpu_lock(self):
        src = (SCRIPTS / "evolution_orchestrator.py").read_text(encoding="utf-8")
        lines = src.splitlines()
        for i, line in enumerate(lines):
            if "run_protenix_inference(" in line and "def " not in line:
                # Look back up to 6 lines for either ``with gpu_lock:``
                # (the typical pattern) or ``with lock_cm:`` (the variant
                # used in _resolve_best_baseline where the caller may
                # legitimately pass None for single-worker mode and we
                # fall back to a no-op _DummyLock).  Both are correct.
                window = "\n".join(lines[max(0, i - 6): i + 1])
                assert ("with gpu_lock:" in window) or ("with lock_cm:" in window), (
                    f"run_protenix_inference at line {i + 1} not protected "
                    f"by a `with gpu_lock:`/`with lock_cm:` block. Surrounding code:\n{window}"
                )
