"""Regression tests for Phase 2 (security) fixes.

Covers:
  * C-1  CORS doesn't combine allow_origins=* with credentials
  * C-2  API key required on mutating endpoints + WS
  * C-3  SSH host-key strictness (StrictHostKeyChecking=accept-new)
  * C-4  Runs API total reflects the real DB count, not the page length
  * C-6  rsync timeout bound is configured
  * C-9  Dockerfile.controller runs as non-root and doesn't blanket-trust XFF
  * C-18 /api/structure-file containment under realpath / commonpath
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
@pytest.fixture
def fresh_app(monkeypatch):
    """Import dashboard_backend.main with a clean module cache so env vars
    set via monkeypatch take effect on module import."""
    def _factory(env: dict[str, str] | None = None):
        if env:
            for k, v in env.items():
                monkeypatch.setenv(k, v)
        # Drop cached modules so module-level code re-runs.
        for mod in [
            "dashboard_backend.main",
            "dashboard_backend.api.runs",
            "dashboard_backend.auth",
            "dashboard_backend",
        ]:
            sys.modules.pop(mod, None)
        import dashboard_backend.main as main_mod  # noqa: PLC0415
        return main_mod
    return _factory


# ===========================================================================
# C-1 CORS
# ===========================================================================
class TestC1_CORS:
    def test_default_origins_explicit_not_wildcard(self, fresh_app, monkeypatch):
        monkeypatch.delenv("CASCADE_ALLOWED_ORIGINS", raising=False)
        main_mod = fresh_app()
        client = TestClient(main_mod.app)
        # An allowed origin (vite default) must get CORS headers
        r = client.options(
            "/api/overview",
            headers={
                "Origin": "http://localhost:5173",
                "Access-Control-Request-Method": "GET",
            },
        )
        assert r.headers.get("access-control-allow-origin") == "http://localhost:5173"
        # A random origin must NOT be echoed back
        r2 = client.options(
            "/api/overview",
            headers={
                "Origin": "https://attacker.example",
                "Access-Control-Request-Method": "GET",
            },
        )
        assert r2.headers.get("access-control-allow-origin") != "*"
        assert r2.headers.get("access-control-allow-origin") != "https://attacker.example"

    def test_wildcard_disables_credentials(self, fresh_app):
        main_mod = fresh_app({"CASCADE_ALLOWED_ORIGINS": "*"})
        client = TestClient(main_mod.app)
        r = client.options(
            "/api/overview",
            headers={
                "Origin": "https://anywhere.example",
                "Access-Control-Request-Method": "GET",
            },
        )
        # When wildcard is explicit, credentials must not be allowed.
        assert r.headers.get("access-control-allow-credentials") != "true"


# ===========================================================================
# C-2 API key gating
# ===========================================================================
class TestC2_AuthMutating:
    def test_unset_key_blocks_mutating_routes(self, fresh_app, monkeypatch):
        monkeypatch.delenv("CASCADE_API_KEY", raising=False)
        main_mod = fresh_app()
        client = TestClient(main_mod.app)
        # POST /api/runs must 401 with no body validation
        r = client.post(
            "/api/runs",
            json={
                "baseline_ids": ["b1"],
                "offer_id": 1,
                "max_generations": 1,
                "variants_per_gen": 1,
            },
        )
        assert r.status_code == 401
        assert "CASCADE_API_KEY" in r.json()["detail"]
        # DELETE /api/runs/{id}
        r2 = client.delete("/api/runs/nonexistent")
        assert r2.status_code == 401
        # GET /api/vast/offers
        r3 = client.get("/api/vast/offers")
        assert r3.status_code == 401

    def test_disabled_token_opens_routes(self, fresh_app):
        main_mod = fresh_app({"CASCADE_API_KEY": "disabled"})
        client = TestClient(main_mod.app)
        # POST hits the body validation now (no 401)
        r = client.post(
            "/api/runs",
            json={
                "baseline_ids": ["b1"],
                "offer_id": 1,
                "max_generations": 1,
                "variants_per_gen": 1,
            },
        )
        assert r.status_code != 401

    def test_wrong_key_rejected(self, fresh_app):
        main_mod = fresh_app({"CASCADE_API_KEY": "real-secret"})
        client = TestClient(main_mod.app)
        r = client.delete("/api/runs/x", headers={"X-Cascade-Api-Key": "WRONG"})
        assert r.status_code == 401
        # Header-based correct key
        r2 = client.delete(
            "/api/runs/x", headers={"X-Cascade-Api-Key": "real-secret"}
        )
        # 401 should NOT recur; it'll 404 for the missing run instead
        assert r2.status_code != 401
        # Query-string fallback
        r3 = client.delete("/api/runs/x?api_key=real-secret")
        assert r3.status_code != 401

    def test_read_only_routes_remain_open(self, fresh_app, monkeypatch):
        monkeypatch.delenv("CASCADE_API_KEY", raising=False)
        main_mod = fresh_app()
        client = TestClient(main_mod.app)
        # Even with no key, /api/overview and /api/variants stay reachable
        # for the React UI's anonymous fetch path.
        r = client.get("/api/overview")
        assert r.status_code == 200, r.text
        r2 = client.get("/api/variants")
        assert r2.status_code == 200


# ===========================================================================
# C-3 SSH host-key strictness
# ===========================================================================
class TestC3_SshHostKeys:
    def test_ssh_args_use_accept_new_not_no(self):
        from dashboard_backend.vast.runner import SshEndpoint
        ep = SshEndpoint(host="x.y", port=22)
        args = ep.ssh_args()
        joined = " ".join(args)
        assert "StrictHostKeyChecking=accept-new" in joined
        assert "StrictHostKeyChecking=no" not in joined
        assert "UserKnownHostsFile=/dev/null" not in joined or os.environ.get(
            "CASCADE_SSH_KNOWN_HOSTS"
        ) == "/dev/null"

    def test_known_hosts_env_override(self, monkeypatch, tmp_path):
        custom = tmp_path / "khs"
        monkeypatch.setenv("CASCADE_SSH_KNOWN_HOSTS", str(custom))
        # Re-import to pick up env on _known_hosts_path
        sys.modules.pop("dashboard_backend.vast.runner", None)
        from dashboard_backend.vast.runner import SshEndpoint
        ep = SshEndpoint(host="x.y", port=22)
        args = " ".join(ep.ssh_args())
        assert f"UserKnownHostsFile={custom}" in args
        assert custom.exists(), "known_hosts file should be created"


# ===========================================================================
# C-4 RunsStore.count + paginated total
# ===========================================================================
class TestC4_PaginationTotal:
    def test_runs_store_count_returns_db_total(self, tmp_path):
        from dashboard_backend.vast.runs_store import (
            Run, RunsStore, RunStatus,
        )
        import time
        db = tmp_path / "runs.db"
        store = RunsStore(db)
        # Create three runs across two statuses.
        for i in range(3):
            r = Run(
                id=f"run-{i}",
                label=f"r{i}",
                status=RunStatus.QUEUED if i < 2 else RunStatus.RUNNING,
                baseline_ids=[],
                crrna_lookup_ids=[],
                max_generations=1,
                variants_per_gen=1,
                workers=1,
                offer_id=None,
                instance_id=None,
                ssh_host=None,
                ssh_port=None,
                ssh_user="root",
                docker_image="",
                dph_usd=None,
                created_at=time.time(),
            )
            store.create(r)
        assert store.count() == 3
        assert store.count(statuses=[RunStatus.QUEUED]) == 2
        assert store.count(statuses=[RunStatus.RUNNING]) == 1

    def test_runs_api_total_uses_real_count(self, fresh_app, monkeypatch, tmp_path):
        # Point CASCADE_ROOT at a clean dir so we get a fresh runs.db
        monkeypatch.setenv("CASCADE_ROOT", str(tmp_path))
        (tmp_path / "metadata").mkdir()
        main_mod = fresh_app({"CASCADE_API_KEY": "disabled"})
        # The runs.py module caches singletons in module-level `_state`; clear
        # it so our fresh CASCADE_ROOT produces a fresh RunsStore.
        from dashboard_backend.api import runs as runs_api
        runs_api._state.clear()
        client = TestClient(main_mod.app)
        r = client.get("/api/runs?limit=10")
        assert r.status_code == 200
        body = r.json()
        assert "total" in body and "rows" in body
        assert "offset" in body and "limit" in body
        # Empty DB -> 0 from store.count(), not from len(page).
        assert body["total"] == 0
        # Now seed one row and confirm total tracks it (proves we're not
        # just measuring len(page)).
        from dashboard_backend.vast.runs_store import Run, RunsStore, RunStatus
        import time
        store = RunsStore(tmp_path / "metadata" / "runs.db")
        store.create(Run(
            id="run-1", label="r1", status=RunStatus.QUEUED, baseline_ids=[],
            crrna_lookup_ids=[], max_generations=1, variants_per_gen=1,
            workers=1, offer_id=None, instance_id=None, ssh_host=None,
            ssh_port=None, ssh_user="root", docker_image="",
            dph_usd=None, created_at=time.time(),
        ))
        # Force the API to re-open the store too
        runs_api._state.clear()
        r2 = client.get("/api/runs?limit=10")
        assert r2.json()["total"] == 1


# ===========================================================================
# C-6 rsync timeout
# ===========================================================================
class TestC6_RsyncTimeout:
    def test_default_timeout_constant_present(self):
        from dashboard_backend.vast import runner
        assert hasattr(runner, "DEFAULT_RSYNC_TIMEOUT_S")
        assert runner.DEFAULT_RSYNC_TIMEOUT_S > 0
        src = (PROJECT_ROOT / "dashboard_backend" / "vast" / "runner.py").read_text()
        assert "asyncio.wait_for(" in src
        assert "DEFAULT_RSYNC_TIMEOUT_S" in src


# ===========================================================================
# C-9 Dockerfile non-root + restricted XFF
# ===========================================================================
class TestC9_DockerfileHardening:
    def test_runs_as_non_root_and_xff_restricted(self):
        df = (PROJECT_ROOT / "Dockerfile.controller").read_text()
        assert "USER cascade" in df, "controller must drop privileges"
        assert "useradd" in df, "must create a non-root user"
        # No blanket trust for X-Forwarded-* headers
        assert '--forwarded-allow-ips "*"' not in df
        assert "CASCADE_FORWARDED_ALLOWED_IPS" in df


# ===========================================================================
# C-18 Symlink path containment
# ===========================================================================
class TestC18_StructureFileContainment:
    def _make_root(self, tmp_path):
        root = tmp_path / "cascade_root"
        (root / "outputs").mkdir(parents=True)
        cif = root / "outputs" / "model.cif"
        cif.write_text("data_dummy\n")
        return root

    def test_rejects_dotdot_traversal(self, fresh_app, monkeypatch, tmp_path):
        root = self._make_root(tmp_path)
        monkeypatch.setenv("CASCADE_ROOT", str(root))
        main_mod = fresh_app({"CASCADE_API_KEY": "disabled"})
        client = TestClient(main_mod.app)
        r = client.get("/api/structure-file?path=../../etc/passwd")
        # Either 400 (rejected) or 404 (path doesn't exist) is acceptable;
        # the security property is that we never serve content outside root.
        assert r.status_code in (400, 404)
        if r.status_code == 200:
            pytest.fail("must not return data for ../.. traversal")

    @pytest.mark.skipif(os.name == "nt", reason="symlink creation needs admin on Windows")
    def test_rejects_symlink_escape(self, fresh_app, monkeypatch, tmp_path):
        root = self._make_root(tmp_path)
        outside = tmp_path / "outside.cif"
        outside.write_text("data_secret\n")
        sym = root / "outputs" / "evil.cif"
        os.symlink(outside, sym)
        monkeypatch.setenv("CASCADE_ROOT", str(root))
        main_mod = fresh_app({"CASCADE_API_KEY": "disabled"})
        client = TestClient(main_mod.app)
        r = client.get("/api/structure-file?path=outputs/evil.cif")
        assert r.status_code == 400, (
            "symlink that resolves outside cascade_root must be rejected"
        )

    def test_allows_legitimate_file(self, fresh_app, monkeypatch, tmp_path):
        root = self._make_root(tmp_path)
        monkeypatch.setenv("CASCADE_ROOT", str(root))
        main_mod = fresh_app({"CASCADE_API_KEY": "disabled"})
        client = TestClient(main_mod.app)
        r = client.get("/api/structure-file?path=outputs/model.cif")
        assert r.status_code == 200, r.text
        assert "data_dummy" in r.text
