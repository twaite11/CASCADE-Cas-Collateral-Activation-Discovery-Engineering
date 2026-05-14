"""Runs API smoke tests with the Vast.ai controller faked out.

Only exercises the router wiring + request/response shapes; full
controller behaviour is tested elsewhere (and by the CLI-mocked
provisioner tests).
"""
import sys
from pathlib import Path

from fastapi.testclient import TestClient

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _build_app(tmp_root: Path):
    """Build a minimal FastAPI app wired to a temp-rooted RunsStore + fakes."""
    import os

    os.environ["CASCADE_ROOT"] = str(tmp_root)
    # Phase-2 (C-2) auth: by default the mutating endpoints require an API key.
    # These tests don't exercise the auth path; opt out explicitly so they
    # continue to test only the underlying behaviour.
    os.environ.setdefault("CASCADE_API_KEY", "disabled")
    # Ensure needed dirs exist for load_config()
    (tmp_root / "metadata").mkdir(parents=True, exist_ok=True)
    (tmp_root / "outputs" / "run_logs").mkdir(parents=True, exist_ok=True)
    (tmp_root / "outputs" / "runs").mkdir(parents=True, exist_ok=True)
    (tmp_root / "outputs" / "phase1_screening").mkdir(parents=True, exist_ok=True)
    (tmp_root / "outputs" / "validated_baseline_ids.txt").write_text(
        "", encoding="utf-8"
    )
    (tmp_root / "metadata" / "variant_domain_metadata.json").write_text(
        "{}", encoding="utf-8"
    )
    (tmp_root / "jsons").mkdir(parents=True, exist_ok=True)

    # Re-import after CASCADE_ROOT is set so load_config() picks it up.
    for mod in list(sys.modules):
        if mod.startswith("dashboard_backend"):
            sys.modules.pop(mod)

    from fastapi import FastAPI

    from dashboard_backend.api import runs as runs_api
    from dashboard_backend.vast.runs_store import Run, RunStatus

    app = FastAPI()
    app.include_router(runs_api.router)

    # Replace controller + provisioner with fakes.
    class _FakeProv:
        async def search_offers(self, **_):
            from dashboard_backend.vast.provisioner import VastOffer

            return [
                VastOffer(
                    id=111,
                    gpu_name="A100_SXM4",
                    num_gpus=1,
                    gpu_ram_gb=80.0,
                    cpu_cores=16,
                    cpu_ram_gb=128.0,
                    disk_space_gb=200.0,
                    dph_total=1.20,
                    inet_down_mbps=500.0,
                    inet_up_mbps=300.0,
                    datacenter="us-east",
                    reliability=0.99,
                    raw={},
                )
            ]

    class _FakeCtrl:
        def __init__(self, store):
            self.store = store

        async def launch(self, **kwargs):
            run = Run(
                id="test1234",
                label=kwargs.get("label") or "run-test1234",
                status=RunStatus.QUEUED,
                baseline_ids=kwargs["baseline_ids"],
                crrna_lookup_ids=kwargs.get("crrna_lookup_ids") or kwargs["baseline_ids"],
                max_generations=kwargs["max_generations"],
                variants_per_gen=kwargs["variants_per_gen"],
                workers=kwargs["workers"],
                offer_id=kwargs["offer_id"],
            )
            self.store.create(run)
            return run

        async def cancel(self, run_id):
            self.store.update(run_id, status=RunStatus.CANCELLED)

    def _fake_prov():
        return _FakeProv()

    def _fake_ctrl():
        return _FakeCtrl(runs_api.get_store(runs_api._config()))

    app.dependency_overrides[runs_api.get_provisioner] = _fake_prov
    app.dependency_overrides[runs_api.get_controller] = _fake_ctrl
    return app


def test_list_offers(tmpdir):
    app = _build_app(Path(tmpdir))
    with TestClient(app) as client:
        resp = client.get("/api/vast/offers?gpu_name=A100_SXM4")
    assert resp.status_code == 200
    data = resp.json()
    assert "query" in data
    assert data["offers"][0]["id"] == 111
    assert data["offers"][0]["dph_total"] == 1.20


def test_launch_and_list_and_get(tmpdir):
    app = _build_app(Path(tmpdir))
    with TestClient(app) as client:
        resp = client.post(
            "/api/runs",
            json={
                "baseline_ids": ["B1", "B2"],
                "offer_id": 111,
                "max_generations": 4,
                "variants_per_gen": 3,
                "workers": 2,
                "label": "alpha",
            },
        )
        assert resp.status_code == 202, resp.text
        run = resp.json()
        assert run["id"] == "test1234"
        assert run["status"] == "queued"
        assert run["baseline_ids"] == ["B1", "B2"]

        lst = client.get("/api/runs").json()
        assert lst["total"] == 1
        assert lst["rows"][0]["id"] == "test1234"

        detail = client.get("/api/runs/test1234").json()
        assert detail["run"]["id"] == "test1234"
        assert detail["log_tail"] == []

        cancel = client.delete("/api/runs/test1234")
        assert cancel.status_code == 202
        after = client.get("/api/runs/test1234").json()
        assert after["run"]["status"] == "cancelled"


def test_launch_validates_crrna_length(tmpdir):
    app = _build_app(Path(tmpdir))
    with TestClient(app) as client:
        resp = client.post(
            "/api/runs",
            json={
                "baseline_ids": ["B1", "B2"],
                "crrna_lookup_ids": ["only_one"],
                "offer_id": 111,
            },
        )
    assert resp.status_code == 400
    assert "crrna_lookup_ids" in resp.json()["detail"]
