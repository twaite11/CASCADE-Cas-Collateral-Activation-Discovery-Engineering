import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dashboard_backend.vast import Run, RunsStore, RunStatus


def _mk_run(**overrides) -> Run:
    base = dict(
        id=RunsStore.new_id(),
        label="test",
        status=RunStatus.QUEUED,
        baseline_ids=["a", "b"],
        crrna_lookup_ids=["a", "b"],
        max_generations=3,
        variants_per_gen=2,
        workers=1,
        offer_id=12345,
    )
    base.update(overrides)
    return Run(**base)


def test_create_get_update_list_roundtrip(tmpdir):
    store = RunsStore(Path(tmpdir) / "runs.db")
    run = _mk_run(label="alpha")
    store.create(run)

    fetched = store.get(run.id)
    assert fetched is not None
    assert fetched.label == "alpha"
    assert fetched.status == RunStatus.QUEUED
    assert fetched.baseline_ids == ["a", "b"]

    store.update(
        run.id,
        status=RunStatus.RUNNING,
        instance_id=42,
        ssh_host="1.2.3.4",
        ssh_port=22022,
    )
    after = store.get(run.id)
    assert after.status == RunStatus.RUNNING
    assert after.instance_id == 42
    assert after.ssh_host == "1.2.3.4"
    assert after.ssh_port == 22022

    # Active-only filter excludes terminal states
    other = _mk_run(label="beta", status=RunStatus.COMPLETED)
    store.create(other)
    active = {r.id for r in store.active_runs()}
    assert run.id in active
    assert other.id not in active


def test_update_list_fields_serialize(tmpdir):
    store = RunsStore(Path(tmpdir) / "runs.db")
    run = _mk_run(baseline_ids=["x"])
    store.create(run)
    store.update(run.id, baseline_ids=["x", "y", "z"])
    got = store.get(run.id)
    assert got.baseline_ids == ["x", "y", "z"]


def test_is_terminal_flag():
    assert _mk_run(status=RunStatus.COMPLETED).is_terminal
    assert _mk_run(status=RunStatus.FAILED).is_terminal
    assert _mk_run(status=RunStatus.CANCELLED).is_terminal
    assert not _mk_run(status=RunStatus.RUNNING).is_terminal
    assert not _mk_run(status=RunStatus.PROVISIONING).is_terminal
