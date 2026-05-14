"""Concurrency regression tests for Option B fan-out (one VPS per baseline).

Three runs finishing close together hammer the controller's shared state
+ filesystem.  These tests cover the surfaces that *could* race:

  * metadata/runs.db -- WAL + busy_timeout = 5s, three rows with distinct ids
  * outputs/rl_gym_data/rl_training_dataset.jsonl -- shared append target,
    now protected by promoter._RL_APPEND_LOCK + O_APPEND + fsync
  * outputs/{optimized_switches,fast_eval,high_fidelity_scoring} -- already
    namespaced by `{run_id}__` prefix, but we want a regression test
  * dashboard_backend.api.runs._state singleton -- check-then-create race
    now serialised by _state_lock
"""
from __future__ import annotations

import concurrent.futures
import json
import os
import sys
import threading
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# ===========================================================================
# runs.db: WAL means 3 concurrent writers with distinct ids don't clash
# ===========================================================================
class TestRunsDbConcurrentWrites:
    def test_three_workers_can_create_runs_in_parallel(self, tmp_path):
        from dashboard_backend.vast.runs_store import Run, RunStatus, RunsStore

        store = RunsStore(tmp_path / "runs.db")

        def make_run(i: int) -> Run:
            return Run(
                id=RunsStore.new_id(),
                label=f"run-{i}",
                status=RunStatus.QUEUED,
                baseline_ids=[f"baseline_{i}"],
                crrna_lookup_ids=[f"baseline_{i}"],
                max_generations=12,
                variants_per_gen=5,
                workers=1,
                offer_id=None,
            )

        runs = [make_run(i) for i in range(3)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
            list(ex.map(store.create, runs))

        assert store.count() == 3
        all_ids = {r.id for r in store.list()}
        assert all_ids == {r.id for r in runs}

    def test_concurrent_status_updates_on_distinct_runs(self, tmp_path):
        from dashboard_backend.vast.runs_store import Run, RunStatus, RunsStore

        store = RunsStore(tmp_path / "runs.db")
        run_ids = []
        for i in range(3):
            r = Run(
                id=RunsStore.new_id(),
                label=f"run-{i}",
                status=RunStatus.QUEUED,
                baseline_ids=[f"baseline_{i}"],
                crrna_lookup_ids=[f"baseline_{i}"],
                max_generations=12,
                variants_per_gen=5,
                workers=1,
                offer_id=None,
            )
            store.create(r)
            run_ids.append(r.id)

        def update(rid: str) -> None:
            store.update(rid, status=RunStatus.RUNNING)

        # Run 50 concurrent updates split across 3 run_ids.
        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            list(ex.map(update, [run_ids[i % 3] for i in range(50)]))

        for rid in run_ids:
            assert store.get(rid).status == RunStatus.RUNNING


# ===========================================================================
# promoter: RL JSONL append is locked + atomic
# ===========================================================================
class TestPromoterRlAppendIsThreadSafe:
    def test_promoter_uses_threading_lock_around_append(self):
        """Source-level audit: the lock and O_APPEND + fsync triple is in place."""
        src = (ROOT / "dashboard_backend" / "vast" / "promoter.py").read_text(
            encoding="utf-8"
        )
        assert "_RL_APPEND_LOCK = threading.Lock()" in src
        assert "with _RL_APPEND_LOCK:" in src
        assert "os.O_APPEND" in src
        assert "os.fsync(fd)" in src
        # Must NOT use the old plain `open("a")` pattern for the RL append.
        assert 'rl_dest.open(\n            "a"' not in src
        assert "rl_dest.open(\"a\"" not in src

    def test_three_concurrent_promotes_produce_no_torn_records(self, tmp_path):
        """End-to-end: spin up 3 promoter calls in parallel and confirm the
        merged dataset is valid JSONL with the expected number of records.
        """
        from dashboard_backend.vast.promoter import promote_run

        cascade_root = tmp_path / "cascade"
        artifacts_dir = tmp_path / "artifacts"
        cascade_root.mkdir()
        artifacts_dir.mkdir()

        # Each run produces 100 records with a unique run-tag, so we can
        # confirm none were dropped or duplicated.
        run_ids = ["run-A", "run-B", "run-C"]
        records_per_run = 100
        for rid in run_ids:
            rl_dir = artifacts_dir / rid / "rl_gym_data"
            rl_dir.mkdir(parents=True)
            with (rl_dir / "rl_training_dataset.jsonl").open("w") as f:
                for i in range(records_per_run):
                    f.write(json.dumps({"variant_id": f"{rid}-{i}", "fitness": 0.5}) + "\n")

        # Fire all three promoters concurrently.
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
            list(ex.map(
                lambda rid: promote_run(rid, artifacts_dir, cascade_root),
                run_ids,
            ))

        merged = cascade_root / "outputs" / "rl_gym_data" / "rl_training_dataset.jsonl"
        assert merged.exists()
        seen = []
        with merged.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Every line must parse cleanly -- no torn records.
                rec = json.loads(line)
                seen.append(rec["variant_id"])

        assert len(seen) == 3 * records_per_run, (
            f"Expected {3 * records_per_run} records, got {len(seen)} -- "
            "concurrent promoters dropped/duplicated rows"
        )
        # Every run's records should appear exactly `records_per_run` times.
        for rid in run_ids:
            count = sum(1 for v in seen if v.startswith(rid + "-"))
            assert count == records_per_run, f"{rid}: expected {records_per_run}, got {count}"


# ===========================================================================
# Artifact dirs: every copied file is `<run_id>__`-prefixed
# ===========================================================================
class TestArtifactDirsAreNamespaced:
    def test_three_runs_writing_same_filename_do_not_clobber(self, tmp_path):
        from dashboard_backend.vast.promoter import promote_run

        cascade_root = tmp_path / "cascade"
        artifacts_dir = tmp_path / "artifacts"
        cascade_root.mkdir()
        artifacts_dir.mkdir()

        # All three runs produce an `optimized_switches/foo.fasta` -- the
        # promoter MUST prefix each with `<run_id>__` so they coexist.
        for rid in ["alpha", "bravo", "charlie"]:
            d = artifacts_dir / rid / "optimized_switches"
            d.mkdir(parents=True)
            (d / "foo.fasta").write_text(f">{rid}\nMKL\n")

            promote_run(rid, artifacts_dir, cascade_root)

        files = sorted(
            p.name for p in (cascade_root / "outputs" / "optimized_switches").iterdir()
        )
        assert files == ["alpha__foo.fasta", "bravo__foo.fasta", "charlie__foo.fasta"]


# ===========================================================================
# _state singleton: concurrent first-init does not produce duplicates
# ===========================================================================
class TestApiRunsStateSingletonLock:
    def test_state_lock_exists(self):
        src = (ROOT / "dashboard_backend" / "api" / "runs.py").read_text(
            encoding="utf-8"
        )
        assert "_state_lock = threading.Lock()" in src
        # Each lazy-init helper must hold the lock.
        for fn in ("def get_store", "def get_provisioner", "def get_controller"):
            idx = src.index(fn)
            tail = src[idx : idx + 600]
            assert "_state_lock" in tail, f"{fn}: missing lock guard"

    def test_concurrent_get_store_returns_same_instance(self, tmp_path, monkeypatch):
        """Ten concurrent get_store() calls must all return identity-equal
        RunsStore instances -- otherwise we'd have three sqlite connections
        all pointing at the same DB with three different WAL views.
        """
        from dashboard_backend.api import runs as runs_api
        from dashboard_backend.config import load_config

        monkeypatch.setenv("CASCADE_ROOT", str(tmp_path))
        (tmp_path / "metadata").mkdir()
        # Clear the singleton cache from a previous test run.
        runs_api._state.clear()

        cfg = load_config()
        results: list[object] = []
        lock = threading.Lock()

        def grab() -> None:
            s = runs_api.get_store(cfg)
            with lock:
                results.append(s)

        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as ex:
            list(ex.map(lambda _: grab(), range(10)))

        # All ten threads must see the same singleton.
        assert len(results) == 10
        assert all(r is results[0] for r in results)
