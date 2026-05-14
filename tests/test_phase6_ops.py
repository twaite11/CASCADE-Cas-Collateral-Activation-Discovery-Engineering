"""Phase 6 operability regression tests.

Covers:
  * C-17 atomic JSONL writes + aggregator skips malformed lines
  * C-21..C-23 LogHub backpressure (slow-consumer sentinel) +
                channel GC for closed runs
"""
from __future__ import annotations

import asyncio
import importlib
import json
import os
import sys
import time
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ===========================================================================
# C-17 atomic JSONL writes
# ===========================================================================
class TestC17_AtomicJsonl:
    def test_save_uses_os_open_with_o_append(self):
        src = (SCRIPTS / "evolution_orchestrator.py").read_text(encoding="utf-8")
        # Make sure we no longer use the plain `with open(dest, "a")` pattern.
        assert 'with open(dest, "a", encoding="utf-8") as f:' not in src
        # Verify the new O_APPEND + fsync pattern is in place.
        assert "os.O_APPEND" in src
        assert "os.fsync(fd)" in src

    def test_aggregator_skips_malformed_lines(self, tmp_path, monkeypatch):
        """Drop in a worker dataset with one valid + one corrupt line and
        confirm the aggregator preserves the valid record and drops the
        malformed one with a warning."""
        import evolution_orchestrator as eo

        gym_dir = tmp_path / "gym"
        worker_dir = gym_dir / "worker_0"
        worker_dir.mkdir(parents=True)
        wf = worker_dir / "rl_training_dataset.jsonl"
        good = json.dumps({"variant_id": "v1", "fitness": 0.42})
        bad = "{this is not, json"
        wf.write_text(good + "\n" + bad + "\n", encoding="utf-8")

        merged_path = tmp_path / "rl_training_dataset.jsonl"
        monkeypatch.setattr(eo, "RL_TRAINING_DATASET", str(merged_path))
        monkeypatch.setattr(eo, "GYM_DIR", str(gym_dir))

        eo._aggregate_worker_results()

        merged = merged_path.read_text(encoding="utf-8").splitlines()
        # Only the valid record survives.
        assert len(merged) == 1
        decoded = json.loads(merged[0])
        assert decoded["variant_id"] == "v1"


# ===========================================================================
# C-21..C-23 LogHub backpressure + GC
# ===========================================================================
class TestLogHub:
    def test_slow_consumer_gets_sentinel(self):
        from dashboard_backend.vast.log_hub import LogHub, SLOW_CONSUMER_NOTICE

        async def scenario():
            hub = LogHub(buffer_size=4, queue_max=2)
            q = await hub.subscribe("run-x")
            # Fill the subscriber's queue to capacity.
            for line in ["a\n", "b\n"]:
                await hub.publish("run-x", line)
            # The 3rd publish should overflow the queue, eject the
            # subscriber, and leave the disconnect sentinel pending.
            await hub.publish("run-x", "c\n")
            # Drain the queue and verify the last item is the sentinel.
            drained: list[str] = []
            while not q.empty():
                drained.append(q.get_nowait())
            assert drained[-1] == SLOW_CONSUMER_NOTICE
            # The subscriber is no longer tracked, so further publishes
            # don't fan out to it.
            await hub.publish("run-x", "d\n")
            assert q.empty()

        asyncio.run(scenario())

    def test_purge_closed_channels(self):
        from dashboard_backend.vast.log_hub import LogHub

        async def scenario():
            hub = LogHub(buffer_size=4)
            await hub.publish("run-y", "hello\n")
            await hub.close("run-y")
            assert hub.channel_count() == 1
            # Not yet old enough -- nothing purged.
            n = await hub.purge_closed_channels(older_than_seconds=600.0)
            assert n == 0
            assert hub.channel_count() == 1
            # Force the closed_at backwards so the channel is "old".
            hub._channels["run-y"].closed_at -= 1000.0
            n = await hub.purge_closed_channels(older_than_seconds=600.0)
            assert n == 1
            assert hub.channel_count() == 0

        asyncio.run(scenario())

    def test_purge_preserves_channels_with_active_subscribers(self):
        """A subscriber still attached to a closed channel keeps it alive
        until they explicitly unsubscribe.  (Otherwise we'd race-cancel
        an in-flight stream replay.)"""
        from dashboard_backend.vast.log_hub import LogHub

        async def scenario():
            hub = LogHub(buffer_size=4)
            await hub.publish("run-z", "hi\n")
            q = await hub.subscribe("run-z")
            await hub.close("run-z")
            hub._channels["run-z"].closed_at -= 1000.0
            # Subscriber present -> no purge.
            assert (await hub.purge_closed_channels(older_than_seconds=1.0)) == 0
            await hub.unsubscribe("run-z", q)
            assert (await hub.purge_closed_channels(older_than_seconds=1.0)) == 1

        asyncio.run(scenario())

    def test_controller_starts_gc_task(self):
        """The controller's _ensure_gc_running stub should spin up a
        background asyncio task on first launch().  We don't exercise a
        full launch (too many fakes needed) -- just confirm the hook is
        wired by inspecting source.
        """
        src = (PROJECT_ROOT / "dashboard_backend" / "vast" / "controller.py").read_text(
            encoding="utf-8"
        )
        assert "_ensure_gc_running" in src
        assert "purge_closed_channels" in src
        assert "asyncio.create_task(self._gc_loop()" in src
