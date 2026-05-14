"""In-process fan-out hub for per-run log lines.

The runner publishes every line to the hub; every WebSocket subscriber gets
its own bounded queue. A ring buffer per run lets late-joiners replay the
last N lines immediately (xterm.js scrollback).

Zero external dependencies (no Redis). Single-process only — fine for the
controller since it runs as one asyncio process.

C-21..C-23 (Phase 6 ops):
  * Slow subscribers used to be silently discarded with no signal -- we now
    push a one-shot ``[disconnected: slow consumer]`` sentinel so the
    client sees *why* its stream went dead.
  * Channels for terminated runs accumulated indefinitely.  ``close()``
    now records a close timestamp; ``purge_closed_channels(older_than)``
    GCs them.  The controller calls this periodically (see
    ``controller.py``).
"""
from __future__ import annotations

import asyncio
import time
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Deque


SLOW_CONSUMER_NOTICE = "\n[disconnected: slow consumer]\n"
STREAM_CLOSED_NOTICE = "\n[stream closed]\n"


@dataclass
class _RunChannel:
    buffer: Deque[str]
    subscribers: set[asyncio.Queue[str]]
    closed: bool = False
    closed_at: float = 0.0
    dropped_subscribers: int = field(default=0)


class LogHub:
    def __init__(self, *, buffer_size: int = 2000, queue_max: int = 5000) -> None:
        self.buffer_size = buffer_size
        self.queue_max = queue_max
        self._channels: dict[str, _RunChannel] = defaultdict(
            lambda: _RunChannel(buffer=deque(maxlen=self.buffer_size), subscribers=set())
        )
        self._lock = asyncio.Lock()

    def snapshot(self, run_id: str) -> list[str]:
        ch = self._channels.get(run_id)
        return list(ch.buffer) if ch else []

    def channel_count(self) -> int:
        """Number of channels currently tracked (closed + open).  Useful
        for tests + the /api/health debug endpoint."""
        return len(self._channels)

    async def publish(self, run_id: str, line: str) -> None:
        async with self._lock:
            ch = self._channels[run_id]
            if ch.closed:
                return
            ch.buffer.append(line)
            dead: list[asyncio.Queue[str]] = []
            for q in ch.subscribers:
                try:
                    q.put_nowait(line)
                except asyncio.QueueFull:
                    dead.append(q)
            for q in dead:
                # C-22: tell the slow consumer why it's being dropped so the
                # UI can show a meaningful disconnect reason instead of
                # silently dying.  The queue is at capacity (that's why we
                # got here) so first evict its oldest item to make room
                # for the sentinel.
                try:
                    q.get_nowait()
                except asyncio.QueueEmpty:
                    pass
                try:
                    q.put_nowait(SLOW_CONSUMER_NOTICE)
                except asyncio.QueueFull:
                    pass
                ch.subscribers.discard(q)
                ch.dropped_subscribers += 1

    async def close(self, run_id: str) -> None:
        async with self._lock:
            ch = self._channels.get(run_id)
            if ch is None:
                return
            ch.closed = True
            ch.closed_at = time.monotonic()  # C-23: record for GC
            for q in list(ch.subscribers):
                try:
                    q.put_nowait(STREAM_CLOSED_NOTICE)
                except asyncio.QueueFull:
                    pass

    async def subscribe(self, run_id: str) -> asyncio.Queue[str]:
        q: asyncio.Queue[str] = asyncio.Queue(maxsize=self.queue_max)
        async with self._lock:
            ch = self._channels[run_id]
            ch.subscribers.add(q)
            for line in ch.buffer:
                try:
                    q.put_nowait(line)
                except asyncio.QueueFull:
                    break
            if ch.closed:
                try:
                    q.put_nowait(STREAM_CLOSED_NOTICE)
                except asyncio.QueueFull:
                    pass
        return q

    async def unsubscribe(self, run_id: str, q: asyncio.Queue[str]) -> None:
        async with self._lock:
            ch = self._channels.get(run_id)
            if ch:
                ch.subscribers.discard(q)

    async def purge_closed_channels(self, older_than_seconds: float = 600.0) -> int:
        """C-23: drop channels that have been ``closed`` for longer than
        ``older_than_seconds`` *and* have no remaining subscribers.

        Returns the number of channels purged.  Safe to call concurrently;
        holds the same internal lock as publish/subscribe.
        """
        now = time.monotonic()
        purged = 0
        async with self._lock:
            for rid in list(self._channels.keys()):
                ch = self._channels[rid]
                if not ch.closed:
                    continue
                if ch.subscribers:
                    continue
                if (now - ch.closed_at) < older_than_seconds:
                    continue
                del self._channels[rid]
                purged += 1
        return purged


hub = LogHub()
