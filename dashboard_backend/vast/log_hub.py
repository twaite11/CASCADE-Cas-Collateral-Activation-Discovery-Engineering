"""In-process fan-out hub for per-run log lines.

The runner publishes every line to the hub; every WebSocket subscriber gets
its own bounded queue. A ring buffer per run lets late-joiners replay the
last N lines immediately (xterm.js scrollback).

Zero external dependencies (no Redis). Single-process only — fine for the
controller since it runs as one asyncio process.
"""
from __future__ import annotations

import asyncio
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Deque


@dataclass
class _RunChannel:
    buffer: Deque[str]
    subscribers: set[asyncio.Queue[str]]
    closed: bool = False


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
                ch.subscribers.discard(q)

    async def close(self, run_id: str) -> None:
        async with self._lock:
            ch = self._channels.get(run_id)
            if ch is None:
                return
            ch.closed = True
            for q in list(ch.subscribers):
                try:
                    q.put_nowait("\n[stream closed]\n")
                except asyncio.QueueFull:
                    pass

    async def subscribe(self, run_id: str) -> asyncio.Queue[str]:
        q: asyncio.Queue[str] = asyncio.Queue(maxsize=self.queue_max)
        async with self._lock:
            ch = self._channels[run_id]
            ch.subscribers.add(q)
            # Pre-fill with the scrollback so a new client sees history.
            for line in ch.buffer:
                try:
                    q.put_nowait(line)
                except asyncio.QueueFull:
                    break
            if ch.closed:
                try:
                    q.put_nowait("\n[stream closed]\n")
                except asyncio.QueueFull:
                    pass
        return q

    async def unsubscribe(self, run_id: str, q: asyncio.Queue[str]) -> None:
        async with self._lock:
            ch = self._channels.get(run_id)
            if ch:
                ch.subscribers.discard(q)


hub = LogHub()
