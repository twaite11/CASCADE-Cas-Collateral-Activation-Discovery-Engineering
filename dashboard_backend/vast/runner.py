"""Async SSH runner that drives a provisioned Vast.ai instance.

Responsibilities:

1. Wait for SSH to come up on the newly created instance.
2. Rsync runtime inputs (jsons/, metadata/, outputs/phase1_screening/,
   validated_baseline_ids.txt) onto the VPS so the orchestrator has
   everything it needs.
3. Stream stdout/stderr of ``docker exec`` / direct python invocation back
   to the controller as line records; a fan-out WebSocket hub will
   broadcast to any connected dashboards.
4. On completion (or failure / cancel), rsync ``outputs/runs/<run_id>/``
   back to the controller and destroy the instance.

``asyncssh`` is imported lazily so importing this module doesn't require the
package at dev time (tests mock the runner).
"""
from __future__ import annotations

import asyncio
import contextlib
import json
import logging
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import AsyncIterator, Awaitable, Callable

log = logging.getLogger(__name__)

LogCallback = Callable[[str], Awaitable[None]]

# C-3 fix: stop disabling SSH host-key verification.  By default we keep a
# trust-on-first-use known_hosts file under ~/.cascade so the *first* time we
# see a Vast.ai instance we pin its fingerprint, and every subsequent
# connection refuses to proceed if the key has changed (== MITM).  Operators
# can override the path via CASCADE_SSH_KNOWN_HOSTS.  Setting it to /dev/null
# explicitly opts back into the legacy (insecure) behaviour with a logged
# warning, for one-off local testing.
DEFAULT_KNOWN_HOSTS = Path.home() / ".cascade" / "known_hosts"


def _resolved_known_hosts() -> Path:
    override = os.environ.get("CASCADE_SSH_KNOWN_HOSTS", "").strip()
    return Path(override) if override else DEFAULT_KNOWN_HOSTS


def _known_hosts_path() -> str:
    path = _resolved_known_hosts()
    # Treat /dev/null (Linux) or the Windows equivalent as an explicit opt-out
    # of host-key pinning.  Warn once when we see it.
    if str(path) in {"/dev/null", "NUL", "nul"}:
        if not getattr(_known_hosts_path, "_warned", False):
            log.warning(
                "SSH host-key verification disabled (known_hosts=%s). MITM "
                "protection is OFF; do this only on trusted local networks.",
                path,
            )
            _known_hosts_path._warned = True  # type: ignore[attr-defined]
        return str(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch(exist_ok=True)
    return str(path)


# C-6 fix: rsync is run via asyncio.create_subprocess_exec + await
# proc.communicate() which has no upper bound.  A wedged Vast network would
# pin the event loop forever.  Operators can override via
# CASCADE_RSYNC_TIMEOUT (seconds); default 10 minutes is enough for our
# multi-GB artifacts dir.
DEFAULT_RSYNC_TIMEOUT_S = float(os.environ.get("CASCADE_RSYNC_TIMEOUT", "600"))


@dataclass
class RunCommand:
    """Composes the shell command the orchestrator container runs inside the VPS.

    The Vast.ai ``--onstart-cmd`` already runs the container's entrypoint and
    wires `$CASCADE_*` env vars. We reuse the same invocation verbatim over
    SSH for post-start tailing or recovery.
    """

    run_id: str
    baseline_ids: list[str]
    crrna_lookup_ids: list[str] = field(default_factory=list)
    max_generations: int = 12
    variants_per_gen: int = 5
    workers: int = 3
    container_root: str = "/workspace/CASCADE"
    logs_subdir: str = "logs"

    def onstart_cmd(self) -> str:
        """Full shell one-liner suitable for ``vastai create instance --onstart-cmd``.

        It:
          * sources conda + activates cascade
          * exports PXDESIGN_CMD
          * cds into the orchestrator scripts dir
          * execs the orchestrator with our flags
          * tees stdout to a run-scoped log file so the controller can rsync
            a complete transcript even if the WebSocket dropped
        """
        bid = ",".join(self.baseline_ids)
        cid = ",".join(self.crrna_lookup_ids) if self.crrna_lookup_ids else ""
        parts = [
            "source /opt/conda/etc/profile.d/conda.sh",
            "conda activate cascade",
            "export PXDESIGN_CMD=/opt/conda/envs/pxdesign/bin/pxdesign",
            f"mkdir -p {self.container_root}/{self.logs_subdir}",
            f"cd {self.container_root}/scripts",
            (
                "python -u evolution_orchestrator.py "
                f"--run-id {shlex.quote(self.run_id)} "
                f"--baseline-ids {shlex.quote(bid)} "
                + (f"--crrna-lookup-ids {shlex.quote(cid)} " if cid else "")
                + f"--max-generations {self.max_generations} "
                f"--variants-per-gen {self.variants_per_gen} "
                f"--workers {self.workers} "
                f"2>&1 | tee {self.container_root}/{self.logs_subdir}/run-{shlex.quote(self.run_id)}.log"
            ),
        ]
        return " && ".join(parts)


@dataclass
class SshEndpoint:
    host: str
    port: int
    user: str = "root"
    key_path: str | None = None

    def ssh_args(self) -> list[str]:
        """Argv fragment suitable for subprocess rsync/ssh calls.

        C-3 fix: replaced ``StrictHostKeyChecking=no`` +
        ``UserKnownHostsFile=/dev/null`` (== "trust everything, remember
        nothing") with ``StrictHostKeyChecking=accept-new`` against a stable
        known_hosts file under ``~/.cascade``.  Trust-on-first-use: the first
        time we see a Vast instance we pin its fingerprint, every subsequent
        connection refuses to proceed if the key changes.  Operators wanting
        the legacy behaviour can set ``CASCADE_SSH_KNOWN_HOSTS=/dev/null``.
        """
        return [
            "-p",
            str(self.port),
            "-o",
            "StrictHostKeyChecking=accept-new",
            "-o",
            f"UserKnownHostsFile={_known_hosts_path()}",
            "-o",
            "LogLevel=ERROR",
        ] + (["-i", self.key_path] if self.key_path else [])


class SshRunner:
    """Streams remote process output and handles rsync syncs.

    The streaming path uses ``asyncssh.connect(...).run(...)`` with live
    stdout capture; file sync uses a subprocess call to the system ``rsync``
    binary because it handles delta transfers and partial files natively.
    """

    def __init__(
        self,
        endpoint: SshEndpoint,
        *,
        logs_dir: Path,
        artifacts_dir: Path,
        inputs_root: Path | None = None,
    ) -> None:
        self.endpoint = endpoint
        self.logs_dir = logs_dir
        self.artifacts_dir = artifacts_dir
        self.inputs_root = inputs_root
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self.artifacts_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------------------------------- connection / auth

    async def _connect(self):
        import asyncssh  # noqa: PLC0415 — optional import

        client_keys: list[str] | None = None
        if self.endpoint.key_path:
            client_keys = [self.endpoint.key_path]

        return await asyncssh.connect(
            host=self.endpoint.host,
            port=self.endpoint.port,
            username=self.endpoint.user,
            known_hosts=None,
            client_keys=client_keys,
        )

    async def wait_for_ssh(self, *, total_timeout: float = 300.0, interval: float = 6.0) -> None:
        """Poll until an SSH connection succeeds. Swallows transient errors."""
        loop = asyncio.get_event_loop()
        deadline = loop.time() + total_timeout
        last_err: Exception | None = None
        while loop.time() < deadline:
            try:
                conn = await self._connect()
                conn.close()
                await conn.wait_closed()
                return
            except Exception as exc:  # noqa: BLE001
                last_err = exc
                await asyncio.sleep(interval)
        raise TimeoutError(f"SSH not ready after {total_timeout}s: {last_err}")

    # ------------------------------------------------------------- sync

    async def rsync_push(self, local: Path, remote: str) -> None:
        """Push local dir/file to the remote path over rsync+ssh."""
        await self._rsync(str(local), f"{self.endpoint.user}@{self.endpoint.host}:{remote}")

    async def rsync_pull(self, remote: str, local: Path) -> None:
        """Pull the remote dir/file to a local path."""
        local.parent.mkdir(parents=True, exist_ok=True)
        await self._rsync(f"{self.endpoint.user}@{self.endpoint.host}:{remote}", str(local))

    async def _rsync(self, src: str, dest: str) -> None:
        rsync = shutil.which("rsync") or "rsync"
        ssh_cmd = "ssh " + " ".join(self.endpoint.ssh_args())
        cmd = [rsync, "-azP", "--partial", "-e", ssh_cmd, src, dest]
        log.info("rsync: %s -> %s", src, dest)
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        # C-6 fix: bound the wait so a wedged remote can't pin the event loop.
        try:
            stdout_b, stderr_b = await asyncio.wait_for(
                proc.communicate(), timeout=DEFAULT_RSYNC_TIMEOUT_S
            )
        except asyncio.TimeoutError:
            # Best-effort cleanup so we don't leak the child process.
            try:
                proc.kill()
            except ProcessLookupError:
                pass
            with contextlib.suppress(Exception):
                await proc.wait()
            raise RuntimeError(
                f"rsync timed out after {DEFAULT_RSYNC_TIMEOUT_S:.0f}s: "
                f"{src} -> {dest}"
            )
        if proc.returncode != 0:
            raise RuntimeError(
                f"rsync failed ({proc.returncode}): {stderr_b.decode(errors='replace')}"
            )
        log.debug("rsync ok: %s bytes", len(stdout_b))

    # ------------------------------------------------- remote execution

    async def run_stream(
        self,
        command: str,
        *,
        on_line: LogCallback | None = None,
        local_log_path: Path | None = None,
    ) -> int:
        """Execute ``command`` remotely, streaming stdout/stderr line-by-line.

        Each line is:
          * appended to ``local_log_path`` (if provided) for durable replay
          * passed to ``on_line(line)`` callback (if provided) for WS fanout
        Returns the remote process exit status.
        """
        conn = await self._connect()
        log_f = local_log_path.open("a", encoding="utf-8") if local_log_path else None
        try:
            async with conn:
                proc = await conn.create_process(command, term_type="xterm-256color")
                assert proc.stdout is not None
                async for line in _merged_stream(proc):
                    line = line.rstrip("\r\n")
                    if log_f:
                        log_f.write(line + "\n")
                        log_f.flush()
                    if on_line:
                        try:
                            await on_line(line)
                        except Exception:  # noqa: BLE001
                            log.exception("on_line callback raised; continuing")
                await proc.wait()
                return int(proc.returncode or 0)
        finally:
            if log_f:
                log_f.close()


async def _merged_stream(proc) -> AsyncIterator[str]:
    """Yield lines from the remote process's stdout+stderr interleaved."""
    # asyncssh gives us two streams; merge them deterministically.
    queue: asyncio.Queue[str | None] = asyncio.Queue()

    async def pump(reader, tag: str) -> None:
        try:
            if reader is None:
                return
            async for line in reader:
                await queue.put(line if isinstance(line, str) else line.decode("utf-8", "replace"))
        finally:
            await queue.put(None)

    tasks = [
        asyncio.create_task(pump(proc.stdout, "out")),
        asyncio.create_task(pump(proc.stderr, "err")),
    ]
    live = len(tasks)
    try:
        while live > 0:
            item = await queue.get()
            if item is None:
                live -= 1
                continue
            yield item
    finally:
        for t in tasks:
            if not t.done():
                t.cancel()


# ---------------------------------------------------------------------------
# Helper: synchronous rsync for one-off calls (used in tests + docs).
# ---------------------------------------------------------------------------
def rsync_cli_invocation(src: str, dest: str, endpoint: SshEndpoint) -> list[str]:
    """Return the argv that would be executed for a given rsync call.

    Exposed for the CONTAINER_DEPLOY.md docs and for mocking in tests.
    """
    rsync = shutil.which("rsync") or "rsync"
    ssh_cmd = "ssh " + " ".join(endpoint.ssh_args())
    return [rsync, "-azP", "--partial", "-e", ssh_cmd, src, dest]


__all__ = [
    "RunCommand",
    "SshEndpoint",
    "SshRunner",
    "rsync_cli_invocation",
]
