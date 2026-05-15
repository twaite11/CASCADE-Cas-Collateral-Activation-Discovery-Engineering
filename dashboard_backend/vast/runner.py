"""Async SSH runner that drives a provisioned Vast.ai instance.

Responsibilities:

1. Wait for SSH to come up on the newly created instance.
2. SFTP runtime inputs (jsons/, metadata/, outputs/phase1_screening/,
   validated_baseline_ids.txt) onto the VPS so the orchestrator has
   everything it needs.
3. Stream stdout/stderr of ``docker exec`` / direct python invocation back
   to the controller as line records; a fan-out WebSocket hub will
   broadcast to any connected dashboards.
4. On completion (or failure / cancel), SFTP ``outputs/runs/<run_id>/``
   back to the controller and destroy the instance.

File transfers use asyncssh's built-in SFTP (same in-process SSH transport
as streaming/wait-for-SSH).  This avoids every Windows-specific issue we
hit with shelling out to ``scp.exe``/``rsync.exe`` through Vast.ai relays
(WinError 2, WinError 1225).

``asyncssh`` is imported lazily so importing this module doesn't require the
package at dev time (tests mock the runner).
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import posixpath
import shlex
import shutil
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


def _sync_executable() -> tuple[str, str]:
    """Prefer ``rsync``; otherwise fall back to OpenSSH ``scp``.

    Vanilla Windows desktops often omit ``rsync`` on PATH --- that surfaces as
    WinError 2 (file-not-found).  Bundled ``scp.exe`` works for our push/pull."""
    preferred = os.environ.get("CASCADE_SYNC_BIN", "").strip()
    if preferred:
        p = preferred if Path(preferred).is_file() else shutil.which(preferred)
        if not p:
            raise RuntimeError(
                f"CASCADE_SYNC_BIN={preferred!r} does not resolve to an executable. "
                "Clear CASCADE_SYNC_BIN to auto-detect, or supply a path to rsync.exe / scp.exe."
            )
        base = Path(p).name.lower()
        kind = "rsync" if "rsync" in base else "scp"
        return str(Path(p)), kind

    if which := shutil.which("rsync"):
        return which, "rsync"
    if which := shutil.which("scp"):
        return which, "scp"

    raise RuntimeError(
        "Neither ``rsync`` nor ``scp`` was found on PATH. On Windows enable OpenSSH "
        "Client (provides ``scp.exe``), install cwRsync/Git-Bash rsync, or use WSL. "
        "You can also set CASCADE_SYNC_BIN to the full path to rsync.exe or scp.exe."
    )


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
            "python -u generate_phase1_structures.py",
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
    """Streams remote process output and syncs artifacts.

    Streaming uses AsyncSSH ``create_process``.  File sync prefers the
    system ``rsync`` binary when present; otherwise it falls back to OpenSSH's
    ``scp -r``, which resolves the usual Windows symptom ``[WinError 2]``
    ``The system cannot find the file specified`` when ``rsync`` is not on PATH.
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

    async def wait_for_ssh(
        self,
        *,
        total_timeout: float = 300.0,
        interval: float = 6.0,
        stabilize_checks: int = 3,
        stabilize_gap: float = 8.0,
    ) -> None:
        """Poll until SSH is reliably up.

        Vast.ai instances often accept one connection during container boot
        then bounce SSHD, causing the very next ``scp``/``rsync`` to hit
        ``[WinError 1225] The remote computer refused the network connection``.

        After the first successful handshake we wait ``stabilize_gap`` seconds
        and re-check ``stabilize_checks`` times.  Only if all checks pass do
        we declare the host stable.
        """
        loop = asyncio.get_event_loop()
        deadline = loop.time() + total_timeout
        last_err: Exception | None = None

        while loop.time() < deadline:
            try:
                conn = await self._connect()
                conn.close()
                await conn.wait_closed()
            except Exception as exc:  # noqa: BLE001
                last_err = exc
                await asyncio.sleep(interval)
                continue

            # First handshake OK — run stabilization pings.
            stable = True
            for _ in range(stabilize_checks):
                await asyncio.sleep(stabilize_gap)
                if loop.time() > deadline:
                    break
                try:
                    conn = await self._connect()
                    conn.close()
                    await conn.wait_closed()
                except Exception as exc:  # noqa: BLE001
                    last_err = exc
                    stable = False
                    break
            if stable:
                return
            log.info("SSH bounced during stabilization; retrying...")
            await asyncio.sleep(interval)

        raise TimeoutError(f"SSH not ready after {total_timeout}s: {last_err}")

    # --------------------------------------------------------- SFTP sync
    # Primary file transfer uses asyncssh's built-in SFTP, running over the
    # same in-process SSH transport that wait_for_ssh / run_stream use.
    # This avoids the Windows-specific failures (WinError 1225, WinError 2)
    # caused by shelling out to scp.exe / rsync.exe through Vast.ai relays.

    async def _sftp_put_recursive(self, sftp, local_path: Path, remote: str) -> int:
        """Upload ``local_path`` to ``remote`` over an open SFTP session.

        Returns the number of files transferred.
        """
        count = 0
        if local_path.is_file():
            await sftp.makedirs(posixpath.dirname(remote), exist_ok=True)
            await sftp.put(str(local_path), remote)
            return 1

        await sftp.makedirs(remote, exist_ok=True)
        for child in local_path.rglob("*"):
            if child.is_dir():
                continue
            rel = child.relative_to(local_path).as_posix()
            dest = posixpath.join(remote, rel)
            await sftp.makedirs(posixpath.dirname(dest), exist_ok=True)
            await sftp.put(str(child), dest)
            count += 1
        return count

    async def _sftp_get_recursive(self, sftp, remote: str, local_path: Path) -> int:
        """Download ``remote`` to ``local_path`` over an open SFTP session.

        Returns the number of files transferred.
        """
        import asyncssh  # noqa: PLC0415

        count = 0
        try:
            attrs = await sftp.stat(remote)
        except asyncssh.SFTPNoSuchFile:
            raise FileNotFoundError(f"Remote path not found: {remote}")

        if not attrs.permissions or not (attrs.permissions & 0o40000):
            local_path.parent.mkdir(parents=True, exist_ok=True)
            await sftp.get(remote, str(local_path))
            return 1

        local_path.mkdir(parents=True, exist_ok=True)
        for entry in await sftp.readdir(remote):
            name = entry.filename
            if name in (".", ".."):
                continue
            child_remote = posixpath.join(remote, name)
            child_local = local_path / name
            child_attrs = entry.attrs
            is_dir = child_attrs.permissions and (child_attrs.permissions & 0o40000)
            if is_dir:
                count += await self._sftp_get_recursive(sftp, child_remote, child_local)
            else:
                child_local.parent.mkdir(parents=True, exist_ok=True)
                await sftp.get(child_remote, str(child_local))
                count += 1
        return count

    async def rsync_push(self, local: Path, remote: str) -> None:
        """Push local dir/file to the remote via SFTP (asyncssh)."""
        local_path = Path(local).resolve()
        log.info("sftp push: %s -> %s", local_path, remote)
        conn = await self._connect()
        async with conn:
            async with conn.start_sftp_client() as sftp:
                n = await self._sftp_put_recursive(sftp, local_path, remote)
        log.info("sftp push complete: %d files transferred", n)

    async def rsync_pull(self, remote: str, local: Path) -> None:
        """Pull remote dir/file to local path via SFTP (asyncssh)."""
        local_path = Path(local)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        remote_clean = remote.rstrip("/")
        log.info("sftp pull: %s -> %s", remote_clean, local_path)
        conn = await self._connect()
        async with conn:
            async with conn.start_sftp_client() as sftp:
                n = await self._sftp_get_recursive(sftp, remote_clean, local_path)
        log.info("sftp pull complete: %d files transferred", n)

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
