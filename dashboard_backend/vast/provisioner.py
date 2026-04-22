"""Thin async wrapper around the ``vastai`` CLI.

The controller shells out to ``vastai`` (pip installable, auth via
``~/.config/vastai/vast_api_key``) instead of driving the REST API directly,
because the CLI is guaranteed to keep up with upstream changes and handles
token refresh transparently. All methods return parsed JSON via ``--raw``.
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import shlex
import shutil
from dataclasses import dataclass, field
from typing import Any, Sequence

log = logging.getLogger(__name__)


class VastCliError(RuntimeError):
    """Raised when ``vastai`` exits non-zero or returns unparseable output."""

    def __init__(self, cmd: Sequence[str], returncode: int, stderr: str) -> None:
        super().__init__(
            f"vastai failed ({returncode}): {' '.join(shlex.quote(c) for c in cmd)}\n{stderr.strip()}"
        )
        self.cmd = list(cmd)
        self.returncode = returncode
        self.stderr = stderr


@dataclass
class VastOffer:
    id: int
    gpu_name: str | None
    num_gpus: int | None
    gpu_ram_gb: float | None
    cpu_cores: int | None
    cpu_ram_gb: float | None
    disk_space_gb: float | None
    dph_total: float | None
    inet_down_mbps: float | None
    inet_up_mbps: float | None
    datacenter: str | None
    reliability: float | None
    raw: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "gpu_name": self.gpu_name,
            "num_gpus": self.num_gpus,
            "gpu_ram_gb": self.gpu_ram_gb,
            "cpu_cores": self.cpu_cores,
            "cpu_ram_gb": self.cpu_ram_gb,
            "disk_space_gb": self.disk_space_gb,
            "dph_total": self.dph_total,
            "inet_down_mbps": self.inet_down_mbps,
            "inet_up_mbps": self.inet_up_mbps,
            "datacenter": self.datacenter,
            "reliability": self.reliability,
        }


@dataclass
class VastInstance:
    id: int
    status: str | None
    actual_status: str | None
    ssh_host: str | None
    ssh_port: int | None
    dph_total: float | None
    image: str | None
    label: str | None
    raw: dict[str, Any] = field(default_factory=dict)

    @property
    def ready(self) -> bool:
        # Vast.ai marks "running" when the container is up and SSH is listening.
        if not (self.ssh_host and self.ssh_port):
            return False
        return (self.actual_status or self.status or "").lower() == "running"

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "status": self.status,
            "actual_status": self.actual_status,
            "ssh_host": self.ssh_host,
            "ssh_port": self.ssh_port,
            "dph_total": self.dph_total,
            "image": self.image,
            "label": self.label,
            "ready": self.ready,
        }


class VastProvisioner:
    """Async facade over the ``vastai`` CLI.

    Each method wraps a single CLI invocation with a timeout and a structured
    error. The ``vastai`` binary path can be overridden via ``bin_path`` for
    tests or multi-env setups.
    """

    def __init__(
        self,
        bin_path: str | None = None,
        api_key: str | None = None,
        default_timeout: float = 60.0,
    ) -> None:
        self.bin_path = bin_path or shutil.which("vastai") or "vastai"
        self.api_key = api_key or os.environ.get("VAST_API_KEY", "").strip() or None
        self.default_timeout = default_timeout

    # ------------------------------------------------------------------ core

    async def _run(
        self,
        args: Sequence[str],
        *,
        timeout: float | None = None,
        parse_json: bool = True,
    ) -> Any:
        cmd = [self.bin_path, *args]
        env = os.environ.copy()
        if self.api_key:
            # CLI accepts --api-key or $VAST_API_KEY; using env avoids leaking
            # the key into process args visible to other local users.
            env["VAST_API_KEY"] = self.api_key

        log.debug("vastai: %s", " ".join(shlex.quote(c) for c in cmd))
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env=env,
        )
        try:
            stdout_b, stderr_b = await asyncio.wait_for(
                proc.communicate(), timeout=timeout or self.default_timeout
            )
        except asyncio.TimeoutError:
            proc.kill()
            await proc.wait()
            raise VastCliError(cmd, -1, f"timeout after {timeout or self.default_timeout}s")

        stdout = stdout_b.decode("utf-8", errors="replace")
        stderr = stderr_b.decode("utf-8", errors="replace")

        if proc.returncode != 0:
            raise VastCliError(cmd, proc.returncode or 1, stderr or stdout)

        if not parse_json:
            return stdout

        stripped = stdout.strip()
        if not stripped:
            return None
        try:
            return json.loads(stripped)
        except json.JSONDecodeError as exc:
            raise VastCliError(cmd, 0, f"non-JSON output: {exc}: {stripped[:500]}")

    # --------------------------------------------------------------- helpers

    @staticmethod
    def _offer_from_raw(row: dict[str, Any]) -> VastOffer:
        gpu_ram_mb = _maybe_float(row.get("gpu_ram"))
        cpu_ram_mb = _maybe_float(row.get("cpu_ram"))
        return VastOffer(
            id=int(row.get("id") or row.get("ask_id") or 0),
            gpu_name=row.get("gpu_name"),
            num_gpus=row.get("num_gpus"),
            gpu_ram_gb=(gpu_ram_mb / 1024) if gpu_ram_mb else None,
            cpu_cores=row.get("cpu_cores"),
            cpu_ram_gb=(cpu_ram_mb / 1024) if cpu_ram_mb else None,
            disk_space_gb=_maybe_float(row.get("disk_space")),
            dph_total=_maybe_float(row.get("dph_total")),
            inet_down_mbps=_maybe_float(row.get("inet_down")),
            inet_up_mbps=_maybe_float(row.get("inet_up")),
            datacenter=row.get("datacenter") or row.get("geolocation"),
            reliability=_maybe_float(row.get("reliability2") or row.get("reliability")),
            raw=row,
        )

    @staticmethod
    def _instance_from_raw(row: dict[str, Any]) -> VastInstance:
        # Vast.ai's "show instance" returns different shapes depending on
        # version; ssh_host/ssh_port live either top-level or inside "ports".
        ssh_host = row.get("ssh_host") or row.get("public_ipaddr")
        ssh_port = row.get("ssh_port")
        if not ssh_port:
            ports = row.get("ports") or {}
            # Vast exposes container port 22 mapped to a host port; find the mapping.
            if isinstance(ports, dict):
                mapping = ports.get("22/tcp") or ports.get("22")
                if isinstance(mapping, list) and mapping:
                    hp = mapping[0]
                    if isinstance(hp, dict):
                        ssh_port = int(hp.get("HostPort") or 0) or None
                        ssh_host = ssh_host or hp.get("HostIp")
        return VastInstance(
            id=int(row.get("id") or 0),
            status=row.get("status") or row.get("cur_state"),
            actual_status=row.get("actual_status"),
            ssh_host=ssh_host,
            ssh_port=int(ssh_port) if ssh_port else None,
            dph_total=_maybe_float(row.get("dph_total")),
            image=row.get("image_uuid") or row.get("image"),
            label=row.get("label"),
            raw=row,
        )

    # --------------------------------------------------------------- commands

    async def search_offers(
        self,
        query: str = "gpu_name=A100_SXM4 num_gpus=1 rentable=true verified=true",
        *,
        on_demand: bool = True,
        limit: int = 10,
        order: str = "dph_total",
    ) -> list[VastOffer]:
        args = [
            "search",
            "offers",
            query,
            "--raw",
            "-o",
            order,
            "--limit",
            str(limit),
        ]
        if on_demand:
            args.append("--on-demand")
        data = await self._run(args, timeout=45.0)
        if not isinstance(data, list):
            return []
        return [self._offer_from_raw(row) for row in data]

    async def create_instance(
        self,
        offer_id: int,
        image: str,
        *,
        disk_gb: int = 80,
        onstart_cmd: str,
        env_vars: dict[str, str] | None = None,
        label: str | None = None,
        docker_login: str | None = None,
    ) -> int:
        """Create an SSH direct-mode instance and return its numeric id."""
        env_str = ""
        if env_vars:
            # vastai --env is a docker-style string. Use -e k=v pairs.
            env_str = " ".join(f"-e {shlex.quote(f'{k}={v}')}" for k, v in env_vars.items())

        args = [
            "create",
            "instance",
            str(offer_id),
            "--image",
            image,
            "--ssh",
            "--direct",
            "--disk",
            str(disk_gb),
            "--onstart-cmd",
            onstart_cmd,
            "--raw",
        ]
        if env_str:
            args.extend(["--env", env_str])
        if label:
            args.extend(["--label", label])
        if docker_login:
            args.extend(["--login", docker_login])

        data = await self._run(args, timeout=90.0)
        # CLI returns {"success": true, "new_contract": <id>} on success.
        if not isinstance(data, dict):
            raise VastCliError(
                [self.bin_path, *args], 0, f"unexpected response: {data!r}"
            )
        if not data.get("success", True):
            raise VastCliError([self.bin_path, *args], 1, json.dumps(data))
        inst = data.get("new_contract") or data.get("instance_id") or data.get("id")
        if inst is None:
            raise VastCliError([self.bin_path, *args], 0, f"no instance id in: {data!r}")
        return int(inst)

    async def show_instance(self, instance_id: int) -> VastInstance | None:
        data = await self._run(
            ["show", "instance", str(instance_id), "--raw"], timeout=30.0
        )
        if not data:
            return None
        if isinstance(data, list):
            data = data[0] if data else None
        if not isinstance(data, dict):
            return None
        return self._instance_from_raw(data)

    async def destroy_instance(self, instance_id: int) -> None:
        await self._run(
            ["destroy", "instance", str(instance_id), "--raw"],
            timeout=45.0,
            parse_json=False,
        )

    async def wait_until_ready(
        self,
        instance_id: int,
        *,
        total_timeout: float = 600.0,
        poll_interval: float = 8.0,
    ) -> VastInstance:
        """Poll ``show instance`` until SSH is up or timeout fires."""
        deadline = asyncio.get_event_loop().time() + total_timeout
        last: VastInstance | None = None
        while True:
            try:
                inst = await self.show_instance(instance_id)
            except VastCliError as exc:
                log.warning("show_instance failed (will retry): %s", exc)
                inst = None
            if inst is not None:
                last = inst
                if inst.ready:
                    return inst
            if asyncio.get_event_loop().time() >= deadline:
                raise TimeoutError(
                    f"instance {instance_id} not ready after {total_timeout}s; last={last}"
                )
            await asyncio.sleep(poll_interval)


def _maybe_float(v: Any) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
