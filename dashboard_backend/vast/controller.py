"""Top-level asyncio orchestration: provision a VPS, rsync inputs, stream
logs, rsync artifacts, destroy. All side effects (DB writes, log fan-out)
flow through injected dependencies so this module is unit-testable with
fakes.
"""
from __future__ import annotations

import asyncio
import logging
import time
from pathlib import Path
from typing import Awaitable, Callable

from .log_hub import LogHub
from .provisioner import VastProvisioner
from .runner import RunCommand, SshEndpoint, SshRunner
from .runs_store import Run, RunStatus, RunsStore

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Artifact promotion — imported lazily in runner so this module stays minimal
# ---------------------------------------------------------------------------
def _default_promoter(run: Run, artifacts_dir: Path, cascade_root: Path) -> None:
    """No-op stub. The artifact-sync todo replaces this with a real promoter."""
    return


PromoteFn = Callable[[Run, Path, Path], None]


class VastController:
    def __init__(
        self,
        *,
        store: RunsStore,
        provisioner: VastProvisioner,
        hub: LogHub,
        orchestrator_image: str,
        logs_dir: Path,
        artifacts_dir: Path,
        cascade_root: Path,
        inputs_to_push: list[Path] | None = None,
        ssh_key_path: str | None = None,
        promote_fn: PromoteFn = _default_promoter,
    ) -> None:
        self.store = store
        self.provisioner = provisioner
        self.hub = hub
        self.orchestrator_image = orchestrator_image
        self.logs_dir = logs_dir
        self.artifacts_dir = artifacts_dir
        self.cascade_root = cascade_root
        self.inputs_to_push = inputs_to_push or []
        self.ssh_key_path = ssh_key_path or None
        self.promote_fn = promote_fn

        self._tasks: dict[str, asyncio.Task[None]] = {}
        # C-23: periodically purge log-hub channels for finished runs so
        # _channels doesn't grow without bound over the controller's
        # lifetime.  Kicked off lazily on the first launch() so unit tests
        # that never await launch() don't have to deal with the task.
        self._gc_task: asyncio.Task[None] | None = None
        self._gc_interval_s: float = 300.0
        self._gc_age_s: float = 600.0

    async def _gc_loop(self) -> None:
        """C-23: background coroutine that drops closed channels older
        than ``self._gc_age_s``.  Cancelled cleanly on controller shutdown.
        """
        while True:
            try:
                await asyncio.sleep(self._gc_interval_s)
                purged = await self.hub.purge_closed_channels(self._gc_age_s)
                if purged:
                    log.info("log_hub GC: purged %d closed channels", purged)
            except asyncio.CancelledError:
                raise
            except Exception as exc:  # noqa: BLE001
                log.warning("log_hub GC iteration failed: %s", exc)

    def _ensure_gc_running(self) -> None:
        if self._gc_task is None or self._gc_task.done():
            self._gc_task = asyncio.create_task(self._gc_loop(), name="log-hub-gc")

    async def shutdown(self) -> None:
        """Cancel the GC task; call from FastAPI lifespan teardown."""
        if self._gc_task and not self._gc_task.done():
            self._gc_task.cancel()
            try:
                await self._gc_task
            except (asyncio.CancelledError, Exception):
                pass
            self._gc_task = None

    # ----------------------------------------------------------- public API

    async def launch(
        self,
        *,
        offer_id: int,
        baseline_ids: list[str],
        crrna_lookup_ids: list[str] | None = None,
        max_generations: int = 12,
        variants_per_gen: int = 5,
        workers: int = 3,
        label: str | None = None,
    ) -> Run:
        run_id = RunsStore.new_id()
        run = Run(
            id=run_id,
            label=label or f"run-{run_id}",
            status=RunStatus.QUEUED,
            baseline_ids=baseline_ids,
            crrna_lookup_ids=crrna_lookup_ids or list(baseline_ids),
            max_generations=max_generations,
            variants_per_gen=variants_per_gen,
            workers=workers,
            offer_id=offer_id,
            docker_image=self.orchestrator_image,
            logs_path=str(self.logs_dir / f"{run_id}.log"),
            artifacts_path=str(self.artifacts_dir / run_id),
        )
        self.store.create(run)

        self._ensure_gc_running()
        task = asyncio.create_task(self._drive(run), name=f"run-{run_id}")
        self._tasks[run_id] = task
        return run

    async def cancel(self, run_id: str) -> None:
        run = self.store.get(run_id)
        if run is None:
            return
        # Best-effort: destroy instance first (stops billing), then mark.
        if run.instance_id:
            try:
                await self.provisioner.destroy_instance(run.instance_id)
            except Exception as exc:  # noqa: BLE001
                log.warning("destroy_instance during cancel failed: %s", exc)
        self.store.update(
            run_id,
            status=RunStatus.CANCELLED,
            finished_at=time.time(),
        )
        await self.hub.close(run_id)
        task = self._tasks.get(run_id)
        if task and not task.done():
            task.cancel()

    # ---------------------------------------------------------- task body

    async def _drive(self, run: Run) -> None:
        hub = self.hub
        store = self.store
        run_id = run.id

        async def emit(line: str) -> None:
            await hub.publish(run_id, line)

        try:
            # 1. PROVISION
            store.update(run_id, status=RunStatus.PROVISIONING)
            await emit(f"[controller] provisioning offer {run.offer_id}...")

            cmd = RunCommand(
                run_id=run_id,
                baseline_ids=run.baseline_ids,
                crrna_lookup_ids=run.crrna_lookup_ids,
                max_generations=run.max_generations,
                variants_per_gen=run.variants_per_gen,
                workers=run.workers,
            )

            instance_id = await self.provisioner.create_instance(
                offer_id=run.offer_id or 0,
                image=self.orchestrator_image,
                disk_gb=80,
                onstart_cmd=cmd.onstart_cmd(),
                env_vars={
                    "CASCADE_RUN_ID": run_id,
                    "CASCADE_BASELINE_IDS": ",".join(run.baseline_ids),
                    "CASCADE_CRRNA_LOOKUP_IDS": ",".join(run.crrna_lookup_ids),
                    "CASCADE_MAX_GENERATIONS": str(run.max_generations),
                    "CASCADE_VARIANTS_PER_GEN": str(run.variants_per_gen),
                    "CASCADE_WORKERS": str(run.workers),
                },
                label=run.label,
            )
            store.update(run_id, instance_id=instance_id)
            await emit(f"[controller] instance {instance_id} requested; waiting for SSH...")

            # 2. WAIT FOR SSH
            store.update(run_id, status=RunStatus.STARTING)
            inst = await self.provisioner.wait_until_ready(instance_id)
            store.update(
                run_id,
                ssh_host=inst.ssh_host,
                ssh_port=inst.ssh_port,
                dph_usd=inst.dph_total,
                started_at=time.time(),
            )
            await emit(f"[controller] ssh up: {inst.ssh_host}:{inst.ssh_port}")
            await emit("[controller] waiting 15s for container to finish booting...")
            await asyncio.sleep(15)

            endpoint = SshEndpoint(
                host=inst.ssh_host or "",
                port=inst.ssh_port or 22,
                user="root",
                key_path=self.ssh_key_path,
            )
            runner = SshRunner(
                endpoint,
                logs_dir=self.logs_dir,
                artifacts_dir=self.artifacts_dir,
            )

            # 3. PUSH INPUTS (retry with backoff — Vast SSH can bounce
            #    during container startup even after wait_for_ssh returns).
            max_upload_attempts = 4
            for p in self.inputs_to_push:
                if not p.exists():
                    await emit(f"[controller] skip input: {p} (missing)")
                    continue
                remote = f"/workspace/CASCADE/{p.relative_to(self.cascade_root).as_posix()}"
                last_upload_err: Exception | None = None
                for attempt in range(1, max_upload_attempts + 1):
                    await emit(
                        f"[controller] upload ({attempt}/{max_upload_attempts}): "
                        f"{p.name} -> {remote}"
                    )
                    try:
                        await runner.rsync_push(p, remote)
                        last_upload_err = None
                        break
                    except Exception as exc:  # noqa: BLE001
                        last_upload_err = exc
                        if attempt < max_upload_attempts:
                            delay = 10 * attempt
                            await emit(
                                f"[controller] upload attempt {attempt} failed "
                                f"({exc}); retrying in {delay}s..."
                            )
                            await asyncio.sleep(delay)
                if last_upload_err is not None:
                    await emit(
                        f"[controller] WARN upload failed after "
                        f"{max_upload_attempts} attempts: {last_upload_err}"
                    )

            # 4. STREAM: the orchestrator is already running via --onstart-cmd,
            #    so we attach to its tee'd log with a follow tail.
            store.update(run_id, status=RunStatus.RUNNING)
            tail_cmd = (
                f"touch /workspace/CASCADE/logs/run-{run_id}.log && "
                f"tail -n +1 -F /workspace/CASCADE/logs/run-{run_id}.log"
            )
            local_log = self.logs_dir / f"{run_id}.log"
            rc = await runner.run_stream(
                tail_cmd, on_line=emit, local_log_path=local_log
            )
            await emit(f"[controller] remote tail exited rc={rc}")

            # 5. PULL ARTIFACTS (with retry — network can be flaky at end of run)
            store.update(run_id, status=RunStatus.SYNCING)
            remote_artifacts = f"/workspace/CASCADE/outputs/runs/{run_id}/"
            local_artifacts = self.artifacts_dir / run_id
            pull_ok = False
            for attempt in range(1, 4):
                await emit(
                    f"[controller] download artifacts ({attempt}/3): "
                    f"{remote_artifacts} -> {local_artifacts}"
                )
                try:
                    await runner.rsync_pull(remote_artifacts, local_artifacts)
                    pull_ok = True
                    break
                except Exception as exc:  # noqa: BLE001
                    if attempt < 3:
                        delay = 10 * attempt
                        await emit(
                            f"[controller] download attempt {attempt} failed "
                            f"({exc}); retrying in {delay}s..."
                        )
                        await asyncio.sleep(delay)
                    else:
                        await emit(f"[controller] WARN download failed after 3 attempts: {exc}")

            # 6. PROMOTE
            try:
                self.promote_fn(run, local_artifacts, self.cascade_root)
                await emit("[controller] artifacts promoted into canonical outputs/")
            except Exception as exc:  # noqa: BLE001
                await emit(f"[controller] WARN promotion failed: {exc}")

            store.update(
                run_id,
                status=RunStatus.COMPLETED,
                finished_at=time.time(),
            )
            await emit("[controller] run completed.")

        except asyncio.CancelledError:
            await emit("[controller] run cancelled.")
            raise
        except Exception as exc:  # noqa: BLE001
            log.exception("run %s failed", run_id)
            self.store.update(
                run_id,
                status=RunStatus.FAILED,
                finished_at=time.time(),
                error=str(exc)[:2000],
            )
            await self.hub.publish(run_id, f"[controller] ERROR: {exc}")
        finally:
            # Always try to destroy the instance so we stop paying.
            current = self.store.get(run_id)
            if current and current.instance_id and current.status in (
                RunStatus.COMPLETED,
                RunStatus.FAILED,
            ):
                try:
                    await self.provisioner.destroy_instance(current.instance_id)
                    await self.hub.publish(run_id, "[controller] instance destroyed.")
                except Exception as exc:  # noqa: BLE001
                    log.warning("destroy on exit failed: %s", exc)
            await self.hub.close(run_id)
            self._tasks.pop(run_id, None)
