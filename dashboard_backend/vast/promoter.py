"""Promote a run's synced artifacts into the canonical outputs/ tree.

After the controller rsyncs ``outputs/runs/<run_id>/`` back, this module
copies the subset of files the dashboard's existing
:class:`dashboard_backend.service.DashboardService` indexes into:

  * ``outputs/optimized_switches/``
  * ``outputs/rl_gym_data/rl_training_dataset.jsonl`` (append-merged)
  * ``outputs/fast_eval/`` (for 3Dmol /api/structure-file lookups)
  * ``outputs/high_fidelity_scoring/`` (ditto, for base-tier tertiary)

Collision strategy: every filename gets a ``<run_id>__`` prefix so parallel
runs that happen to optimize the same variant label don't clobber each
other's outputs. The RL dataset is append-only so we just concatenate.

Use :func:`promote_run` from the controller's ``promote_fn`` hook and from
CLI tools (``python -m dashboard_backend.vast.promoter <run_id>``).
"""
from __future__ import annotations

import argparse
import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

log = logging.getLogger(__name__)


@dataclass
class PromotionReport:
    run_id: str
    copied: list[str] = field(default_factory=list)
    appended_records: int = 0
    skipped: list[str] = field(default_factory=list)

    @property
    def total(self) -> int:
        return len(self.copied) + self.appended_records


# Relative sub-paths inside the run tree to mirror into the canonical tree.
# Only directories the dashboard needs to see end up here; everything else
# stays under outputs/runs/<run_id>/ for provenance/debugging.
_DIR_MIRRORS = [
    ("optimized_switches", "outputs/optimized_switches"),
    ("fast_eval", "outputs/fast_eval"),
    ("high_fidelity_scoring", "outputs/high_fidelity_scoring"),
]
_RL_SUBPATH = ("rl_gym_data/rl_training_dataset.jsonl", "outputs/rl_gym_data/rl_training_dataset.jsonl")


def _iter_copyable(src_dir: Path) -> Iterable[Path]:
    if not src_dir.exists():
        return []
    return (p for p in src_dir.rglob("*") if p.is_file())


def _safe_copy(src: Path, dest: Path, prefix: str) -> Path:
    """Copy ``src`` into ``dest`` directory with a collision-safe name.

    Preserves relative sub-tree structure under ``dest`` but prefixes the
    top-level basename with ``{prefix}__`` so parallel runs don't overwrite
    each other. Uses ``shutil.copy2`` to keep mtime for the dashboard's
    "latest first" ordering.
    """
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / f"{prefix}__{src.name}"
    shutil.copy2(src, target)
    return target


def promote_run(
    run_id: str,
    artifacts_dir: Path,
    cascade_root: Path,
    *,
    service_cache_bump: "object | None" = None,
) -> PromotionReport:
    """Copy artifacts from ``artifacts_dir/<run_id>`` into the canonical tree.

    ``service_cache_bump`` — if provided and has ``invalidate()`` — will be
    invoked after promotion so the dashboard re-reads outputs. Typed as
    object so we don't hard-couple to DashboardService here.
    """
    report = PromotionReport(run_id=run_id)

    run_root = artifacts_dir / run_id if artifacts_dir.name != run_id else artifacts_dir
    if not run_root.exists():
        log.warning("promote_run: run root does not exist: %s", run_root)
        report.skipped.append(str(run_root))
        return report

    # 1. Mirror whole directories into canonical tree (collision-prefixed).
    for sub, dest_rel in _DIR_MIRRORS:
        src_dir = run_root / sub
        if not src_dir.exists():
            continue
        dest_dir = cascade_root / dest_rel
        for src in _iter_copyable(src_dir):
            rel_parent = src.parent.relative_to(src_dir)
            target_dir = dest_dir / rel_parent
            try:
                target = _safe_copy(src, target_dir, run_id)
                report.copied.append(str(target.relative_to(cascade_root)))
            except Exception as exc:  # noqa: BLE001
                log.warning("promote copy failed: %s -> %s: %s", src, target_dir, exc)
                report.skipped.append(str(src))

    # 2. Append the RL dataset records to the canonical JSONL so the
    #    dashboard's DashboardService.load_variants() picks them up.
    rl_src = run_root / _RL_SUBPATH[0]
    rl_dest = cascade_root / _RL_SUBPATH[1]
    if rl_src.exists():
        rl_dest.parent.mkdir(parents=True, exist_ok=True)
        with rl_src.open("r", encoding="utf-8") as fin, rl_dest.open(
            "a", encoding="utf-8"
        ) as fout:
            for line in fin:
                if line.strip():
                    fout.write(line.rstrip("\n") + "\n")
                    report.appended_records += 1

    # 3. Best-effort cache invalidate on the DashboardService singleton.
    if service_cache_bump is not None and hasattr(service_cache_bump, "invalidate"):
        try:
            service_cache_bump.invalidate()
        except Exception as exc:  # noqa: BLE001
            log.warning("service cache invalidate failed: %s", exc)

    log.info(
        "promoted run %s: copied=%d appended=%d skipped=%d",
        run_id,
        len(report.copied),
        report.appended_records,
        len(report.skipped),
    )
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _main() -> int:
    parser = argparse.ArgumentParser(
        description="Promote a synced CASCADE run into canonical outputs/."
    )
    parser.add_argument("run_id", help="Run identifier (the subdir under outputs/runs/)")
    parser.add_argument(
        "--artifacts-dir",
        default="outputs/runs",
        help="Parent dir holding run-scoped synced artifacts (default: outputs/runs).",
    )
    parser.add_argument(
        "--cascade-root",
        default=".",
        help="CASCADE repo root where canonical outputs/ lives (default: cwd).",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    report = promote_run(
        args.run_id,
        Path(args.artifacts_dir).resolve(),
        Path(args.cascade_root).resolve(),
    )
    print(
        f"run_id={report.run_id} copied={len(report.copied)} "
        f"appended={report.appended_records} skipped={len(report.skipped)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
