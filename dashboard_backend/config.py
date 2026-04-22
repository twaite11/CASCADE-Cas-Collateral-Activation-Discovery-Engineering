from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DashboardConfig:
    cascade_root: Path
    sqlite_db_path: Path
    rl_dataset_path: Path
    domain_metadata_path: Path
    optimized_dir: Path
    fast_eval_dir: Path
    hf_eval_dir: Path
    min_iptm: float
    min_af2_ig: float
    max_on_distance: float
    # Vast.ai / runs controller paths (optional; only used by /api/runs)
    runs_db_path: Path = Path("metadata/runs.db")
    runs_logs_dir: Path = Path("outputs/run_logs")
    runs_artifacts_dir: Path = Path("outputs/runs")
    orchestrator_image: str = "ghcr.io/twaite11/cascade-orchestrator:latest"
    vast_ssh_key_path: str = ""


def load_config() -> DashboardConfig:
    root_override = os.environ.get("CASCADE_ROOT", "").strip()
    cascade_root = (
        Path(root_override).resolve()
        if root_override
        else Path(__file__).resolve().parents[1]
    )

    return DashboardConfig(
        cascade_root=cascade_root,
        sqlite_db_path=cascade_root / "metadata" / "cas13_variants.db",
        rl_dataset_path=cascade_root / "outputs" / "rl_gym_data" / "rl_training_dataset.jsonl",
        domain_metadata_path=cascade_root / "metadata" / "variant_domain_metadata.json",
        optimized_dir=cascade_root / "outputs" / "optimized_switches",
        fast_eval_dir=cascade_root / "outputs" / "fast_eval",
        hf_eval_dir=cascade_root / "outputs" / "high_fidelity_scoring",
        min_iptm=float(os.environ.get("DASH_MIN_IPTM", "0.85")),
        min_af2_ig=float(os.environ.get("DASH_MIN_AF2_IG", "0.80")),
        max_on_distance=float(os.environ.get("DASH_MAX_ON_DISTANCE", "12.0")),
        runs_db_path=cascade_root / "metadata" / "runs.db",
        runs_logs_dir=cascade_root / "outputs" / "run_logs",
        runs_artifacts_dir=cascade_root / "outputs" / "runs",
        orchestrator_image=os.environ.get(
            "CASCADE_ORCH_IMAGE",
            "ghcr.io/twaite11/cascade-orchestrator:latest",
        ),
        vast_ssh_key_path=os.environ.get("VAST_SSH_KEY_PATH", "").strip(),
    )
