import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dashboard_backend.config import DashboardConfig
from dashboard_backend.service import DashboardService
from dashboard_backend.storage.base import VariantCatalogStore


class DummyStore(VariantCatalogStore):
    def fetch_all_variants(self):
        return []

    def ping(self):
        return {"ok": True, "backend": "dummy"}


def _mk_cfg(root: Path) -> DashboardConfig:
    return DashboardConfig(
        cascade_root=root,
        sqlite_db_path=root / "metadata" / "cas13_variants.db",
        rl_dataset_path=root / "outputs" / "rl_gym_data" / "rl_training_dataset.jsonl",
        domain_metadata_path=root / "metadata" / "variant_domain_metadata.json",
        optimized_dir=root / "outputs" / "optimized_switches",
        fast_eval_dir=root / "outputs" / "fast_eval",
        hf_eval_dir=root / "outputs" / "high_fidelity_scoring",
        min_iptm=0.85,
        min_af2_ig=0.80,
        max_on_distance=12.0,
    )


def test_hybrid_optimized_classification(tmpdir):
    root = Path(tmpdir)
    (root / "outputs" / "rl_gym_data").mkdir(parents=True)
    (root / "metadata").mkdir(parents=True)
    (root / "outputs" / "optimized_switches").mkdir(parents=True)
    (root / "outputs" / "fast_eval").mkdir(parents=True)
    (root / "outputs" / "high_fidelity_scoring").mkdir(parents=True)

    record = {
        "variant_id": "Laaaaaa_g01_v00",
        "generation": 1,
        "baseline_id": "baseline_id",
        "crrna_lookup_id": "baseline_id",
        "fitness": 10.0,
        "off_dist_A": 30.0,
        "on_dist_A": 8.0,
        "iptm": 0.9,
        "af2_ig": 0.85,
        "is_elite": False,
    }
    (root / "outputs" / "rl_gym_data" / "rl_training_dataset.jsonl").write_text(
        json.dumps(record) + "\n", encoding="utf-8"
    )
    (root / "metadata" / "variant_domain_metadata.json").write_text(
        json.dumps({"baseline_id": {"domains": {"HEPN1": {"start": 1, "end": 10}, "HEPN2": {"start": 20, "end": 30}}}}),
        encoding="utf-8",
    )

    svc = DashboardService(_mk_cfg(root), DummyStore())
    rows = svc.load_variants()
    assert len(rows) == 1
    assert rows[0]["optimized_switch"] is True
    assert "threshold_pass" in rows[0]["optimized_reasons"]


def test_overview_schema_warning(tmpdir):
    root = Path(tmpdir)
    (root / "outputs" / "rl_gym_data").mkdir(parents=True)
    (root / "metadata").mkdir(parents=True)
    (root / "outputs" / "optimized_switches").mkdir(parents=True)
    (root / "outputs" / "fast_eval").mkdir(parents=True)
    (root / "outputs" / "high_fidelity_scoring").mkdir(parents=True)

    # Missing many required fields on purpose
    (root / "outputs" / "rl_gym_data" / "rl_training_dataset.jsonl").write_text(
        json.dumps({"variant_id": "only_id"}) + "\n", encoding="utf-8"
    )
    (root / "metadata" / "variant_domain_metadata.json").write_text("{}", encoding="utf-8")

    svc = DashboardService(_mk_cfg(root), DummyStore())
    overview = svc.get_overview()
    warnings = overview.get("warnings", [])
    assert any(str(w).startswith("schema_drift_missing_required_fields") for w in warnings)
