import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dashboard_backend.api.baselines import build_baselines
from dashboard_backend.config import DashboardConfig
from dashboard_backend.storage.base import VariantCatalogRecord, VariantCatalogStore


class DummyCatalog(VariantCatalogStore):
    def __init__(self, records=None):
        self._records = records or []

    def fetch_all_variants(self):
        return list(self._records)

    def ping(self):
        return {"ok": True}


def _cfg(root: Path) -> DashboardConfig:
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


def _write_metadata(root: Path, baselines: dict) -> None:
    (root / "metadata").mkdir(parents=True, exist_ok=True)
    (root / "metadata" / "variant_domain_metadata.json").write_text(
        json.dumps(baselines), encoding="utf-8"
    )


def _write_base_json(root: Path, baseline_id: str, repeat: str, spacer: str) -> None:
    (root / "jsons").mkdir(parents=True, exist_ok=True)
    payload = [
        {
            "name": baseline_id,
            "sequences": [
                {"proteinChain": {"sequence": "M" * 50, "count": 1}},
                {"rnaSequence": {"sequence": f"{repeat}{spacer}", "count": 1}},
            ],
        }
    ]
    (root / "jsons" / f"{baseline_id}.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )


def test_build_baselines_merges_all_sources(tmpdir):
    root = Path(tmpdir)
    _write_metadata(
        root,
        {
            "bl_alpha": {
                "sequence_length": 700,
                "subtype": "cas13a",
                "domains": {"HEPN1": {"start": 10, "end": 100}, "HEPN2": {"start": 400, "end": 500}},
                "crRNA_repeat_used": "CGGGUGUAGCUCAGUUGGUUAGAGCGCC",
            },
            "bl_beta": {
                "sequence_length": 800,
                "subtype": "cas13b",
                "domains": {"HEPN1": {"start": 5, "end": 95}, "HEPN2": {"start": 600, "end": 700}},
                "crRNA_repeat_used": "GAGACG",
            },
        },
    )
    _write_base_json(root, "bl_alpha", "CGGGUGUAGCUCAGUUGGUUAGAGCGCC", "ACGUACGU")

    (root / "outputs" / "phase1_screening" / "bl_alpha_pred").mkdir(parents=True)
    (root / "outputs").mkdir(exist_ok=True)
    (root / "outputs" / "validated_baseline_ids.txt").write_text("bl_alpha\n", encoding="utf-8")

    catalog = DummyCatalog(
        [
            VariantCatalogRecord(
                sequence_id="bl_alpha",
                sra_accession="SRR000001",
                score=0.92,
                hepn1_start=10,
                hepn1_end=100,
                hepn2_start=400,
                hepn2_end=500,
                status="kept",
                reason=None,
            )
        ]
    )

    rows = build_baselines(_cfg(root), catalog)
    assert len(rows) == 2

    alpha = next(r for r in rows if r.baseline_id == "bl_alpha")
    assert alpha.validated is True
    assert alpha.has_phase1_structure is True
    assert alpha.crrna_repeat == "CGGGUGUAGCUCAGUUGGUUAGAGCGCC"
    assert alpha.crrna_spacer == "ACGUACGU"
    assert alpha.sra_accession == "SRR000001"
    assert alpha.hepn2_end == 500

    beta = next(r for r in rows if r.baseline_id == "bl_beta")
    assert beta.validated is False
    assert beta.has_phase1_structure is False
    assert beta.crrna_spacer is None  # no base json -> no spacer known

    # Sort order: validated + phase1-structured first
    assert rows[0].baseline_id == "bl_alpha"
