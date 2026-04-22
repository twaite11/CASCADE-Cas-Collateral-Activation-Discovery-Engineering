import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dashboard_backend.vast.promoter import promote_run


def _populate_run_tree(root: Path, run_id: str) -> Path:
    """Build a realistic synced run tree the controller would produce."""
    run_root = root / run_id
    # Optimized switches
    opt = run_root / "optimized_switches"
    opt.mkdir(parents=True)
    (opt / "Laaa_g01_v00_optimal.fasta").write_text(">x\nMA\n", encoding="utf-8")
    (opt / "Laaa_g01_v00_ternary_complex.pdb").write_text("REMARK dummy\n", encoding="utf-8")
    # Fast eval JSON
    fe = run_root / "fast_eval"
    fe.mkdir(parents=True)
    (fe / "Laaa_g01_v00_ON.json").write_text(json.dumps({"metric": 1}), encoding="utf-8")
    # High fidelity PDB
    hf = run_root / "high_fidelity_scoring"
    hf.mkdir(parents=True)
    (hf / "Laaa_g01_v00_ternary_complex.pdb").write_text("HEADER\n", encoding="utf-8")
    # RL gym dataset
    gym = run_root / "rl_gym_data"
    gym.mkdir(parents=True)
    (gym / "rl_training_dataset.jsonl").write_text(
        json.dumps({"variant_id": "Laaa_g01_v00", "fitness": 5.0}) + "\n",
        encoding="utf-8",
    )
    return run_root


def test_promote_run_mirrors_dirs_and_appends_rl(tmpdir):
    cascade_root = Path(tmpdir) / "cascade"
    (cascade_root / "outputs" / "rl_gym_data").mkdir(parents=True)
    (cascade_root / "outputs" / "rl_gym_data" / "rl_training_dataset.jsonl").write_text(
        json.dumps({"variant_id": "pre_existing", "fitness": 1.0}) + "\n",
        encoding="utf-8",
    )

    artifacts_dir = cascade_root / "outputs" / "runs"
    _populate_run_tree(artifacts_dir, "runABC")

    report = promote_run("runABC", artifacts_dir, cascade_root)

    assert report.appended_records == 1
    # Pre-existing record preserved, new one appended
    rl = (cascade_root / "outputs" / "rl_gym_data" / "rl_training_dataset.jsonl").read_text(
        encoding="utf-8"
    )
    assert "pre_existing" in rl
    assert "Laaa_g01_v00" in rl

    # Optimized switches mirrored with collision-safe prefix
    fa = cascade_root / "outputs" / "optimized_switches" / "runABC__Laaa_g01_v00_optimal.fasta"
    assert fa.exists()
    pdb = cascade_root / "outputs" / "optimized_switches" / "runABC__Laaa_g01_v00_ternary_complex.pdb"
    assert pdb.exists()

    # Fast eval + HF mirrored
    assert (
        cascade_root / "outputs" / "fast_eval" / "runABC__Laaa_g01_v00_ON.json"
    ).exists()
    assert (
        cascade_root
        / "outputs"
        / "high_fidelity_scoring"
        / "runABC__Laaa_g01_v00_ternary_complex.pdb"
    ).exists()

    # Total counts reflect all copied files + appended records
    assert len(report.copied) >= 4
    assert report.total >= 5


def test_promote_run_missing_tree_returns_skipped(tmpdir):
    cascade_root = Path(tmpdir) / "cascade"
    artifacts_dir = cascade_root / "outputs" / "runs"
    artifacts_dir.mkdir(parents=True)
    report = promote_run("nonexistent", artifacts_dir, cascade_root)
    assert report.copied == []
    assert report.appended_records == 0
    assert len(report.skipped) == 1


def test_promote_run_invalidates_service_cache(tmpdir):
    cascade_root = Path(tmpdir) / "cascade"
    artifacts_dir = cascade_root / "outputs" / "runs"
    _populate_run_tree(artifacts_dir, "runXYZ")
    (cascade_root / "outputs" / "rl_gym_data").mkdir(parents=True, exist_ok=True)

    class _FakeService:
        def __init__(self):
            self.invalidated = False

        def invalidate(self):
            self.invalidated = True

    svc = _FakeService()
    promote_run("runXYZ", artifacts_dir, cascade_root, service_cache_bump=svc)
    assert svc.invalidated is True
