"""Finish steps 5/6/9/10 of the quarantine that were skipped/missed.

- Move per-ORF subdirectories under outputs/phase1_screening (numeric Campaign-1 IDs).
- Move per-baseline files/subdirs under outputs/optimized_switches that are not from
  the three confirmed Cas13a lineages.
- Quarantine metadata/variant_domain_metadata.json (it indexes the wiped catalog).
- Move temporary audit scripts to quarantine for traceability.
"""
from __future__ import annotations

import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAMP = "20260513"
QDIR = ROOT / "quarantine" / f"non_cas13_audit_{STAMP}"
QDIR.mkdir(parents=True, exist_ok=True)

CONFIRMED_PREFIXES = (
    "NZ_JAASWF010000012.1_ORF_f1_6252",
    "NZ_JAAROR010000001.1_ORF_f1_65475",
    "NZ_JAARYF010000017.1_ORF_f1_4771",
)


def quarantine(src: Path, reason: str) -> None:
    if not src.exists():
        print(f"  SKIP (not present): {src}")
        return
    rel = src.relative_to(ROOT)
    dst = QDIR / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    shutil.move(str(src), str(dst))
    kind = "dir " if src.is_dir() else "file"
    print(f"  MOVE {kind}: {rel}  [{reason}]")


def main() -> None:
    print(f"# CASCADE quarantine finish -- {datetime.now(timezone.utc).isoformat()}")
    print(f"# Quarantine dir: {QDIR.relative_to(ROOT)}")

    # Step 5: phase1_screening per-ORF subdirectories (mostly numeric Campaign 1 IDs)
    print("\n## phase1_screening per-ORF subdirs")
    p1 = ROOT / "outputs" / "phase1_screening"
    if p1.exists():
        for entry in sorted(p1.iterdir()):
            keep = any(entry.name.startswith(p) for p in CONFIRMED_PREFIXES)
            if not keep:
                quarantine(entry, "non-confirmed-Cas13a phase1 screening output")

    # Step 6: optimized_switches per-baseline subdirs/files
    print("\n## optimized_switches non-confirmed lineages")
    opt = ROOT / "outputs" / "optimized_switches"
    if opt.exists():
        for entry in sorted(opt.iterdir()):
            if entry.name == ".ipynb_checkpoints":
                quarantine(entry, "ipynb checkpoints noise")
                continue
            keep = any(entry.name.startswith(p) for p in CONFIRMED_PREFIXES)
            if not keep:
                quarantine(entry, "evolved variant from unverified lineage")

    # Step 9: variant_domain_metadata.json
    print("\n## variant_domain_metadata.json")
    vdm = ROOT / "metadata" / "variant_domain_metadata.json"
    if vdm.exists():
        quarantine(vdm, "indexes wiped Campaign 1 catalog; will be regenerated")

    # Step 10: temporary audit scripts (kept for provenance)
    print("\n## temporary audit scripts -> quarantine")
    for scratch in (
        "verify_candidates.py",
        "inspect_db.py",
        "audit_variants_db.py",
        "audit_variants_classification.csv",
    ):
        p = ROOT / scratch
        if p.exists():
            quarantine(p, "temporary audit script (kept for traceability)")

    print("\nDone.")


if __name__ == "__main__":
    main()
