#!/usr/bin/env python3
"""Build a Zenodo-ready results bundle (mining + curated GPU artifacts).

Does not upload. Creates dist/zenodo/CASCADE_results_<date>/ and a zip.
Upload the zip at https://zenodo.org/deposit/new (CC-BY-4.0 recommended for data).

Usage:
    python scripts/build_zenodo_bundle.py
    python scripts/build_zenodo_bundle.py --include-cif-samples
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _copy(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dest, dirs_exist_ok=True)
    else:
        shutil.copy2(src, dest)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--include-cif-samples",
        action="store_true",
        help="Include one sample_0.cif per confirmed Phase-1 baseline (larger zip)",
    )
    ap.add_argument(
        "--out-root",
        type=Path,
        default=ROOT / "dist" / "zenodo",
        help="Output directory for the bundle folder + zip",
    )
    args = ap.parse_args()

    stamp = date.today().isoformat()
    bundle = args.out_root / f"CASCADE_results_{stamp}"
    if bundle.exists():
        shutil.rmtree(bundle)
    bundle.mkdir(parents=True)

    # --- Always include: mining + docs + confirmed baselines ---
    copies: list[tuple[str, str]] = [
        ("outputs/mining_v3_campaign1", "mining/campaign1"),
        ("outputs/mining_v3_campaign2", "mining/campaign2"),
        ("outputs/validated_baseline_ids.txt", "baselines/validated_baseline_ids.txt"),
        ("data/mined_hits/confirmed_novel_cas13.fasta", "baselines/confirmed_novel_cas13.fasta"),
        ("data/mined_hits/confirmed_novel_cas13_metadata.csv", "baselines/confirmed_novel_cas13_metadata.csv"),
        ("docs/AUDIT_2026-05-13.md", "docs/AUDIT_2026-05-13.md"),
        ("docs/REMINE_2026-05-13.md", "docs/REMINE_2026-05-13.md"),
        ("docs/DUAL_USE.md", "docs/DUAL_USE.md"),
        ("THIRD_PARTY_NOTICES.md", "docs/THIRD_PARTY_NOTICES.md"),
        ("LICENSE", "LICENSE"),
        ("NOTICE", "NOTICE"),
    ]

    # --- GPU / evolution artifacts if present locally ---
    opt = ROOT / "outputs" / "optimized_switches"
    if opt.is_dir():
        for p in opt.glob("*.fasta"):
            copies.append((str(p.relative_to(ROOT)), f"gpu/optimized_switches/{p.name}"))

    rl = ROOT / "outputs" / "rl_gym_data"
    if rl.is_dir():
        for p in rl.rglob("*"):
            if p.is_file() and p.suffix in {".jsonl", ".json"}:
                rel = p.relative_to(ROOT)
                copies.append((str(rel), f"gpu/rl_gym_data/{p.relative_to(rl).as_posix()}"))

    # Phase-1 confidence summaries (small); optional one CIF per baseline
    phase1 = ROOT / "outputs" / "phase1_screening"
    if phase1.is_dir():
        for p in phase1.rglob("*summary_confidence*.json"):
            if ".ipynb_checkpoints" in str(p):
                continue
            rel = p.relative_to(phase1)
            copies.append(
                (str(p.relative_to(ROOT)), f"gpu/phase1_confidence/{rel.as_posix()}")
            )
        if args.include_cif_samples:
            for p in phase1.rglob("*_sample_0.cif"):
                if "positive_control" in str(p):
                    continue
                rel = p.relative_to(phase1)
                copies.append(
                    (str(p.relative_to(ROOT)), f"gpu/phase1_cif_samples/{rel.as_posix()}")
                )

    manifest_files: list[dict] = []
    missing: list[str] = []
    for src_rel, dest_rel in copies:
        src = ROOT / src_rel
        if not src.exists():
            missing.append(src_rel)
            continue
        dest = bundle / dest_rel
        _copy(src, dest)
        if src.is_file():
            manifest_files.append(
                {"path": dest_rel, "sha256": _sha256(dest), "bytes": dest.stat().st_size}
            )
        else:
            for f in dest.rglob("*"):
                if f.is_file():
                    rel = f.relative_to(bundle).as_posix()
                    manifest_files.append(
                        {"path": rel, "sha256": _sha256(f), "bytes": f.stat().st_size}
                    )

    meta = {
        "title": "CASCADE frozen mining and curated GPU results",
        "date": stamp,
        "software": "https://github.com/twaite11/CASCADE-Cas-Collateral-Activation-Discovery-Engineering",
        "license_suggestion": "CC-BY-4.0 for this data deposit; Apache-2.0 for software",
        "notes": [
            "Mining campaigns use mining_v3 (strict HEPN + anti-signatures + reciprocal check).",
            "GPU artifacts are computational rankings only — not wet-lab validated.",
            "Model weights are not included.",
        ],
        "missing_sources": missing,
        "files": sorted(manifest_files, key=lambda x: x["path"]),
    }
    (bundle / "MANIFEST.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    (bundle / "README.md").write_text(
        "\n".join(
            [
                "# CASCADE results bundle",
                "",
                f"Built: {stamp}",
                "",
                "Contents:",
                "- `mining/` — mining_v3 Campaign 1/2 summaries, hits, rejections",
                "- `baselines/` — three confirmed L. booriae Cas13a + validated IDs",
                "- `gpu/` — curated evolution / Phase-1 confidence artifacts (if present on build host)",
                "- `docs/` — audit, remine, dual-use, third-party notices",
                "- `MANIFEST.json` — SHA-256 checksums",
                "",
                "## Zenodo upload",
                "",
                "1. Create a new upload at https://zenodo.org/deposit/new",
                "2. Upload the sibling `.zip` of this folder",
                "3. Title: CASCADE: mining_v3 remine and curated in silico evolution artifacts",
                "4. License: Creative Commons Attribution 4.0 International",
                "5. Link related identifier to the GitHub repo (and JOSS DOI when available)",
                "6. Paste the Zenodo DOI into `CITATION.cff` and `paper/paper.md`",
                "",
            ]
        ),
        encoding="utf-8",
    )

    zip_path = args.out_root / f"CASCADE_results_{stamp}.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for f in bundle.rglob("*"):
            if f.is_file():
                zf.write(f, arcname=f"{bundle.name}/{f.relative_to(bundle).as_posix()}")

    print(f"Bundle: {bundle}")
    print(f"Zip:    {zip_path} ({zip_path.stat().st_size / 1e6:.1f} MB)")
    print(f"Files:  {len(manifest_files)}  missing sources: {len(missing)}")
    if missing:
        for m in missing:
            print(f"  missing: {m}")


if __name__ == "__main__":
    main()
