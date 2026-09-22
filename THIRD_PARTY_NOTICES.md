# Third-party software and weights

CASCADE (Apache-2.0) is an **orchestration** layer. It does **not** redistribute
foundation-model weights. Install generators and structure oracles yourself and
point CASCADE at them via environment variables (`EVAL_CMD`, `PXDESIGN_CMD`,
`PROTEINMPNN_DIR`, checkpoint paths).

| Component | Role in CASCADE | Upstream license (code) | Weights |
|:----------|:----------------|:------------------------|:--------|
| **[Protenix](https://github.com/bytedance/Protenix)** | Structure prediction oracle (Python) | Apache-2.0 | Download separately from ByteDance / Protenix releases; not shipped here |
| **[PXDesign](https://github.com/bytedance/PXDesign)** | Linker / binder generation (Gen 1) | Apache-2.0 | Download separately; not shipped here |
| **[Cattle-Prod](https://github.com/twaite11/cattle-prod)** | Optional Rust Protenix-compatible eval CLI | Apache-2.0 (fork of Protenix lineage) | Convert / place checkpoints yourself; not shipped here |
| **[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)** | Inverse folding on cached backbones (Gen 2+) | MIT (Baker lab / Dauparas et al.) | Use upstream weights; not shipped here |
| Biopython, NumPy, FastAPI, etc. | Orchestration / dashboard | See package metadata | N/A |

## What you may redistribute from *this* repo

- CASCADE source under Apache-2.0 (see `LICENSE` and `NOTICE`)
- Frozen mining tables and curated result bundles you build with
  `scripts/build_zenodo_bundle.py` (your own run outputs)

## What you must **not** assume is redistributable via CASCADE

- Protenix / PXDesign / Cattle-Prod / ProteinMPNN **checkpoint files**
- Any weights placed under `/workspace/models/` or similar local paths
- Docker images that bake third-party weights (build those yourself; do not
  push proprietary or restricted weights to public registries without a
  license check)

Always verify the current upstream LICENSE and model-use terms before
commercial redistribution of combined stacks.
