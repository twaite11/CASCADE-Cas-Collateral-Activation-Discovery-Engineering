# CASCADE Tests

Unit tests for pipeline components. **No GPU, Protenix, or PXDesign required.**

## Setup

```bash
pip install -r requirements-dev.txt
# or: pip install pytest
```

## Run

From project root:

```bash
pytest tests/ -v
```

## What's Tested

| Module | Tests |
|--------|-------|
| **pdb_kinematics** | HEPN distance calculation, Protenix score extraction |
| **protenix_eval** | Mismatch sequences, OFF/ON JSON generation, off-target JSON |
| **01_parse_and_annotate** | DB init, FASTA/CSV load, HEPN detection, Protenix JSON output |
| **evolution_orchestrator** | Fitness, EvolutionGym, mutation extraction, HEPN indices |
| **03_pxdesign_wrapper** | Freeze config generation |
| **run_protenix (mocked)** | CLI command construction |

## Fixtures

- `tmpdir`: Temporary directory
- `sample_fasta`, `sample_csv`: FASTA + metadata with HEPN motifs
- `sample_metadata_json`: Variant domain metadata
- `minimal_pdb`: PDB with residues 10 and 50, 30Å apart
- `protenix_summary_json`: Mock Protenix output
- `variant_fasta`, `baseline_json`: For mutation and evaluation tests
