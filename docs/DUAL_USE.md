# Dual-use and intended use

## Intended use

CASCADE is **computational research software** for:

1. Mining Cas13-like effectors from metagenomic contigs with strict QC
2. Bootstrapping baselines with native CRISPR repeats and structure screens
3. Constrained *in silico* design of inter-domain linkers scored by predicted
   OFF vs ON catalytic geometry (and mismatch leak)

The reference trigger set in `data/fusion_targets.json` is **tumor-restricted
fusion RNA junctions** (e.g. BCR-ABL1, EWSR1-FLI1). That choice reflects an
oncology *hypothesis*, not a validated therapy.

## Explicit non-claims

- CASCADE does **not** produce a clinical product or therapeutic.
- CASCADE does **not** demonstrate cell killing, collateral RNase activity in
  lysates, or animal efficacy.
- Structure-prediction distances (e.g. His NE2–NE2 thresholds) are **design
  proxies**, not experimentally calibrated catalytic measurements.
- Outputs labeled “elite” or “optimal” are **computational rankings** only.

## Dual-use note

Cas13 collateral cleavage is widely studied for diagnostics and research.
CASCADE’s design loop (freeze recognition/catalysis, evolve mechanical
coupling, score two conformational states) could in principle be pointed at
other RNA triggers. Users are responsible for complying with institutional
biosafety rules, dual-use review, and applicable law. Do not use this software
to design systems intended to cause harm to people, animals, or crops outside
authorized, reviewed research programs.

We publish the pipeline to support **transparent methods**, **false-positive
mining audits**, and **reproducible in silico search**—not to release a
ready-to-deploy biological weapon or drug.
