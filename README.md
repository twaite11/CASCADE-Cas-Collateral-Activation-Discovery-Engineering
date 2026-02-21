CASCADE

Weaponizing Cas13 Collateral Cleavage: An AI-Driven Structural Pipeline for Engineered Suicide Switches in Targeted Oncology

The Purpose & Guiding Tenet

Traditional gene therapies utilizing Cas13 view its "collateral" trans-cleavage (the indiscriminate shredding of bystander RNA upon activation) as a catastrophic flaw that causes cellular toxicity. CASCADE reverses this paradigm. Our goal is to discover, validate, and engineer novel Cas13e-like proteins from metagenomic data to act as highly specific biological suicide switches. These engineered ribonucleoproteins (RNPs) will remain completely dormant in healthy cells. However, when the engineered Cas13 detects and binds a highly specific tumor fusion RNA (e.g., BCR-ABL, EWS-FLI1), it will trigger a massive conformational shift. This aligns the catalytic HEPN domains on the exterior of the protein, deliberately initiating collateral trans-cleavage to induce rapid, programmed cell death strictly within the tumor environment.

The Core Engineering Philosophy

Nature Provides the Targeting (Do Not Touch): The wild-type crRNA binding pocket (the Recognition lobe: NTD and Helical-1 domains) has been optimized by evolution to secure its specific CRISPR repeat. During our directed evolution and engineering phases, we strictly freeze the REC lobe. We utilize the native crRNA discovered from the metagenomic contigs to ensure maximum binding affinity without needing to reinvent the guide RNA mechanics.

Physics Provides the Switch (Engineer This): The thermodynamic transition between the inactive (OFF) and active (ON) states is controlled by the Inter-Domain Linkers (IDLs) and the Helical-2 domain. We engineer these specific regions to hyper-stabilize the dormant state, creating a high energetic barrier that ensures zero leakiness. Crucially, we optimize these domains to maximize trans-cleavage post-activation, ensuring the enzyme stays persistently locked in its highly lethal, active state only after the specific tumor fusion RNA physically forces the switch open.

Expected Inputs

The pipeline begins with raw deep-learning hits mined from Sequence Read Archive (SRA) contigs. It expects two highly specific files located in data/mined_hits/:

The Sequence Data: deep_hits_*.fasta

Contains the predicted Cas13e-like protein sequences (ORFs).

The Metagenomic Metadata: deep_hits_*_metadata.csv

Contains the sequence IDs, SRA accessions, and crucially, the overlapping k-mers representing the native CRISPR direct repeat domains associated with each Cas13 hit.

Pipeline Architecture

Phase 1: High-Throughput Validation & Database Anchoring

Ingestion: Parses the FASTA and CSV data into a memory-safe local SQLite database (cas13_variants.db), allowing the pipeline to handle hundreds of thousands of sequences without RAM overload.

HEPN Anchoring: Utilizes the highly conserved $R\phi X_3H$ sequence signature to map the exact 3D structural boundaries of the catalytic HEPN1 and HEPN2 domains.

True Cas Validation: Runs a rapid structural prediction using protenix-mini to filter out hits that fail to form a bilobed structure or fail to bind the native crRNA.

Phase 2: The Active Learning Directed Evolution Loop

Unlike traditional static pipelines, CASCADE utilizes an autonomous Reinforcement Learning (RL) Gym to continuously mutate, evaluate, and learn from structural predictions based on strict biophysics.

Generation & Constraint (pxdesign_wrapper.py): Generates structural variants by aggressively mutating the IDLs while mathematically freezing the spatial coordinates of the REC lobe based on the SQLite metadata.

The Fast Filter (protenix_eval.py & pdb_kinematics.py): Variants are subjected to rapid OFF/ON state predictions. The spatial math engine calculates the precise Euclidean distance between the $R\phi X_3H$ motifs. We select variants where the motifs are >25Å apart in the OFF state, but snap to <12Å apart in the ON state.

Full Ternary Oracle Scoring: Variants that exhibit massive shifts are promoted to the high-fidelity protenix_base model. They are rigorously scored (ipTM, AF2-IG) exclusively in their full Ternary complex (Protein + crRNA + Tumor Fusion RNA).

The Gym Feedback Matrix: Mutations from failed variants are penalized, while mutations from elite ternary structures are highly rewarded. The Gym exports these weights as a physical bias matrix (mpnn_bias_gen_X.json), continuously steering the PXDesign diffusion engine toward the optimal thermodynamic switch.

Project Structure

/workspace/CASCADE/
│
├── data/
│   └── mined_hits/               # Drop your FASTAs and CSVs here
│
├── databases/                    # Localized NT-RNA, Rfam, RNAcentral (IP Protection)
│
├── metadata/
│   ├── cas13_variants.db         # Local SQLite database
│   └── variant_domain_metadata.json # HEPN mapping boundaries
│
├── jsons/                        # Protenix/AlphaFold3 compatible input payloads
│
├── scripts/
│   ├── 01_parse_and_annotate.py  # Ingests hits, maps domains, builds DB & JSONs
│   ├── 02_run_screening.sh       # GPU-optimized Phase 1 Protenix-Mini filtering
│   ├── pxdesign_wrapper.py       # Invokes PXDesign diffusion generation
│   ├── evolution_orchestrator.py # Active Learning Gym (Master Controller)
│   └── utils/
│       ├── protenix_eval.py      # Generates ON/OFF payloads & runs inference
│       └── pdb_kinematics.py     # 3D Math Engine & Confidence JSON parser
│
└── outputs/
    ├── phase1_screening/         # Initial baseline PDBs
    ├── generation_queue/         # FASTA outputs from PXDesign
    ├── fast_eval/                # Mini-model OFF/ON structural checks
    ├── high_fidelity_scoring/    # Base-model Ternary predictions
    ├── rl_gym_data/              # Dynamic MPNN Bias matrices
    └── optimized_switches/       # Elite sequences & PDBs ready for synthesis
