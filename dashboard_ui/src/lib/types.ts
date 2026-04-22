export type Baseline = {
  baseline_id: string;
  subtype: string | null;
  sequence_length: number | null;
  hepn1_start: number | null;
  hepn1_end: number | null;
  hepn2_start: number | null;
  hepn2_end: number | null;
  crrna_repeat: string | null;
  crrna_spacer: string | null;
  crrna_lookup_id: string;
  has_phase1_structure: boolean;
  validated: boolean;
  sra_accession: string | null;
  score: number | null;
  status: string | null;
  reason: string | null;
  base_json: string | null;
};

export type BaselinesResponse = {
  total: number;
  offset: number;
  limit: number;
  rows: Baseline[];
};

export type VastOffer = {
  id: number;
  gpu_name: string | null;
  num_gpus: number | null;
  gpu_ram_gb: number | null;
  cpu_cores: number | null;
  cpu_ram_gb: number | null;
  disk_space_gb: number | null;
  dph_total: number | null;
  inet_down_mbps: number | null;
  inet_up_mbps: number | null;
  datacenter: string | null;
  reliability: number | null;
};

export type OffersResponse = {
  query: string;
  offers: VastOffer[];
};

export type RunStatus =
  | "queued"
  | "provisioning"
  | "starting"
  | "running"
  | "syncing"
  | "completed"
  | "failed"
  | "cancelled";

export type Run = {
  id: string;
  label: string;
  status: RunStatus;
  created_at: number;
  started_at: number | null;
  finished_at: number | null;
  cost_usd: number | null;
  dph_usd: number | null;
  baseline_ids: string[];
  max_generations: number;
  variants_per_gen: number;
  instance_id: number | null;
  ssh_host: string | null;
  ssh_port: number | null;
  error: string | null;
};

export type RunsResponse = {
  total: number;
  rows: Run[];
};

export type RunDetail = {
  run: Run;
  log_tail: string[];
};

export type OptimizedArtifacts = {
  fasta?: string;
  structure?: string;
  crrna?: string;
};

export type EvalArtifacts = {
  on_structure?: string;
  off_structure?: string;
  offtarget_structure?: string;
  on_summary?: string;
  off_summary?: string;
  offtarget_summary?: string;
};

export type DomainMetadata = {
  hepn1_start?: number;
  hepn1_end?: number;
  hepn2_start?: number;
  hepn2_end?: number;
  [key: string]: unknown;
};

export type OptimizedVariant = {
  variant_id: string;
  baseline_id: string | null;
  generation: number;
  fitness: number;
  iptm?: number | null;
  af2_ig?: number | null;
  off_dist_A?: number | null;
  on_dist_A?: number | null;
  hepn_shift_A?: number | null;
  optimized_reasons?: string[];
  optimized_artifacts?: OptimizedArtifacts;
  eval_artifacts?: EvalArtifacts;
  domain_metadata?: DomainMetadata;
};

export type OptimizedResponse = {
  total: number;
  rows: OptimizedVariant[];
};
