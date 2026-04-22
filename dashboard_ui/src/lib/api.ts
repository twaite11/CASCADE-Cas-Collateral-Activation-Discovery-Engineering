import type {
  Baseline,
  BaselinesResponse,
  OffersResponse,
  Run,
  RunDetail,
  RunsResponse,
} from "./types";

const API_BASE = ""; // same-origin: FastAPI serves the SPA under /dashboard

async function http<T>(path: string, init?: RequestInit): Promise<T> {
  const resp = await fetch(API_BASE + path, {
    headers: { "Content-Type": "application/json" },
    ...init,
  });
  if (!resp.ok) {
    const text = await resp.text().catch(() => "");
    throw new Error(`${resp.status} ${resp.statusText}: ${text}`);
  }
  if (resp.status === 204) return undefined as T;
  return resp.json() as Promise<T>;
}

export const api = {
  // ---------------------------- baselines
  listBaselines: (params: {
    search?: string;
    validated_only?: boolean;
    with_structure_only?: boolean;
    subtype?: string;
    limit?: number;
    offset?: number;
  }) => {
    const qs = new URLSearchParams();
    if (params.search) qs.set("search", params.search);
    if (params.validated_only) qs.set("validated_only", "true");
    if (params.with_structure_only) qs.set("with_structure_only", "true");
    if (params.subtype) qs.set("subtype", params.subtype);
    if (params.limit != null) qs.set("limit", String(params.limit));
    if (params.offset != null) qs.set("offset", String(params.offset));
    return http<BaselinesResponse>(`/api/baselines?${qs.toString()}`);
  },
  getBaseline: (id: string) => http<Baseline>(`/api/baselines/${id}`),

  // ---------------------------- vast offers
  listOffers: (params: {
    gpu_name?: string;
    num_gpus?: number;
    min_vram_gb?: number;
    max_dph_total?: number;
    limit?: number;
  }) => {
    const qs = new URLSearchParams();
    if (params.gpu_name) qs.set("gpu_name", params.gpu_name);
    if (params.num_gpus) qs.set("num_gpus", String(params.num_gpus));
    if (params.min_vram_gb) qs.set("min_vram_gb", String(params.min_vram_gb));
    if (params.max_dph_total) qs.set("max_dph_total", String(params.max_dph_total));
    if (params.limit) qs.set("limit", String(params.limit));
    return http<OffersResponse>(`/api/vast/offers?${qs.toString()}`);
  },

  // ---------------------------- runs
  launchRun: (body: {
    baseline_ids: string[];
    crrna_lookup_ids?: string[];
    offer_id: number;
    max_generations?: number;
    variants_per_gen?: number;
    workers?: number;
    label?: string;
  }) =>
    http<Run>("/api/runs", {
      method: "POST",
      body: JSON.stringify(body),
    }),
  listRuns: (params?: { status?: string[]; limit?: number; offset?: number }) => {
    const qs = new URLSearchParams();
    params?.status?.forEach((s) => qs.append("status", s));
    if (params?.limit) qs.set("limit", String(params.limit));
    if (params?.offset) qs.set("offset", String(params.offset));
    return http<RunsResponse>(`/api/runs?${qs.toString()}`);
  },
  getRun: (id: string, tail_lines = 500) =>
    http<RunDetail>(`/api/runs/${id}?tail_lines=${tail_lines}`),
  cancelRun: (id: string) =>
    http<{ ok: boolean; run_id: string }>(`/api/runs/${id}`, {
      method: "DELETE",
    }),

  // ---------------------------- existing dashboard endpoints
  overview: () => http<unknown>("/api/overview"),
  variants: () => http<unknown>("/api/variants"),
  optimized: () => http<unknown>("/api/optimized-switches"),
  production: () => http<unknown>("/api/production"),
  variantDetail: (id: string) => http<unknown>(`/api/variant/${id}`),
  compare: (a: string, b: string) =>
    http<unknown>(`/api/compare?a=${encodeURIComponent(a)}&b=${encodeURIComponent(b)}`),
};

// ---------------------------- WS helper
export function openLogStream(
  runId: string,
  onLine: (line: string) => void,
  onClose?: () => void,
): () => void {
  const proto = window.location.protocol === "https:" ? "wss:" : "ws:";
  const ws = new WebSocket(`${proto}//${window.location.host}/api/runs/${runId}/logs`);
  ws.onmessage = (ev) => onLine(String(ev.data));
  ws.onclose = () => onClose?.();
  ws.onerror = () => onClose?.();
  return () => {
    try {
      ws.close();
    } catch {
      /* noop */
    }
  };
}
