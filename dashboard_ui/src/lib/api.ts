import type {
  Baseline,
  BaselinesResponse,
  OffersResponse,
  OptimizedResponse,
  Run,
  RunDetail,
  RunsResponse,
} from "./types";

const API_BASE = ""; // same-origin: FastAPI serves the SPA under /dashboard

/**
 * Phase-2/C-2 client-side API-key support.
 *
 * The backend now requires `X-Cascade-Api-Key` (or `?api_key=`) on every
 * mutating endpoint and on the WebSocket.  We surface the same key to
 * read-only endpoints too so a wrong key is detected immediately on first
 * page load (better failure mode than getting 200s on /api/variants and
 * 401s later when the user clicks "New run").
 *
 * Storage order:
 *   1. `window.CASCADE_API_KEY` injected by the controller's index.html
 *      (preferred for prod deploys -- never lives in localStorage).
 *   2. `localStorage["cascade_api_key"]` for developer convenience.
 *
 * If neither is set we send no header at all, which the backend interprets
 * as "anonymous" and only allows on read-only routes.
 */
declare global {
  interface Window {
    CASCADE_API_KEY?: string;
  }
}

function getApiKey(): string | null {
  if (typeof window === "undefined") return null;
  const fromWindow = window.CASCADE_API_KEY?.trim();
  if (fromWindow) return fromWindow;
  try {
    const fromLS = window.localStorage?.getItem("cascade_api_key")?.trim();
    if (fromLS) return fromLS;
  } catch {
    /* localStorage may be blocked in some sandboxes */
  }
  return null;
}

export function setApiKey(key: string | null): void {
  if (typeof window === "undefined") return;
  try {
    if (key) window.localStorage.setItem("cascade_api_key", key);
    else window.localStorage.removeItem("cascade_api_key");
  } catch {
    /* noop */
  }
}

/** Merge the API key into a RequestInit.  Exported so non-`api.ts` callers
 *  like `Structure3D` can also use it for their own fetches. */
export function withAuth(init: RequestInit = {}): RequestInit {
  const headers = new Headers(init.headers ?? {});
  const key = getApiKey();
  if (key) headers.set("X-Cascade-Api-Key", key);
  return { ...init, headers };
}

async function http<T>(path: string, init?: RequestInit): Promise<T> {
  const headers = new Headers(init?.headers ?? {});
  if (!headers.has("Content-Type")) headers.set("Content-Type", "application/json");
  const key = getApiKey();
  if (key) headers.set("X-Cascade-Api-Key", key);

  const resp = await fetch(API_BASE + path, { ...init, headers });
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
  optimized: (limit = 500) =>
    http<OptimizedResponse>(`/api/optimized-switches?limit=${limit}`),
  production: () => http<unknown>("/api/production"),
  variantDetail: (id: string, signal?: AbortSignal) =>
    http<Record<string, unknown>>(
      `/api/variant/${encodeURIComponent(id)}`,
      { signal },
    ),
  compare: (left: string, right: string) =>
    http<unknown>(
      `/api/compare?left=${encodeURIComponent(left)}&right=${encodeURIComponent(right)}`,
    ),
};

export function structureFileUrl(relPath: string): string {
  return `/api/structure-file?path=${encodeURIComponent(relPath)}`;
}

// ---------------------------- WS helper
export function openLogStream(
  runId: string,
  onLine: (line: string) => void,
  onClose?: (reason?: string) => void,
): () => void {
  const proto = window.location.protocol === "https:" ? "wss:" : "ws:";
  // C-2: pass the API key as a query param -- browsers can't set custom
  // headers on a WS handshake, so the backend accepts either path.
  const key = getApiKey();
  const auth = key ? `?api_key=${encodeURIComponent(key)}` : "";
  const ws = new WebSocket(
    `${proto}//${window.location.host}/api/runs/${runId}/logs${auth}`,
  );
  ws.onmessage = (ev) => onLine(String(ev.data));
  ws.onclose = (ev) => onClose?.(ev.reason || `closed (code ${ev.code})`);
  ws.onerror = () => onClose?.("websocket error");
  return () => {
    try {
      ws.close();
    } catch {
      /* noop */
    }
  };
}
