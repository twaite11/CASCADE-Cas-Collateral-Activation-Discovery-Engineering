import { useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { Search, Filter as FilterIcon } from "lucide-react";

import { withAuth } from "@/lib/api";
import { cn } from "@/lib/utils";

import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { useVariantDrawer } from "@/stores/variantDrawer";

interface VariantsResponse {
  total: number;
  offset: number;
  limit: number;
  rows: Array<{
    variant_id: string;
    baseline_id?: string | null;
    generation?: number | null;
    fitness?: number | null;
    iptm?: number | null;
    af2_ig?: number | null;
    on_dist_A?: number | null;
    off_dist_A?: number | null;
    is_elite?: boolean;
    optimized_switch?: boolean;
  }>;
}

const PAGE_SIZE = 200;

/**
 * F-1 fix: VariantsPage used to be a "coming soon" stub.  This is the real
 * implementation, mirroring the patterns of OptimizedSidebar / BaselinesPage:
 *   * server-side pagination (limit + offset),
 *   * client-side text filter on variant_id,
 *   * filter chips for `elite` / `optimized_switch`,
 *   * click-row to open the VariantDetailDrawer (shared global store).
 */
export function VariantsPage() {
  const [offset, setOffset] = useState(0);
  const [search, setSearch] = useState("");
  const [eliteOnly, setEliteOnly] = useState(false);
  const [optimizedOnly, setOptimizedOnly] = useState(false);

  const q = useQuery<VariantsResponse, Error>({
    queryKey: ["variants", { offset, search, eliteOnly, optimizedOnly }],
    queryFn: async ({ signal }) => {
      const qs = new URLSearchParams();
      qs.set("limit", String(PAGE_SIZE));
      qs.set("offset", String(offset));
      if (search) qs.set("search", search);
      if (eliteOnly) qs.set("elite_only", "true");
      if (optimizedOnly) qs.set("optimized_only", "true");
      const resp = await fetch(`/api/variants?${qs}`, withAuth({ signal }));
      if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText}`);
      return (await resp.json()) as VariantsResponse;
    },
    placeholderData: (prev) => prev,
    staleTime: 10_000,
  });

  const rows = q.data?.rows ?? [];
  const total = q.data?.total ?? 0;
  const openDrawer = useVariantDrawer((s) => s.open);

  const pageEnd = Math.min(offset + PAGE_SIZE, total);

  return (
    <div className="space-y-4">
      <header>
        <h1 className="text-2xl font-semibold tracking-tight">Variants</h1>
        <p className="text-sm text-muted-foreground">
          {q.isLoading
            ? "loading…"
            : `${rows.length.toLocaleString()} shown · ${total.toLocaleString()} matching`}
        </p>
      </header>

      <Card>
        <CardContent className="flex flex-wrap items-end gap-3 py-3">
          <div className="flex flex-1 items-center gap-2 min-w-[240px]">
            <Search className="h-4 w-4 text-muted-foreground" />
            <Input
              value={search}
              onChange={(e) => {
                setOffset(0);
                setSearch(e.target.value);
              }}
              placeholder="Variant ID contains…"
              className="h-8 text-sm"
            />
          </div>
          <Button
            size="sm"
            variant={eliteOnly ? "default" : "outline"}
            onClick={() => {
              setOffset(0);
              setEliteOnly((v) => !v);
            }}
          >
            <FilterIcon className="mr-1 h-3.5 w-3.5" />
            elite
          </Button>
          <Button
            size="sm"
            variant={optimizedOnly ? "default" : "outline"}
            onClick={() => {
              setOffset(0);
              setOptimizedOnly((v) => !v);
            }}
          >
            <FilterIcon className="mr-1 h-3.5 w-3.5" />
            optimized
          </Button>
        </CardContent>
      </Card>

      {q.error && (
        <Card>
          <CardContent className="py-4 text-sm text-red-400">
            Failed to load variants: {(q.error as Error).message}
          </CardContent>
        </Card>
      )}

      <Card>
        <CardContent className="p-0">
          <div className="grid grid-cols-[1fr_120px_80px_80px_80px_80px_80px] gap-2 border-b bg-muted/30 px-3 py-2 text-[11px] font-medium uppercase text-muted-foreground">
            <div>variant</div>
            <div>baseline</div>
            <div className="text-right">gen</div>
            <div className="text-right">fitness</div>
            <div className="text-right">ipTM</div>
            <div className="text-right">AF2-IG</div>
            <div className="text-right">flags</div>
          </div>
          {rows.map((r) => (
            <button
              key={r.variant_id}
              onClick={() => openDrawer(r.variant_id)}
              className="grid w-full grid-cols-[1fr_120px_80px_80px_80px_80px_80px] gap-2 border-b px-3 py-2 text-left text-xs hover:bg-accent/30"
            >
              <span className="truncate font-mono">{r.variant_id}</span>
              <span className="truncate text-muted-foreground">
                {r.baseline_id ?? "—"}
              </span>
              <span className="text-right font-mono">{r.generation ?? "—"}</span>
              <span
                className={cn(
                  "text-right font-mono",
                  r.is_elite && "text-amber-300",
                )}
              >
                {fmt(r.fitness, 3)}
              </span>
              <span className="text-right font-mono">{fmt(r.iptm, 3)}</span>
              <span className="text-right font-mono">{fmt(r.af2_ig, 3)}</span>
              <span className="flex justify-end gap-1">
                {r.is_elite && <Badge variant="default" className="text-[9px]">elite</Badge>}
                {r.optimized_switch && (
                  <Badge variant="outline" className="text-[9px]">opt</Badge>
                )}
              </span>
            </button>
          ))}
          {!q.isLoading && rows.length === 0 && (
            <div className="px-3 py-6 text-center text-xs text-muted-foreground">
              No variants match the current filter.
            </div>
          )}
        </CardContent>
      </Card>

      <footer className="flex items-center justify-between text-xs text-muted-foreground">
        <span>
          {total === 0
            ? "0 results"
            : `showing ${offset + 1}–${pageEnd} of ${total.toLocaleString()}`}
        </span>
        <div className="flex gap-2">
          <Button
            size="sm"
            variant="outline"
            disabled={offset === 0 || q.isLoading}
            onClick={() => setOffset(Math.max(0, offset - PAGE_SIZE))}
          >
            Prev
          </Button>
          <Button
            size="sm"
            variant="outline"
            disabled={pageEnd >= total || q.isLoading}
            onClick={() => setOffset(offset + PAGE_SIZE)}
          >
            Next
          </Button>
        </div>
      </footer>
    </div>
  );
}

function fmt(n: number | null | undefined, digits = 2): string {
  if (n === null || n === undefined || Number.isNaN(n)) return "—";
  return Number(n).toFixed(digits);
}
