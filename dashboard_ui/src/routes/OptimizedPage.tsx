import { useMemo, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { Search, Sparkles } from "lucide-react";

import { api } from "@/lib/api";
import { Badge } from "@/components/ui/badge";
import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { useVariantDrawer } from "@/stores/variantDrawer";

/**
 * F-1 fix: full Optimized Switches view (previously a stub).  Shares its
 * data source with `OptimizedSidebar` but adds table layout, server-side
 * sort by fitness, and pagination so the operator can survey >25 hits.
 */
export function OptimizedPage() {
  const [search, setSearch] = useState("");
  const { data, isLoading, error } = useQuery({
    queryKey: ["optimized-switches", "page"],
    queryFn: () => api.optimized(2000),
    staleTime: 15_000,
  });
  const openDrawer = useVariantDrawer((s) => s.open);

  const rows = data?.rows ?? [];
  const filtered = useMemo(() => {
    if (!search) return rows;
    const q = search.toLowerCase();
    return rows.filter(
      (r) =>
        r.variant_id.toLowerCase().includes(q) ||
        (r.baseline_id ?? "").toLowerCase().includes(q),
    );
  }, [rows, search]);

  return (
    <div className="space-y-4">
      <header>
        <h1 className="flex items-center gap-2 text-2xl font-semibold tracking-tight">
          <Sparkles className="h-5 w-5 text-amber-400" />
          Optimized Switches
        </h1>
        <p className="text-sm text-muted-foreground">
          {isLoading
            ? "loading…"
            : `${filtered.length.toLocaleString()} / ${rows.length.toLocaleString()} matching`}
        </p>
      </header>

      <Card>
        <CardContent className="flex items-center gap-2 py-3">
          <Search className="h-4 w-4 text-muted-foreground" />
          <Input
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            placeholder="variant_id or baseline_id contains…"
            className="h-8 text-sm"
          />
        </CardContent>
      </Card>

      {error && (
        <Card>
          <CardContent className="py-4 text-sm text-red-400">
            Failed to load: {(error as Error).message}
          </CardContent>
        </Card>
      )}

      <Card>
        <CardContent className="p-0">
          <div className="grid grid-cols-[1fr_120px_70px_80px_80px_80px_80px_80px] gap-2 border-b bg-muted/30 px-3 py-2 text-[11px] font-medium uppercase text-muted-foreground">
            <div>variant</div>
            <div>baseline</div>
            <div className="text-right">gen</div>
            <div className="text-right">fitness</div>
            <div className="text-right">ipTM</div>
            <div className="text-right">AF2-IG</div>
            <div className="text-right">ON Å</div>
            <div className="text-right">OFF Å</div>
          </div>
          {filtered.map((r) => (
            <button
              key={r.variant_id}
              onClick={() => openDrawer(r.variant_id)}
              className="grid w-full grid-cols-[1fr_120px_70px_80px_80px_80px_80px_80px] gap-2 border-b px-3 py-2 text-left text-xs hover:bg-accent/30"
            >
              <span className="flex items-center gap-2 truncate font-mono">
                <span className="truncate">{r.variant_id}</span>
                {r.optimized_reasons?.includes("elite") && (
                  <Badge variant="default" className="text-[9px]">elite</Badge>
                )}
              </span>
              <span className="truncate text-muted-foreground">
                {r.baseline_id ?? "—"}
              </span>
              <span className="text-right font-mono">{r.generation}</span>
              <span className="text-right font-mono text-amber-300">
                {fmt(r.fitness, 3)}
              </span>
              <span className="text-right font-mono">{fmt(r.iptm, 3)}</span>
              <span className="text-right font-mono">{fmt(r.af2_ig, 3)}</span>
              <span className="text-right font-mono">{fmt(r.on_dist_A, 2)}</span>
              <span className="text-right font-mono">{fmt(r.off_dist_A, 2)}</span>
            </button>
          ))}
          {!isLoading && filtered.length === 0 && (
            <div className="px-3 py-6 text-center text-xs text-muted-foreground">
              No optimized switches match the current filter.
            </div>
          )}
        </CardContent>
      </Card>
    </div>
  );
}

function fmt(n: number | null | undefined, digits = 2): string {
  if (n === null || n === undefined || Number.isNaN(n)) return "—";
  return Number(n).toFixed(digits);
}
