import { useQuery } from "@tanstack/react-query";
import { Factory, TrendingUp } from "lucide-react";

import { withAuth } from "@/lib/api";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";

interface ProductionSummary {
  totals: {
    generated: number;
    passed_filter: number;
    optimized: number;
    elite: number;
  };
  funnel: {
    generated_to_passed: number;
    passed_to_optimized: number;
    optimized_to_elite: number;
  };
  lineages: Array<{
    lineage_id: string;
    generated: number;
    optimized: number;
    elite: number;
    optimized_yield: number;
    elite_yield: number;
  }>;
}

/**
 * F-1 fix: real Production view (replaces the "coming soon" stub). Reads
 * `/api/production`, which the backend already computes.  Displays the
 * lineage yield funnel + a leaderboard of top lineages by optimized
 * yield.
 */
export function ProductionPage() {
  const { data, isLoading, error } = useQuery({
    queryKey: ["production-summary"],
    queryFn: async ({ signal }) => {
      const resp = await fetch("/api/production", withAuth({ signal }));
      if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText}`);
      return (await resp.json()) as ProductionSummary;
    },
    refetchInterval: 30_000,
  });

  return (
    <div className="space-y-4">
      <header>
        <h1 className="flex items-center gap-2 text-2xl font-semibold tracking-tight">
          <Factory className="h-5 w-5" />
          Production
        </h1>
        <p className="text-sm text-muted-foreground">
          Funnel + per-lineage yield across the full RL evolution dataset.
        </p>
      </header>

      {error && (
        <Card>
          <CardContent className="py-4 text-sm text-red-400">
            Failed to load: {(error as Error).message}
          </CardContent>
        </Card>
      )}

      {isLoading && (
        <Card>
          <CardContent className="py-4 text-sm text-muted-foreground">
            loading…
          </CardContent>
        </Card>
      )}

      {data && (
        <>
          <div className="grid grid-cols-2 gap-3 md:grid-cols-4">
            <Stat
              label="Generated"
              value={data.totals.generated.toLocaleString()}
            />
            <Stat
              label="Passed filter"
              value={data.totals.passed_filter.toLocaleString()}
              sub={`${pct(data.funnel.generated_to_passed)} of generated`}
            />
            <Stat
              label="Optimized"
              value={data.totals.optimized.toLocaleString()}
              sub={`${pct(data.funnel.passed_to_optimized)} of passed`}
            />
            <Stat
              label="Elite"
              value={data.totals.elite.toLocaleString()}
              sub={`${pct(data.funnel.optimized_to_elite)} of optimized`}
              accent="text-amber-300"
            />
          </div>

          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="flex items-center gap-2 text-sm">
                <TrendingUp className="h-4 w-4" />
                Lineages (top 50 by optimized yield)
              </CardTitle>
            </CardHeader>
            <CardContent className="p-0">
              <div className="grid grid-cols-[1fr_80px_80px_80px_80px_80px] gap-2 border-b bg-muted/30 px-3 py-2 text-[11px] font-medium uppercase text-muted-foreground">
                <div>lineage</div>
                <div className="text-right">gen</div>
                <div className="text-right">opt</div>
                <div className="text-right">elite</div>
                <div className="text-right">opt yield</div>
                <div className="text-right">elite yield</div>
              </div>
              {data.lineages.map((l) => (
                <div
                  key={l.lineage_id || "(no-lineage)"}
                  className="grid grid-cols-[1fr_80px_80px_80px_80px_80px] gap-2 border-b px-3 py-2 text-xs"
                >
                  <span className="truncate font-mono">
                    {l.lineage_id || "—"}
                  </span>
                  <span className="text-right font-mono">{l.generated}</span>
                  <span className="text-right font-mono">{l.optimized}</span>
                  <span className="text-right font-mono">{l.elite}</span>
                  <span className="text-right font-mono">
                    {pct(l.optimized_yield)}
                  </span>
                  <span className="text-right font-mono text-amber-300">
                    {pct(l.elite_yield)}
                  </span>
                </div>
              ))}
            </CardContent>
          </Card>
        </>
      )}
    </div>
  );
}

function Stat({
  label,
  value,
  sub,
  accent,
}: {
  label: string;
  value: string;
  sub?: string;
  accent?: string;
}) {
  return (
    <Card>
      <CardContent className="py-3">
        <div className="text-[11px] uppercase tracking-wide text-muted-foreground">
          {label}
        </div>
        <div className={`text-2xl font-semibold ${accent ?? ""}`}>{value}</div>
        {sub && (
          <div className="text-[10px] text-muted-foreground">{sub}</div>
        )}
      </CardContent>
    </Card>
  );
}

function pct(n: number): string {
  if (!Number.isFinite(n)) return "—";
  return `${(n * 100).toFixed(1)}%`;
}
