import { useEffect, useState, useRef } from "react";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import {
  CircleDot,
  Square,
  Timer,
  DollarSign,
  ChevronDown,
  ChevronRight,
  Dna,
} from "lucide-react";

import { api } from "@/lib/api";
import type { Run, RunStatus } from "@/lib/types";
import { cn, classifyStatus, formatDph, formatRelative } from "@/lib/utils";
import { toast } from "@/components/ui/toast";

import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { LogTerminal } from "@/components/LogTerminal";

interface Props {
  run: Run;
}

const ACTIVE: Set<RunStatus> = new Set([
  "queued",
  "provisioning",
  "starting",
  "running",
  "syncing",
]);

export function RunCard({ run }: Props) {
  const qc = useQueryClient();
  const [expanded, setExpanded] = useState(ACTIVE.has(run.status));
  const lastStatus = useRef<RunStatus>(run.status);

  const detailQ = useQuery({
    enabled: expanded,
    queryKey: ["run-detail", run.id],
    queryFn: () => api.getRun(run.id, 500),
    refetchInterval: ACTIVE.has(run.status) ? 5000 : false,
  });

  useEffect(() => {
    if (lastStatus.current !== run.status) {
      if (run.status === "completed") {
        toast({
          kind: "success",
          title: `Run ${run.label || run.id.slice(0, 8)} completed`,
          description: `${run.baseline_ids.length} baselines, cost ~$${(run.cost_usd ?? 0).toFixed(2)}`,
        });
      } else if (run.status === "failed") {
        toast({
          kind: "error",
          title: `Run ${run.label || run.id.slice(0, 8)} failed`,
          description: run.error ?? "see logs for details",
        });
      } else if (run.status === "cancelled") {
        toast({
          kind: "warn",
          title: `Run ${run.label || run.id.slice(0, 8)} cancelled`,
        });
      }
      lastStatus.current = run.status;
    }
  }, [run.status, run.baseline_ids.length, run.cost_usd, run.error, run.id, run.label]);

  const active = ACTIVE.has(run.status);
  const statusVariant = classifyStatus(run.status);
  const elapsedSec =
    (run.finished_at ?? Date.now() / 1000) - (run.started_at ?? run.created_at);
  const elapsed = formatDuration(elapsedSec);

  async function handleCancel() {
    if (!confirm(`Cancel run ${run.id.slice(0, 8)} and destroy its instance?`)) {
      return;
    }
    try {
      await api.cancelRun(run.id);
      toast({ kind: "info", title: "Cancel requested" });
      qc.invalidateQueries({ queryKey: ["runs"] });
    } catch (e: any) {
      toast({
        kind: "error",
        title: "Cancel failed",
        description: String(e?.message ?? e),
      });
    }
  }

  return (
    <Card className={cn(active && "ring-1 ring-primary/30")}>
      <CardHeader className="flex-row items-center justify-between space-y-0 pb-3">
        <button
          onClick={() => setExpanded((v) => !v)}
          className="flex flex-1 items-center gap-2 text-left"
        >
          {expanded ? (
            <ChevronDown className="h-4 w-4 text-muted-foreground" />
          ) : (
            <ChevronRight className="h-4 w-4 text-muted-foreground" />
          )}
          <CardTitle className="flex items-center gap-2 text-base">
            <CircleDot
              className={cn(
                "h-3 w-3",
                active ? "animate-pulse text-emerald-400" : "text-muted-foreground",
              )}
            />
            {run.label || `run ${run.id.slice(0, 8)}`}
          </CardTitle>
          <Badge variant={statusVariant}>{run.status}</Badge>
        </button>

        <div className="flex items-center gap-3 text-xs text-muted-foreground">
          <span className="inline-flex items-center gap-1">
            <Dna className="h-3 w-3" /> {run.baseline_ids.length}
          </span>
          <span className="inline-flex items-center gap-1">
            <Timer className="h-3 w-3" /> {elapsed}
          </span>
          <span className="inline-flex items-center gap-1">
            <DollarSign className="h-3 w-3" />
            {run.cost_usd != null ? `$${run.cost_usd.toFixed(2)}` : "—"}
            {run.dph_usd != null && (
              <span className="opacity-60">({formatDph(run.dph_usd)})</span>
            )}
          </span>
          {active && (
            <Button size="sm" variant="destructive" onClick={handleCancel}>
              <Square className="mr-1 h-3 w-3" /> cancel
            </Button>
          )}
        </div>
      </CardHeader>

      {expanded && (
        <CardContent className="space-y-3">
          <div className="flex flex-wrap gap-x-6 gap-y-1 text-xs text-muted-foreground">
            <span>
              <b>id</b> <code className="font-mono">{run.id}</code>
            </span>
            <span>
              <b>instance</b> {run.instance_id ?? "—"}
            </span>
            {run.ssh_host && (
              <span>
                <b>ssh</b>{" "}
                <code className="font-mono">
                  {run.ssh_host}:{run.ssh_port}
                </code>
              </span>
            )}
            <span>
              <b>gen/var</b> {run.max_generations}×{run.variants_per_gen}
            </span>
            <span>
              <b>started</b> {formatRelative(run.started_at ?? run.created_at)}
            </span>
          </div>

          <div className="flex flex-wrap gap-1">
            {run.baseline_ids.map((id) => (
              <Badge
                variant="outline"
                key={id}
                className="font-mono text-[11px]"
              >
                {id}
              </Badge>
            ))}
          </div>

          <LogTerminal
            runId={run.id}
            initialLines={detailQ.data?.log_tail ?? []}
          />
        </CardContent>
      )}
    </Card>
  );
}

function formatDuration(seconds: number): string {
  if (!isFinite(seconds) || seconds < 0) return "—";
  const h = Math.floor(seconds / 3600);
  const m = Math.floor((seconds % 3600) / 60);
  const s = Math.floor(seconds % 60);
  if (h > 0) return `${h}h ${m}m`;
  if (m > 0) return `${m}m ${s}s`;
  return `${s}s`;
}
