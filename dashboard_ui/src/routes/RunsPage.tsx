import { useMemo } from "react";
import { useQuery } from "@tanstack/react-query";
import { Link } from "react-router-dom";
import { Rocket, Inbox } from "lucide-react";

import { api } from "@/lib/api";
import type { Run, RunStatus } from "@/lib/types";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { RunCard } from "@/components/RunCard";

const ACTIVE: ReadonlySet<RunStatus> = new Set([
  "queued",
  "provisioning",
  "starting",
  "running",
  "syncing",
]);

export function RunsPage() {
  const { data, isLoading, error } = useQuery({
    queryKey: ["runs"],
    queryFn: () => api.listRuns({ limit: 200 }),
    refetchInterval: 4000,
  });

  const runs = data?.rows ?? [];
  const { active, finished } = useMemo(() => {
    const a: Run[] = [];
    const f: Run[] = [];
    for (const r of runs) {
      if (ACTIVE.has(r.status)) a.push(r);
      else f.push(r);
    }
    return { active: a, finished: f };
  }, [runs]);

  return (
    <div className="space-y-6">
      <header className="flex items-end justify-between">
        <div>
          <h1 className="text-2xl font-semibold tracking-tight">Runs</h1>
          <p className="text-sm text-muted-foreground">
            {isLoading
              ? "loading…"
              : `${active.length} active • ${finished.length} finished`}
          </p>
        </div>
        <Link to="/baselines">
          <Button>
            <Rocket className="mr-2 h-4 w-4" />
            New run
          </Button>
        </Link>
      </header>

      {error && (
        <Card>
          <CardContent className="py-4 text-sm text-red-400">
            Failed to load runs: {(error as Error).message}
          </CardContent>
        </Card>
      )}

      <section className="space-y-3">
        <h2 className="text-sm font-medium text-muted-foreground">
          Active runs
        </h2>
        {active.length === 0 ? (
          <EmptyState
            title="No active runs"
            subtitle="Launch one from the Baselines tab."
          />
        ) : (
          <div className="space-y-3">
            {active.map((r) => (
              <RunCard key={r.id} run={r} />
            ))}
          </div>
        )}
      </section>

      {finished.length > 0 && (
        <section className="space-y-3">
          <h2 className="text-sm font-medium text-muted-foreground">
            Recent runs
          </h2>
          <div className="space-y-3">
            {finished.slice(0, 25).map((r) => (
              <RunCard key={r.id} run={r} />
            ))}
          </div>
        </section>
      )}
    </div>
  );
}

function EmptyState({
  title,
  subtitle,
}: {
  title: string;
  subtitle: string;
}) {
  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2 text-base">
          <Inbox className="h-4 w-4" />
          {title}
        </CardTitle>
      </CardHeader>
      <CardContent className="text-sm text-muted-foreground">
        {subtitle}
      </CardContent>
    </Card>
  );
}
