import { BaselinesTable } from "@/components/BaselinesTable";

export function BaselinesPage() {
  return (
    <div className="space-y-6">
      <header>
        <h1 className="text-2xl font-semibold tracking-tight">Baselines</h1>
        <p className="text-sm text-muted-foreground">
          Pick enzyme + crRNA IDs, then launch parallel evolution runs on
          Vast.ai.
        </p>
      </header>
      <BaselinesTable />
    </div>
  );
}
