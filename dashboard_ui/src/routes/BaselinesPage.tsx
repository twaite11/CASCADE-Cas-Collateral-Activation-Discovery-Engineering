import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";

export function BaselinesPage() {
  return (
    <div className="space-y-4">
      <h1 className="text-2xl font-semibold tracking-tight">Baselines</h1>
      <p className="text-sm text-muted-foreground">
        Pick enzyme + crRNA IDs to launch parallel evolution runs.
      </p>
      <Card>
        <CardHeader>
          <CardTitle>Coming soon</CardTitle>
        </CardHeader>
        <CardContent className="text-sm text-muted-foreground">
          The Baselines tab wires up in the next commit.
        </CardContent>
      </Card>
    </div>
  );
}
