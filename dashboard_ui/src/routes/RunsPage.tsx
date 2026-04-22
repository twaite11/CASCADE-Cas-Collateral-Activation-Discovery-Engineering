import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";

export function RunsPage() {
  return (
    <div className="space-y-4">
      <h1 className="text-2xl font-semibold tracking-tight">Runs</h1>
      <p className="text-sm text-muted-foreground">
        Live parallel evolution runs on Vast.ai. Launch, watch, and cancel here.
      </p>
      <Card>
        <CardHeader>
          <CardTitle>Coming soon</CardTitle>
        </CardHeader>
        <CardContent className="text-sm text-muted-foreground">
          The Runs tab wires up in the next commit. In the meantime,{" "}
          <code className="rounded bg-muted px-1">POST /api/runs</code> is
          already live from the backend.
        </CardContent>
      </Card>
    </div>
  );
}
