import { Card, CardHeader, CardTitle, CardContent } from "@/components/ui/card";

/**
 * Right-side collapsible "Optimized Switches" leaderboard. Stubbed here; full
 * implementation (3Dmol toggle, HEPN coloring, filters) arrives in the
 * frontend-sidebar-3d commit.
 */
export function OptimizedSidebar() {
  return (
    <aside className="hidden w-96 shrink-0 overflow-y-auto border-l bg-card/30 p-4 xl:block">
      <Card>
        <CardHeader>
          <CardTitle className="text-base">Optimized Switches</CardTitle>
        </CardHeader>
        <CardContent className="text-xs text-muted-foreground">
          Leaderboard with 3D toggle arrives in the next commit.
        </CardContent>
      </Card>
    </aside>
  );
}
