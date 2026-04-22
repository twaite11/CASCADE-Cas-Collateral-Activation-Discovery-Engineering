import { create } from "zustand";
import { Keyboard } from "lucide-react";

import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
} from "@/components/ui/dialog";

interface HelpState {
  open: boolean;
  set: (open: boolean) => void;
  toggle: () => void;
}

export const useHelp = create<HelpState>((set, get) => ({
  open: false,
  set: (open) => set({ open }),
  toggle: () => set({ open: !get().open }),
}));

const SHORTCUTS: Array<{ keys: string; description: string }> = [
  { keys: "?", description: "Toggle this help panel" },
  { keys: "Esc", description: "Close dialogs / drawers / help" },
  { keys: "g r", description: "Go to Runs" },
  { keys: "g b", description: "Go to Baselines" },
  { keys: "g v", description: "Go to Variants" },
  { keys: "g o", description: "Go to Optimized Switches" },
  { keys: "g p", description: "Go to Production" },
  { keys: "n", description: "From Baselines: launch new run with selection" },
  { keys: "/", description: "Focus the nearest search input" },
];

export function HelpPanel() {
  const open = useHelp((s) => s.open);
  const set = useHelp((s) => s.set);

  return (
    <Dialog open={open} onOpenChange={set}>
      <DialogContent>
        <DialogHeader>
          <DialogTitle className="flex items-center gap-2">
            <Keyboard className="h-5 w-5" />
            Keyboard shortcuts
          </DialogTitle>
          <DialogDescription>
            Speed up navigation around the CASCADE dashboard.
          </DialogDescription>
        </DialogHeader>
        <div className="space-y-1">
          {SHORTCUTS.map((s) => (
            <div
              key={s.keys}
              className="flex items-center justify-between rounded-md border-b border-border/40 py-2 last:border-0"
            >
              <span className="text-sm text-muted-foreground">
                {s.description}
              </span>
              <span className="font-mono text-xs">
                {s.keys.split(" ").map((k) => (
                  <kbd
                    key={k}
                    className="ml-1 rounded border bg-muted px-1.5 py-0.5"
                  >
                    {k}
                  </kbd>
                ))}
              </span>
            </div>
          ))}
        </div>
      </DialogContent>
    </Dialog>
  );
}
