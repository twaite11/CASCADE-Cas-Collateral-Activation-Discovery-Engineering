import { create } from "zustand";
import { X, CheckCircle2, AlertTriangle, AlertCircle, Info } from "lucide-react";
import { cn } from "@/lib/utils";

export type ToastKind = "success" | "error" | "warn" | "info";

export interface ToastItem {
  id: string;
  kind: ToastKind;
  title: string;
  description?: string;
}

interface ToastState {
  items: ToastItem[];
  push: (t: Omit<ToastItem, "id">) => void;
  dismiss: (id: string) => void;
}

export const useToasts = create<ToastState>((set) => ({
  items: [],
  push: (t) => {
    const id = crypto.randomUUID();
    set((s) => ({ items: [...s.items, { id, ...t }] }));
    setTimeout(() => {
      set((s) => ({ items: s.items.filter((i) => i.id !== id) }));
    }, 6000);
  },
  dismiss: (id) =>
    set((s) => ({ items: s.items.filter((i) => i.id !== id) })),
}));

export function toast(item: Omit<ToastItem, "id">) {
  useToasts.getState().push(item);
}

const ICONS: Record<ToastKind, typeof CheckCircle2> = {
  success: CheckCircle2,
  error: AlertCircle,
  warn: AlertTriangle,
  info: Info,
};

const COLORS: Record<ToastKind, string> = {
  success: "border-emerald-600/40 bg-emerald-950/80 text-emerald-100",
  error: "border-red-600/40 bg-red-950/80 text-red-100",
  warn: "border-amber-500/40 bg-amber-950/80 text-amber-100",
  info: "border-sky-600/40 bg-sky-950/80 text-sky-100",
};

export function Toaster() {
  const items = useToasts((s) => s.items);
  const dismiss = useToasts((s) => s.dismiss);
  return (
    <div className="pointer-events-none fixed bottom-4 right-4 z-[100] flex w-80 flex-col gap-2">
      {items.map((t) => {
        const Icon = ICONS[t.kind];
        return (
          <div
            key={t.id}
            className={cn(
              "pointer-events-auto rounded-md border px-4 py-3 shadow-lg backdrop-blur",
              COLORS[t.kind],
            )}
          >
            <div className="flex items-start gap-2">
              <Icon className="mt-0.5 h-4 w-4 shrink-0" />
              <div className="flex-1 text-sm">
                <div className="font-medium">{t.title}</div>
                {t.description && (
                  <div className="mt-0.5 text-xs opacity-80">
                    {t.description}
                  </div>
                )}
              </div>
              <button
                onClick={() => dismiss(t.id)}
                className="text-current opacity-60 hover:opacity-100"
                aria-label="Dismiss"
              >
                <X className="h-4 w-4" />
              </button>
            </div>
          </div>
        );
      })}
    </div>
  );
}
