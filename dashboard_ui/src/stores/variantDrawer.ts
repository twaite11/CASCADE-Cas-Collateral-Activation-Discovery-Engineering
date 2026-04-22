import { create } from "zustand";

interface State {
  variantId: string | null;
  compareWith: string | null;
  open: (variantId: string, compareWith?: string | null) => void;
  close: () => void;
  setCompareWith: (id: string | null) => void;
}

export const useVariantDrawer = create<State>((set) => ({
  variantId: null,
  compareWith: null,
  open: (variantId, compareWith = null) => set({ variantId, compareWith }),
  close: () => set({ variantId: null, compareWith: null }),
  setCompareWith: (compareWith) => set({ compareWith }),
}));
