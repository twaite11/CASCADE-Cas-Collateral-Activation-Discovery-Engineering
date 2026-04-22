import { Link, NavLink, Route, Routes, Navigate } from "react-router-dom";
import { FlaskConical, Play, Dna, Sparkles, Factory } from "lucide-react";

import { cn } from "@/lib/utils";
import { Toaster } from "@/components/ui/toast";
import { OptimizedSidebar } from "@/components/OptimizedSidebar";
import { RunsPage } from "@/routes/RunsPage";
import { BaselinesPage } from "@/routes/BaselinesPage";
import { VariantsPage } from "@/routes/VariantsPage";
import { OptimizedPage } from "@/routes/OptimizedPage";
import { ProductionPage } from "@/routes/ProductionPage";

const NAV = [
  { to: "/runs", label: "Runs", icon: Play },
  { to: "/baselines", label: "Baselines", icon: FlaskConical },
  { to: "/variants", label: "Variants", icon: Dna },
  { to: "/optimized", label: "Optimized", icon: Sparkles },
  { to: "/production", label: "Production", icon: Factory },
] as const;

function TopNav() {
  return (
    <header className="sticky top-0 z-40 flex h-14 items-center border-b bg-background/80 backdrop-blur">
      <div className="flex w-full items-center gap-6 px-6">
        <Link to="/runs" className="flex items-center gap-2 font-semibold">
          <span className="inline-flex h-7 w-7 items-center justify-center rounded-md bg-gradient-to-br from-sky-500 to-indigo-500 text-xs font-bold text-white">
            C
          </span>
          <span>CASCADE</span>
          <span className="text-xs font-normal text-muted-foreground">
            evolution dashboard
          </span>
        </Link>
        <nav className="flex items-center gap-1">
          {NAV.map(({ to, label, icon: Icon }) => (
            <NavLink
              key={to}
              to={to}
              className={({ isActive }) =>
                cn(
                  "inline-flex items-center gap-1.5 rounded-md px-3 py-1.5 text-sm font-medium transition-colors",
                  isActive
                    ? "bg-secondary text-secondary-foreground"
                    : "text-muted-foreground hover:bg-accent hover:text-accent-foreground",
                )
              }
            >
              <Icon className="h-4 w-4" />
              {label}
            </NavLink>
          ))}
        </nav>
      </div>
    </header>
  );
}

export default function App() {
  return (
    <div className="flex h-screen flex-col">
      <TopNav />
      <div className="flex min-h-0 flex-1">
        <main className="min-w-0 flex-1 overflow-y-auto px-6 py-6">
          <Routes>
            <Route path="/" element={<Navigate to="/runs" replace />} />
            <Route path="/runs" element={<RunsPage />} />
            <Route path="/baselines" element={<BaselinesPage />} />
            <Route path="/variants" element={<VariantsPage />} />
            <Route path="/optimized" element={<OptimizedPage />} />
            <Route path="/production" element={<ProductionPage />} />
          </Routes>
        </main>
        <OptimizedSidebar />
      </div>
      <Toaster />
    </div>
  );
}
