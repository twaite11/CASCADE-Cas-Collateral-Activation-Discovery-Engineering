function toFixed(v, digits = 3) {
  if (v === null || v === undefined || Number.isNaN(Number(v))) return "-";
  return Number(v).toFixed(digits);
}

function metricCell(label, value) {
  return `<div class="metric"><div class="label">${label}</div><div class="value">${value}</div></div>`;
}

async function fetchJson(url) {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`HTTP ${resp.status} for ${url}`);
  return resp.json();
}

function esc(v) {
  if (v === null || v === undefined) return "";
  return String(v)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
}

function renderOverview(data) {
  const warnings = (data.warnings || []).join(", ");
  const el = document.getElementById("overviewGrid");
  el.innerHTML = [
    metricCell("Total Variants", data.total_variants ?? 0),
    metricCell("Max Generation", data.max_generation ?? 0),
    metricCell("Optimized Switches", data.optimized_switches ?? 0),
    metricCell("Elite Count", data.elite_count ?? 0),
    metricCell("Validated Baseline Records", data.validated_baseline_records ?? 0),
    metricCell("Mean Fitness", toFixed(data.mean_fitness, 2)),
    metricCell("Mean ipTM", toFixed(data.mean_iptm, 3)),
    metricCell("Mean AF2-IG", toFixed(data.mean_af2_ig, 3)),
    metricCell("Last Update", data.dataset_last_updated ?? "N/A"),
    metricCell("Warnings", warnings || "none"),
  ].join("");
}

function renderHealth(data) {
  const el = document.getElementById("healthGrid");
  const stage = data.stage_stats || {};
  el.innerHTML = [
    metricCell("Records", data.total_records ?? 0),
    metricCell("Structure Success Rate", `${toFixed((data.structure_success_rate ?? 0) * 100, 1)}%`),
    metricCell("Full Ternary Rate", `${toFixed((data.full_ternary_rate ?? 0) * 100, 1)}%`),
    metricCell("Failure Count", data.failure_count ?? 0),
    metricCell("Best Variant", data.best_variant_id ?? "N/A"),
    metricCell("Best Fitness", toFixed(data.best_fitness, 2)),
    metricCell("Mini Eval Success", stage.mini_eval_success ?? 0),
    metricCell("Base Eval Success", stage.base_eval_success ?? 0),
    metricCell("MSA Reuse Seen", stage.msa_reuse_observed ?? 0),
    metricCell("Variants/hour", toFixed(data.throughput_variants_per_hour ?? 0, 2)),
  ].join("");
}

function renderProduction(data) {
  const totals = data.totals || {};
  const funnel = data.funnel || {};
  const el = document.getElementById("productionGrid");
  el.innerHTML = [
    metricCell("Generated", totals.generated ?? 0),
    metricCell("Passed Filter", totals.passed_filter ?? 0),
    metricCell("Optimized", totals.optimized ?? 0),
    metricCell("Elite", totals.elite ?? 0),
    metricCell("Generated->Passed", `${toFixed((funnel.generated_to_passed ?? 0) * 100, 1)}%`),
    metricCell("Passed->Optimized", `${toFixed((funnel.passed_to_optimized ?? 0) * 100, 1)}%`),
    metricCell("Optimized->Elite", `${toFixed((funnel.optimized_to_elite ?? 0) * 100, 1)}%`),
  ].join("");

  const rows = data.lineages || [];
  renderRows("lineageTable", rows, [
    (r) => esc(r.lineage_id ?? "-"),
    (r) => r.generated ?? "-",
    (r) => r.optimized ?? "-",
    (r) => r.elite ?? "-",
    (r) => `${toFixed((r.optimized_yield ?? 0) * 100, 1)}%`,
    (r) => `${toFixed((r.elite_yield ?? 0) * 100, 1)}%`,
  ], false);
}

function renderRows(tableId, rows, columns, clickable = true) {
  const tbody = document.querySelector(`#${tableId} tbody`);
  tbody.innerHTML = rows
    .map((r) => {
      const tds = columns.map((col) => `<td>${col(r)}</td>`).join("");
      if (!clickable) return `<tr>${tds}</tr>`;
      return `<tr data-variant-id="${esc(r.variant_id ?? "")}">${tds}</tr>`;
    })
    .join("");
}

function attachRowClicks(tableId) {
  const rows = document.querySelectorAll(`#${tableId} tbody tr[data-variant-id]`);
  for (const row of rows) {
    row.addEventListener("click", async () => {
      const id = row.getAttribute("data-variant-id");
      if (id) {
        await renderDetail(id);
      }
    });
  }
}

function renderOptimizedCards(rows) {
  const el = document.getElementById("optimizedCards");
  const top = (rows || []).slice(0, 6);
  if (!top.length) {
    el.innerHTML = `<div class="muted">No optimized variants yet.</div>`;
    return;
  }
  el.innerHTML = top
    .map(
      (r) => `
      <article class="optCard" data-variant-id="${esc(r.variant_id ?? "")}">
        <div class="optTitle">${esc(r.variant_id ?? "-")}</div>
        <div class="optMeta">Gen ${esc(r.generation ?? "-")} | Fitness ${toFixed(r.fitness, 2)}</div>
        <div class="optMeta">ipTM ${toFixed(r.iptm, 3)} | AF2-IG ${toFixed(r.af2_ig, 3)} | ON ${toFixed(r.on_dist_A, 2)}A</div>
        <div class="optMeta">Shift ${toFixed(r.hepn_shift_A, 2)}A</div>
        <div class="optReasons">${esc((r.optimized_reasons || []).join(", "))}</div>
      </article>`
    )
    .join("");
  const cards = document.querySelectorAll(".optCard[data-variant-id]");
  for (const card of cards) {
    card.addEventListener("click", async () => {
      const id = card.getAttribute("data-variant-id");
      if (id) {
        await renderDetail(id);
      }
    });
  }
}

function artifactLinks(artifacts) {
  const links = [];
  if (artifacts?.structure)
    links.push(
      `<a target="_blank" href="/api/structure-file?path=${encodeURIComponent(
        artifacts.structure
      )}">structure</a>`
    );
  if (artifacts?.fasta) links.push(`<code>${esc(artifacts.fasta)}</code>`);
  if (artifacts?.crrna) links.push(`<code>${esc(artifacts.crrna)}</code>`);
  return links.length ? links.join("<br/>") : "<span class='muted'>None</span>";
}

function collectStructureOptions(row) {
  const options = [];
  if (row.structure_path) options.push({ label: "primary", path: row.structure_path });
  const evalArtifacts = row.eval_artifacts || {};
  const optimized = row.optimized_artifacts || {};
  for (const [k, v] of Object.entries(evalArtifacts)) {
    if (k.endsWith("_structure")) options.push({ label: k, path: v });
  }
  if (optimized.structure) options.push({ label: "optimized_structure", path: optimized.structure });
  const uniq = [];
  const seen = new Set();
  for (const o of options) {
    if (!o.path || seen.has(o.path)) continue;
    seen.add(o.path);
    uniq.push(o);
  }
  return uniq;
}

function updateStructureSelect(options) {
  const sel = document.getElementById("structureSelect");
  sel.innerHTML = "";
  for (const o of options) {
    const opt = document.createElement("option");
    opt.value = o.path;
    opt.textContent = `${o.label}: ${o.path}`;
    sel.appendChild(opt);
  }
}

function getDomainRanges(row) {
  const domains = row?.domain_metadata?.domains || {};
  const h1 = domains.HEPN1 || {};
  const h2 = domains.HEPN2 || {};
  return {
    h1Start: Number(h1.start || 0),
    h1End: Number(h1.end || 0),
    h2Start: Number(h2.start || 0),
    h2End: Number(h2.end || 0),
  };
}

function applyDomainColoring(viewer, ranges) {
  viewer.setStyle({}, { cartoon: { color: "lightgray" } });
  if (ranges.h1Start > 0 && ranges.h1End >= ranges.h1Start) {
    viewer.setStyle(
      { resi: `${ranges.h1Start}-${ranges.h1End}` },
      { cartoon: { color: "orange" }, stick: { colorscheme: "orangeCarbon", radius: 0.2 } }
    );
  }
  if (ranges.h2Start > 0 && ranges.h2End >= ranges.h2Start) {
    viewer.setStyle(
      { resi: `${ranges.h2Start}-${ranges.h2End}` },
      { cartoon: { color: "cyan" }, stick: { colorscheme: "cyanCarbon", radius: 0.2 } }
    );
  }
  viewer.zoomTo();
  viewer.render();
}

async function loadIntoViewer(targetId, structurePath, ranges) {
  if (!structurePath) return;
  const viewer = $3Dmol.createViewer(targetId, { backgroundColor: "black" });
  const ext = structurePath.toLowerCase().endsWith(".pdb") ? "pdb" : "cif";
  const url = `/api/structure-file?path=${encodeURIComponent(structurePath)}`;
  const text = await fetch(url).then((r) => r.text());
  viewer.addModel(text, ext);
  applyDomainColoring(viewer, ranges);
}

async function renderDetail(variantId) {
  const panel = document.getElementById("detailPanel");
  panel.innerHTML = `<div class="muted">Loading ${esc(variantId)}...</div>`;
  try {
    const row = await fetchJson(`/api/variant/${encodeURIComponent(variantId)}`);
    const domain = row.domain_metadata || {};
    const catalog = row.catalog_metadata || {};
    const lineage = row.lineage_summary || {};
    const recent = row.lineage_recent || [];
    const evalArtifacts = row.eval_artifacts || {};
    const summaryLinks = Object.entries(evalArtifacts)
      .filter(([k]) => k.endsWith("_summary"))
      .map(([k, v]) => `<div>${esc(k)}: <code>${esc(v)}</code></div>`)
      .join("");
    panel.innerHTML = `
      <div class="detailHeader">
        <h3>${esc(row.variant_id || "-")}</h3>
        <div class="muted">Baseline: ${esc(row.baseline_id || "-")} | Gen ${esc(row.generation || "-")}</div>
      </div>
      <div class="grid">
        ${metricCell("Fitness", toFixed(row.fitness, 2))}
        ${metricCell("OFF Dist (A)", toFixed(row.off_dist_A, 2))}
        ${metricCell("ON Dist (A)", toFixed(row.on_dist_A, 2))}
        ${metricCell("HEPN Shift (A)", toFixed(row.hepn_shift_A, 2))}
        ${metricCell("ipTM", toFixed(row.iptm, 3))}
        ${metricCell("AF2-IG", toFixed(row.af2_ig, 3))}
      </div>
      <h4>Optimized classification</h4>
      <div>${row.optimized_switch ? "yes" : "no"} (${esc((row.optimized_reasons || []).join(", ") || "n/a")})</div>
      <h4>Optimized artifacts</h4>
      <div>${artifactLinks(row.optimized_artifacts || {})}</div>
      <h4>Eval summary paths</h4>
      <div>${summaryLinks || "<span class='muted'>None</span>"}</div>
      <h4>Validation metadata</h4>
      <pre>${esc(JSON.stringify(row.validation_metadata || {}, null, 2))}</pre>
      <h4>Domain metadata</h4>
      <pre>${esc(JSON.stringify(domain, null, 2))}</pre>
      <h4>Catalog metadata</h4>
      <pre>${esc(JSON.stringify(catalog, null, 2))}</pre>
      <h4>Lineage summary</h4>
      <pre>${esc(JSON.stringify(lineage, null, 2))}</pre>
      <h4>Recent lineage variants</h4>
      <pre>${esc(JSON.stringify(recent, null, 2))}</pre>
    `;
    document.getElementById("viewerVariant").value = row.variant_id || "";
    const opts = collectStructureOptions(row);
    updateStructureSelect(opts);
  } catch (err) {
    panel.innerHTML = `<div class="muted">Failed to load detail for ${esc(variantId)}.</div>`;
    console.error(err);
  }
}

function buildVariantQuery() {
  const q = new URLSearchParams();
  q.set("limit", "500");
  const search = document.getElementById("variantSearch").value.trim();
  const generation = document.getElementById("generationFilter").value.trim();
  const lineage = document.getElementById("lineageFilter").value.trim();
  const optimizedOnly = document.getElementById("optimizedOnly").checked;
  const eliteOnly = document.getElementById("eliteOnly").checked;
  const minFitness = document.getElementById("minFitness").value.trim();
  const minIptm = document.getElementById("minIptm").value.trim();
  const minAf2ig = document.getElementById("minAf2ig").value.trim();
  if (search) q.set("search", search);
  if (generation) q.set("generation", generation);
  if (lineage) q.set("lineage", lineage);
  if (optimizedOnly) q.set("optimized_only", "true");
  if (eliteOnly) q.set("elite_only", "true");
  if (minFitness) q.set("min_fitness", minFitness);
  if (minIptm) q.set("min_iptm", minIptm);
  if (minAf2ig) q.set("min_af2_ig", minAf2ig);
  return q.toString();
}

async function refreshAll() {
  try {
    const query = buildVariantQuery();
    const [overview, health, production, variants, optimized] = await Promise.all([
      fetchJson("/api/overview"),
      fetchJson("/api/pipeline-health"),
      fetchJson("/api/production"),
      fetchJson(`/api/variants?${query}`),
      fetchJson("/api/optimized-switches?limit=300"),
    ]);
    renderOverview(overview);
    renderHealth(health);
    renderProduction(production);
    renderOptimizedCards(optimized.rows || []);
    renderRows("optimizedTable", optimized.rows || [], [
      (r) => esc(r.variant_id ?? "-"),
      (r) => r.generation ?? "-",
      (r) => toFixed(r.fitness, 2),
      (r) => toFixed(r.iptm, 3),
      (r) => toFixed(r.af2_ig, 3),
      (r) => toFixed(r.on_dist_A, 2),
      (r) => esc((r.optimized_reasons || []).join(", ")),
    ]);
    attachRowClicks("optimizedTable");

    renderRows("variantsTable", variants.rows || [], [
      (r) => esc(r.variant_id ?? "-"),
      (r) => r.generation ?? "-",
      (r) => esc(r.baseline_id ?? "-"),
      (r) => toFixed(r.fitness, 2),
      (r) => toFixed(r.off_dist_A, 2),
      (r) => toFixed(r.on_dist_A, 2),
      (r) => toFixed(r.iptm, 3),
      (r) => toFixed(r.af2_ig, 3),
      (r) => (r.optimized_switch ? "yes" : "no"),
    ]);
    attachRowClicks("variantsTable");
  } catch (err) {
    console.error(err);
  }
}

let timer = null;
function startAutoRefresh() {
  const sec = Number(document.getElementById("refreshSeconds").value || 10);
  if (timer) clearInterval(timer);
  timer = setInterval(refreshAll, sec * 1000);
}

async function loadPrimaryViewer() {
  const variant = document.getElementById("viewerVariant").value.trim();
  if (!variant) return;
  const row = await fetchJson(`/api/variant/${encodeURIComponent(variant)}`);
  const ranges = getDomainRanges(row);
  const selected = document.getElementById("structureSelect").value;
  if (!selected) return;
  await loadIntoViewer("molViewerPrimary", selected, ranges);
}

async function loadCompareViewer() {
  const left = document.getElementById("viewerVariant").value.trim();
  const right = document.getElementById("compareVariant").value.trim();
  if (!left || !right) return;
  const cmp = await fetchJson(`/api/compare?left=${encodeURIComponent(left)}&right=${encodeURIComponent(right)}`);
  const rightRow = cmp.right;
  const rightOpts = collectStructureOptions(rightRow);
  if (!rightOpts.length) return;
  const ranges = getDomainRanges(rightRow);
  await loadIntoViewer("molViewerCompare", rightOpts[0].path, ranges);
}

document.getElementById("refreshNow").addEventListener("click", refreshAll);
document.getElementById("refreshSeconds").addEventListener("change", startAutoRefresh);
document.getElementById("applyFilters").addEventListener("click", refreshAll);
document.getElementById("loadStructure").addEventListener("click", loadPrimaryViewer);
document.getElementById("loadCompare").addEventListener("click", loadCompareViewer);
document.getElementById("variantSearch").addEventListener("keydown", (e) => {
  if (e.key === "Enter") refreshAll();
});
document.getElementById("generationFilter").addEventListener("keydown", (e) => {
  if (e.key === "Enter") refreshAll();
});
document.getElementById("lineageFilter").addEventListener("keydown", (e) => {
  if (e.key === "Enter") refreshAll();
});
refreshAll();
startAutoRefresh();
