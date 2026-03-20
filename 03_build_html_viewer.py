#!/usr/bin/env python3
"""Step 03: Build self-contained HTML VenomsBase viewer.

Creates a single HTML file with embedded Plotly.js scatter plot + protein data.
No server needed — Todd opens it in Chrome and explores.

Features:
- Interactive 2D scatter (UMAP or PCA)
- Color by any tag dimension (dropdown)
- Click point → side panel with all tags
- Search by name/identifier
- Genomic context panel (placeholder → populated if GFF data available)

Inputs:
  - data/venombase_unified.csv
  - data/VenomsBase_Prototype.parquetbundle (for UMAP/PCA coordinates)

Output:
  - data/VenomsBase_Explorer.html (~2-5 MB, self-contained)
"""

import csv
import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

REPO_ROOT = Path(__file__).resolve().parent
DATA_DIR = REPO_ROOT / "data"


def load_bundle_projections(bundle_path):
    """Extract UMAP/PCA coordinates from parquetbundle."""
    delimiter = b"---PARQUET_DELIMITER---"
    raw = bundle_path.read_bytes()
    parts = raw.split(delimiter)

    import io
    ann_df = pd.read_parquet(io.BytesIO(parts[0]))
    proj_df = pd.read_parquet(io.BytesIO(parts[2]))

    # Get UMAP coordinates
    umap = proj_df[proj_df["projection_name"] == "UMAP_2"].set_index("identifier")
    pca = proj_df[proj_df["projection_name"] == "PCA_2"].set_index("identifier")

    return ann_df, umap, pca


def load_genomic_context():
    """Load genomic windows from GFF data via src/genomic_context.py."""
    from src.genomic_context import load_all_genomic_context
    return load_all_genomic_context(config.SP_DIR)


def build_html(annotations_df, umap_coords, pca_coords, genomic_data, output_path):
    """Build self-contained HTML with Plotly scatter + tag panel + genomic viewer."""

    # Prepare protein data as JSON
    proteins = []
    for _, row in annotations_df.iterrows():
        pid = row["identifier"]
        p = dict(row)
        # Add coordinates
        if pid in umap_coords.index:
            p["umap_x"] = float(umap_coords.loc[pid, "x"])
            p["umap_y"] = float(umap_coords.loc[pid, "y"])
        if pid in pca_coords.index:
            p["pca_x"] = float(pca_coords.loc[pid, "x"])
            p["pca_y"] = float(pca_coords.loc[pid, "y"])
        proteins.append(p)

    # Determine color-by columns (non-empty tag columns with <50 unique values)
    color_columns = []
    for col in annotations_df.columns:
        if col in ("identifier", "name", "description", "sequence.full", "sequence.mature"):
            continue
        nunique = annotations_df[col].nunique()
        if 2 <= nunique <= 50:
            color_columns.append(col)

    # Bake genomic windows into protein records
    # Each window is a scaffold region with venom genes + flanking
    genomic_windows = []
    for win in genomic_data:
        genomic_windows.append({
            "species": win["species"],
            "scaffold": win["scaffold"],
            "n_venom": win["n_venom"],
            "loci": [{"n": l["name"], "s": l["start"], "e": l["end"],
                       "d": l["strand"], "v": l["is_venom"]} for l in win["loci"]],
        })

    proteins_json = json.dumps(proteins, default=str)
    color_columns_json = json.dumps(color_columns)
    genomic_json_str = json.dumps(genomic_windows)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>VenomsBase Explorer</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
:root {{
  --bg: #0a0a0a; --bg2: #111; --bg3: #1a1a2e; --bg4: #0d0d1a;
  --fg: #e0e0e0; --fg2: #aaa; --fg3: #888; --fg4: #555;
  --border: #333; --accent: #7c9eff; --input-bg: #2a2a3e;
  --gene: #3a5a8a; --gene-stroke: #5a7aaa; --gene-hl: #F38400;
  --seq-color: #6a9; --cat-border: #222;
}}
body.light {{
  --bg: #f5f5f5; --bg2: #fff; --bg3: #e8eaf0; --bg4: #eef0f8;
  --fg: #222; --fg2: #555; --fg3: #777; --fg4: #aaa;
  --border: #ccc; --accent: #3355aa; --input-bg: #fff;
  --gene: #6b8ec0; --gene-stroke: #4a6e9a; --gene-hl: #d06000;
  --seq-color: #2a7a55; --cat-border: #ddd;
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; background: var(--bg); color: var(--fg); height: 100vh; overflow: hidden; }}
#app {{ display: grid; grid-template-columns: 1fr 360px; grid-template-rows: 44px 1fr 220px; height: 100vh; }}
#toolbar {{ grid-column: 1 / -1; background: var(--bg3); padding: 6px 16px; display: flex; align-items: center; gap: 12px; border-bottom: 1px solid var(--border); }}
#toolbar h1 {{ font-size: 15px; color: var(--accent); font-weight: 600; white-space: nowrap; }}
#toolbar select, #toolbar input {{ background: var(--input-bg); border: 1px solid var(--border); color: var(--fg); padding: 3px 6px; border-radius: 3px; font-size: 12px; }}
#toolbar input {{ width: 180px; }}
#toolbar label {{ font-size: 12px; color: var(--fg3); white-space: nowrap; }}
#theme-btn {{ background: none; border: 1px solid var(--border); color: var(--fg3); padding: 3px 8px; border-radius: 3px; cursor: pointer; font-size: 12px; }}
#scatter {{ background: var(--bg); }}
#panel {{ background: var(--bg2); border-left: 1px solid var(--border); overflow-y: auto; padding: 10px; }}
#panel h2 {{ font-size: 13px; color: var(--accent); margin-bottom: 6px; }}
.tag {{ margin-bottom: 3px; font-size: 11px; line-height: 1.4; }}
.tag-key {{ color: var(--fg3); }}
.tag-val {{ color: var(--fg); word-break: break-all; }}
.tag-seq {{ font-family: 'SF Mono', Menlo, monospace; font-size: 10px; color: var(--seq-color); max-height: 50px; overflow-y: auto; word-break: break-all; }}
.cat-section {{ margin-top: 6px; padding-top: 5px; border-top: 1px solid var(--cat-border); }}
.cat-label {{ font-size: 10px; color: var(--accent); margin-bottom: 3px; text-transform: uppercase; letter-spacing: 0.5px; }}
#genomic {{ grid-column: 1 / -1; background: var(--bg4); border-top: 1px solid var(--border); padding: 8px 12px; }}
#genomic-header {{ display: flex; align-items: center; gap: 12px; margin-bottom: 4px; }}
#genomic-header h3 {{ font-size: 12px; color: var(--accent); }}
#genomic-info {{ font-size: 11px; color: var(--fg3); }}
#genomic-canvas {{ width: 100%; height: 175px; display: block; }}
.no-selection {{ color: var(--fg4); font-size: 12px; padding: 16px; text-align: center; }}
#stats {{ font-size: 11px; color: var(--fg3); margin-left: auto; white-space: nowrap; }}
</style>
</head>
<body>
<div id="app">
  <div id="toolbar">
    <h1>VenomsBase</h1>
    <label>Color: <select id="colorBy"></select></label>
    <label>Proj: <select id="projection"><option value="umap">UMAP</option><option value="pca">PCA</option></select></label>
    <input type="text" id="search" placeholder="Search...">
    <button id="theme-btn">Light</button>
    <span id="stats"></span>
  </div>
  <div id="scatter"></div>
  <div id="panel"><div class="no-selection">Click a protein to inspect</div></div>
  <div id="genomic">
    <div id="genomic-header"><h3>Genomic Context</h3><span id="genomic-info"></span></div>
    <canvas id="genomic-canvas"></canvas>
  </div>
</div>
<script>
const PROTEINS = {proteins_json};
const COLOR_COLUMNS = {color_columns_json};
const GWINDOWS = {genomic_json_str};
const PALETTE = ['#F3C300','#875692','#F38400','#A1CAF1','#BE0032','#C2B280','#848482','#008856','#E68FAC','#0067A5','#F99379','#604E97','#F6A600','#B3446C','#DCD300','#882D17','#8DB600','#654522','#E25822','#2B3D26','#F2F3F4','#555555'];

let currentColorCol = '', currentProj = 'umap', isDark = true;
const $ = id => document.getElementById(id);
const css = getComputedStyle(document.body);
const getVar = v => getComputedStyle(document.body).getPropertyValue(v).trim();

// Theme toggle
$('theme-btn').addEventListener('click', () => {{
  isDark = !isDark;
  document.body.classList.toggle('light', !isDark);
  $('theme-btn').textContent = isDark ? 'Light' : 'Dark';
  renderScatter();
  if (lastProtein) drawGenomicContext(lastProtein);
}});

// Populate color dropdown
const colorSel = $('colorBy');
COLOR_COLUMNS.forEach(col => {{ const o = document.createElement('option'); o.value = col; o.textContent = col; colorSel.appendChild(o); }});
const defaultCol = COLOR_COLUMNS.find(c => c === 'venom.group') || COLOR_COLUMNS.find(c => c === 'venom.family') || COLOR_COLUMNS[0] || '';
colorSel.value = defaultCol; currentColorCol = defaultCol;

function getCoords(proj) {{
  return PROTEINS.filter(p => p[proj+'_x'] !== undefined).map(p => ({{ x:p[proj+'_x'], y:p[proj+'_y'], id:p.identifier, name:p.name||p.identifier, protein:p }}));
}}

function buildTraces(colorCol, proj) {{
  const pts = getCoords(proj);
  if (!colorCol) return [{{ x:pts.map(p=>p.x), y:pts.map(p=>p.y), text:pts.map(p=>p.name), customdata:pts.map(p=>p.protein), mode:'markers', type:'scattergl', marker:{{size:5, color:'#7c9eff', opacity:0.8}}, hovertemplate:'%{{text}}<extra></extra>' }}];
  const groups = {{}};
  pts.forEach(p => {{ const v = p.protein[colorCol] || '(empty)'; (groups[v] = groups[v]||[]).push(p); }});
  return Object.entries(groups).sort((a,b) => b[1].length - a[1].length).map(([val, pts], i) => ({{
    x:pts.map(p=>p.x), y:pts.map(p=>p.y), text:pts.map(p=>p.name), customdata:pts.map(p=>p.protein),
    mode:'markers', type:'scattergl',
    name: (val.length>25 ? val.substring(0,25)+'..' : val) + ' ('+pts.length+')',
    marker:{{ size:5, color:PALETTE[i%PALETTE.length], opacity:0.8 }},
    hovertemplate:'%{{text}}<extra>'+val.substring(0,40)+'</extra>'
  }}));
}}

function renderScatter() {{
  const bg = getVar('--bg'), fg4 = getVar('--fg4'), fg3 = getVar('--fg3');
  const traces = buildTraces(currentColorCol, currentProj);
  $('stats').textContent = traces.reduce((s,t) => s+t.x.length, 0) + ' proteins';
  Plotly.react('scatter', traces, {{
    paper_bgcolor: bg, plot_bgcolor: bg,
    margin: {{ t:8, b:28, l:36, r:8 }},
    xaxis: {{ showgrid:false, zeroline:false, color:fg4 }},
    yaxis: {{ showgrid:false, zeroline:false, color:fg4 }},
    legend: {{ font:{{size:10, color:fg3}}, bgcolor:'rgba(0,0,0,0)', x:1, xanchor:'right', y:1 }},
    dragmode:'pan', hovermode:'closest'
  }}, {{ responsive:true, scrollZoom:true }});
}}

let lastProtein = null;
function showProteinPanel(protein) {{
  lastProtein = protein;
  const panel = $('panel');
  let h = '<h2>' + (protein.name || protein.identifier) + '</h2>';
  if (protein.description) h += '<div style="color:var(--fg2);font-size:11px;margin-bottom:6px">' + protein.description + '</div>';
  const cats = {{}};
  Object.entries(protein).forEach(([k,v]) => {{
    if (!v || v==='' || k.endsWith('_x') || k.endsWith('_y')) return;
    const c = k.includes('.') ? k.split('.')[0] : 'core';
    (cats[c] = cats[c]||[]).push([k,v]);
  }});
  const order = ['core','taxonomy','venom','sequence','structure','prediction','identity','source'];
  const labels = {{ core:'Identity', taxonomy:'Taxonomy', venom:'Venom', sequence:'Sequence', structure:'Structure', prediction:'Predictions', identity:'Cross-references', source:'Source' }};
  order.forEach(c => {{
    if (!cats[c]) return;
    h += '<div class="cat-section"><div class="cat-label">'+(labels[c]||c)+'</div>';
    cats[c].forEach(([k,v]) => {{
      const dk = k.includes('.') ? k.split('.').slice(1).join('.') : k;
      const isSeq = k.includes('sequence.full') || k.includes('sequence.mature');
      h += isSeq && v.length>20
        ? '<div class="tag"><span class="tag-key">'+dk+':</span><div class="tag-seq">'+v+'</div></div>'
        : '<div class="tag"><span class="tag-key">'+dk+':</span> <span class="tag-val">'+v+'</span></div>';
    }});
    h += '</div>';
  }});
  panel.innerHTML = h;
}}

function drawGenomicContext(protein) {{
  const canvas = $('genomic-canvas');
  const ctx = canvas.getContext('2d');
  const dpr = window.devicePixelRatio || 1;
  const W = canvas.parentElement.clientWidth - 24;
  const H = 175;
  canvas.width = W * dpr; canvas.height = H * dpr;
  canvas.style.width = W + 'px'; canvas.style.height = H + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);

  const fgCol = getVar('--fg'), fg3Col = getVar('--fg3'), fg4Col = getVar('--fg4');
  const borderCol = getVar('--border'), hlCol = getVar('--gene-hl');
  const venomCol = getVar('--gene'), venomStroke = getVar('--gene-stroke');
  const greyFill = isDark ? '#2a2a2a' : '#ccc';
  const greyStroke = isDark ? '#444' : '#aaa';

  // Find best matching genomic window for this protein
  const pSpecies = (protein['taxonomy.species']||'').toLowerCase();
  const pName = (protein.name||'').toLowerCase();
  const pFamily = (protein['venom.family']||'').toLowerCase();
  let bestWin = null, bestGeneIdx = -1;

  for (const win of GWINDOWS) {{
    const wSpecies = win.species.toLowerCase();
    // Score: species match + gene name match
    const speciesHit = pSpecies && (wSpecies.includes(pSpecies.split(' ')[0]) || pSpecies.includes(wSpecies.split(' ')[0]));
    if (!speciesHit && !pFamily) continue;

    // Look for specific gene name match in this window
    for (let i = 0; i < win.loci.length; i++) {{
      const gn = win.loci[i].n.toLowerCase();
      const nameParts = pName.replace(/[^a-z0-9]/g,' ').split(/\\s+/).filter(p => p.length > 2);
      if (nameParts.some(part => gn.includes(part))) {{
        bestWin = win; bestGeneIdx = i; break;
      }}
    }}
    if (bestWin) break;

    // Fallback: species match with venom family match
    if (speciesHit && pFamily) {{
      const fam3 = pFamily.substring(0,3);
      for (let i = 0; i < win.loci.length; i++) {{
        if (win.loci[i].v && win.loci[i].n.toLowerCase().includes(fam3)) {{
          bestWin = win; bestGeneIdx = i; break;
        }}
      }}
    }}
    // Fallback: any species match to a window with venom genes
    if (!bestWin && speciesHit && win.n_venom > 0) {{
      bestWin = win;
    }}
    if (bestWin) break;
  }}

  // Final fallback: show any window with most venom genes (for demo)
  if (!bestWin && GWINDOWS.length > 0) {{
    bestWin = GWINDOWS.reduce((a, b) => a.n_venom > b.n_venom ? a : b);
  }}

  if (!bestWin || bestWin.loci.length === 0) {{
    ctx.fillStyle = fg4Col; ctx.font = '12px sans-serif';
    ctx.fillText('No genomic context available', 20, H/2);
    $('genomic-info').textContent = '';
    return;
  }}

  const loci = bestWin.loci;
  const minPos = Math.min(...loci.map(l => l.s));
  const maxPos = Math.max(...loci.map(l => l.e));
  const range = maxPos - minPos || 1;
  const mx = 40, trackY = Math.round(H * 0.52), gH = 20, arrowW = 7;
  const scale = (W - mx * 2) / range;

  // Info
  const nVenom = loci.filter(l => l.v).length;
  $('genomic-info').textContent = bestWin.species + '  ·  ' + bestWin.scaffold + '  ·  ' + (range/1000).toFixed(0) + ' kb  ·  ' + nVenom + ' venom / ' + loci.length + ' total genes';

  // Chromosome backbone
  ctx.strokeStyle = borderCol; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(mx - 8, trackY); ctx.lineTo(W - mx + 8, trackY); ctx.stroke();

  // Scale bar
  const niceStep = Math.pow(10, Math.floor(Math.log10(range / 4)));
  const sbPx = niceStep * scale;
  ctx.strokeStyle = fg4Col; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(W - mx - sbPx, H - 10); ctx.lineTo(W - mx, H - 10); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(W - mx - sbPx, H - 14); ctx.lineTo(W - mx - sbPx, H - 6); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(W - mx, H - 14); ctx.lineTo(W - mx, H - 6); ctx.stroke();
  ctx.fillStyle = fg4Col; ctx.font = '9px sans-serif'; ctx.textAlign = 'center';
  ctx.fillText(niceStep >= 1000 ? (niceStep/1000) + ' kb' : niceStep + ' bp', W - mx - sbPx/2, H - 1);
  ctx.textAlign = 'left';

  // Draw each gene
  loci.forEach((gene, i) => {{
    const x1 = mx + (gene.s - minPos) * scale;
    const x2 = mx + (gene.e - minPos) * scale;
    const w = Math.max(x2 - x1, 8);
    const isHL = i === bestGeneIdx;
    const isVenom = gene.v;
    const above = gene.d === '+';
    const y = above ? trackY - gH - 4 : trackY + 4;

    // Colors: venom = blue/accent, non-venom = grey, highlighted = orange
    let fill, stroke;
    if (isHL) {{ fill = hlCol; stroke = hlCol; }}
    else if (isVenom) {{ fill = venomCol; stroke = venomStroke; }}
    else {{ fill = greyFill; stroke = greyStroke; }}

    ctx.fillStyle = fill; ctx.strokeStyle = stroke; ctx.lineWidth = 1;
    const aw = Math.min(arrowW, w * 0.35);
    ctx.beginPath();
    if (gene.d === '+') {{
      ctx.moveTo(x1, y); ctx.lineTo(x1+w-aw, y); ctx.lineTo(x1+w, y+gH/2);
      ctx.lineTo(x1+w-aw, y+gH); ctx.lineTo(x1, y+gH);
    }} else {{
      ctx.moveTo(x1+aw, y); ctx.lineTo(x1+w, y); ctx.lineTo(x1+w, y+gH);
      ctx.lineTo(x1+aw, y+gH); ctx.lineTo(x1, y+gH/2);
    }}
    ctx.closePath(); ctx.fill(); ctx.stroke();

    // Stem
    ctx.strokeStyle = isHL ? hlCol+'66' : (isVenom ? venomStroke+'44' : greyStroke+'44');
    ctx.lineWidth = 0.5;
    ctx.beginPath();
    ctx.moveTo(x1+w/2, above ? y+gH : y);
    ctx.lineTo(x1+w/2, trackY);
    ctx.stroke();

    // Label — always for venom genes, only if space for non-venom
    if (isVenom || isHL || w > 40) {{
      ctx.fillStyle = isHL ? '#fff' : (isVenom ? fgCol : fg3Col);
      ctx.font = (isHL ? 'bold ' : '') + (isVenom ? '10px' : '9px') + ' sans-serif';
      const lbl = gene.n.length > 16 ? gene.n.substring(0,15)+'..' : gene.n;
      const ty = above ? y - 2 : y + gH + 11;
      ctx.fillText(lbl, x1, ty);
    }}
  }});

  // Highlight box around matched gene
  if (bestGeneIdx >= 0 && bestGeneIdx < loci.length) {{
    const g = loci[bestGeneIdx];
    const x1 = mx + (g.s - minPos) * scale;
    const x2 = mx + (g.e - minPos) * scale;
    const w = Math.max(x2-x1, 8);
    const above = g.d === '+';
    const y = above ? trackY - gH - 4 : trackY + 4;
    ctx.strokeStyle = hlCol; ctx.lineWidth = 2; ctx.setLineDash([3,2]);
    ctx.strokeRect(x1 - 4, y - 4, w + 8, gH + 8);
    ctx.setLineDash([]);
  }}
}}

// Events
colorSel.addEventListener('change', () => {{ currentColorCol = colorSel.value; renderScatter(); }});
$('projection').addEventListener('change', e => {{ currentProj = e.target.value; renderScatter(); }});

$('search').addEventListener('input', e => {{
  const q = e.target.value.toLowerCase();
  if (q.length < 2) return;
  const m = PROTEINS.filter(p => (p.name||'').toLowerCase().includes(q) || (p.identifier||'').toLowerCase().includes(q) || (p.description||'').toLowerCase().includes(q));
  if (m.length === 1) {{ showProteinPanel(m[0]); drawGenomicContext(m[0]); }}
}});

renderScatter();

setTimeout(() => {{
  $('scatter').on('plotly_click', data => {{
    if (data.points && data.points[0]) {{
      const p = data.points[0].customdata;
      showProteinPanel(p); drawGenomicContext(p);
    }}
  }});
}}, 300);
</script>
</body>
</html>"""

    output_path.write_text(html)
    size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"  Written: {output_path} ({size_mb:.1f} MB)")


def main():
    DATA_DIR.mkdir(exist_ok=True)

    print("Building VenomsBase HTML Explorer...")

    # Load bundle for projections
    print("\n1. Loading projections from bundle...")
    bundle_path = DATA_DIR / "VenomsBase_Prototype.parquetbundle"
    ann_df, umap, pca = load_bundle_projections(bundle_path)
    print(f"  {len(ann_df)} proteins with projections")

    # Load genomic context (baked windows)
    print("\n2. Loading genomic context...")
    genomic = load_genomic_context()
    total_genes = sum(len(w["loci"]) for w in genomic)
    total_venom = sum(w["n_venom"] for w in genomic)
    print(f"  {len(genomic)} windows, {total_genes} genes ({total_venom} venom)")

    # Build HTML
    print("\n3. Building HTML...")
    output = DATA_DIR / "VenomsBase_Explorer.html"
    build_html(ann_df, umap, pca, genomic, output)

    print("\nDone! Open in browser:")
    print(f"  open {output}")


if __name__ == "__main__":
    main()
