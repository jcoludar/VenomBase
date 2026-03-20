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
    """Parse GFF data for genomic context panel.

    Extracts gene positions from SP and 3FTx GFF annotations.
    Returns dict: scaffold -> list of {name, start, end, strand, evalue}.
    """
    genomic_data = {}

    # Parse SP GFFs (exon-level BLAST hits → group by gene name → get gene bounds)
    for gff_name in ["TS.gff", "Crvi_mi2.gff"]:
        gff_path = config.SP_DIR / "genomic_annotations" / gff_name
        if not gff_path.exists():
            continue

        species = "Thamnophis_sirtalis" if "TS" in gff_name else "Crotalus_viridis"
        genes = {}  # gene_name -> {scaffold, min_start, max_end, strand, best_evalue}

        with open(gff_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                scaffold = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]

                # Extract gene name
                name = ""
                for attr in attrs.split(";"):
                    if attr.startswith("Name="):
                        name = attr.split("=")[1].split("_e")[0]  # Remove exon suffix
                        break

                if not name:
                    continue

                key = f"{species}|{scaffold}|{name}"
                if key not in genes:
                    genes[key] = {
                        "species": species,
                        "scaffold": scaffold,
                        "name": name,
                        "start": min(start, end),
                        "end": max(start, end),
                        "strand": strand,
                    }
                else:
                    genes[key]["start"] = min(genes[key]["start"], start, end)
                    genes[key]["end"] = max(genes[key]["end"], start, end)

        # Group by scaffold
        for key, gene in genes.items():
            scaffold_key = f"{gene['species']}|{gene['scaffold']}"
            if scaffold_key not in genomic_data:
                genomic_data[scaffold_key] = []
            genomic_data[scaffold_key].append(gene)

    # Sort genes by position within each scaffold
    for key in genomic_data:
        genomic_data[key].sort(key=lambda g: g["start"])

    return genomic_data


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

    # Build genomic context JSON (simplified for embedding)
    genomic_json = {}
    for scaffold_key, genes in genomic_data.items():
        genomic_json[scaffold_key] = [
            {"name": g["name"], "start": g["start"], "end": g["end"], "strand": g["strand"]}
            for g in genes
        ]

    proteins_json = json.dumps(proteins, default=str)
    color_columns_json = json.dumps(color_columns)
    genomic_json_str = json.dumps(genomic_json)

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
const GENOMIC = {genomic_json_str};
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

  const bgCol = getVar('--bg4'), fgCol = getVar('--fg'), fg3Col = getVar('--fg3'), fg4Col = getVar('--fg4');
  const borderCol = getVar('--border'), geneCol = getVar('--gene'), geneStroke = getVar('--gene-stroke'), hlCol = getVar('--gene-hl');

  // Match protein to genomic region
  const species = (protein['taxonomy.species']||'').replace(/ /g,'_');
  const pname = (protein.name||'').toLowerCase();
  const pfam = (protein['venom.family']||'').toLowerCase();
  let matchedRegion = null, matchedGene = null;

  // Strategy: match by species prefix in scaffold key, then by gene name similarity
  for (const [key, genes] of Object.entries(GENOMIC)) {{
    const keySpecies = key.split('|')[0].toLowerCase();
    // Try species match first
    const speciesMatch = species && keySpecies.includes(species.split('_')[0].toLowerCase());
    if (!speciesMatch) continue;

    for (const gene of genes) {{
      const gn = gene.name.toLowerCase();
      // Match by protein name parts
      const nameParts = pname.replace(/[^a-z0-9]/g,' ').split(/\\s+/).filter(p => p.length > 2);
      const hit = nameParts.some(part => gn.includes(part));
      if (hit) {{
        matchedRegion = {{ key, genes }};
        matchedGene = gene;
        break;
      }}
    }}
    // If we found a species match but no specific gene, show the region anyway
    if (!matchedRegion && genes.length > 0) {{
      matchedRegion = {{ key, genes }};
    }}
    if (matchedRegion) break;
  }}

  // Fallback: match by venom family to any region
  if (!matchedRegion && pfam) {{
    for (const [key, genes] of Object.entries(GENOMIC)) {{
      for (const gene of genes) {{
        if (gene.name.toLowerCase().includes(pfam.substring(0,3))) {{
          matchedRegion = {{ key, genes }};
          matchedGene = gene;
          break;
        }}
      }}
      if (matchedRegion) break;
    }}
  }}

  if (!matchedRegion) {{
    ctx.fillStyle = fg4Col;
    ctx.font = '12px sans-serif';
    ctx.fillText('No genomic context available for this protein', 20, H/2 - 8);
    ctx.font = '11px sans-serif';
    ctx.fillText('Requires GFF annotation linking protein to genomic coordinates', 20, H/2 + 10);
    $('genomic-info').textContent = '';
    return;
  }}

  const genes = matchedRegion.genes;
  const minPos = Math.min(...genes.map(g => g.start));
  const maxPos = Math.max(...genes.map(g => g.end));
  const range = maxPos - minPos || 1;
  const mx = 50, trackY = H * 0.55, gH = 22, arrowW = 8;
  const scale = (W - mx*2) / range;

  // Info label
  const parts = matchedRegion.key.split('|');
  $('genomic-info').textContent = parts[0].replace(/_/g,' ') + '  ·  ' + parts[1] + '  ·  ' + (range/1000).toFixed(0) + ' kb region  ·  ' + genes.length + ' genes';

  // Scale bar
  const niceStep = Math.pow(10, Math.floor(Math.log10(range/4)));
  const scaleBarBp = niceStep;
  const scaleBarPx = scaleBarBp * scale;
  ctx.strokeStyle = fg4Col; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(mx, H-12); ctx.lineTo(mx + scaleBarPx, H-12); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(mx, H-16); ctx.lineTo(mx, H-8); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(mx+scaleBarPx, H-16); ctx.lineTo(mx+scaleBarPx, H-8); ctx.stroke();
  ctx.fillStyle = fg4Col; ctx.font = '9px sans-serif';
  ctx.fillText((scaleBarBp >= 1000 ? (scaleBarBp/1000)+'kb' : scaleBarBp+'bp'), mx + scaleBarPx/2 - 10, H-2);

  // Chromosome line
  ctx.strokeStyle = borderCol; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(mx-10, trackY); ctx.lineTo(W-mx+10, trackY); ctx.stroke();

  // Tick marks
  ctx.strokeStyle = fg4Col; ctx.lineWidth = 0.5;
  for (let p = Math.ceil(minPos/niceStep)*niceStep; p <= maxPos; p += niceStep) {{
    const x = mx + (p - minPos) * scale;
    ctx.beginPath(); ctx.moveTo(x, trackY-3); ctx.lineTo(x, trackY+3); ctx.stroke();
  }}

  // Draw genes
  genes.forEach((gene, i) => {{
    const x1 = mx + (gene.start - minPos) * scale;
    const x2 = mx + (gene.end - minPos) * scale;
    const w = Math.max(x2 - x1, 6);
    const isHL = gene === matchedGene;
    const above = gene.strand === '+';
    const y = above ? trackY - gH - 6 : trackY + 6;

    // Arrow body
    ctx.fillStyle = isHL ? hlCol : geneCol;
    ctx.strokeStyle = isHL ? hlCol : geneStroke;
    ctx.lineWidth = 1;
    ctx.beginPath();
    if (gene.strand === '+') {{
      const aw = Math.min(arrowW, w * 0.3);
      ctx.moveTo(x1, y); ctx.lineTo(x1+w-aw, y); ctx.lineTo(x1+w, y+gH/2);
      ctx.lineTo(x1+w-aw, y+gH); ctx.lineTo(x1, y+gH);
    }} else {{
      const aw = Math.min(arrowW, w * 0.3);
      ctx.moveTo(x1+aw, y); ctx.lineTo(x1+w, y); ctx.lineTo(x1+w, y+gH);
      ctx.lineTo(x1+aw, y+gH); ctx.lineTo(x1, y+gH/2);
    }}
    ctx.closePath(); ctx.fill(); ctx.stroke();

    // Stem to track
    ctx.strokeStyle = isHL ? hlCol+'88' : borderCol;
    ctx.lineWidth = 0.5;
    ctx.beginPath();
    ctx.moveTo(x1+w/2, above ? y+gH : y);
    ctx.lineTo(x1+w/2, trackY);
    ctx.stroke();

    // Label (only if gene wide enough or highlighted)
    if (w > 30 || isHL) {{
      ctx.fillStyle = isHL ? '#fff' : fgCol;
      ctx.font = (isHL ? 'bold ' : '') + '9px sans-serif';
      const lbl = gene.name.length > 14 ? gene.name.substring(0,13)+'..' : gene.name;
      const textY = above ? y - 3 : y + gH + 10;
      ctx.fillText(lbl, x1 + 2, textY);
    }}
  }});

  // Highlight ring
  if (matchedGene) {{
    const x1 = mx + (matchedGene.start - minPos) * scale;
    const x2 = mx + (matchedGene.end - minPos) * scale;
    const w = Math.max(x2-x1, 6);
    const above = matchedGene.strand === '+';
    const y = above ? trackY - gH - 6 : trackY + 6;
    ctx.strokeStyle = hlCol; ctx.lineWidth = 2;
    ctx.strokeRect(x1-3, y-3, w+6, gH+6);
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

    # Load genomic context
    print("\n2. Loading genomic context...")
    genomic = load_genomic_context()
    total_genes = sum(len(g) for g in genomic.values())
    print(f"  {len(genomic)} scaffolds, {total_genes} genes")

    # Build HTML
    print("\n3. Building HTML...")
    output = DATA_DIR / "VenomsBase_Explorer.html"
    build_html(ann_df, umap, pca, genomic, output)

    print("\nDone! Open in browser:")
    print(f"  open {output}")


if __name__ == "__main__":
    main()
