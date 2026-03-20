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
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; background: #0a0a0a; color: #e0e0e0; height: 100vh; overflow: hidden; }}
#app {{ display: grid; grid-template-columns: 1fr 380px; grid-template-rows: 48px 1fr 240px; height: 100vh; }}
#toolbar {{ grid-column: 1 / -1; background: #1a1a2e; padding: 8px 16px; display: flex; align-items: center; gap: 16px; border-bottom: 1px solid #333; }}
#toolbar h1 {{ font-size: 16px; color: #7c9eff; font-weight: 600; }}
#toolbar select, #toolbar input {{ background: #2a2a3e; border: 1px solid #444; color: #e0e0e0; padding: 4px 8px; border-radius: 4px; font-size: 13px; }}
#toolbar input {{ width: 200px; }}
#scatter {{ background: #0a0a0a; }}
#panel {{ background: #111; border-left: 1px solid #333; overflow-y: auto; padding: 12px; }}
#panel h2 {{ font-size: 14px; color: #7c9eff; margin-bottom: 8px; }}
#panel .tag {{ margin-bottom: 4px; font-size: 12px; }}
#panel .tag-key {{ color: #888; }}
#panel .tag-val {{ color: #e0e0e0; word-break: break-all; }}
#panel .tag-seq {{ font-family: monospace; font-size: 10px; color: #6a9; max-height: 60px; overflow-y: auto; }}
#genomic {{ grid-column: 1 / -1; background: #0d0d1a; border-top: 1px solid #333; padding: 8px; overflow-x: auto; }}
#genomic h3 {{ font-size: 12px; color: #7c9eff; margin-bottom: 4px; }}
#genomic-canvas {{ width: 100%; height: 190px; }}
.no-selection {{ color: #555; font-size: 13px; padding: 20px; text-align: center; }}
#stats {{ font-size: 12px; color: #888; }}
</style>
</head>
<body>
<div id="app">
  <div id="toolbar">
    <h1>VenomsBase Explorer</h1>
    <label>Color by: <select id="colorBy"></select></label>
    <label>Projection: <select id="projection">
      <option value="umap">UMAP</option>
      <option value="pca">PCA</option>
    </select></label>
    <input type="text" id="search" placeholder="Search protein name...">
    <span id="stats"></span>
  </div>
  <div id="scatter"></div>
  <div id="panel">
    <div class="no-selection">Click a protein point to see its tags</div>
  </div>
  <div id="genomic">
    <h3>Genomic Context</h3>
    <canvas id="genomic-canvas"></canvas>
  </div>
</div>
<script>
const PROTEINS = {proteins_json};
const COLOR_COLUMNS = {color_columns_json};
const GENOMIC = {genomic_json_str};

// Kelly's 22 maximally distinct colors
const PALETTE = [
  '#F3C300','#875692','#F38400','#A1CAF1','#BE0032','#C2B280',
  '#848482','#008856','#E68FAC','#0067A5','#F99379','#604E97',
  '#F6A600','#B3446C','#DCD300','#882D17','#8DB600','#654522',
  '#E25822','#2B3D26','#F2F3F4','#222222'
];

let currentColorCol = '';
let currentProj = 'umap';

// Populate color dropdown
const colorSel = document.getElementById('colorBy');
COLOR_COLUMNS.forEach(col => {{
  const opt = document.createElement('option');
  opt.value = col;
  opt.textContent = col;
  colorSel.appendChild(opt);
}});
if (COLOR_COLUMNS.length > 0) {{
  // Default to venom.family or first available
  const defaultCol = COLOR_COLUMNS.find(c => c === 'venom.family') || COLOR_COLUMNS[0];
  colorSel.value = defaultCol;
  currentColorCol = defaultCol;
}}

function getCoords(proj) {{
  return PROTEINS.filter(p => p[proj + '_x'] !== undefined).map(p => ({{
    x: p[proj + '_x'],
    y: p[proj + '_y'],
    id: p.identifier,
    name: p.name || p.identifier,
    protein: p
  }}));
}}

function buildTraces(colorCol, proj) {{
  const pts = getCoords(proj);
  if (!colorCol) {{
    return [{{
      x: pts.map(p => p.x),
      y: pts.map(p => p.y),
      text: pts.map(p => p.name),
      customdata: pts.map(p => p.protein),
      mode: 'markers',
      type: 'scattergl',
      marker: {{ size: 4, color: '#7c9eff', opacity: 0.7 }},
      hovertemplate: '%{{text}}<extra></extra>'
    }}];
  }}

  // Group by color value
  const groups = {{}};
  pts.forEach(p => {{
    const val = p.protein[colorCol] || '(empty)';
    if (!groups[val]) groups[val] = [];
    groups[val].push(p);
  }});

  const sortedGroups = Object.entries(groups).sort((a, b) => b[1].length - a[1].length);
  return sortedGroups.map(([val, pts], i) => ({{
    x: pts.map(p => p.x),
    y: pts.map(p => p.y),
    text: pts.map(p => p.name),
    customdata: pts.map(p => p.protein),
    mode: 'markers',
    type: 'scattergl',
    name: val.substring(0, 30) + (val.length > 30 ? '...' : '') + ' (' + pts.length + ')',
    marker: {{ size: 4, color: PALETTE[i % PALETTE.length], opacity: 0.7 }},
    hovertemplate: '%{{text}}<extra>' + val.substring(0, 40) + '</extra>'
  }}));
}}

function renderScatter() {{
  const traces = buildTraces(currentColorCol, currentProj);
  const total = traces.reduce((s, t) => s + t.x.length, 0);
  document.getElementById('stats').textContent = total + ' proteins';

  Plotly.react('scatter', traces, {{
    paper_bgcolor: '#0a0a0a',
    plot_bgcolor: '#0a0a0a',
    margin: {{ t: 10, b: 30, l: 40, r: 10 }},
    xaxis: {{ showgrid: false, zeroline: false, color: '#555' }},
    yaxis: {{ showgrid: false, zeroline: false, color: '#555' }},
    legend: {{ font: {{ size: 10, color: '#aaa' }}, bgcolor: 'rgba(0,0,0,0.5)', x: 1, xanchor: 'right' }},
    dragmode: 'pan',
    hovermode: 'closest'
  }}, {{ responsive: true, scrollZoom: true }});
}}

function showProteinPanel(protein) {{
  const panel = document.getElementById('panel');
  let html = '<h2>' + (protein.name || protein.identifier) + '</h2>';
  if (protein.description) {{
    html += '<div style="color:#aaa;font-size:12px;margin-bottom:8px">' + protein.description + '</div>';
  }}

  // Group tags by category
  const categories = {{}};
  Object.entries(protein).forEach(([key, val]) => {{
    if (!val || val === '' || key.endsWith('_x') || key.endsWith('_y')) return;
    const cat = key.includes('.') ? key.split('.')[0] : 'core';
    if (!categories[cat]) categories[cat] = [];
    categories[cat].push([key, val]);
  }});

  const catOrder = ['core', 'taxonomy', 'venom', 'sequence', 'structure', 'prediction', 'identity', 'source'];
  const catLabels = {{
    core: 'Identity', taxonomy: 'Taxonomy', venom: 'Venom', sequence: 'Sequence',
    structure: 'Structure', prediction: 'Predictions', identity: 'Cross-references', source: 'Source'
  }};

  catOrder.forEach(cat => {{
    if (!categories[cat]) return;
    html += '<div style="margin-top:8px;padding-top:6px;border-top:1px solid #222">';
    html += '<div style="font-size:11px;color:#7c9eff;margin-bottom:4px;text-transform:uppercase">' + (catLabels[cat] || cat) + '</div>';
    categories[cat].forEach(([key, val]) => {{
      const displayKey = key.includes('.') ? key.split('.').slice(1).join('.') : key;
      const isSeq = key.includes('sequence.full') || key.includes('sequence.mature');
      if (isSeq && val.length > 20) {{
        html += '<div class="tag"><span class="tag-key">' + displayKey + ':</span> <div class="tag-seq">' + val + '</div></div>';
      }} else {{
        html += '<div class="tag"><span class="tag-key">' + displayKey + ':</span> <span class="tag-val"> ' + val + '</span></div>';
      }}
    }});
    html += '</div>';
  }});
  panel.innerHTML = html;
}}

function drawGenomicContext(protein) {{
  const canvas = document.getElementById('genomic-canvas');
  const ctx = canvas.getContext('2d');
  const W = canvas.width = canvas.parentElement.clientWidth - 16;
  const H = canvas.height = 190;
  ctx.clearRect(0, 0, W, H);

  // Find matching scaffold for this protein
  const species = (protein['taxonomy.species'] || '').replace(/ /g, '_');
  let matchedRegion = null;
  let matchedGene = null;

  for (const [key, genes] of Object.entries(GENOMIC)) {{
    for (const gene of genes) {{
      const nameMatch = protein.name && gene.name.includes(protein.name.split('_')[0]);
      const familyMatch = protein['venom.family'] && gene.name.toLowerCase().includes(protein['venom.family'].toLowerCase().substring(0, 3));
      if (nameMatch || familyMatch) {{
        matchedRegion = {{ key, genes }};
        matchedGene = gene;
        break;
      }}
    }}
    if (matchedRegion) break;
  }}

  if (!matchedRegion) {{
    ctx.fillStyle = '#555';
    ctx.font = '13px sans-serif';
    ctx.fillText('No genomic context available for this protein', 20, H / 2);
    ctx.fillStyle = '#444';
    ctx.font = '11px sans-serif';
    ctx.fillText('Genomic coordinates needed (GFF annotation)', 20, H / 2 + 20);
    return;
  }}

  // Draw gene arrows
  const genes = matchedRegion.genes;
  const minPos = Math.min(...genes.map(g => g.start));
  const maxPos = Math.max(...genes.map(g => g.end));
  const range = maxPos - minPos || 1;
  const margin = 60;
  const trackY = H / 2;
  const geneH = 24;
  const scale = (W - margin * 2) / range;

  // Scaffold line
  ctx.strokeStyle = '#333';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(margin, trackY);
  ctx.lineTo(W - margin, trackY);
  ctx.stroke();

  // Label
  ctx.fillStyle = '#666';
  ctx.font = '10px sans-serif';
  const parts = matchedRegion.key.split('|');
  ctx.fillText(parts[0] + ' : ' + parts[1], margin, 14);
  ctx.fillText(minPos.toLocaleString() + ' – ' + maxPos.toLocaleString() + ' bp', margin, 26);

  // Draw genes
  genes.forEach((gene, i) => {{
    const x1 = margin + (gene.start - minPos) * scale;
    const x2 = margin + (gene.end - minPos) * scale;
    const w = Math.max(x2 - x1, 3);
    const isHighlight = gene === matchedGene;
    const row = i % 2;  // stagger
    const y = trackY - geneH / 2 - (row ? geneH + 4 : 0);

    // Gene body
    ctx.fillStyle = isHighlight ? '#F38400' : '#3a5a8a';
    ctx.strokeStyle = isHighlight ? '#F38400' : '#5a7aaa';
    ctx.lineWidth = 1;

    // Arrow shape
    ctx.beginPath();
    if (gene.strand === '+') {{
      ctx.moveTo(x1, y);
      ctx.lineTo(x1 + w - 6, y);
      ctx.lineTo(x1 + w, y + geneH / 2);
      ctx.lineTo(x1 + w - 6, y + geneH);
      ctx.lineTo(x1, y + geneH);
    }} else {{
      ctx.moveTo(x1 + 6, y);
      ctx.lineTo(x1 + w, y);
      ctx.lineTo(x1 + w, y + geneH);
      ctx.lineTo(x1 + 6, y + geneH);
      ctx.lineTo(x1, y + geneH / 2);
    }}
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    // Gene name
    ctx.fillStyle = isHighlight ? '#fff' : '#ccc';
    ctx.font = (isHighlight ? 'bold ' : '') + '10px sans-serif';
    const label = gene.name.length > 12 ? gene.name.substring(0, 12) + '..' : gene.name;
    ctx.fillText(label, x1 + 3, y + geneH / 2 + 3);

    // Connector to track
    if (row) {{
      ctx.strokeStyle = '#333';
      ctx.beginPath();
      ctx.moveTo(x1 + w / 2, y + geneH);
      ctx.lineTo(x1 + w / 2, trackY);
      ctx.stroke();
    }}
  }});
}}

// Event handlers
colorSel.addEventListener('change', () => {{ currentColorCol = colorSel.value; renderScatter(); }});
document.getElementById('projection').addEventListener('change', (e) => {{ currentProj = e.target.value; renderScatter(); }});

document.getElementById('scatter').addEventListener('plotly_click', (e) => {{
  if (e.detail && e.detail.points && e.detail.points[0]) {{
    const protein = e.detail.points[0].customdata;
    showProteinPanel(protein);
    drawGenomicContext(protein);
  }}
}});

// Plotly click via internal API
document.getElementById('scatter').on && document.getElementById('scatter').on('plotly_click', (data) => {{
  if (data.points && data.points[0]) {{
    const protein = data.points[0].customdata;
    showProteinPanel(protein);
    drawGenomicContext(protein);
  }}
}});

// Search
document.getElementById('search').addEventListener('input', (e) => {{
  const q = e.target.value.toLowerCase();
  if (q.length < 2) {{ renderScatter(); return; }}
  // Highlight matching points (just re-render with filter info)
  const matches = PROTEINS.filter(p =>
    (p.name || '').toLowerCase().includes(q) ||
    (p.identifier || '').toLowerCase().includes(q) ||
    (p.description || '').toLowerCase().includes(q)
  );
  if (matches.length === 1) {{
    showProteinPanel(matches[0]);
    drawGenomicContext(matches[0]);
  }}
}});

// Initial render
renderScatter();

// Fix Plotly click for scattergl
setTimeout(() => {{
  const scatterDiv = document.getElementById('scatter');
  scatterDiv.on('plotly_click', (data) => {{
    if (data.points && data.points[0]) {{
      const protein = data.points[0].customdata;
      showProteinPanel(protein);
      drawGenomicContext(protein);
    }}
  }});
}}, 500);
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
