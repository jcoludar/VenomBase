#!/usr/bin/env python3
"""Step 05: Build unified snake venom VenomsBase viewer.

Merges three venom protein families into one interactive HTML explorer:
  - 3FTx (Nature Comms 2023, Koludarov et al.) — toxin groups only (no Ly6)
  - KLK/SP serine proteases (BMC Biology 2021, Barua & Koludarov)
  - PLA2 phospholipases (MBE 2020)

Computes UMAP + PCA on the merged ProtT5 embedding matrix, parses genomic
context from SP GFFs + PLA2 SM5 GFF, and builds a self-contained HTML file.

Output: data/VenomsBase_Snakes.html
"""

import csv
import importlib.util
import json
import re
import sys
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "tools"))
import config
from src.genomic_context import load_all_genomic_context

DATA_DIR = Path(__file__).resolve().parent / "data"

# ── Harmonized color palettes ────────────────────────────────────────────────

# Family-level colors (for "color by venom_family")
FAMILY_COLORS = {
    "3FTx": "#1565C0",
    "Serine protease": "#7B1FA2",
    "PLA2": "#E64B35",
}

# 3FTx subgroup colors (Nature Comms 2023, Fig 3) — prefixed
THRFTX_COLORS = {
    "3FTx Short-chain": "#64B5F6",
    "3FTx Long-chain": "#1565C0",
    "3FTx Plesiotypic": "#4CAF50",
    "3FTx Non-standard": "#827717",
}

# KLK/SP colors
KLK_COLORS = {
    "KLK": "#7B1FA2",
}

# PLA2 gene group colors (bioRxiv 2020, Fig 1) — prefixed
PLA2_COLORS = {
    "Pla2g2 E": "#E91E63",
    "Pla2g2 F": "#E64A19",
    "Pla2g2 C": "#FF9800",
    "Pla2g2 D1": "#4CAF50",
    "Pla2g2 D2": "#2E7D32",
    "Pla2g2 D3": "#66BB6A",
    "Pla2g2 A": "#7B1FA2",
    "Pla2g2 V1": "#283593",
    "Pla2g2 V2": "#3F51B5",
    "Pla2g2 B": "#455A64",
    "Pla2g2 Ga": "#6D4C41",
    "Pla2g2 Gb": "#827717",
    "Pla2g2 Gc": "#F9A825",
    "Pla2g2 Gk": "#CDDC39",
    "Pla2g2 G0": "#FFEB3B",
    "Pla2g2 preGk": "#FFF176",
    "Pla2g2 pre-g2": "#BDBDBD",
    "Pla2g2 V": "#3949AB",
    "Pla2g2 pre-V": "#5C6BC0",
    "Pla2g2 preg2": "#BDBDBD",
    "Pla2g2": "#E64B35",  # fallback for unknown group
}

# Merge all subgroup colors into one dict for the HTML color override
ALL_GROUP_COLORS = {}
ALL_GROUP_COLORS.update(THRFTX_COLORS)
ALL_GROUP_COLORS.update(KLK_COLORS)
ALL_GROUP_COLORS.update(PLA2_COLORS)

# 3FTx groups to keep (drop all Ly6 groups)
TOXIN_GROUPS = {"Short-chain", "Long-chain", "Plesiotypic", "Non-standard"}

# ── PLA2 SM5 GFF parser ─────────────────────────────────────────────────────

PLA2_VENOM_PREFIX = "Pla2g2"
PLA2_FLANKING = 2  # non-venom genes to show on each side


def parse_pla2_sm5_gff(gff_path):
    """Parse PLA2 SM5 GFF (Geneious manual annotations) into genomic windows.

    Structure:
    - Multiple ##sequence-region headers delimit scaffolds
    - Feature types: gene, pseudogene, CDS, exon, source, Assembly_gap
    - 'source' features contain species name: Name=source <Species>
    - Gene names starting with Pla2g2 are venom; others (OTUD3, UBXN10, etc.) are flanking
    """
    # Pass 1: collect scaffold regions + species + gene features
    scaffolds = {}  # scaffold_id -> {start, end, species, genes: []}
    current_scaffold = None

    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Sequence-region header
            m = re.match(r"^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)", line)
            if m:
                current_scaffold = m.group(1)
                if current_scaffold not in scaffolds:
                    scaffolds[current_scaffold] = {
                        "start": int(m.group(2)),
                        "end": int(m.group(3)),
                        "species": "",
                        "genes": [],
                    }
                continue

            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            scaffold_id = parts[0]
            feat_type = parts[2]

            # Extract Name attribute
            name = ""
            for attr in parts[8].split(";"):
                attr = attr.strip()
                if attr.startswith("Name="):
                    name = attr.split("=", 1)[1]

            # Source features → species name
            if feat_type == "source" and name.startswith("source "):
                species_raw = name[len("source "):]
                # Clean up breed/source prefixes
                if "source " in species_raw:
                    species_raw = species_raw.split("source ")[-1]
                if scaffold_id in scaffolds:
                    scaffolds[scaffold_id]["species"] = species_raw.strip()
                continue

            # Only keep gene and pseudogene features
            if feat_type not in ("gene", "pseudogene"):
                continue

            if not name or scaffold_id not in scaffolds:
                continue

            s, e = int(parts[3]), int(parts[4])
            start, end = min(s, e), max(s, e)
            strand = parts[6] if parts[6] in ("+", "-") else "+"

            # Determine if this is a venom gene
            is_venom = PLA2_VENOM_PREFIX in name

            # Clean gene name for display
            display_name = name
            # Strip parenthetical notes for display
            if "(" in display_name:
                display_name = display_name.split("(")[0].strip()

            scaffolds[scaffold_id]["genes"].append({
                "name": display_name,
                "full_name": name,
                "start": start,
                "end": end,
                "strand": strand,
                "is_venom": is_venom,
                "is_pseudo": feat_type == "pseudogene" or "pseudogene" in name.lower(),
            })

    # Build windows: for each scaffold with venom genes, create a window
    windows = []
    for scaffold_id, data in scaffolds.items():
        genes = data["genes"]
        if not genes:
            continue

        # Sort by position
        genes.sort(key=lambda g: g["start"])

        # Deduplicate overlapping genes (keep first occurrence by name)
        seen_names = set()
        deduped = []
        for g in genes:
            key = (g["name"], g["start"])
            if key not in seen_names:
                seen_names.add(key)
                deduped.append(g)
        genes = deduped

        venom_indices = [i for i, g in enumerate(genes) if g["is_venom"]]
        if not venom_indices:
            continue

        # Window: venom genes + flanking
        first_v = min(venom_indices)
        last_v = max(venom_indices)
        win_start = max(0, first_v - PLA2_FLANKING)
        win_end = min(len(genes) - 1, last_v + PLA2_FLANKING)
        window_genes = genes[win_start:win_end + 1]

        species = data["species"] or scaffold_id
        n_venom = sum(1 for g in window_genes if g["is_venom"])

        windows.append({
            "species": species,
            "scaffold": scaffold_id,
            "n_venom": n_venom,
            "loci": [{
                "n": g["name"],
                "s": g["start"],
                "e": g["end"],
                "d": g["strand"],
                "v": g["is_venom"],
            } for g in window_genes],
        })

    return windows


# ── Data loaders ─────────────────────────────────────────────────────────────

def load_3ftx_data():
    """Load 3FTx proteins + ProtT5 embeddings. Filter to toxin groups only."""
    # Annotations
    proteins = []
    with open(config.THRFTX_DATASET, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            identifier = row.get("identifier", "").strip()
            if not identifier:
                continue
            group = row.get("group", "").strip()
            if group not in TOXIN_GROUPS:
                continue  # Drop all Ly6 groups

            species_name = row.get("species", "").strip()
            genus = row.get("genus", "").strip()
            family = row.get("family", "").strip()

            taxon = row.get("taxon_of_interest", "").strip()
            proteins.append({
                "identifier": identifier,
                "name": identifier.split("|")[-1] if "|" in identifier else identifier,
                "description": f"3FTx {group}",
                "species": species_name,
                "taxon_group": taxon if taxon else family,
                "venom_family": "3FTx",
                "venom_group": f"3FTx {group}",
                "source_dataset": "3FTx (Koludarov 2023, Nat Comms)",
            })

    # Embeddings from zip
    sd4_zip = config.THRFTX_DIR / "embeddings" / "SD4_Protein_embeddings.zip"
    embeddings = {}
    if sd4_zip.exists():
        with zipfile.ZipFile(sd4_zip) as z:
            with z.open("SM8_Protein_embeddings/embs_3ftx.h5") as zf:
                with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
                    tmp.write(zf.read())
                    tmp.flush()
                    with h5py.File(tmp.name, "r") as f:
                        for key in f.keys():
                            arr = np.array(f[key])
                            # Shape is (1, 1024) — squeeze to (1024,)
                            if arr.ndim == 2:
                                arr = arr.mean(axis=0)
                            embeddings[key] = arr.astype(np.float32)

    print(f"  3FTx: {len(proteins)} toxin proteins (Ly6 filtered out), "
          f"{len(embeddings)} embeddings available")
    return proteins, embeddings


def load_sp_data():
    """Load KLK/SP serine protease proteins + ProtT5 embeddings."""
    # Load UniProt/prefix metadata
    sp_meta = {}
    meta_tsv = config.SP_DIR / "sp_uniprot_metadata.tsv"
    if meta_tsv.exists():
        with open(meta_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                sp_meta[row["accession"]] = {
                    "protein_name": row.get("protein_name", ""),
                    "species": row.get("species", ""),
                }

    # Parse FASTA for identifiers
    proteins = []
    fasta = config.SP_DIR / "sequences" / "SP_complete.fasta"
    if not fasta.exists():
        print(f"  WARNING: SP FASTA not found: {fasta}")
        return proteins, {}

    current_id = None
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    meta = sp_meta.get(current_id, {})
                    species = meta.get("species", "")
                    pname = meta.get("protein_name", "")
                    # Clean species: remove parenthetical synonyms
                    if "(" in species:
                        species = species.split("(")[0].strip()
                    proteins.append({
                        "identifier": current_id,
                        "name": pname if pname and pname != current_id else current_id,
                        "description": "Serine protease (KLK family)",
                        "species": species,
                        "taxon_group": "Reptilia" if not species else "",
                        "venom_family": "Serine protease",
                        "venom_group": "KLK",
                        "source_dataset": "KLK/SP (Barua & Koludarov 2021, BMC Biology)",
                    })
                current_id = line[1:].split()[0]
    if current_id:
        meta = sp_meta.get(current_id, {})
        species = meta.get("species", "")
        pname = meta.get("protein_name", "")
        if "(" in species:
            species = species.split("(")[0].strip()
        proteins.append({
            "identifier": current_id,
            "name": pname if pname and pname != current_id else current_id,
            "description": "Serine protease (KLK family)",
            "species": species,
            "taxon_group": "Reptilia" if not species else "",
            "venom_family": "Serine protease",
            "venom_group": "KLK",
            "source_dataset": "KLK/SP (Barua & Koludarov 2021, BMC Biology)",
        })

    # ProtT5 embeddings (per-protein, 1024-dim)
    emb_path = config.SP_DIR / "sp_prott5.h5"
    embeddings = {}
    if emb_path.exists():
        with h5py.File(emb_path, "r") as f:
            for key in f.keys():
                arr = np.array(f[key])
                if arr.ndim == 2:
                    arr = arr.mean(axis=0)
                embeddings[key] = arr.astype(np.float32)
        print(f"  SP: {len(proteins)} proteins, {len(embeddings)} ProtT5 embeddings")
    else:
        print(f"  SP: {len(proteins)} proteins, NO embeddings (file missing)")

    return proteins, embeddings


def load_pla2_data():
    """Load PLA2 proteins + ProtT5 embeddings from annotations.csv + h5."""
    se_root = Path(__file__).resolve().parent.parent.parent
    ann_path = se_root / "data" / "Pla2" / "annotations.csv"
    emb_path = se_root / "data" / "Pla2" / "pla2_prott5.h5"

    # Annotations
    proteins = []
    with open(ann_path, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            identifier = row.get("identifier", "").strip()
            if not identifier:
                continue
            gene = row.get("Gene", "").strip()
            gene_group = row.get("Gene Group", "").strip()
            clade = row.get("Clade", "").strip()
            species = row.get("Species", "").strip()

            proteins.append({
                "identifier": identifier,
                "name": gene if gene else identifier,
                "description": f"Pla2g2 {gene_group}" if gene_group else "Pla2g2",
                "species": species,
                "taxon_group": clade if clade else "",
                "venom_family": "PLA2",
                "venom_group": f"Pla2g2 {gene_group}" if gene_group else "Pla2g2",
                "source_dataset": "Pla2g2 (Koludarov 2020, bioRxiv)",
            })

    # Embeddings (per-protein, 1024-dim — already mean-pooled)
    embeddings = {}
    if emb_path.exists():
        with h5py.File(emb_path, "r") as f:
            for key in f.keys():
                arr = np.array(f[key])
                if arr.ndim == 2:
                    arr = arr.mean(axis=0)
                embeddings[key] = arr.astype(np.float32)
        print(f"  PLA2: {len(proteins)} proteins, {len(embeddings)} ProtT5 embeddings")
    else:
        print(f"  PLA2: {len(proteins)} proteins, NO embeddings (file missing)")

    return proteins, embeddings


# ── Main pipeline ────────────────────────────────────────────────────────────

def main():
    DATA_DIR.mkdir(exist_ok=True)
    print("=" * 60)
    print("Building VenomsBase — Snake Venoms (unified viewer)")
    print("=" * 60)

    # 1. Load all three datasets
    print("\n1. Loading datasets...")
    print("   3FTx (toxin groups only)...")
    ftx_proteins, ftx_emb = load_3ftx_data()

    print("   Serine proteases (KLK)...")
    sp_proteins, sp_emb = load_sp_data()

    print("   PLA2 phospholipases...")
    pla2_proteins, pla2_emb = load_pla2_data()

    # 2. Merge all proteins
    print("\n2. Merging datasets...")
    all_proteins = ftx_proteins + sp_proteins + pla2_proteins
    all_embeddings = {}
    all_embeddings.update(ftx_emb)
    all_embeddings.update(sp_emb)
    all_embeddings.update(pla2_emb)
    print(f"  Total: {len(all_proteins)} proteins, {len(all_embeddings)} embeddings")

    # Match proteins to embeddings
    matched = [(p, all_embeddings[p["identifier"]])
               for p in all_proteins if p["identifier"] in all_embeddings]
    unmatched = [p for p in all_proteins if p["identifier"] not in all_embeddings]
    print(f"  Matched: {len(matched)} proteins with embeddings")
    if unmatched:
        families = defaultdict(int)
        for p in unmatched:
            families[p["venom_family"]] += 1
        print(f"  Unmatched: {len(unmatched)} — " +
              ", ".join(f"{k}: {v}" for k, v in families.items()))

    # 3. Compute projections (UMAP + PCA)
    print("\n3. Computing projections on merged embedding matrix...")
    emb_matrix = np.array([m[1] for m in matched])
    print(f"  Matrix shape: {emb_matrix.shape}")

    from sklearn.decomposition import PCA
    pca = PCA(n_components=2, random_state=42)
    pca_2d = pca.fit_transform(emb_matrix)
    print(f"  PCA done (var explained: {pca.explained_variance_ratio_.sum():.2%})")

    import umap
    reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1,
                        random_state=42, metric="cosine")
    umap_2d = reducer.fit_transform(emb_matrix)
    print(f"  UMAP done")

    # 4. Build annotation DataFrame
    print("\n4. Building annotation table...")
    prot_list = [m[0] for m in matched]
    records = []
    for p in prot_list:
        rec = {k: str(v)[:253] if v else "" for k, v in p.items()}
        rec.pop("sequence", None)
        records.append(rec)

    ann_df = pd.DataFrame(records)
    ann_df = ann_df.fillna("")  # clean NaN for JSON serialization

    umap_df = pd.DataFrame({
        "identifier": [p["identifier"] for p in prot_list],
        "x": umap_2d[:, 0].astype(float),
        "y": umap_2d[:, 1].astype(float),
    }).set_index("identifier")

    pca_df = pd.DataFrame({
        "identifier": [p["identifier"] for p in prot_list],
        "x": pca_2d[:, 0].astype(float),
        "y": pca_2d[:, 1].astype(float),
    }).set_index("identifier")

    # Summary per family
    for fam in ["3FTx", "Serine protease", "PLA2"]:
        n = sum(1 for p in prot_list if p["venom_family"] == fam)
        print(f"  {fam}: {n} proteins in scatter")

    # 5. Load genomic context
    print("\n5. Loading genomic context...")

    # SP genomic context (Thamnophis + Crotalus GFFs)
    sp_windows = load_all_genomic_context(config.SP_DIR)
    sp_venom = sum(w["n_venom"] for w in sp_windows)
    sp_total = sum(len(w["loci"]) for w in sp_windows)
    print(f"  SP: {len(sp_windows)} windows, {sp_total} genes ({sp_venom} venom)")

    # PLA2 genomic context (SM5 GFF)
    se_root = Path(__file__).resolve().parent.parent.parent
    sm5_path = se_root / "projects" / "Pla2" / "SM5_AnnotationsPla2g2.txt"
    pla2_windows = []
    if sm5_path.exists():
        pla2_windows = parse_pla2_sm5_gff(sm5_path)
        pla2_venom = sum(w["n_venom"] for w in pla2_windows)
        pla2_total = sum(len(w["loci"]) for w in pla2_windows)
        print(f"  PLA2: {len(pla2_windows)} windows, {pla2_total} genes ({pla2_venom} venom)")
    else:
        print(f"  PLA2 SM5 GFF not found: {sm5_path}")

    # Tag windows with family for genomic context matching
    for w in sp_windows:
        w["family"] = "Serine protease"
    for w in pla2_windows:
        w["family"] = "PLA2"

    # Merge all genomic windows
    all_windows = sp_windows + pla2_windows
    print(f"  Total: {len(all_windows)} genomic windows")

    # Serialize windows to compact JSON format
    gw_json = []
    for win in all_windows:
        if win["loci"] and isinstance(win["loci"][0], dict) and "n" in win["loci"][0]:
            # Already in compact format (from PLA2 parser)
            gw_json.append(win)
        else:
            gw_json.append({
                "species": win["species"],
                "scaffold": win["scaffold"],
                "n_venom": win["n_venom"],
                "family": win.get("family", ""),
                "loci": [{"n": l["name"], "s": l["start"], "e": l["end"],
                           "d": l["strand"], "v": l["is_venom"]} for l in win["loci"]],
            })

    # 6. Build HTML via 03's build_html function
    print("\n6. Building HTML...")
    output = DATA_DIR / "VenomsBase_Snakes.html"

    spec = importlib.util.spec_from_file_location(
        "viewer", str(Path(__file__).resolve().parent / "03_build_html_viewer.py"))
    viewer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(viewer)

    # Build color map from our harmonized palettes
    color_map = {}
    color_map.update(FAMILY_COLORS)
    color_map.update(ALL_GROUP_COLORS)

    viewer.build_html(ann_df, umap_df, pca_df, gw_json, output,
                      color_map=color_map,
                      title="VenomsBase \u2014 Snake Venoms")

    # Summary
    print("\n" + "=" * 60)
    print("Done! Viewer ready:")
    print(f"  open {output}")
    print(f"\n  Proteins: {len(prot_list)}")
    for fam in ["3FTx", "Serine protease", "PLA2"]:
        n = sum(1 for p in prot_list if p["venom_family"] == fam)
        print(f"    {fam}: {n}")
    print(f"  Genomic windows: {len(all_windows)}")
    print(f"  Color-by columns: venom_family ({len(FAMILY_COLORS)} colors), "
          f"venom_group ({len(ALL_GROUP_COLORS)} colors)")
    print("=" * 60)


if __name__ == "__main__":
    main()
