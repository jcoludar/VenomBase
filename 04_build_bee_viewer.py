#!/usr/bin/env python3
"""Step 04: Build bee-only VenomsBase viewer.

Extracts bee venom protein sequences from aligned FASTAs + ToxProt reference,
generates ProtT5 embeddings via BiocentralAPI, and builds a self-contained HTML.

Output: data/VenomsBase_Bees.html
"""

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "tools"))
import config

DATA_DIR = Path(__file__).resolve().parent / "data"


def extract_bee_sequences():
    """Get sequences for bee venom proteins from aligned FASTAs + ToxProt."""
    seqs = {}  # id -> {name, sequence, family, source}

    # 1. ToxProt reference (clean, curated, has UniProt IDs)
    toxprot = config.BEE_VENOM_DIR / "uniprot-reviewed_Apoidea_291020.fasta"
    if toxprot.exists():
        current_id = current_header = None
        current_seq = []
        with open(toxprot) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        seq = "".join(current_seq).replace("-", "")
                        # Parse sp|P01500|APAM_APIME Apamin OS=...
                        parts = current_header.split("|")
                        uid = parts[1] if len(parts) > 1 else current_id
                        name_part = current_header.split(" ", 1)
                        desc = name_part[1].split(" OS=")[0] if len(name_part) > 1 else ""
                        seqs[uid] = {"name": uid, "sequence": seq, "description": desc,
                                     "family": "", "source": "ToxProt"}
                    current_header = line[1:]
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id:
                seq = "".join(current_seq).replace("-", "")
                parts = current_header.split("|")
                uid = parts[1] if len(parts) > 1 else current_id
                desc = current_header.split(" ", 1)[1].split(" OS=")[0] if " " in current_header else ""
                seqs[uid] = {"name": uid, "sequence": seq, "description": desc,
                             "family": "", "source": "ToxProt"}

    # 2. Aligned FASTAs (one per family, contains cross-species sequences)
    aligned_dir = config.BEE_ALIGNED
    if aligned_dir.exists():
        for fasta_file in sorted(aligned_dir.glob("*.fasta")):
            family = fasta_file.stem.replace("_cleaned", "")
            current_id = None
            current_seq = []
            with open(fasta_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if current_id:
                            seq = "".join(current_seq).replace("-", "").replace(".", "")
                            if len(seq) > 10 and current_id not in seqs:
                                seqs[current_id] = {"name": current_id, "sequence": seq,
                                                    "description": family, "family": family,
                                                    "source": "aligned"}
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_id:
                    seq = "".join(current_seq).replace("-", "").replace(".", "")
                    if len(seq) > 10 and current_id not in seqs:
                        seqs[current_id] = {"name": current_id, "sequence": seq,
                                            "description": family, "family": family,
                                            "source": "aligned"}

    return seqs


def match_sequences_to_annotations(seqs, annotations_csv):
    """Match extracted sequences to the master annotations CSV."""
    annotations = {}
    with open(annotations_csv, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            entry = (row.get("Entry") or row.get("\ufeffEntry", "")).strip()
            if entry:
                annotations[entry] = row

    # Build unified protein list
    proteins = []
    matched = 0
    for sid, sdata in seqs.items():
        p = {
            "identifier": sid,
            "name": sdata["name"],
            "description": sdata["description"],
            "sequence": sdata["sequence"],
            "sequence_length": len(sdata["sequence"]),
            "venom_family": sdata["family"],
            "source": sdata["source"],
        }
        # Try to match to annotations
        ann = annotations.get(sid)
        if ann:
            matched += 1
            p["species"] = ann.get("Species", "")
            p["hymenoptera_group"] = ann.get("Hymenoptera group", "")
            p["protein_names"] = ann.get("Protein names", "")
            p["venom_family"] = ann.get("Protein family/Gene name", "") or sdata["family"]
        proteins.append(p)

    print(f"  {len(proteins)} proteins, {matched} matched to annotations")
    return proteins


def generate_embeddings(proteins):
    """Generate ProtT5 embeddings via BiocentralAPI."""
    try:
        from pipelines.protspace_pipeline import _generate_embeddings_api
    except ImportError:
        pass

    # Write temp FASTA
    import tempfile
    fasta_path = Path(tempfile.mktemp(suffix=".fasta"))
    with open(fasta_path, "w") as f:
        for p in proteins:
            if p.get("sequence"):
                f.write(f">{p['identifier']}\n{p['sequence']}\n")

    # Try BiocentralAPI
    try:
        from biocentral_api import BiocentralAPI
        api = BiocentralAPI()
        print("  Generating embeddings via BiocentralAPI...")
        result = api.embed(str(fasta_path), model="prott5")
        import h5py
        embeddings = {}
        with h5py.File(result, "r") as hf:
            for key in hf.keys():
                arr = np.array(hf[key])
                if arr.ndim == 2:
                    arr = arr.mean(axis=0)
                embeddings[key] = arr.astype(np.float32)
        fasta_path.unlink(missing_ok=True)
        return embeddings
    except Exception as e:
        print(f"  BiocentralAPI failed: {e}")
        print("  Using sequence-length PCA as fallback (replace with real embeddings later)")
        fasta_path.unlink(missing_ok=True)

    # Fallback: simple sequence features for projection
    from sklearn.decomposition import PCA
    # Compute amino acid composition as feature vector
    AA = "ACDEFGHIKLMNPQRSTVWY"
    features = []
    ids = []
    for p in proteins:
        seq = p.get("sequence", "")
        if not seq:
            continue
        comp = [seq.count(aa) / max(len(seq), 1) for aa in AA]
        comp.append(len(seq) / 1000.0)  # normalized length
        features.append(comp)
        ids.append(p["identifier"])

    features = np.array(features, dtype=np.float32)
    return {pid: features[i] for i, pid in enumerate(ids)}


def build_viewer_html(proteins, embeddings, genomic_windows, output_path, title="VenomsBase"):
    """Build self-contained HTML viewer."""
    from sklearn.decomposition import PCA

    # Match proteins to embeddings, compute projections
    matched = [(p, embeddings[p["identifier"]]) for p in proteins if p["identifier"] in embeddings]
    print(f"  {len(matched)} proteins with embeddings")

    if len(matched) < 5:
        print("  ERROR: Too few proteins with embeddings")
        return

    emb_matrix = np.array([m[1] for m in matched])
    prot_list = [m[0] for m in matched]

    pca = PCA(n_components=2)
    pca_2d = pca.fit_transform(emb_matrix)

    try:
        import umap
        reducer = umap.UMAP(n_components=2, n_neighbors=min(15, len(matched) - 1),
                            min_dist=0.1, random_state=42)
        umap_2d = reducer.fit_transform(emb_matrix)
    except Exception:
        umap_2d = pca_2d

    # Build JSON data
    json_proteins = []
    for i, p in enumerate(prot_list):
        jp = {k: str(v)[:253] for k, v in p.items() if v and k != "sequence"}
        jp["umap_x"] = float(umap_2d[i, 0])
        jp["umap_y"] = float(umap_2d[i, 1])
        jp["pca_x"] = float(pca_2d[i, 0])
        jp["pca_y"] = float(pca_2d[i, 1])
        # Truncate sequence for display
        seq = p.get("sequence", "")
        if seq:
            jp["sequence"] = seq[:253] + "..." if len(seq) > 256 else seq
        json_proteins.append(jp)

    # Color columns
    color_cols = []
    for col in ["venom_family", "hymenoptera_group", "species", "source"]:
        vals = set(p.get(col, "") for p in json_proteins if p.get(col))
        if 2 <= len(vals) <= 50:
            color_cols.append(col)

    # Read the HTML template from 03 and inject our data
    # (reuse the same template structure)
    from importlib.util import spec_from_file_location, module_from_spec
    viewer_script = Path(__file__).resolve().parent / "03_build_html_viewer.py"

    # Just generate inline — simpler than importing
    proteins_json = json.dumps(json_proteins, default=str)
    color_columns_json = json.dumps(color_cols)
    genomic_json_str = json.dumps(genomic_windows)

    # Read the template from 03's output format
    template_path = Path(__file__).resolve().parent / "03_build_html_viewer.py"
    # Actually, let's just read the last generated HTML and replace the data
    # Simpler: re-import the build function
    # Simplest: just call 03's build_html with our data

    # Build a DataFrame matching what build_html expects
    ann_df = pd.DataFrame(json_proteins)
    ann_df = ann_df.drop(columns=["umap_x", "umap_y", "pca_x", "pca_y"], errors="ignore")

    umap_df = pd.DataFrame({
        "identifier": [p["identifier"] for p in json_proteins],
        "x": [p["umap_x"] for p in json_proteins],
        "y": [p["umap_y"] for p in json_proteins],
    }).set_index("identifier")

    pca_df = pd.DataFrame({
        "identifier": [p["identifier"] for p in json_proteins],
        "x": [p["pca_x"] for p in json_proteins],
        "y": [p["pca_y"] for p in json_proteins],
    }).set_index("identifier")

    # Import and call the HTML builder from 03
    import importlib.util
    spec = importlib.util.spec_from_file_location("viewer", str(viewer_script))
    viewer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(viewer)
    viewer.build_html(ann_df, umap_df, pca_df, genomic_windows, output_path)


def main():
    DATA_DIR.mkdir(exist_ok=True)

    print("Building VenomsBase BEE viewer...")

    # Extract sequences
    print("\n1. Extracting bee venom sequences...")
    seqs = extract_bee_sequences()
    print(f"  {len(seqs)} unique sequences")

    # Match to annotations
    print("\n2. Matching to annotations...")
    proteins = match_sequences_to_annotations(seqs, config.BEE_ANNOTATIONS)

    # Generate embeddings
    print("\n3. Generating embeddings...")
    embeddings = generate_embeddings(proteins)
    print(f"  {len(embeddings)} embeddings")

    # No bee genomic context yet (would need Apis GFF)
    genomic_windows = []

    # Build HTML
    print("\n4. Building HTML...")
    output = DATA_DIR / "VenomsBase_Bees.html"
    build_viewer_html(proteins, embeddings, genomic_windows, output, "VenomsBase — Bee Venoms")

    print(f"\nDone! open {output}")


if __name__ == "__main__":
    main()
