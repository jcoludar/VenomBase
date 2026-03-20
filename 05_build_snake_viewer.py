#!/usr/bin/env python3
"""Step 05: Build snake-only VenomsBase viewer.

Uses 3FTx/Ly6 dataset (1,427 proteins with ProtT5 embeddings) + serine protease
dataset, with snake genomic context from GFF annotations.

Output: data/VenomsBase_Snakes.html
"""

import csv
import json
import sys
import zipfile
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "tools"))
import config
from src.genomic_context import load_all_genomic_context

DATA_DIR = Path(__file__).resolve().parent / "data"


def load_3ftx_data():
    """Load 3FTx/Ly6 proteins + embeddings."""
    # Annotations
    proteins = []
    with open(config.THRFTX_DATASET, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            identifier = row.get("identifier", "").strip()
            if not identifier:
                continue
            group = row.get("group", "").strip()
            proteins.append({
                "identifier": identifier,
                "name": identifier.split("|")[-1] if "|" in identifier else identifier,
                "description": group,
                "species": row.get("species", "").strip(),
                "family": row.get("family", "").strip(),
                "genus": row.get("genus", "").strip(),
                "taxon_group": row.get("taxon_of_interest", "").strip(),
                "venom_family": "3FTx" if "3FTx" in group else ("Ly6" if "Ly6" in group else group),
                "venom_group": group,
                "venom_status": "curated" if "3FTx" in group or "toxin" in group.lower() else (
                    "non-venom_paralog" if "Ly6" in group else "candidate"),
                "membrane": row.get("membran_prediction", "").strip(),
                "uniprot_id": row.get("uniprot_id", "").strip(),
                "genbank_id": row.get("genbank_id", "").strip(),
                "pdb_id": row.get("pdb_id", "").strip(),
                "lddt": row.get("LDDT", "").strip(),
                "source_dataset": "3FTx (Koludarov 2023, Nat Comms)",
            })

    # Embeddings
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
                            if arr.ndim == 2:
                                arr = arr.mean(axis=0)
                            embeddings[key] = arr.astype(np.float32)

    print(f"  3FTx: {len(proteins)} proteins, {len(embeddings)} embeddings")
    return proteins, embeddings


def load_sp_data():
    """Load serine protease proteins from FASTA."""
    proteins = []
    fasta = config.SP_DIR / "sequences" / "SP_complete.fasta"
    if not fasta.exists():
        return proteins, {}

    current_id = None
    current_seq = []
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    seq = "".join(current_seq).replace("-", "")
                    proteins.append({
                        "identifier": current_id,
                        "name": current_id,
                        "description": "Serine protease (KLK/SVL)",
                        "venom_family": "Serine protease",
                        "venom_group": "KLK/SVL",
                        "venom_status": "curated",
                        "sequence": seq,
                        "source_dataset": "SP (Barua & Koludarov 2021, BMC Biology)",
                    })
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id and current_seq:
        seq = "".join(current_seq).replace("-", "")
        proteins.append({
            "identifier": current_id,
            "name": current_id,
            "description": "Serine protease (KLK/SVL)",
            "venom_family": "Serine protease",
            "venom_group": "KLK/SVL",
            "venom_status": "curated",
            "sequence": seq,
            "source_dataset": "SP (Barua & Koludarov 2021, BMC Biology)",
        })

    # AA composition embeddings as fallback for SP
    AA = "ACDEFGHIKLMNPQRSTVWY"
    embeddings = {}
    for p in proteins:
        seq = p.get("sequence", "")
        if seq:
            comp = [seq.count(aa) / max(len(seq), 1) for aa in AA]
            comp.append(len(seq) / 1000.0)
            embeddings[p["identifier"]] = np.array(comp, dtype=np.float32)

    print(f"  SP: {len(proteins)} proteins (AA-composition embeddings)")
    return proteins, embeddings


def main():
    DATA_DIR.mkdir(exist_ok=True)
    print("Building VenomsBase SNAKE viewer...")

    # Load data
    print("\n1. Loading 3FTx/Ly6...")
    ftx_proteins, ftx_emb = load_3ftx_data()

    print("   Loading serine proteases...")
    sp_proteins, sp_emb = load_sp_data()

    # Only use 3FTx for scatter (they have real ProtT5 embeddings)
    # SP proteins included in data but won't have scatter positions
    # unless we compute real embeddings
    all_proteins = ftx_proteins  # + sp_proteins  # SP lacks real embeddings

    # Load genomic context
    print("\n2. Loading snake genomic context...")
    genomic_windows = load_all_genomic_context(config.SP_DIR)
    n_venom = sum(w["n_venom"] for w in genomic_windows)
    n_total = sum(len(w["loci"]) for w in genomic_windows)
    print(f"  {len(genomic_windows)} windows, {n_total} genes ({n_venom} venom)")

    # Compute projections
    print("\n3. Computing projections...")
    from sklearn.decomposition import PCA

    matched = [(p, ftx_emb[p["identifier"]]) for p in all_proteins if p["identifier"] in ftx_emb]
    print(f"  {len(matched)} proteins with ProtT5 embeddings")

    emb_matrix = np.array([m[1] for m in matched])
    pca = PCA(n_components=2)
    pca_2d = pca.fit_transform(emb_matrix)

    import umap
    reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
    umap_2d = reducer.fit_transform(emb_matrix)

    # Build annotation DataFrame
    prot_list = [m[0] for m in matched]
    records = []
    for i, p in enumerate(prot_list):
        rec = {k: str(v)[:253] if v else "" for k, v in p.items()}
        rec.pop("sequence", None)  # don't include raw sequence in scatter data
        records.append(rec)

    ann_df = pd.DataFrame(records)

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

    # Serialize windows
    gw_json = []
    for win in genomic_windows:
        gw_json.append({
            "species": win["species"],
            "scaffold": win["scaffold"],
            "n_venom": win["n_venom"],
            "loci": [{"n": l["name"], "s": l["start"], "e": l["end"],
                       "d": l["strand"], "v": l["is_venom"]} for l in win["loci"]],
        })

    # Build HTML
    print("\n4. Building HTML...")
    output = DATA_DIR / "VenomsBase_Snakes.html"

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "viewer", str(Path(__file__).resolve().parent / "03_build_html_viewer.py"))
    viewer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(viewer)
    viewer.build_html(ann_df, umap_df, pca_df, gw_json, output)

    print(f"\nDone! open {output}")


if __name__ == "__main__":
    main()
