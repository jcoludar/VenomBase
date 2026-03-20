#!/usr/bin/env python3
"""Step 02: Build VenomsBase prototype ProtSpace bundle.

Creates a .parquetbundle from the unified dataset that can be viewed
in ProtSpace. Uses existing ProtT5 embeddings where available,
generates UMAP projections.

The bundle demonstrates the embedding-centric architecture:
- Each protein is a point in embedding space
- All metadata are tag columns (annotation dimensions)
- Click any point → see all its tags

Inputs:
  - data/venombase_unified.csv (from step 01)
  - 3FTx embeddings (SD4: embs_3ftx.h5)

Output:
  - data/VenomsBase_Prototype.parquetbundle
"""

import csv
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

OUTPUT_DIR = Path(__file__).resolve().parent / "data"


def load_3ftx_embeddings():
    """Extract and load 3FTx ProtT5 embeddings from SD4 zip."""
    sd4_zip = config.THRFTX_DIR / "embeddings" / "SD4_Protein_embeddings.zip"
    if not sd4_zip.exists():
        print("  Warning: SD4 embeddings not found")
        return {}

    embeddings = {}
    with zipfile.ZipFile(sd4_zip) as z:
        # Extract the H5 file
        h5_name = "SM8_Protein_embeddings/embs_3ftx.h5"
        with z.open(h5_name) as zf:
            with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
                tmp.write(zf.read())
                tmp.flush()
                with h5py.File(tmp.name, "r") as f:
                    for key in f.keys():
                        arr = np.array(f[key])
                        # Per-residue embeddings (seq_len, features) → mean pool to (features,)
                        if arr.ndim == 2:
                            arr = arr.mean(axis=0)
                        embeddings[key] = arr.astype(np.float32)

    sample = next(iter(embeddings.values())) if embeddings else None
    print(f"  Loaded {len(embeddings)} 3FTx embeddings (dim={sample.shape[0] if sample is not None else '?'})")
    return embeddings


def build_bundle(annotations_df, embeddings_dict, output_path):
    """Build a .parquetbundle file from annotations + embeddings.

    Format: 3 Parquet tables separated by delimiter.
    """
    from sklearn.decomposition import PCA

    delimiter = b"---PARQUET_DELIMITER---"

    # Match annotations to embeddings
    matched_ids = [idx for idx in annotations_df["identifier"] if idx in embeddings_dict]
    print(f"  Matched {len(matched_ids)} proteins with embeddings out of {len(annotations_df)}")

    if len(matched_ids) < 10:
        print("  Too few matched embeddings. Building annotation-only bundle with random projections for demo.")
        # Use random projections for demo
        n = len(annotations_df)
        np.random.seed(42)
        umap_2d = np.random.randn(n, 2) * 5
        pca_2d = np.random.randn(n, 2) * 3
        matched_ids = list(annotations_df["identifier"])
    else:
        # Build embedding matrix for matched IDs
        emb_matrix = np.array([embeddings_dict[idx] for idx in matched_ids])
        annotations_df = annotations_df[annotations_df["identifier"].isin(matched_ids)].copy()

        # PCA projection
        print("  Computing PCA...")
        pca = PCA(n_components=2)
        pca_2d = pca.fit_transform(emb_matrix)

        # UMAP projection
        try:
            import umap
            print("  Computing UMAP...")
            reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
            umap_2d = reducer.fit_transform(emb_matrix)
        except ImportError:
            print("  UMAP not available, using PCA only")
            umap_2d = pca_2d

    # Truncate strings to 256 chars (ProtSpace constraint)
    for col in annotations_df.columns:
        if annotations_df[col].dtype == object:
            annotations_df[col] = annotations_df[col].fillna("").astype(str).apply(
                lambda x: x[:253] + "..." if len(x) > 256 else x
            )

    # Table 1: Annotations
    ann_table = annotations_df.reset_index(drop=True)

    # Table 2: Projection metadata
    proj_meta = pd.DataFrame([
        {"projection_name": "UMAP_2", "dimensions": 2, "info_json": "{}"},
        {"projection_name": "PCA_2", "dimensions": 2, "info_json": "{}"},
    ])

    # Table 3: Projection data
    ids = list(annotations_df["identifier"])
    proj_rows = []
    for i, pid in enumerate(ids):
        proj_rows.append({
            "projection_name": "UMAP_2",
            "identifier": pid,
            "x": float(umap_2d[i, 0]),
            "y": float(umap_2d[i, 1]),
            "z": 0.0,
        })
        proj_rows.append({
            "projection_name": "PCA_2",
            "identifier": pid,
            "x": float(pca_2d[i, 0]),
            "y": float(pca_2d[i, 1]),
            "z": 0.0,
        })
    proj_data = pd.DataFrame(proj_rows)

    # Write bundle
    import io

    parts = []
    for table in [ann_table, proj_meta, proj_data]:
        buf = io.BytesIO()
        table.to_parquet(buf, index=False)
        parts.append(buf.getvalue())

    with open(output_path, "wb") as f:
        f.write(parts[0])
        f.write(delimiter)
        f.write(parts[1])
        f.write(delimiter)
        f.write(parts[2])

    size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"  Written: {output_path} ({size_mb:.1f} MB)")
    print(f"  Proteins: {len(annotations_df)}")
    print(f"  Tag columns: {len(annotations_df.columns)}")


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("Building VenomsBase prototype bundle...")

    # Load unified dataset
    print("\n1. Loading unified dataset...")
    df = pd.read_csv(OUTPUT_DIR / "venombase_unified.csv", dtype=str).fillna("")
    print(f"  {len(df)} proteins, {len(df.columns)} columns")

    # Load existing embeddings
    print("\n2. Loading embeddings...")
    embeddings = load_3ftx_embeddings()

    # Build bundle
    print("\n3. Building ProtSpace bundle...")
    output_path = OUTPUT_DIR / "VenomsBase_Prototype.parquetbundle"
    build_bundle(df, embeddings, output_path)

    print("\nDone! Open in ProtSpace to explore.")
    print(f"  protspace {output_path}")


if __name__ == "__main__":
    main()
