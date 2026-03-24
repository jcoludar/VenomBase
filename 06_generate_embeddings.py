#!/usr/bin/env python3
"""Step 06: Generate ProtT5 embeddings for VenomBase datasets.

Generates per-protein ProtT5 embeddings for:
  1. Serine protease (KLK/SVL) proteins — degapped from SP_complete.fasta
  2. Bee venom proteins — degapped from aligned FASTAs + ToxProt

Run with ProteEmbedExplorations venv (has torch + transformers):
  /Users/jcoludar/CascadeProjects/ProteEmbedExplorations/.venv/bin/python 06_generate_embeddings.py

Output:
  data/serine_protease/sp_prott5.h5     — ~491 proteins, (1024,) per protein
  data/bee_venom/bee_prott5.h5          — ~4600 proteins, (1024,) per protein
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

import numpy as np
from embeddings import read_fasta, degap, filter_by_length, save_embeddings, mean_pool
from embeddings.prott5 import extract_prott5


def generate_sp_embeddings():
    """Generate ProtT5 embeddings for serine protease proteins."""
    sp_fasta = REPO_ROOT / "data" / "serine_protease" / "sequences" / "SP_complete.fasta"
    output_h5 = REPO_ROOT / "data" / "serine_protease" / "sp_prott5.h5"

    if not sp_fasta.exists():
        print(f"SKIP: {sp_fasta} not found")
        return

    print("=" * 60)
    print("Serine Protease ProtT5 Embeddings")
    print("=" * 60)

    # Load and degap
    seqs = read_fasta(sp_fasta)
    print(f"  Raw: {len(seqs)} sequences")

    seqs = degap(seqs)
    seqs = filter_by_length(seqs, min_length=30)
    print(f"  After degap+filter: {len(seqs)} sequences")

    lengths = [len(s) for s in seqs.values()]
    print(f"  Length range: {min(lengths)}-{max(lengths)} (mean {np.mean(lengths):.0f})")

    # Extract per-residue, then mean-pool
    per_residue = extract_prott5(seqs, batch_size=4)
    per_protein = mean_pool(per_residue)

    save_embeddings(per_protein, output_h5)
    print(f"  Output: {output_h5}")
    return per_protein


def generate_bee_embeddings():
    """Generate ProtT5 embeddings for bee venom proteins."""
    bee_dir = REPO_ROOT / "data" / "bee_venom"
    aligned_dir = bee_dir / "aligned_sequences"
    toxprot_fasta = bee_dir / "uniprot-reviewed_Apoidea_291020.fasta"
    output_h5 = bee_dir / "bee_prott5.h5"

    print("\n" + "=" * 60)
    print("Bee Venom ProtT5 Embeddings")
    print("=" * 60)

    all_seqs: dict[str, str] = {}

    # 1. ToxProt reference (curated, has UniProt IDs)
    if toxprot_fasta.exists():
        toxprot = read_fasta(toxprot_fasta)
        # Extract clean UniProt IDs from sp|P01500|APAM_APIME format
        clean_toxprot = {}
        for sid, seq in toxprot.items():
            parts = sid.split("|")
            uid = parts[1] if len(parts) > 1 else sid
            clean_toxprot[uid] = seq
        all_seqs.update(clean_toxprot)
        print(f"  ToxProt: {len(clean_toxprot)} sequences")

    # 2. Aligned FASTAs (one per family)
    if aligned_dir.exists():
        for fasta_file in sorted(aligned_dir.glob("*.fasta")):
            family_seqs = read_fasta(fasta_file)
            family_seqs = degap(family_seqs)
            # Only add new IDs (ToxProt takes priority)
            n_new = 0
            for sid, seq in family_seqs.items():
                if sid not in all_seqs and len(seq) > 10:
                    all_seqs[sid] = seq
                    n_new += 1
            if n_new:
                print(f"  {fasta_file.stem}: +{n_new} new")

    print(f"  Total unique: {len(all_seqs)} sequences")

    # Filter
    all_seqs = filter_by_length(all_seqs, min_length=15)
    print(f"  After length filter: {len(all_seqs)} sequences")

    lengths = [len(s) for s in all_seqs.values()]
    print(f"  Length range: {min(lengths)}-{max(lengths)} (mean {np.mean(lengths):.0f})")

    # Extract per-residue, then mean-pool
    per_residue = extract_prott5(all_seqs, batch_size=4)
    per_protein = mean_pool(per_residue)

    save_embeddings(per_protein, output_h5)
    print(f"  Output: {output_h5}")
    return per_protein


def main():
    print("VenomBase — ProtT5 Embedding Generation")
    print(f"Using tools from: {REPO_ROOT / 'tools' / 'embeddings'}\n")

    sp_embs = generate_sp_embeddings()
    bee_embs = generate_bee_embeddings()

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    if sp_embs:
        print(f"  SP:  {len(sp_embs)} proteins × {next(iter(sp_embs.values())).shape}")
    if bee_embs:
        print(f"  Bee: {len(bee_embs)} proteins × {next(iter(bee_embs.values())).shape}")
    print("Done!")


if __name__ == "__main__":
    main()
