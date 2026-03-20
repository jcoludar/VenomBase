#!/usr/bin/env python3
"""Step 01: Build unified VenomsBase protein dataset.

Merges published datasets into a single annotations table with standardized
tag columns. Each row = one protein. Each column = one tag dimension.

Inputs:
  - Bee venom: Hymenoptera_Venoms_CVP_Feb_2022_372.csv (372 proteins)
  - 3FTx/Ly6: SD1_Dataset_and_information.csv (1,426 proteins)
  - Serine proteases: aligned sequence headers (~500 proteins)
  - Bee venom activities: SupplementaryTable4.xlsx

Output:
  - data/venombase_unified.csv — master protein table
"""

import csv
import sys
from pathlib import Path

# Project imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

REPO_ROOT = Path(__file__).resolve().parent
OUTPUT = REPO_ROOT / "data"


def load_bee_venom():
    """Load bee venom evolution dataset (Koludarov et al. 2023, BMC Biology)."""
    proteins = []
    with open(config.BEE_ANNOTATIONS, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            entry = row.get("Entry") or row.get("\ufeffEntry", "")
            if not entry:
                continue
            proteins.append({
                "identifier": entry.strip(),
                "name": row.get("Entry name", "").strip(),
                "description": row.get("Protein names", "").strip(),
                "taxonomy.species": row.get("Species", "").strip(),
                "taxonomy.group": row.get("Hymenoptera group", "").strip(),
                "venom.family": row.get("Protein family/Gene name", "").strip(),
                "venom.mvp_group": row.get("mvp members", "").strip(),
                "sequence.length": row.get("Length", "").strip(),
                "venom.status": "curated",
                "source.dataset": "bee_venom_evolution",
                "source.doi": "10.1186/s12915-023-01656-5",
                "source.db": "UniProtKB/ToxProt",
                "taxonomy.clade": "Hymenoptera",
            })
    return proteins


def load_3ftx():
    """Load 3FTx/Ly6 dataset (Koludarov et al. 2023, Nature Comms)."""
    proteins = []
    with open(config.THRFTX_DATASET, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            identifier = row.get("identifier", "").strip()
            if not identifier:
                continue

            # Determine venom status from group
            group = row.get("group", "").strip()
            if "3FTx" in group or "toxin" in group.lower():
                venom_status = "curated"
            elif "Ly6" in group:
                venom_status = "non-venom_paralog"
            else:
                venom_status = "candidate"

            proteins.append({
                "identifier": identifier,
                "name": identifier.split("|")[-1] if "|" in identifier else identifier,
                "description": group,
                "taxonomy.species": row.get("species", "").strip(),
                "taxonomy.family": row.get("family", "").strip(),
                "taxonomy.genus": row.get("genus", "").strip(),
                "taxonomy.taxon_id": row.get("taxon_id", "").strip(),
                "taxonomy.group": row.get("taxon_of_interest", "").strip(),
                "venom.family": "3FTx" if "3FTx" in group else "Ly6",
                "venom.group": group,
                "venom.status": venom_status,
                "sequence.full": row.get("full_seq", "").strip(),
                "sequence.mature": row.get("mature_seq", "").strip(),
                "sequence.length": str(len(row.get("full_seq", "").strip())),
                "structure.pdb_id": row.get("pdb_id", "").strip(),
                "structure.lddt": row.get("LDDT", "").strip(),
                "structure.tm_score": row.get("TM_score", "").strip(),
                "structure.rmsd": row.get("RMSD", "").strip(),
                "prediction.membrane": row.get("membran_prediction", "").strip(),
                "identity.uniprot_id": row.get("uniprot_id", "").strip(),
                "identity.genbank_id": row.get("genbank_id", "").strip(),
                "identity.refseq_id": row.get("refseq_id", "").strip(),
                "source.dataset": "3ftx_evolution",
                "source.doi": "10.1038/s41467-023-36SEQ",  # placeholder
                "source.db": row.get("db", "").strip(),
                "taxonomy.clade": "Squamata" if "Serpentes" in row.get("taxon_of_interest", "") or "snake" in group.lower() else "Vertebrata",
            })
    return proteins


def load_serine_protease():
    """Load serine protease dataset (Barua*, Koludarov* 2021, BMC Biology).

    Parse from FASTA headers since we don't have a clean CSV.
    """
    proteins = []
    fasta_path = config.SP_DIR / "sequences" / "SP_complete.fasta"
    if not fasta_path.exists():
        print(f"  Warning: {fasta_path} not found, skipping SP dataset")
        return proteins

    current_id = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    seq = "".join(current_seq)
                    proteins.append({
                        "identifier": f"SP_{current_id}",
                        "name": current_id,
                        "description": "Serine protease (KLK/SVL family)",
                        "venom.family": "Serine protease",
                        "venom.group": "KLK/SVL",
                        "venom.status": "curated",
                        "sequence.full": seq,
                        "sequence.length": str(len(seq)),
                        "source.dataset": "serine_protease_evolution",
                        "source.doi": "10.1186/s12915-021-01191-1",
                        "taxonomy.clade": "Vertebrata",
                    })
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    # Last entry
    if current_id and current_seq:
        seq = "".join(current_seq)
        proteins.append({
            "identifier": f"SP_{current_id}",
            "name": current_id,
            "description": "Serine protease (KLK/SVL family)",
            "venom.family": "Serine protease",
            "venom.group": "KLK/SVL",
            "venom.status": "curated",
            "sequence.full": seq,
            "sequence.length": str(len(seq)),
            "source.dataset": "serine_protease_evolution",
            "source.doi": "10.1186/s12915-021-01191-1",
            "taxonomy.clade": "Vertebrata",
        })

    return proteins


def main():
    OUTPUT.mkdir(exist_ok=True)

    print("Building unified VenomsBase dataset...")

    # Load all sources
    print("  Loading bee venom evolution...")
    bee = load_bee_venom()
    print(f"    {len(bee)} proteins")

    print("  Loading 3FTx/Ly6 evolution...")
    ftx = load_3ftx()
    print(f"    {len(ftx)} proteins")

    print("  Loading serine protease evolution...")
    sp = load_serine_protease()
    print(f"    {len(sp)} proteins")

    # Merge
    all_proteins = bee + ftx + sp
    print(f"\n  Total: {len(all_proteins)} proteins")

    # Collect all tag columns
    all_columns = set()
    for p in all_proteins:
        all_columns.update(p.keys())

    # Define column order: identity first, then tags grouped by category
    identity_cols = ["identifier", "name", "description"]
    taxonomy_cols = sorted(c for c in all_columns if c.startswith("taxonomy."))
    venom_cols = sorted(c for c in all_columns if c.startswith("venom."))
    sequence_cols = sorted(c for c in all_columns if c.startswith("sequence."))
    structure_cols = sorted(c for c in all_columns if c.startswith("structure."))
    prediction_cols = sorted(c for c in all_columns if c.startswith("prediction."))
    identity_extra = sorted(c for c in all_columns if c.startswith("identity."))
    source_cols = sorted(c for c in all_columns if c.startswith("source."))

    columns = (
        identity_cols + taxonomy_cols + venom_cols + sequence_cols +
        structure_cols + prediction_cols + identity_extra + source_cols
    )

    # Write unified CSV
    output_path = OUTPUT / "venombase_unified.csv"
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for p in all_proteins:
            writer.writerow(p)

    print(f"\n  Written: {output_path}")
    print(f"  Columns: {len(columns)}")
    print(f"  Rows: {len(all_proteins)}")

    # Stats
    datasets = {}
    for p in all_proteins:
        ds = p.get("source.dataset", "unknown")
        datasets[ds] = datasets.get(ds, 0) + 1
    print("\n  By dataset:")
    for ds, count in sorted(datasets.items()):
        print(f"    {ds}: {count}")

    families = {}
    for p in all_proteins:
        fam = p.get("venom.family", "unknown")
        families[fam] = families.get(fam, 0) + 1
    print(f"\n  Venom families: {len(families)}")
    for fam, count in sorted(families.items(), key=lambda x: -x[1])[:15]:
        print(f"    {fam}: {count}")


if __name__ == "__main__":
    main()
