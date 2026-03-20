"""Build genomic context annotations for VenomsBase proteins.

Parses GFF BLAST-hit data into clean gene loci, identifies venom genes,
and creates windowed genomic neighborhoods to bake into protein records.
"""

from collections import defaultdict
from pathlib import Path


# Gene name prefixes that indicate venom serine proteases
VENOM_PREFIXES = [
    "Prmu", "SVL", "SVsp", "Crvi", "Thsi", "Helo",
    "Blbr", "Sopa", "Nosc", "Pste", "Geja", "Anca",
]

FLANKING_NON_VENOM = 2  # how many non-venom genes to show on each side


def parse_gff_to_loci(gff_path, species):
    """Parse BLAST-hit GFF into distinct gene loci by clustering overlapping hits."""
    hits = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            scaffold = parts[0]
            s, e = int(parts[3]), int(parts[4])
            start, end = min(s, e), max(s, e)
            strand = parts[6]
            name = ""
            evalue = 1.0
            for attr in parts[8].split(";"):
                if attr.startswith("Name="):
                    name = attr.split("=")[1]
                elif attr.startswith("Note="):
                    try:
                        evalue = float(attr.split("=")[1])
                    except ValueError:
                        evalue = 1.0
            gene_name = name.rsplit("_e", 1)[0] if "_e" in name else name
            if not gene_name:
                continue
            hits.append({
                "scaffold": scaffold, "start": start, "end": end,
                "strand": strand, "query_gene": gene_name, "evalue": evalue,
            })

    # Group by scaffold, cluster overlapping hits into loci
    by_scaffold = defaultdict(list)
    for h in hits:
        by_scaffold[h["scaffold"]].append(h)

    result = {}
    for scaffold, scaffold_hits in by_scaffold.items():
        scaffold_hits.sort(key=lambda h: h["start"])
        loci = []
        for hit in scaffold_hits:
            merged = False
            for locus in loci:
                if hit["start"] <= locus["end"] + 5000:
                    locus["end"] = max(locus["end"], hit["end"])
                    locus["hits"].append(hit)
                    merged = True
                    break
            if not merged:
                loci.append({"start": hit["start"], "end": hit["end"], "hits": [hit]})

        # Pick best gene name per locus
        clean_loci = []
        for locus in loci:
            best = min(locus["hits"], key=lambda h: h["evalue"])
            gn = best["query_gene"]
            is_venom = any(gn.startswith(p) for p in VENOM_PREFIXES)
            clean_loci.append({
                "name": gn,
                "start": locus["start"],
                "end": locus["end"],
                "strand": best["strand"],
                "is_venom": is_venom,
            })

        clean_loci.sort(key=lambda l: l["start"])
        result[scaffold] = {"species": species, "loci": clean_loci}

    return result


def build_windows(scaffolds_data):
    """Build venom-focused genomic windows.

    For each scaffold with venom genes:
    - Include ALL venom genes (even if far apart)
    - Include up to FLANKING_NON_VENOM non-venom genes on each side of the cluster
    - Non-venom genes between venom genes are included too
    """
    windows = []

    for scaffold, data in scaffolds_data.items():
        loci = data["loci"]
        species = data["species"]

        venom_indices = [i for i, l in enumerate(loci) if l["is_venom"]]
        if not venom_indices:
            continue

        # Window: from (first_venom - FLANKING) to (last_venom + FLANKING)
        first_v = min(venom_indices)
        last_v = max(venom_indices)
        win_start = max(0, first_v - FLANKING_NON_VENOM)
        win_end = min(len(loci) - 1, last_v + FLANKING_NON_VENOM)

        window_loci = loci[win_start:win_end + 1]

        windows.append({
            "species": species,
            "scaffold": scaffold,
            "loci": window_loci,
            "n_venom": len(venom_indices),
            "n_total": len(window_loci),
            "region_start": window_loci[0]["start"],
            "region_end": window_loci[-1]["end"],
        })

    return windows


def load_all_genomic_context(sp_data_dir):
    """Load and process all GFF data into baked genomic context."""
    all_scaffolds = {}

    gff_dir = sp_data_dir / "genomic_annotations"
    for gff_name, species in [("TS.gff", "Thamnophis sirtalis"), ("Crvi_mi2.gff", "Crotalus viridis")]:
        gff_path = gff_dir / gff_name
        if gff_path.exists():
            scaffolds = parse_gff_to_loci(gff_path, species)
            all_scaffolds.update(scaffolds)

    windows = build_windows(all_scaffolds)
    return windows
