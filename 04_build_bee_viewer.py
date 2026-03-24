#!/usr/bin/env python3
"""Step 04: Build bee venom VenomsBase viewer.

Loads pre-computed ProtT5 embeddings (4,606 proteins), matches to the
Hymenoptera_Venoms_CVP annotation CSV, assigns PBVP-family colors,
parses 739 GFF3 genomic annotations for gene-neighborhood context,
and builds a self-contained HTML viewer via 03_build_html_viewer.py.

Output: data/VenomsBase_Bees.html
"""

import csv
import importlib.util
import re
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "tools"))
import config

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = Path(__file__).resolve().parent / "data"

# Published PBVP family palette (Koludarov et al. 2023 BMC Biology)
BEE_PBVP_COLORS = {
    "Melittin": "#D32F2F",
    "Apamin": "#FBC02D",
    "Hyaluronidase": "#F9A825",
    "Icarapin": "#7CB342",
    "MCDP": "#1A237E",
    "Phospholipase A2": "#00838F",
    "Secapin": "#1976D2",
    "Venom Allergen": "#7B1FA2",
    "Venom acid phosphatase": "#E65100",
    "Dipeptidyl peptidase": "#BF360C",
    "Venom serine protease": "#C62828",
    "Tertiapin": "#F57F17",
    "Kunitz": "#E6A000",
    "Mastoparan": "#795548",
    "Poneritoxin": "#546E7A",
    "Phospholipase A1": "#00695C",
    "MRJP/Yellow": "#FFB300",
    "Serpin": "#8D6E63",
    "Venom peptide": "#78909C",
    "Chitinase": "#A1887F",
    "C-type lectin": "#AB47BC",
    "Transglutaminase": "#5D4037",
    "Ubiquitin": "#607D8B",
    "14-3-3": "#455A64",
    "Ferredoxin": "#66BB6A",
    "Glycosyltransferase": "#A1887F",
    "IGF-like": "#7986CB",
    "NUDK": "#B0BEC5",
    "Phosphodiesterase": "#4DB6AC",
    "Peroxiredoxin": "#81C784",
    "Peroxidase": "#AED581",
    "Propolyphenoloxidase": "#DCE775",
    "SLIT3": "#90A4AE",
    "Unclassified": "#BDBDBD",
}

# GFF3 filename family code -> canonical family name
FAMILY_CODE_MAP = {
    "MEL": "Melittin",
    "APA": "Apamin",
    "APH": "Venom acid phosphatase",
    "HYAL": "Hyaluronidase",
    "ICA": "Icarapin",
    "PLA2": "Phospholipase A2",
    "SECA": "Secapin",
    "SP": "Venom serine protease",
    "VA": "Venom Allergen",
    "DPP4": "Dipeptidyl peptidase",
}

# FASTA filename (aligned_sequences/) -> canonical family name
FASTA_FAMILY_MAP = {
    "APAM": "Apamin",
    "ATHR": "Venom peptide",
    "BOMB": "Mastoparan",
    "CAES": "Venom peptide",
    "CALG": "Venom peptide",
    "CAPE": "Venom peptide",
    "CARE": "Venom peptide",
    "CHIT": "Chitinase",
    "CLEC": "C-type lectin",
    "CYPR": "Venom peptide",
    "DPP4": "Dipeptidyl peptidase",
    "FERR": "Ferredoxin",
    "GT5": "Glycosyltransferase",
    "HYAL": "Hyaluronidase",
    "ICAR": "Icarapin",
    "IGFL": "IGF-like",
    "KAZA": "Kunitz",
    "MRJP": "MRJP/Yellow",
    "NEW1": "Venom peptide",
    "NEW3": "Venom peptide",
    "NUDK": "NUDK",
    "PAPH": "Venom acid phosphatase",
    "PDES": "Phosphodiesterase",
    "PERN": "Peroxiredoxin",
    "PERS": "Peroxidase",
    "PLA1": "Phospholipase A1",
    "PLA2": "Phospholipase A2",
    "PRPO": "Propolyphenoloxidase",
    "SERP": "Serpin",
    "SLIT3": "SLIT3",
    "SP": "Venom serine protease",
    "SPIN": "Venom peptide",
    "SRDP": "Venom peptide",
    "SUDI": "Venom peptide",
    "THDP": "Venom peptide",
    "TRGE": "Transglutaminase",
    "UBIQ": "Ubiquitin",
    "VA3": "Venom Allergen",
    "VAPH": "Venom acid phosphatase",
    "ACPH": "Venom acid phosphatase",
    "MEL": "Melittin",
    "SECA": "Secapin",
    "1433": "14-3-3",
    # Additional prefixes found in H5 keys
    "CAGL": "Venom peptide",
    "GLT5": "Glycosyltransferase",
    "NEW2": "Venom peptide",
    "SEPR": "Venom serine protease",
}

# Known scaffold-prefix -> species mapping for the bee study
# NC_01576* = Apis mellifera (Amel_4.5), NC_03763*/NC_03764*/NC_03765* = Apis mellifera (HAv3.1)
# NC_03950* = Bombus terrestris, NC_04575* = Bombus impatiens
# NC_05266* = Megachile rotundata, NW_003565*/NW_003789*/NW_003790*/NW_003791* = Apis mellifera (Amel_4.5)
# NW_006263* = Apis mellifera (HAv3.1), NW_012026* = Bombus terrestris
# NW_014332*/NW_014333* = Habropoda laboriosa, NW_014569* = Dufourea novaeangliae
# NW_015148*/NW_015149* = Eufriesea mexicana, NW_015373*/NW_015374* = Megachile rotundata
# NW_016017*/NW_016018*/NW_016019* = Apis dorsata
# NW_016907*/NW_016908*/NW_016909*/NW_016910*/NW_016914*/NW_016915*/NW_016916*/NW_016918*/NW_016926* = Ceratina calcarata
# NW_017100* = Polistes dominula, NW_017131*/NW_017132*/..NW_017179* = Polistes canadensis
# NW_020229*/NW_020230*/NW_020311* = Apis cerana, NW_020521* = Nomia melanderi
# NW_021683* = Cephus cinctus, NW_022279* = Orussus abietinus
# NW_022339*/NW_022340*/NW_022341*/NW_022343* = Trichogramma pretiosum
# NW_022639* = Nasonia vitripennis, NW_022882*/NW_022883*/..NW_022900* = Vespa mandarinia
# NW_023009* = Athalia rosae
# KB640601/KQ43* = Apis mellifera (alt assembly)
# NIJG0* = Euglossa dilemma, WNWW0* = Colletes gigas
# WUUM0* = Xylocopa violacea
# R*_*, contig_*, scaffold_*, ctg*_* = Halictus scabiosae
SCAFFOLD_SPECIES = {
    # NC_ chromosomes (9-char prefixes for precision)
    "NC_015762": "Apis mellifera",
    "NC_015764": "Apis mellifera",
    "NC_015765": "Apis mellifera",
    "NC_015766": "Apis mellifera",
    "NC_015767": "Apis mellifera",
    "NC_015769": "Apis mellifera",
    "NC_015771": "Apis mellifera",
    "NC_015772": "Apis mellifera",
    "NC_015773": "Apis mellifera",
    "NC_015774": "Apis mellifera",
    "NC_015775": "Apis mellifera",
    "NC_015776": "Apis mellifera",
    "NC_015777": "Apis mellifera",
    "NC_037638": "Apis mellifera",
    "NC_037639": "Apis mellifera",
    "NC_037641": "Apis mellifera",
    "NC_037642": "Apis mellifera",
    "NC_037644": "Apis mellifera",
    "NC_037645": "Apis mellifera",
    "NC_037648": "Apis mellifera",
    "NC_037649": "Apis mellifera",
    "NC_037650": "Apis mellifera",
    "NC_037651": "Apis mellifera",
    "NC_037652": "Apis mellifera",
    "NC_037653": "Apis mellifera",
    "NC_039506": "Bombus terrestris",
    "NC_039507": "Bombus terrestris",
    "NC_039508": "Bombus terrestris",
    "NC_039509": "Bombus terrestris",
    "NC_039510": "Bombus terrestris",
    "NC_039511": "Bombus terrestris",
    "NC_039512": "Bombus terrestris",
    "NC_039515": "Bombus terrestris",
    "NC_039519": "Bombus terrestris",
    "NC_045757": "Bombus impatiens",
    "NC_045758": "Bombus impatiens",
    "NC_045759": "Bombus impatiens",
    "NC_045760": "Bombus impatiens",
    "NC_045761": "Bombus impatiens",
    "NC_052665": "Megachile rotundata",
    "NC_052667": "Megachile rotundata",
    "NC_052673": "Megachile rotundata",
    "NC_052678": "Megachile rotundata",
    # NW_ scaffolds (9-char prefixes — NCBI assembly-specific ranges)
    "NW_003565": "Apis mellifera",
    "NW_003789": "Apis mellifera",
    "NW_003790": "Apis mellifera",
    "NW_003791": "Apis mellifera",
    "NW_003797": "Apis mellifera",
    "NW_003798": "Apis mellifera",
    "NW_003800": "Apis mellifera",
    "NW_006263": "Apis mellifera",
    "NW_006264": "Apis mellifera",
    "NW_012026": "Bombus terrestris",
    "NW_012027": "Bombus terrestris",
    "NW_012028": "Bombus terrestris",
    "NW_012033": "Bombus terrestris",
    "NW_012075": "Bombus terrestris",
    "NW_012077": "Bombus terrestris",
    "NW_012078": "Bombus terrestris",
    "NW_012160": "Bombus terrestris",
    "NW_014332": "Habropoda laboriosa",
    "NW_014333": "Habropoda laboriosa",
    "NW_014569": "Dufourea novaeangliae",
    "NW_015148": "Eufriesea mexicana",
    "NW_015149": "Eufriesea mexicana",
    "NW_015373": "Megachile rotundata",
    "NW_015374": "Megachile rotundata",
    "NW_016017": "Apis dorsata",
    "NW_016018": "Apis dorsata",
    "NW_016019": "Apis dorsata",
    "NW_016907": "Ceratina calcarata",
    "NW_016908": "Ceratina calcarata",
    "NW_016909": "Ceratina calcarata",
    "NW_016910": "Ceratina calcarata",
    "NW_016914": "Ceratina calcarata",
    "NW_016915": "Ceratina calcarata",
    "NW_016916": "Ceratina calcarata",
    "NW_016918": "Ceratina calcarata",
    "NW_016926": "Ceratina calcarata",
    "NW_017100": "Polistes dominula",
    "NW_017101": "Polistes dominula",
    "NW_017131": "Polistes canadensis",
    "NW_017132": "Polistes canadensis",
    "NW_017134": "Polistes canadensis",
    "NW_017135": "Polistes canadensis",
    "NW_017138": "Polistes canadensis",
    "NW_017139": "Polistes canadensis",
    "NW_017140": "Polistes canadensis",
    "NW_017143": "Polistes canadensis",
    "NW_017145": "Polistes canadensis",
    "NW_017147": "Polistes canadensis",
    "NW_017149": "Polistes canadensis",
    "NW_017150": "Polistes canadensis",
    "NW_017152": "Polistes canadensis",
    "NW_017156": "Polistes canadensis",
    "NW_017157": "Polistes canadensis",
    "NW_017162": "Polistes canadensis",
    "NW_017163": "Polistes canadensis",
    "NW_017172": "Polistes canadensis",
    "NW_017176": "Polistes canadensis",
    "NW_017177": "Polistes canadensis",
    "NW_017179": "Polistes canadensis",
    "NW_020229": "Apis cerana",
    "NW_020230": "Apis cerana",
    "NW_020311": "Apis cerana",
    "NW_020521": "Nomia melanderi",
    "NW_020522": "Nomia melanderi",
    "NW_021683": "Cephus cinctus",
    "NW_022279": "Orussus abietinus",
    "NW_022339": "Trichogramma pretiosum",
    "NW_022340": "Trichogramma pretiosum",
    "NW_022341": "Trichogramma pretiosum",
    "NW_022343": "Trichogramma pretiosum",
    "NW_022639": "Nasonia vitripennis",
    "NW_022882": "Vespa mandarinia",
    "NW_022883": "Vespa mandarinia",
    "NW_022884": "Vespa mandarinia",
    "NW_022887": "Vespa mandarinia",
    "NW_022888": "Vespa mandarinia",
    "NW_022889": "Vespa mandarinia",
    "NW_022890": "Vespa mandarinia",
    "NW_022891": "Vespa mandarinia",
    "NW_022900": "Vespa mandarinia",
    "NW_023009": "Athalia rosae",
    # Short-prefix assemblies
    "KB640601": "Apis mellifera",
    "KQ435": "Apis mellifera",
    "KQ436": "Apis mellifera",
    "NIJG0": "Euglossa dilemma",
    "WNWW0": "Colletes gigas",
    "WUUM0": "Xylocopa violacea",
}

# Catch-all prefixes for Halictus scabiosae (contig_*, scaffold_*, R*_*, ctg*_*)
HALICTUS_PATTERNS = re.compile(r"^(contig_|scaffold_|R\d+_|ctg\d+_)")


# ---------------------------------------------------------------------------
# 1. Load embeddings
# ---------------------------------------------------------------------------

def load_embeddings(h5_path):
    """Load pre-computed ProtT5 per-protein embeddings from H5."""
    embeddings = {}
    with h5py.File(h5_path, "r") as f:
        for key in f.keys():
            arr = np.array(f[key])
            # Already per-protein (1024,) — no pooling needed
            if arr.ndim == 2:
                arr = arr.mean(axis=0)
            embeddings[key] = arr.astype(np.float32)
    return embeddings


# ---------------------------------------------------------------------------
# 2. Load annotations and match to embeddings
# ---------------------------------------------------------------------------

def fuzzy_family(raw_family):
    """Map a raw protein family string to a canonical BEE_PBVP_COLORS key."""
    if not raw_family:
        return "Unclassified"
    low = raw_family.lower().strip()
    if not low or low == "nan":
        return "Unclassified"
    # Direct substring matches (order matters — check specific before general)
    mapping = [
        ("melittin", "Melittin"),
        ("apamin", "Apamin"),
        ("hyaluronidase", "Hyaluronidase"),
        ("icarapin", "Icarapin"),
        ("mcdp", "MCDP"),
        ("mast cell degranulating", "MCDP"),
        ("phospholipase a2", "Phospholipase A2"),
        ("phospholipase a1", "Phospholipase A1"),
        ("secapin", "Secapin"),
        ("venom allergen", "Venom Allergen"),
        ("allergen", "Venom Allergen"),
        ("venom acid phosphatase", "Venom acid phosphatase"),
        ("acid phosphatase", "Venom acid phosphatase"),
        ("dipeptidyl peptidase", "Dipeptidyl peptidase"),
        ("dpp4", "Dipeptidyl peptidase"),
        ("venom serine protease", "Venom serine protease"),
        ("serine protease", "Venom serine protease"),
        ("serine carboxypeptidase", "Venom serine protease"),
        ("tertiapin", "Tertiapin"),
        ("kunitz", "Kunitz"),
        ("mastoparan", "Mastoparan"),
        ("mastoporan", "Mastoparan"),
        ("bombolittin", "Mastoparan"),
        ("poneritoxin", "Poneritoxin"),
        ("ectatoxin", "Poneritoxin"),
        ("myrmecitoxin", "Poneritoxin"),
        ("myrmeciitoxin", "Poneritoxin"),
        ("myrmicitoxin", "Poneritoxin"),
        ("pseudomyrmecitoxin", "Poneritoxin"),
        ("mrjp", "MRJP/Yellow"),
        ("yellow", "MRJP/Yellow"),
        ("serpin", "Serpin"),
        ("chitinase", "Chitinase"),
        ("c-type lectin", "C-type lectin"),
        ("lectin", "C-type lectin"),
        ("transglutaminase", "Transglutaminase"),
        ("ubiquitin", "Ubiquitin"),
        ("trehalase", "Venom enzyme"),
        ("aspartylglucosaminidase", "Venom enzyme"),
        ("metalloproteinase", "Venom enzyme"),
        ("carboxylesterase", "Venom enzyme"),
        ("methyltransferase", "Venom enzyme"),
        # Named venom peptides → Venom peptide
        ("ampulexin", "Venom peptide"),
        ("anoplin", "Venom peptide"),
        ("bradykinin", "Venom peptide"),
        ("chemotactic", "Venom peptide"),
        ("crabrolin", "Venom peptide"),
        ("decoralin", "Venom peptide"),
        ("dominulin", "Venom peptide"),
        ("fulvonin", "Venom peptide"),
        ("lasioglossin", "Venom peptide"),
        ("macropin", "Venom peptide"),
        ("orientotoxin", "Venom peptide"),
        ("osmin", "Venom peptide"),
        ("panurgin", "Venom peptide"),
        ("paulistine", "Venom peptide"),
        ("polybine", "Venom peptide"),
        ("pompilidotoxin", "Venom peptide"),
        ("poneratoxin", "Poneritoxin"),
        ("protonectin", "Venom peptide"),
        ("protopolybiakinin", "Venom peptide"),
        ("sylverin", "Venom peptide"),
        ("tachikinin", "Venom peptide"),
        ("vespakin", "Venom peptide"),
        ("vespulakinin", "Venom peptide"),
        ("vespin", "Venom peptide"),
        ("waspkinin", "Venom peptide"),
        ("xylopin", "Venom peptide"),
        ("codesan", "Venom peptide"),
        ("halictin", "Venom peptide"),
        ("histamin releasing", "Venom peptide"),
        ("venom peptide", "Venom peptide"),
        ("venom protein", "Venom peptide"),
        ("venom protease", "Venom serine protease"),
        ("small venom", "Venom peptide"),
        ("cystein-rich", "Venom peptide"),
    ]
    for pat, canon in mapping:
        if pat in low:
            return canon
    # Return the raw string — preserves original info, never "Other"
    return raw_family.strip()


def family_from_h5_key(key):
    """Infer venom family from an H5 key that came from aligned FASTAs.

    Keys like 'APAM_g1', 'PLA2_g4', 'SP_g3' have the FASTA filename prefix
    before the first underscore followed by '_g'.
    """
    # Pattern: PREFIX_gN or PREFIX_gN_pN
    m = re.match(r"^([A-Z0-9]+)_g\d+", key)
    if m:
        prefix = m.group(1)
        return FASTA_FAMILY_MAP.get(prefix)  # None if not found, not "Other"
    return None


def build_fasta_id_to_family():
    """Build {sequence_id: family} from aligned_sequences/*.fasta filenames."""
    mapping = {}
    aligned_dir = config.BEE_ALIGNED
    if not aligned_dir.exists():
        return mapping
    for fasta_file in aligned_dir.glob("*.fasta"):
        stem = fasta_file.stem.replace("_cleaned", "")
        family = FASTA_FAMILY_MAP.get(stem, stem)
        with open(fasta_file) as f:
            for line in f:
                if line.startswith(">"):
                    sid = line[1:].split()[0]
                    if sid not in mapping:
                        mapping[sid] = family
    print(f"  Built FASTA→family mapping: {len(mapping)} IDs across {len(set(mapping.values()))} families")
    return mapping


def load_species_map():
    """Load protein ID → species mapping from pre-built TSV."""
    species_map = {}
    tsv = config.BEE_VENOM_DIR / "protein_species_map.tsv"
    if tsv.exists():
        with open(tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row.get("species"):
                    species_map[row["protein_id"]] = row["species"]
        print(f"  Species map: {len(species_map)} proteins mapped")
    return species_map


# Complete species → Hymenoptera group mapping.
# Covers all 23 genome species from Koludarov et al. 2023 + outgroups + ToxProt refs.
SPECIES_TO_GROUP = {
    # 2_Honeybees_bumblebees (Apidae sensu lato + close relatives)
    "Apis mellifera": "Honeybees & bumblebees",
    "Apis mellifera ligustica": "Honeybees & bumblebees",
    "Apis cerana": "Honeybees & bumblebees",
    "Apis dorsata": "Honeybees & bumblebees",
    "Apis florea": "Honeybees & bumblebees",
    "Bombus terrestris": "Honeybees & bumblebees",
    "Bombus impatiens": "Honeybees & bumblebees",
    "Bombus vosnesenski": "Honeybees & bumblebees",
    "Bombus pensylvanicus": "Honeybees & bumblebees",
    "Bombus ignitus": "Honeybees & bumblebees",
    "Bombus lapidarius": "Honeybees & bumblebees",
    "Eufriesea mexicana": "Honeybees & bumblebees",
    "Euglossa dilemma": "Honeybees & bumblebees",
    "Ceratina calcarata": "Honeybees & bumblebees",
    "Habropoda laboriosa": "Honeybees & bumblebees",
    "Xylocopa violacea": "Honeybees & bumblebees",
    # 1_Solitary_bees (non-Apidae bees)
    "Colletes gigas": "Solitary bees",
    "Dufourea novaeangliae": "Solitary bees",
    "Megachile rotundata": "Solitary bees",
    "Megalopta genalis": "Solitary bees",
    "Nomia melanderi": "Solitary bees",
    "Lasioglossum albipes": "Solitary bees",
    "Halictus scabiosae": "Solitary bees",
    "Osmia rufa": "Solitary bees",
    "Panurgus calcaratus": "Solitary bees",
    "Macropis fulvipes": "Solitary bees",
    "Melecta albifrons": "Solitary bees",
    # 3_Ants (Formicidae)
    "Solenopsis invicta": "Ants",
    "Solenopsis geminata": "Ants",
    "Solenopsis richteri": "Ants",
    "Solenopsis saevissima": "Ants",
    "Harpegnathos saltator": "Ants",
    "Vollenhovia emeryi": "Ants",
    "Paraponera clavata": "Ants",
    "Dinoponera australis": "Ants",
    "Dinoponera quadriceps": "Ants",
    "Neoponera apicalis": "Ants",
    "Neoponera goeldii": "Ants",
    "Neoponera inversa": "Ants",
    "Brachyponera chinensis": "Ants",
    "Ectatomma brunneum": "Ants",
    "Ectatomma tuberculatum": "Ants",
    "Odontomachus monticola": "Ants",
    "Myrmecia banksi": "Ants",
    "Myrmecia gulosa": "Ants",
    "Myrmecia pilosula": "Ants",
    "Myrmica rubra": "Ants",
    "Pseudomyrmex penetrator": "Ants",
    "Pseudomyrmex triplarinus": "Ants",
    "Tetramorium bicarinatum": "Ants",
    "Manica rubida": "Ants",
    "Anochetus emarginatus": "Ants",
    "Dolichovespula arenaria": "Wasps",  # actually a wasp, mislabeled in CSV
    # 4_TrueWasps (Vespidae)
    "Vespa mandarinia": "Wasps",
    "Vespa crabro": "Wasps",
    "Vespa affinis": "Wasps",
    "Vespa analis": "Wasps",
    "Vespa basalis": "Wasps",
    "Vespa bicolor": "Wasps",
    "Vespa magnifica": "Wasps",
    "Vespa orientalis": "Wasps",
    "Vespa tropica": "Wasps",
    "Vespa velutina": "Wasps",
    "Vespa velutina nigrithorax": "Wasps",
    "Vespa xanthoptera": "Wasps",
    "Polistes dominula": "Wasps",
    "Polistes canadensis": "Wasps",
    "Polistes annularis": "Wasps",
    "Polistes exclamans": "Wasps",
    "Polistes fuscatus": "Wasps",
    "Polistes gallicus": "Wasps",
    "Polistes hebraeus": "Wasps",
    "Polistes jokahamae": "Wasps",
    "Polistes lanio": "Wasps",
    "Polistes major": "Wasps",
    "Vespula germanica": "Wasps",
    "Vespula lewisii": "Wasps",
    "Vespula maculifrons": "Wasps",
    "Vespula pensylvanica": "Wasps",
    "Vespula squamosa": "Wasps",
    "Vespula vidua": "Wasps",
    "Vespula vulgaris": "Wasps",
    "Vespula flavopilosa": "Wasps",
    "Dolichovespula maculata": "Wasps",
    "Agelaia pallipes pallipes": "Wasps",
    "Polybia occidentalis": "Wasps",
    "Polybia paulista": "Wasps",
    "Polybia scutellaris rioplatensis": "Wasps",
    "Protonectarina sylveirae": "Wasps",
    "Protopolybia exigua": "Wasps",
    "Parapolybia indica": "Wasps",
    "Icaria sp.": "Wasps",
    "Eumenes fraterculus": "Wasps",
    "Eumenes pomiformis": "Wasps",
    "Eumenes rubrofemoratus": "Wasps",
    "Eumenes rubronotatus": "Wasps",
    "Orancistrocerus drewseni": "Wasps",
    "Oreumenes decoratus": "Wasps",
    "Rhynchium brunneum": "Wasps",
    "Megascolia flavifrons": "Wasps",
    "Sphecius speciosus": "Wasps",
    # 5_Parasitoids
    "Nasonia vitripennis": "Parasitoids",
    "Trichogramma pretiosum": "Parasitoids",
    "Ampulex compressa": "Parasitoids",
    "Asobara tabida": "Parasitoids",
    "Chelonus sp. nr. curvimaculatus": "Parasitoids",
    "Cotesia rubecula": "Parasitoids",
    "Eulophus pennicornis": "Parasitoids",
    "Microctonus hyperodae": "Parasitoids",
    "Pimpla hypochondriaca": "Parasitoids",
    "Anoplius samariensis": "Parasitoids",
    "Batozonellus maculifrons": "Parasitoids",
    "Cyphononyx dorsalis": "Parasitoids",
    "Isodontia harmandi": "Parasitoids",
    "Sphex argentatus argentatus": "Parasitoids",
    "Anterhynchium flavomarginatum micado": "Parasitoids",
    # Symphyta outgroups
    "Athalia rosae": "Sawflies",
    "Cephus cinctus": "Sawflies",
    "Orussus abietinus": "Sawflies",
}

# Known PDB structures from bee venom (all Apis mellifera)
PDB_SPECIES = {
    "1FCQ_A": "Apis mellifera",   # hyaluronidase
    "1POC_A": "Apis mellifera",   # phospholipase A2
    "1QNX_A": "Apis mellifera",   # venom allergen 5
    "2ATM_A": "Apis mellifera",   # hyaluronidase
    "2VZN_A": "Apis mellifera",   # venom allergen 3
}


def load_reference_fasta():
    """Parse uniprot-reviewed_Apoidea reference FASTA for species + descriptions."""
    ref_fasta = config.BEE_VENOM_DIR / "uniprot-reviewed_Apoidea_291020.fasta"
    ref_data = {}  # uid -> {species, description}
    if not ref_fasta.exists():
        return ref_data
    with open(ref_fasta) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip()
                # >sp|P01501|MEL_APIME Melittin OS=Apis mellifera OX=7460 ...
                if "|" in parts:
                    uid = parts.split("|")[1]
                else:
                    uid = parts.split()[0]
                # Extract OS=species
                m = re.search(r"OS=(.+?)(?:\s+OX=|\s+GN=|\s+PE=|$)", parts)
                species = m.group(1).strip() if m else ""
                # Description: text between entry name and OS=
                desc_m = re.search(r"\|[A-Z0-9_]+\s+(.+?)\s+OS=", parts)
                desc = desc_m.group(1).strip() if desc_m else ""
                ref_data[uid] = {"species": species, "description": desc}
    print(f"  Reference FASTA: {len(ref_data)} entries for species/description fallback")
    return ref_data


def resolve_species(key, species_map, ref_data=None):
    """Get species for a protein ID from all available sources."""
    # 1. Pre-built species map (covers 4,391 / 4,606)
    if key in species_map:
        return species_map[key]
    # 2. PDB structures
    if key in PDB_SPECIES:
        return PDB_SPECIES[key]
    # 3. Reference FASTA
    if ref_data and key in ref_data:
        return ref_data[key]["species"]
    return ""


def resolve_group(species):
    """Map species name to Hymenoptera group. Uses exact match, then genus match."""
    if not species:
        return ""
    species = species.strip()
    # Exact match
    if species in SPECIES_TO_GROUP:
        return SPECIES_TO_GROUP[species]
    # Strip parenthetical annotations: "Apis mellifera (Honeybee)" → "Apis mellifera"
    clean = re.sub(r'\s*\(.*?\)', '', species).strip()
    if clean in SPECIES_TO_GROUP:
        return SPECIES_TO_GROUP[clean]
    # Genus match
    genus = species.split()[0]
    for sp, group in SPECIES_TO_GROUP.items():
        if sp.startswith(genus + " "):
            return group
    return ""


def load_annotations_and_match(h5_keys):
    """Load CSV annotations and build a unified annotation DataFrame keyed by H5 key."""
    # Build FASTA ID → family mapping from aligned sequences
    fasta_family = build_fasta_id_to_family()

    # Load species mapping + reference FASTA fallback
    species_map = load_species_map()
    ref_data = load_reference_fasta()

    # Read the master CSV
    ann = pd.read_csv(config.BEE_ANNOTATIONS, encoding="utf-8-sig")
    ann.columns = ann.columns.str.strip()

    # Build lookup dicts from CSV: by Entry, by Entry name
    by_entry = {}
    by_entry_name = {}
    for _, row in ann.iterrows():
        entry = str(row.get("Entry", "")).strip()
        entry_name = str(row.get("Entry name", "")).strip()
        if entry:
            by_entry[entry] = row
        if entry_name:
            by_entry_name[entry_name] = row

    records = []
    matched = 0

    for key in h5_keys:
        rec = {"identifier": key}

        # Try matching to CSV
        csv_row = None

        # Direct Entry match (e.g., key = "P01500")
        if key in by_entry:
            csv_row = by_entry[key]
        else:
            # Try stripping version from key (e.g., "A0A1W6EVM7" from "A0A1W6EVM7")
            bare = key.split(".")[0] if "." in key else key
            if bare in by_entry:
                csv_row = by_entry[bare]
            # Try Entry name match (e.g., key = "APAM_APIME")
            elif key in by_entry_name:
                csv_row = by_entry_name[key]
            else:
                # sp|P01501.1|MEL_APIME -> try P01501 and MEL_APIME
                if key.startswith("sp|"):
                    parts = key.split("|")
                    if len(parts) >= 2:
                        uid = parts[1].split(".")[0]
                        if uid in by_entry:
                            csv_row = by_entry[uid]
                    if csv_row is None and len(parts) >= 3:
                        ename = parts[2]
                        if ename in by_entry_name:
                            csv_row = by_entry_name[ename]

        if csv_row is not None:
            matched += 1
            rec["name"] = str(csv_row.get("Entry name", key))
            rec["description"] = str(csv_row.get("Protein names", ""))[:253]
            raw_fam = str(csv_row.get("Protein family/Gene name", ""))
            rec["venom_family"] = fuzzy_family(raw_fam)
            rec["species"] = str(csv_row.get("Species", ""))
            rec["source"] = "UniProt"
        else:
            # Not in CSV — infer family from H5 key pattern or FASTA mapping
            fam = family_from_h5_key(key)
            if not fam and key in fasta_family:
                fam = fasta_family[key]
            rec["name"] = key
            rec["description"] = ""
            # Use reference FASTA description to classify if no family found
            if not fam and key in ref_data:
                desc = ref_data[key].get("description", "").lower()
                fam = fuzzy_family(desc)
                if fam == desc:  # fuzzy_family returned raw string, meaning no match
                    fam = "Venom peptide"  # safe default for ToxProt Apoidea entries
            rec["venom_family"] = fam if fam else "Venom peptide"
            rec["source"] = "aligned_FASTA" if fam else "reference"
            rec["species"] = resolve_species(key, species_map, ref_data)

        # Derive taxon_group from species for ALL proteins (not just CSV-matched)
        rec["taxon_group"] = resolve_group(rec.get("species", ""))

        records.append(rec)

    df = pd.DataFrame(records)

    # Report coverage
    n_species = df["species"].replace("", pd.NA).dropna().nunique()
    n_no_species = (df["species"].fillna("") == "").sum()
    n_no_group = (df["taxon_group"].fillna("") == "").sum()
    n_no_family = df["venom_family"].isin(["Unclassified", ""]).sum()
    print(f"  {len(records)} proteins, {matched} matched to annotations CSV")
    print(f"  Species: {n_species} unique, {n_no_species} missing")
    print(f"  Taxon group: {n_no_group} missing")
    print(f"  Venom family: {n_no_family} unclassified")
    return df


# ---------------------------------------------------------------------------
# 3. Compute projections (UMAP + PCA)
# ---------------------------------------------------------------------------

def compute_projections(embeddings, identifiers):
    """Compute UMAP and PCA 2D projections from embedding matrix."""
    from sklearn.decomposition import PCA

    # Build matrix in identifier order
    emb_matrix = np.array([embeddings[pid] for pid in identifiers])

    # PCA
    pca = PCA(n_components=2, random_state=42)
    pca_2d = pca.fit_transform(emb_matrix)
    pca_df = pd.DataFrame(
        {"x": pca_2d[:, 0], "y": pca_2d[:, 1]},
        index=identifiers,
    )
    pca_df.index.name = "identifier"
    print(f"  PCA variance explained: {pca.explained_variance_ratio_.sum():.1%}")

    # UMAP
    try:
        import umap
        n_neighbors = min(15, len(identifiers) - 1)
        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=n_neighbors,
            min_dist=0.1,
            metric="cosine",
            random_state=42,
        )
        umap_2d = reducer.fit_transform(emb_matrix)
        umap_df = pd.DataFrame(
            {"x": umap_2d[:, 0], "y": umap_2d[:, 1]},
            index=identifiers,
        )
        umap_df.index.name = "identifier"
        print(f"  UMAP computed ({len(identifiers)} points)")
    except ImportError:
        print("  WARNING: umap-learn not installed, falling back to PCA for UMAP view")
        umap_df = pca_df.copy()

    return umap_df, pca_df


# ---------------------------------------------------------------------------
# 4. Parse bee GFF3 genomic context
# ---------------------------------------------------------------------------

def scaffold_to_species(scaffold_id):
    """Map a scaffold accession to its species using prefix lookup."""
    # Check explicit prefixes (longest match first)
    for prefix_len in [9, 8, 7, 6, 5]:
        if prefix_len > len(scaffold_id):
            continue
        prefix = scaffold_id[:prefix_len]
        if prefix in SCAFFOLD_SPECIES:
            return SCAFFOLD_SPECIES[prefix]
    # Halictus catch-all
    if HALICTUS_PATTERNS.match(scaffold_id):
        return "Halictus scabiosae"
    return "Unknown"


def extract_family_code(filename):
    """Extract the venom family code from a GFF3 filename.

    e.g., 'NC_015765.1._MEL.gff_modified_...' -> 'MEL'
    """
    # Pattern: <scaffold>._<CODE>.gff_modified...
    m = re.search(r"\._([A-Z0-9]+)\.gff", filename)
    if m:
        return m.group(1)
    return None


def parse_bee_gff3(gff_path):
    """Parse a single bee GFF3 file into venom gene + flanking gene records.

    Returns list of gene dicts with: name, start, end, strand, is_venom, source
    """
    filename = gff_path.name
    family_code = extract_family_code(filename)
    family_name = FAMILY_CODE_MAP.get(family_code, "Other") if family_code else "Other"

    # Extract scaffold ID from the ##sequence-region pragma or first data line
    scaffold_id = None
    venom_genes = []  # Koludarov_et_al_2021 gene features
    gnomon_genes = {}  # Gnomon gene/mRNA features (deduplicate by position)

    with open(gff_path, errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 2:
                    scaffold_id = parts[1]
                continue
            if line.startswith("#") or not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            if scaffold_id is None:
                scaffold_id = cols[0]

            source = cols[1]
            feature_type = cols[2]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6] if cols[6] in ("+", "-") else "+"

            # Parse Name from attributes
            name = ""
            attrs = cols[8]
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("Name="):
                    name = attr[5:]
                    break

            # Venom genes: Koludarov_et_al_2021 source, gene or CDS type
            if source == "Koludarov_et_al_2021" and feature_type == "gene":
                if name and name != "GAP":
                    venom_genes.append({
                        "name": name,
                        "start": min(start, end),
                        "end": max(start, end),
                        "strand": strand,
                        "is_venom": True,
                        "family": family_name,
                    })

            # Flanking context: Gnomon gene or mRNA features
            elif source == "Gnomon" and feature_type in ("gene", "mRNA"):
                if name:
                    # Use midpoint as dedup key to avoid duplicate gene/mRNA
                    mid = (start + end) // 2
                    key = f"{name}_{mid // 10000}"
                    if key not in gnomon_genes:
                        # Clean up Gnomon name: use LOC* for genes, XM_* for mRNA
                        display_name = name
                        if display_name.startswith("LOC"):
                            display_name = display_name  # keep as-is
                        elif display_name.startswith("XM_") or display_name.startswith("XR_"):
                            # Shorten: XM_003394788.3 -> XM_..4788
                            display_name = display_name[:3] + ".." + display_name[-4:]
                        gnomon_genes[key] = {
                            "name": display_name,
                            "start": min(start, end),
                            "end": max(start, end),
                            "strand": strand,
                            "is_venom": False,
                            "family": "",
                        }

    return scaffold_id, venom_genes, list(gnomon_genes.values())


def build_genomic_windows(gff_dir, max_flanking=3):
    """Parse all bee GFF3 files and build genomic context windows.

    For each file:
    - Collect Koludarov_et_al_2021 gene features (venom genes)
    - Collect nearby Gnomon gene/mRNA features (flanking context)
    - Cluster venom genes on the same scaffold
    - Include up to max_flanking Gnomon genes on each side of the cluster
    """
    gff_files = sorted(gff_dir.glob("*.gff3"))
    # Skip bees_all.gff3 (it's a concatenation)
    gff_files = [f for f in gff_files if f.name != "bees_all.gff3"]

    # Accumulate per-scaffold data across files (same scaffold may appear
    # in multiple files, e.g., NC_015765.1._MEL and NC_015765.1._SP)
    scaffold_venom = defaultdict(list)
    scaffold_flanking = defaultdict(list)
    scaffold_species = {}
    scaffold_families = defaultdict(set)

    for gff_path in gff_files:
        scaffold_id, venom, flanking = parse_bee_gff3(gff_path)
        if not scaffold_id or not venom:
            continue

        species = scaffold_to_species(scaffold_id)
        scaffold_species[scaffold_id] = species
        scaffold_venom[scaffold_id].extend(venom)
        scaffold_flanking[scaffold_id].extend(flanking)
        for v in venom:
            scaffold_families[scaffold_id].add(v.get("family", ""))

    # Build windows per scaffold
    windows = []
    for scaffold_id in sorted(scaffold_venom.keys()):
        venom = scaffold_venom[scaffold_id]
        flanking = scaffold_flanking[scaffold_id]
        species = scaffold_species.get(scaffold_id, "Unknown")

        # Deduplicate venom genes by position
        seen_pos = set()
        unique_venom = []
        for v in sorted(venom, key=lambda g: g["start"]):
            key = (v["start"], v["end"])
            if key not in seen_pos:
                seen_pos.add(key)
                unique_venom.append(v)
        venom = unique_venom

        if not venom:
            continue

        # Determine the venom region boundaries
        venom_start = min(v["start"] for v in venom)
        venom_end = max(v["end"] for v in venom)

        # Find flanking Gnomon genes near the venom cluster
        # Deduplicate flanking by name+position
        seen_flanking = set()
        unique_flanking = []
        for g in flanking:
            key = (g["name"], g["start"] // 5000)
            if key not in seen_flanking:
                seen_flanking.add(key)
                unique_flanking.append(g)

        # Sort all flanking by start position
        unique_flanking.sort(key=lambda g: g["start"])

        # Select flanking genes: those within 200kb of the venom cluster
        # or the closest N on each side
        margin = 200_000
        nearby = [g for g in unique_flanking
                   if g["end"] >= venom_start - margin and g["start"] <= venom_end + margin]

        # Separate into upstream and downstream of the venom cluster
        upstream = [g for g in nearby if g["end"] < venom_start]
        downstream = [g for g in nearby if g["start"] > venom_end]
        # Between venom genes (included regardless)
        between = [g for g in nearby
                   if g["start"] >= venom_start and g["end"] <= venom_end
                   and not g["is_venom"]]

        # Take max_flanking from each side (closest to the cluster)
        upstream_sel = upstream[-max_flanking:] if len(upstream) > max_flanking else upstream
        downstream_sel = downstream[:max_flanking] if len(downstream) > max_flanking else downstream

        # Combine and sort all loci
        all_loci = upstream_sel + venom + between + downstream_sel
        all_loci.sort(key=lambda g: g["start"])

        # Deduplicate (venom may overlap with between)
        final_loci = []
        seen = set()
        for g in all_loci:
            key = (g["start"], g["end"])
            if key not in seen:
                seen.add(key)
                final_loci.append(g)

        families = scaffold_families.get(scaffold_id, set())
        family_str = ", ".join(sorted(f for f in families if f))

        windows.append({
            "species": species,
            "scaffold": scaffold_id,
            "n_venom": len(venom),
            "families": family_str,
            "loci": [
                {
                    "n": g["name"],
                    "s": g["start"],
                    "e": g["end"],
                    "d": g["strand"],
                    "v": g["is_venom"],
                }
                for g in final_loci
            ],
        })

    # Sort windows: most venom genes first
    windows.sort(key=lambda w: (-w["n_venom"], w["species"], w["scaffold"]))
    return windows


# ---------------------------------------------------------------------------
# 5. Build the HTML viewer via 03
# ---------------------------------------------------------------------------

def build_viewer(ann_df, umap_df, pca_df, genomic_windows, output_path):
    """Import build_html from 03_build_html_viewer.py and generate the HTML."""
    viewer_script = Path(__file__).resolve().parent / "03_build_html_viewer.py"
    spec = importlib.util.spec_from_file_location("viewer", str(viewer_script))
    viewer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(viewer)

    viewer.build_html(ann_df, umap_df, pca_df, genomic_windows, output_path,
                      color_map=BEE_PBVP_COLORS,
                      title="VenomsBase \u2014 Bee Venoms")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    DATA_DIR.mkdir(exist_ok=True)
    print("=" * 60)
    print("Building VenomsBase — Bee Venoms viewer")
    print("=" * 60)

    # 1. Load embeddings
    h5_path = REPO_ROOT / "data" / "bee_venom" / "bee_prott5.h5"
    print(f"\n1. Loading embeddings from {h5_path.name}...")
    embeddings = load_embeddings(h5_path)
    print(f"   {len(embeddings)} proteins with 1024-dim ProtT5 embeddings")

    # 2. Load annotations and match
    print("\n2. Matching to annotations CSV...")
    h5_keys = list(embeddings.keys())
    ann_df = load_annotations_and_match(h5_keys)

    # Report family distribution
    family_counts = ann_df["venom_family"].value_counts()
    print("\n   Family distribution:")
    for fam, count in family_counts.items():
        color = BEE_PBVP_COLORS.get(fam, "#9E9E9E")
        print(f"     {fam:30s} {count:5d}  {color}")

    # 3. Compute projections
    print("\n3. Computing UMAP + PCA projections...")
    identifiers = ann_df["identifier"].tolist()
    umap_df, pca_df = compute_projections(embeddings, identifiers)

    # 4. Parse genomic context
    gff_dir = REPO_ROOT / "data" / "bee_venom" / "genomic_annotations"
    print(f"\n4. Parsing genomic context from {gff_dir.name}/...")
    genomic_windows = build_genomic_windows(gff_dir)
    total_genes = sum(len(w["loci"]) for w in genomic_windows)
    total_venom = sum(w["n_venom"] for w in genomic_windows)
    n_species = len(set(w["species"] for w in genomic_windows))
    print(f"   {len(genomic_windows)} windows across {n_species} species")
    print(f"   {total_genes} total loci ({total_venom} venom genes)")

    # Show species breakdown
    species_counts = defaultdict(int)
    for w in genomic_windows:
        species_counts[w["species"]] += 1
    print("\n   Species breakdown:")
    for sp, cnt in sorted(species_counts.items(), key=lambda x: -x[1]):
        print(f"     {sp:30s} {cnt:4d} windows")

    # 5. Build HTML
    print("\n5. Building HTML viewer...")
    output = DATA_DIR / "VenomsBase_Bees.html"
    build_viewer(ann_df, umap_df, pca_df, genomic_windows, output)

    print(f"\nDone! Open in browser:")
    print(f"  open {output}")


if __name__ == "__main__":
    main()
