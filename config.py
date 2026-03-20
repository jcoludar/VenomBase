"""VenomsBase project configuration.

Paths, parameters, and constants for the VenomsBase data architecture prototype.
"""

from pathlib import Path

# --- Paths ---
REPO_ROOT = Path(__file__).resolve().parent
SE_ROOT = REPO_ROOT.parent.parent  # Parent repository root

# Project directories
DATA_DIR = REPO_ROOT / "data"
DOCS_DIR = REPO_ROOT / "docs"
SRC_DIR = REPO_ROOT / "src"
SCRIPTS_DIR = REPO_ROOT / "scripts"
TESTS_DIR = REPO_ROOT / "tests"

# Bee venom evolution study (BMC Biology 2023, Koludarov et al.)
# Data lives in data/bee_venom/ at SE root, symlinked into VenomBase/data/
BEE_VENOM_DIR = DATA_DIR / "bee_venom"
BEE_ANNOTATIONS = BEE_VENOM_DIR / "Hymenoptera_Venoms_CVP_Feb_2022_372.csv"
BEE_SUPP_TABLES = BEE_VENOM_DIR / "supplementary_tables"
BEE_PROTEOMICS = BEE_VENOM_DIR / "proteomics"
BEE_TRANSCRIPTOMICS = BEE_VENOM_DIR / "transcriptomics"
BEE_ALIGNED = BEE_VENOM_DIR / "aligned_sequences"
BEE_TREES = BEE_VENOM_DIR / "trees"
BEE_EMBEDDINGS = BEE_VENOM_DIR / "embeddings"
BEE_BLAST = BEE_VENOM_DIR / "blast_results"

# Serine protease evolution (BMC Biology 2021, Barua*, Koludarov* & Mikheyev)
SP_DIR = DATA_DIR / "serine_protease"

# 3FTx evolution (Nature Communications 2023, Koludarov, Senoner et al.)
THRFTX_DIR = DATA_DIR / "3ftx_evolution"
THRFTX_MASTER_CSV = THRFTX_DIR / "3and6.csv"
THRFTX_DATASET = THRFTX_DIR / "supplementary" / "SD1_Dataset_and_information.csv"

# Symlinked data from parent repository (local dev only)
SNAKE_VENOM_DIR = DATA_DIR / "snake_venom"
CONOTOXIN_DIR = DATA_DIR / "conotoxin"
KUNITZ_DIR = DATA_DIR / "kunitz"
THREE_FTX_DIR = DATA_DIR / "3ftx"
PLA2_DIR = DATA_DIR / "pla2"
ANT_VENOMS_DIR = DATA_DIR / "ant_venoms"
NEMERTEA_DIR = DATA_DIR / "nemertea_toxprot"
BUNDLES_DIR = DATA_DIR / "bundles"

# --- Tools (from parent repository) ---
TOOLS_DIR = SE_ROOT / "tools"
PROTSPACE_DIR = TOOLS_DIR / "protspace" / "protspace"

# --- Parameters ---
# Embedding models
EMBEDDING_MODELS = ["prott5", "esm2"]

# Species in bee venom study
BEE_SPECIES = {
    "apis": "Apis mellifera",
    "halictus": "Halictus scabiosae",
    "xylocopa": "Xylocopa violacea",
}

# Venom protein families (from aligned sequences)
VENOM_FAMILIES = [
    "APAM", "ATHR", "BOMB", "CAES", "CALG", "CAPE", "CARE", "CHIT",
    "CLEC", "CYPR", "DPP4", "FERR", "GT5", "HYAL", "ICAR", "IGFL",
    "KAZA", "MRJP", "NEW1", "NEW3", "NUDK", "PAPH", "PDES", "PERN",
    "PERS", "PLA1", "PLA2", "PRPO", "SERP", "SLIT3", "SP", "SPIN",
    "SRDP", "SUDI", "THDP", "TRGE", "UBIQ", "VA3", "VAPH",
]

# --- Publication Info ---
PUBLICATION = {
    "title": "Prevalent bee venom genes evolved before the aculeate stinger and eusociality",
    "authors": "Koludarov et al.",
    "journal": "BMC Biology",
    "year": 2023,
    "volume": "21:229",
    "doi": "10.1186/s12915-023-01656-5",
    "status": "published",
}
