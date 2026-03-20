#!/usr/bin/env python3
"""Set up symlinks to local SpeciesEmbedding data for development.

This script is for local development only — it symlinks venom-related
data from the parent SpeciesEmbedding repository into data/.
On GitHub, users should use download_data.py instead.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
SE_ROOT = REPO_ROOT.parent.parent  # SpeciesEmbedding root

# Mapping: local name → path relative to SpeciesEmbedding root
SYMLINKS = {
    "bee_venom": "data/bee_venom",
    "serine_protease": "data/serine_protease",
    "3ftx_evolution": "data/3ftx_evolution",
    "snake_venom": "data/snake_venom",
    "conotoxin": "data/conotoxin",
    "kunitz": "data/Kunitz",
    "3ftx": "data/3FTx",
    "pla2": "data/Pla2",
    "ant_venoms": "data/ant_venoms",
    "nemertea_toxprot": "data/Nemertea_ToxProt",
    "bundles": "data/bundles",
    "reference": "data/reference",
    "genomes": "data/genomes",
}


def main():
    DATA_DIR.mkdir(exist_ok=True)

    for name, rel_path in SYMLINKS.items():
        source = SE_ROOT / rel_path
        target = DATA_DIR / name

        if target.exists() or target.is_symlink():
            print(f"  skip {name} (already exists)")
            continue

        if not source.exists():
            print(f"  skip {name} (source not found: {source})")
            continue

        target.symlink_to(source)
        print(f"  link {name} -> {source}")


if __name__ == "__main__":
    main()
