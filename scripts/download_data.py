#!/usr/bin/env python3
"""Download reference datasets for VenomsBase development.

Downloads publicly available venom protein datasets from UniProt/ToxProt
and other sources. For local development data, use setup_local_symlinks.py.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"

# TODO: Add download sources as architecture stabilizes
# - ToxProt (UniProt venom subset)
# - ConoServer exports
# - ArachnoServer exports
# - VenomZone reference data


def main():
    DATA_DIR.mkdir(exist_ok=True)
    print("VenomsBase data downloader")
    print("=" * 40)
    print("Not yet implemented — architecture proposal phase.")
    print("Use scripts/setup_local_symlinks.py for local development.")


if __name__ == "__main__":
    main()
