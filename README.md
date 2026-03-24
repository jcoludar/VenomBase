# VenomsBase / VenomZone

Protein-centric data architecture and web platform for venom research.

## Data Model

Three object types:

- **Protein Record** — The primary data container. Each protein carries its embedding (ProtT5), sequences, structure pointers, taxonomy, classification, expression, proteomics, activity, and literature as typed data slots.
- **Genome Record** — Standalone. Gene annotations and scaffolds from GFF3. Linked to proteins via gene IDs, but not hierarchically above them.
- **Raw Data** — Stored as files (reads, spectra, matrices). Features extracted into protein/genome records.

See `docs/architecture_proposal.md` for the full data model.

## Prototypes

Interactive HTML viewers demonstrating protein-centric exploration:

- **VenomsBase_Snakes.html** — 2,235 snake venom proteins (3FTx + KLK/SP + PLA2) with real ProtT5 embeddings, publication colors, genomic context
- **VenomsBase_Bees.html** — 4,606 bee venom proteins (33 families) with real ProtT5 embeddings, 739 GFF3 genomic regions across 23 species

Open in Chrome. Click proteins to inspect. Color by any annotation dimension.

## Data

Data files are not tracked in this repository. Use the setup scripts:

```bash
python scripts/download_data.py          # Download reference datasets
python scripts/setup_local_symlinks.py   # Symlink to parent repository data (dev only)
```

## Pipeline

```bash
python 04_build_bee_viewer.py     # Bee venom viewer
python 05_build_snake_viewer.py   # Snake venom viewer
```

## Authors

- Ivan Koludarov (TUM)
- Todd A. Castoe (UT Arlington)

## References

Castoe TA, Daly M, Jungo F, Kirchhoff KN, Koludarov I, et al. (2025).
A Vision for VenomsBase: An Integrated Knowledgebase for the Study of Venoms and Their Applications.
*Integrative Organismal Biology*. doi:10.1093/iob/obaf026
