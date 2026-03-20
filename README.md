# VenomsBase Data Architecture

Embedding-centric data architecture for VenomsBase — a community knowledgebase for venom research.

This repository contains the architecture proposal, data model, and reference implementation for organizing venom data as protein-centric embeddings with rich metadata tags.

## Vision

Everything becomes an embedding. Proteins, genes, species — each is a point in a learned manifold. Sequence, structure, taxonomy, expression, bioactivity, genomic context — these are all *tags* on the embedding, not the data itself. The embedding IS the data.

VenomsBase proves this pattern for venomics.

## Architecture

See `docs/architecture_proposal.md` for the full proposal.

## Data

Data files are not tracked in this repository. Use the download scripts:

```bash
python scripts/download_data.py          # Download reference datasets
python scripts/setup_local_symlinks.py   # Symlink to local parent repository data (dev only)
```

## Authors

- Ivan Koludarov (TUM)
- Todd A. Castoe (UT Arlington)

## References

Castoe TA, Daly M, Jungo F, Kirchhoff KN, Koludarov I, et al. (2025).
A Vision for VenomsBase: An Integrated Knowledgebase for the Study of Venoms and Their Applications.
*Integrative Organismal Biology*. doi:10.1093/iob/obaf026
