# VenomsBase: Data Architecture Proposal

> Response to: "Species-Centered Database Architecture for VenomsBase" (Castoe, March 2026)
>
> Ivan Koludarov — March 2026

## Summary

This document proposes a concrete data architecture for VenomZone, the web platform for VenomsBase. The central idea: **the protein is the primary data object**. Each protein carries its own embedding, sequences, annotations, and links to external data. Genomes are parallel first-class objects with their own viewer. All other data types — transcriptomics, proteomics, structural, functional — attach as features on proteins rather than living in separate relational tables.

Two working prototypes are attached (interactive HTML viewers for snake and bee venom proteins) to demonstrate how the protein-centric model works in practice.

---

## 1. Core Concept: Protein as Data Container

A protein language model (ProtT5) encodes each protein sequence as a 1024-dimensional vector. Proteins with similar sequences, structures, and functions land near each other in this space — regardless of species, naming convention, or database of origin. This gives us a continuous similarity landscape that no relational schema can provide.

Each protein is a self-contained record with typed data slots:

```
Protein Record
├── identity
│   ├── accession, gene name, aliases
│   └── cross-references (UniProt, GenBank, PDB)
├── sequences
│   ├── full-length (canonical)
│   ├── mature peptide
│   └── isoforms / splice variants
├── embeddings
│   ├── ProtT5 (1024-dim)
│   └── ESM2, ESM-C (optional)
├── structure
│   ├── AlphaFold prediction
│   └── experimental (PDB)
├── taxonomy
│   ├── species, TaxID, common name
│   └── clade, order, family
├── classification
│   ├── venom family + evidence level
│   ├── venom group / subgroup
│   └── historical classifications (e.g., Dowell 2016)
├── expression
│   ├── tissue-specific TPM/FPKM values
│   └── sample metadata (tissue, stage, method)
├── proteomics
│   ├── detection evidence (PSMs, confidence)
│   └── PTMs and proteoforms
├── activity
│   ├── molecular targets
│   ├── IC50 / EC50 values
│   └── assay metadata
├── genomic context  →  link to Genome Record
└── literature
    └── DOIs, PubMed IDs
```

New data types (e.g., spatial transcriptomics, cryo-EM) become new slots on the protein record — no schema migration required.

## 2. Three Data Object Types

The architecture uses three object types, not one hierarchy:

| Object | Role | Storage | Viewer |
|--------|------|---------|--------|
| **Protein Record** | Primary data container. Carries embedding, sequences, all annotations. | Structured file (Parquet) + embedding vectors (HDF5) | Embedding scatter (UMAP/PCA) + info panel |
| **Genome Record** | Standalone. Gene annotations, scaffolds, regulatory elements. Linked to proteins via gene IDs. | GFF3 + FASTA references | Genomic context viewer (gene neighborhood canvas) |
| **Raw Data** | Stored separately, referenced by protein/genome records. Sequencing reads, mass spectra, expression matrices. | Object storage (files) | Not directly visualized — features extracted into protein/genome records |

The key design choice: **genomes and proteins are parallel entities**, not hierarchical. A genome region contains multiple genes; each gene may correspond to a protein in the embedding space. They link to each other, but neither is subordinate.

This differs from the species-centered hierarchy (Species → Genome → Gene → Transcript → Protein) in an important way: the protein is independently queryable. You can find similar proteins across all species by nearest-neighbor search in embedding space, without traversing a join chain.

## 3. How -Omics Data Fits

All -omics data types follow the same pattern: **raw data stored as files, extracted features become slots on protein records**.

| Data Type | Raw Storage | Feature on Protein | How It's Used |
|-----------|-------------|-------------------|---------------|
| **Transcriptomics** | FASTQ/BAM files | `expression.*` (TPM, tissue, stage) | Color by expression level; filter by tissue |
| **Proteomics** | MS raw files | `proteomics.*` (PSMs, PTMs) | Evidence tags; proteoform annotations |
| **Structure** | PDB/mmCIF files | `structure.*` (AF confidence, method) | 3D viewer; structural clustering |
| **Functional assays** | Assay reports | `activity.*` (IC50, target, method) | Activity overlay; target enrichment |
| **Genomics** | GFF3 + FASTA | `genomic_context` → Genome Record link | Gene neighborhood viewer |
| **Literature** | — | `literature.*` (DOIs) | Cross-references |

This means every module in the species-centered schema has a concrete home:

| Species-Centered Module | Protein-Centric Equivalent |
|------------------------|---------------------------|
| Species & Taxonomy (2A) | `taxonomy.*` slots on protein record |
| Genome Assembly (2B) | Genome Record (standalone object) |
| Gene-Transcript-Protein (2C) | `identity.*` + `sequences.*` slots |
| Venom Classification (2C) | `classification.*` slots with evidence level |
| Functional Genomics (3) | Genome Record regulatory features |
| Transcriptomics (4) | `expression.*` slots |
| Proteomics (5) | `proteomics.*` slots |
| Structure (6) | `structure.*` slots |
| Gene Family (7) | `classification.*` slots |
| Functional Activity (8) | `activity.*` slots |

## 4. Why This Matters: The Reclassification Problem

Dowell et al. (2016) proposed a venom protein family classification that is widely used but increasingly recognized as inaccurate. Our studies (Koludarov et al. 2023) provide updated classifications, but the field still references Dowell's scheme.

**In a relational schema:** updating family assignments means modifying GeneFamily and FamilyMembership tables, cascading through OrthogroupMembership, and potentially breaking saved queries and cross-references.

**In the protein-centric model:** update the `classification.venom_family` slot. The old value is preserved as `classification.venom_family_dowell2016`. No schema changes, no cascading joins, both classifications coexist on the same record.

This generalizes to any reclassification event — nomenclature changes, family reassignments, new evidence. All are slot operations on protein records.

## 5. What the Prototypes Show

The attached HTML viewers demonstrate three interface components:

- **Protein Space:** ProtT5 embedding scatter (UMAP/PCA) as the navigational anchor. Click any protein to inspect it.
- **Information Panel:** All slots displayed on click — species, family, group, source dataset, gene name.
- **Genomic Context:** Gene neighborhood canvas showing venom genes (blue) and flanking genes (grey) from published GFF annotations, matched to the clicked protein.

**Snake viewer** (2,235 proteins): three families — Pla2g2 (452, Koludarov et al. 2020), KLK serine proteases (491, Barua & Koludarov 2021 BMC Biology), 3FTx (1,292 toxin proteins, Koludarov et al. 2023 Nature Comms). Genomic context from published manual GFF annotations.

**Bee viewer** (4,606 proteins): 33 protein families from the Hymenoptera venom arsenal (Koludarov et al. 2023 BMC Biology). Genomic context from 739 manually annotated GFF3 regions across 23 species.

The viewers are deliberately simple — no server, no database, just data and the idea. The point is that protein space works as a navigational interface, and all metadata naturally becomes "color by" dimensions.

## 6. Data Completeness Tiers

Following the platinum/gold/silver/bronze concept, defined by slot coverage:

| Tier | Required Slots | Example |
|------|---------------|---------|
| **Platinum** | embedding + sequence + taxonomy + structure + expression + proteomics + activity | Apamin (*Apis mellifera*) |
| **Gold** | embedding + sequence + taxonomy + structure + classification | Most ToxProt entries |
| **Silver** | embedding + sequence + taxonomy | Computationally predicted venom proteins |
| **Bronze** | sequence + taxonomy only | Unannotated transcriptome hits |

## 7. Implementation Path

**Phase 1 — Protein records + embedding viewer.** Load existing datasets (ToxProt, our published data) as protein records with embedding, sequence, taxonomy, and classification slots. Web interface: embedding scatter with color-by and info panel. This is what the prototypes already demonstrate.

**Phase 2 — Genome viewer + links.** Add genome records from GFF3 annotations. Link proteins to their genomic context. Interactive gene neighborhood canvas (already prototyped in the HTML viewers).

**Phase 3 — Expression + proteomics + activity.** Add quantitative data slots. Import from published supplementary tables, VenomZone, and community submissions.

**Phase 4 — Community features.** User accounts, protein-level comments and annotations, evidence curation, data submission workflows.

## 8. Next Steps

- Review the attached viewers and this data model
- Agree on the protein record format and slot ontology
- Define the initial scope: which species, which data types, what tier to target first
- Discuss infrastructure: hosting, data submission pipeline, community curation model

Happy to jump on a call to walk through it.
