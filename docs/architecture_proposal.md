# VenomsBase Data Architecture Proposal

> Response to: "Species-Centered Database Architecture for VenomsBase" (Castoe, March 2026)
>
> Authors: Ivan Koludarov, Todd A. Castoe
> Date: March 2026

## Executive Summary

We propose an **embedding-centric, protein-first** architecture for VenomsBase that
fulfills all requirements of the species-centered relational schema while adding a
computational layer that enables modern machine learning, visual exploration, and
cross-species comparison at scale.

The core insight: **the protein embedding IS the data object**. Sequence, structure,
taxonomy, expression, bioactivity, genomic context — these are all *tags* (metadata
dimensions) on the embedding, not the data itself. This inverts the traditional
database-first approach: instead of querying tables to find proteins, you navigate a
continuous embedding space where similar proteins cluster together regardless of
species or nomenclature.

This architecture is demonstrated with a working prototype using 2,300+ venom
proteins from three published studies spanning Hymenoptera and Squamata.

## Why Embedding-Centric?

Traditional knowledgebases organize data in relational tables. You search by name,
filter by species, join across tables. This works for known queries but fails for
the core challenge of venom research: **discovery of unknown relationships**.

Embeddings solve this because:
1. **Similar proteins cluster together** — even across species, even without shared nomenclature
2. **Every data type can be projected** — sequence, structure, and function all map to the same space
3. **Visual exploration is native** — researchers can literally see protein relationships
4. **New data integrates seamlessly** — drop a new protein in, it finds its neighbors automatically
5. **Cross-species comparison is implicit** — not a separate analysis step

## Architecture Overview

```
┌─────────────────────────────────────────────────────────┐
│                    VISUALIZATION LAYER                    │
│         ProtSpace / ProtXplore (WebGL scatter)           │
│    Click a point → see all tags → comment / annotate     │
└────────────────────────┬────────────────────────────────┘
                         │
┌────────────────────────┴────────────────────────────────┐
│                    EMBEDDING LAYER                        │
│           The protein IS a point in this space            │
│                                                          │
│   Per protein: embedding vector (ProtT5 / ESM-2)        │
│   Projections: UMAP 2D/3D, PCA, t-SNE                   │
│   Distance = similarity (sequence, structure, function)  │
└────────────────────────┬────────────────────────────────┘
                         │
┌────────────────────────┴────────────────────────────────┐
│                      TAG LAYER                            │
│        Every piece of metadata is a tag on the protein    │
│                                                          │
│   ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐  │
│   │ Identity │ │ Taxonomy │ │ Sequence │ │Structure │  │
│   │ uniprot  │ │ species  │ │ full_seq │ │ pdb_id   │  │
│   │ genbank  │ │ taxon_id │ │ mature   │ │ LDDT     │  │
│   │ refseq   │ │ family   │ │ signal_p │ │ TM_score │  │
│   │ name     │ │ genus    │ │ domains  │ │ RMSD     │  │
│   └──────────┘ └──────────┘ └──────────┘ └──────────┘  │
│   ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐  │
│   │  Venom   │ │ Genomics │ │Expression│ │ Activity │  │
│   │ family   │ │ gene_id  │ │ tissue   │ │ targets  │  │
│   │ group    │ │ gff_loc  │ │ TPM/FPKM │ │ assays   │  │
│   │ status   │ │ synteny  │ │ proteome │ │ toxicity │  │
│   │ evidence │ │ exons    │ │ enriched │ │ IC50     │  │
│   └──────────┘ └──────────┘ └──────────┘ └──────────┘  │
│   ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐  │
│   │Prediction│ │Phylogeny │ │Literature│ │  User    │  │
│   │ GO terms │ │ tree_id  │ │ doi      │ │ comments │  │
│   │ CATH     │ │ orthogrp │ │ pubmed   │ │ ratings  │  │
│   │ membrane │ │ clade    │ │ curated  │ │ flags    │  │
│   │ disorder │ │ bootstrap│ │ evidence │ │ notes    │  │
│   └──────────┘ └──────────┘ └──────────┘ └──────────┘  │
└─────────────────────────────────────────────────────────┘
```

## How This Maps to Todd's Schema

Every module in the species-centered schema maps directly to a tag category:

| Todd's Module | Our Tag Category | Storage |
|---------------|-----------------|---------|
| Species & Taxonomy (2A) | `taxonomy.*` | Inline tags |
| Genome Assembly (2B) | `genomics.assembly_*` | Reference pointer |
| Gene–Transcript–Protein (2C) | `identity.*` + `sequence.*` | Inline + FASTA |
| Venom Classification (2C) | `venom.*` | Inline tags (evidence-based, as Todd specifies) |
| Functional Genomics (3) | `genomics.regulatory_*` | Reference pointer to tracks |
| Transcriptomics (4) | `expression.*` | Inline values + reference pointer |
| Proteomics (5) | `proteomics.*` | Inline detection + reference pointer |
| Structure (6) | `structure.*` | Inline scores + PDB pointer |
| Gene Family / Comparative (7) | `phylogeny.*` | Inline assignments + tree pointer |
| Functional Activity (8) | `activity.*` | Inline values |
| Literature (8) | `literature.*` | DOI/PubMed pointers |

**Key compatibility:** Todd's evidence separation principle is preserved. Venom status
is a tag with an evidence level, not a binary flag. Tags can be updated, versioned,
and traced to their source.

**What we add:** The embedding layer that makes all of this *navigable*. You don't just
query "show me all PLA2 from snakes" — you see PLA2 clustering with PLA2 from bees
and discover shared evolutionary origins visually.

## The Protein Record

```python
@dataclass
class VenomProtein:
    """Core data object. Everything is a tag on the embedding."""

    # Identity (required)
    identifier: str            # Primary key (UniProt, GenBank, or custom)
    name: str                  # Human-readable name

    # Embedding (the actual data object)
    embedding: np.ndarray      # ProtT5 or ESM-2 vector (float16)
    projection_2d: Tuple[float, float]  # UMAP x, y
    projection_3d: Tuple[float, float, float]  # UMAP x, y, z

    # Tags (all optional — grow as data becomes available)
    tags: Dict[str, Any]
    # tags["taxonomy.species"] = "Naja naja"
    # tags["taxonomy.taxon_id"] = 8639
    # tags["taxonomy.family"] = "Elapidae"
    # tags["sequence.full"] = "MKTLLLTLVVVTIVCLDLG..."
    # tags["sequence.mature"] = "LECHNQQSSQPPTTK..."
    # tags["sequence.length"] = 74
    # tags["venom.family"] = "3FTx"
    # tags["venom.group"] = "short-chain"
    # tags["venom.status"] = "curated"
    # tags["venom.evidence"] = ["proteomics", "expression"]
    # tags["structure.pdb_id"] = "1IQ9"
    # tags["structure.lddt"] = 92.3
    # tags["activity.target"] = "nAChR"
    # tags["activity.ic50"] = "2.3 nM"
    # tags["expression.venom_gland_tpm"] = 15234.5
    # tags["phylogeny.orthogroup"] = "OG0001234"
    # tags["genomics.chromosome"] = "NW_020769389"
    # tags["genomics.start"] = 1234567
    # tags["prediction.go_mfo"] = ["GO:0090729"]  # toxin activity
    # tags["prediction.membrane"] = "Secreted"
    # tags["literature.doi"] = ["10.1038/s41467-023-xxxxx"]
```

## Data Completeness Tiers

Following Todd's platinum/gold/silver/bronze concept, but defined by tag coverage:

| Tier | Required Tags | Example |
|------|---------------|---------|
| **Platinum** | embedding + sequence + taxonomy + structure + expression + proteomics + activity | Apamin (Apis mellifera) |
| **Gold** | embedding + sequence + taxonomy + structure + venom classification | Most ToxProt entries |
| **Silver** | embedding + sequence + taxonomy | Computationally predicted venom proteins |
| **Bronze** | sequence + taxonomy only | Unannotated transcriptome hits |

## Storage Format: ProtSpace Bundle

The prototype ships as a `.parquetbundle` — the same format already used by ProtSpace
and ProtXplain. This means:

1. **Immediately visualizable** — drag into ProtSpace web app, explore
2. **Lightweight** — a few MB for thousands of proteins
3. **Self-contained** — annotations + projections + settings in one file
4. **Compatible** — works with existing ProtXplain/ProtSpace infrastructure

For the full database (Phase 2+), the bundle is the **export/visualization format**.
The system-of-record can be PostgreSQL (as Todd proposes) or SQLite, with tags stored
as JSON columns or a tag table. The bundle is generated from the database for
exploration.

## Prototype Datasets

| Dataset | Proteins | Source | Tags Available |
|---------|----------|--------|----------------|
| Bee venom (Koludarov 2023) | 372 | BMC Biology | taxonomy, family, group, length, proteomics, transcriptomics, activity |
| 3FTx/Ly6 (Koludarov 2023) | 1,426 | Nature Comms | taxonomy, family, group, membrane, GO, CATH, structure scores, sequences |
| Serine proteases (Barua & Koludarov 2021) | ~500 | BMC Biology | taxonomy, family, synteny, tissue expression, genomic annotations |
| + ToxProt reference | 8,055 | UniProtKB | taxonomy, family, function, sequence, cross-references |

**Total prototype: ~10,000 venom proteins** across all major venomous lineages.

## Implementation Phases

### Phase 1: Bundle Prototype (Current)
- Merge published datasets into unified protein table
- Generate ProtT5 embeddings for all proteins
- Build `.parquetbundle` with tag columns
- Ship to Todd as attachment — "here, explore this in ProtSpace"

### Phase 2: VenomsBase Backend
- PostgreSQL system-of-record (Todd's relational schema for provenance)
- REST API serving protein records with tags
- Bundle export endpoint (generate ProtSpace bundles on demand)
- ProtXplain integration: VenomsBase as a data source module

### Phase 3: Full Platform
- ProtXplain frontend: genomics viewer, species selector, protein space
- User comments/discussion per protein (ProtColab)
- Community curation workflows
- Species Embedding integration: venom proteins as first vertical

## Relationship to Species Embedding

VenomsBase is the first domain-specific implementation of a broader framework:

```
VenomsBase (venom proteins)
    ↓ generalizes to
ProtXplain (all proteins)
    ↓ generalizes to
Species Embedding (all biological data)
    H(s) = (h_hyp, h_euc) — product manifold per species
```

The protein embedding in VenomsBase is a component of the full species embedding.
As the framework matures, VenomsBase proteins gain additional context: their species'
position in taxonomy space (hyperbolic), their genomic neighborhood (synteny
embeddings), their expression profile (transcriptomic embeddings).

## Why Tags Beat Relational Joins: The Dowell Problem

A concrete example of why the tag architecture matters: Dowell et al. (2016) proposed
a venom protein family classification that is widely used but increasingly recognized
as inaccurate — even by Dowell himself. Our studies (Koludarov et al. 2023) provide
updated classifications, but the field still references Dowell's scheme.

In a relational database, changing the family classification means:
- Updating a `GeneFamily` table
- Cascading changes through `FamilyMembership`, `OrthogroupMembership`
- Rebuilding cross-references, breaking saved queries
- Risk of inconsistency during migration

In the tag architecture:
- Update `venom.family` tag on affected proteins
- Old value preserved in `venom.family_dowell2016` if needed for backward compatibility
- No schema changes, no cascading joins, no downtime
- Both classifications can coexist as parallel tag dimensions

This generalizes: **any reclassification is just a tag update**. Nomenclature changes
(which plague venom research), family reassignments, new evidence for venom status —
all are tag operations, not schema migrations.

## FAIR Compliance

| Principle | Implementation |
|-----------|---------------|
| **Findable** | Every protein has a persistent identifier; metadata searchable via tags |
| **Accessible** | Open-access bundles downloadable; REST API with authentication |
| **Interoperable** | Standard formats (Parquet, FASTA, GFF); UniProt/NCBI cross-references |
| **Reproducible** | Embedding model versions tracked; provenance tags on all derived data |

## Next Steps

1. **Immediate:** Build prototype bundle from 3 published datasets → ship to Todd
2. **Short-term:** Set up GitHub repo, define tag ontology, write REST API spec
3. **Medium-term:** ProtXplain module with VenomsBase backend
4. **Long-term:** Community platform with Species Embedding integration
