# VenomsBase — Next Steps

## Monday Session (2026-03-24)

### Snake Viewer
- [ ] Add KLK/serine protease proteins with real ProtT5 embeddings (not just 3FTx)
- [ ] Add PLA2 data from `data/Pla2/` (452 ProtT5 embeddings already in `pla2_prott5.h5`)
- [ ] Merge 3FTx + KLK/SP + PLA2 into unified snake venom viewer
- [ ] GFF genomic context for all three families (SP GFFs exist, need PLA2 + 3FTx gene annotations)
- [ ] Generate ProtT5 embeddings for SP proteins (currently AA-composition fallback)

### Bee Viewer
- [ ] Generate real ProtT5 embeddings via BiocentralAPI (fix API call — `model` param changed)
- [ ] Better matching: map aligned FASTA headers → UniProt IDs from master CSV
- [ ] Add Apis mellifera venom gene genomic context (need GFF with venom loci from our data)
- [ ] Integrate proteomics/transcriptomics tags (from ODS tables → convert to CSV first)

### Prettify
- [ ] Title bar: "VenomsBase — Snake Venoms" / "VenomsBase — Bee Venoms"
- [ ] Legend improvements (collapsible, counts)
- [ ] Genomic panel: tooltip on gene hover
- [ ] Side panel: clickable UniProt/GenBank links
- [ ] Sequence display: color by residue property

### Send to Todd
- [ ] Write response email/document referencing architecture_proposal.md
- [ ] Attach both HTML viewers + architecture doc
- [ ] Mention Dowell reclassification argument
- [ ] Earmark: protein-level comments/discussion for logged-in users (future feature)
- [ ] Earmark: Workflow platform module integration (Step 3)

## Data Sources for Monday

| Dataset | Location | Proteins | Embeddings | Genomic |
|---------|----------|----------|------------|---------|
| 3FTx/Ly6 | data/3ftx_evolution/ | 1,427 | ProtT5 (SD4) | scaffold FASTAs |
| KLK/SP | data/serine_protease/ | ~500 | needs generation | GFF (TS.gff, Crvi_mi2.gff) |
| PLA2 | data/Pla2/ | 452 | ProtT5 (pla2_prott5.h5) | SyntenyDB + exon_homology |
| Bee venom | data/bee_venom/ | 4,634 | needs generation | needs Apis GFF |
| ToxProt | via UniProt/ToxProt | 8,055 | available | — |
