# METHODS — GSE252285 UCell program localization

- Input object: `/path/to/cochlear_analysis/GSE252285_stage1_20260814/06_objects/no_doublet_robust_qc_tsne_sensitivity/02_no_doublet_annotated_tsne.rds` (read only).
- Fixed object: 32,285 RNA features and 36,111 cells.
- Annotation field: `broad_annotation_no_doublet_sensitivity`; values and ordering were not modified.
- No QC, normalization, dimensional reduction, clustering, doublet filtering, or annotation was rerun.
- Existing t-SNE coordinates were reused only for visualization.
- Species: Mus musculus.
- Gene sets: frozen Mouse MSigDB 2026.1.Mm Hallmark and Reactome symbol GMT files.
- Programs were declared before scoring as exact official sets or deterministic unions listed in `PROGRAM_DEFINITIONS.tsv`.
- No hand-written genes were added; APOLD1/ERRFI1 were not used to define or adapt any signature.
- UCell input: complete sparse RNA counts matrix; no DEG or gene prefiltering before cell ranking.
- UCell parameters: maxRank=1500, ties.method=average, chunk.size=500, ncores=4.
- Seed: 20260823. UCell itself is deterministic for these settings.
- Broad summaries: median, first quartile (Q1), third quartile (Q3), and IQR per fixed broad cell type.
- Gene coverage: object-level signature matches; per cell type, genes detected in at least one cell and median per-cell matched-gene detection fraction.
- Top cell types were ranked by raw median UCell score without excluding low-cell-count or UNRESOLVED groups; low cell counts are flagged.
- Control/HHL/NIHL mapping follows the user-resolved GEO barcode suffix mapping (-3/-1/-2). Conditions are descriptive facets only.
- No condition-level P values, matched-null analysis, or cell-as-replicate inference was performed.

## Versions

- R: 4.4.3
- UCell: 2.10.1
- Seurat: 5.5.1
- SeuratObject: 5.4.0
- Matrix: 1.7.6
- data.table: 1.18.4
- ggplot2: 4.0.3
- patchwork: 1.3.2

## Frozen gene-set files

- `/path/to/cochlear_analysis/GSE252285_UCELL_GSE160639_PROGRAMS_20260823/gene_sets/mh.all.v2026.1.Mm.symbols.gmt`
- `/path/to/cochlear_analysis/GSE252285_UCELL_GSE160639_PROGRAMS_20260823/gene_sets/m2.cp.reactome.v2026.1.Mm.symbols.gmt`
