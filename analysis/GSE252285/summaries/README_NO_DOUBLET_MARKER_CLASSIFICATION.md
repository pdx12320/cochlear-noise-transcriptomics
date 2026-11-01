# GSE252285 all-cluster marker classification

This local bundle summarizes all 28 clusters in the 36,111-cell robust-QC, no-doublet-filter, t-SNE sensitivity branch at resolution 0.6.

## Current marker-revised recommendation

- Clusters 6 and 9: `Other / Erythroid-like`, HIGH confidence for erythroid identity. They may instead use a new broad label, `Erythroid_hematopoietic`, if the atlas taxonomy is expanded.
- Cluster 13: `Immune_hematopoietic / hematopoietic progenitor-like`, MODERATE. Endothelial-associated expression remains a conflict; it is not forced to myeloid progenitor or hemogenic endothelial.
- Cluster 19: `Fibrocyte_mesenchymal / Type-I-OMC-like, perivascular`, HIGH.
- SGN/neuronal cells were not recovered because no cluster has a coherent canonical neuronal program.

The recommended classification is recorded in `ALL_CLUSTERS_MARKER_CLASSIFICATION_REVISED.tsv`. Broad counts are in `BROAD_CLASSIFICATION_COUNTS_REVISED.tsv`.

## Source material

- `source_tables/` contains the complete existing annotation, marker, Top20, composition, resolution and focused-audit tables from the remote analysis branch.
- `source_figures/` contains the complete PDF and PNG figure set from that branch, including the focused clusters 6/9/13/19 audit.
- `source_release/` contains the original branch and focused-audit summaries.

## Important status

This is a marker-audit summary, not an RDS writeback. The existing Seurat object still contains the earlier sensitivity labels in which clusters 6, 9, 13 and 19 are `UNRESOLVED`. No RDS is included in this bundle and no upstream analysis was rerun.

`CONDITION_MAPPING_RESOLVED=YES` (`GEO_-3/LIB3=Control`, `GEO_-1/LIB1=HHL`, `GEO_-2/LIB2=NIHL`); no condition-level differential expression was performed.
