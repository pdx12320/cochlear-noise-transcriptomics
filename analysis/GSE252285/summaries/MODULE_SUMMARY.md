# GSE252285 robust-QC / no-doublet-filter / t-SNE sensitivity

- STATUS: SUCCESS
- ANALYSIS_STATUS: NO_DOUBLETFILTER_SENSITIVITY_NOT_FROZEN
- INPUT_DIMENSION: 32285 genes x 39092 cells
- ROBUST_QC_RETAINED: 36111
- DOUBLETS_REMOVED: 0 (filter not run or used in this branch)
- FINAL_CELLS: 36111
- NUMBER_OF_CLUSTERS_RES0.6: 28
- DISPLAY_EMBEDDING: t-SNE only; UMAP not run
- BROAD_CELLTYPE_COUNTS: Endothelial_vascular=275; Fibrocyte_mesenchymal=2997; Glial=264; Hair_sensory=24; Immune_myeloid=27676; Supporting_epithelial=149; UNRESOLVED=4726
- UNRESOLVED_CLUSTERS: 6,9,13,19
- TOP20_FOUND: 20/20
- APOLD1_CONTEXT: Endothelial_vascular (66.91% expressing); second=Fibrocyte_mesenchymal
- ERRFI1_CONTEXT: Fibrocyte_mesenchymal (51.75% expressing); second=Supporting_epithelial
- CONDITION_MAPPING_STATUS: CONDITION_MAPPING_RESOLVED
- CONDITION_MAPPING: GEO_-3/LIB3=Control; GEO_-1/LIB1=HHL; GEO_-2/LIB2=NIHL
- CONDITION_LEVEL_DE: NOT_PERFORMED
- ANNOTATION_FREEZE_COMPLETE_FOR_THIS_BRANCH: NO

## Interpretation warning

All 36,111 robust-QC cells were retained without a doublet filter. Clusters with mixed canonical lineage panels or unstable evidence were conservatively marked UNRESOLVED. t-SNE geometry was not used as annotation evidence.

The barcode-level reference labels were required annotation priors, combined with canonical-marker support, purity and competing-lineage checks. The minimal input is supplied in `inputs/annotation_reference.tsv`; all 28 cluster decisions were reproduced by the adapted reader. No scDblFinder score/class was imported or used for inclusion, clustering, markers or label assignment. Resolved conditions were added later for descriptive summaries only.

The original execution status strings (`NOT_FROZEN` and `ANNOTATION_FREEZE_COMPLETE_FOR_THIS_BRANCH: NO`) record that the broad biological annotation remained provisional. The exact selected object and current result tables are now versioned in the publication release; these status strings do not mean that the release object is unavailable.
