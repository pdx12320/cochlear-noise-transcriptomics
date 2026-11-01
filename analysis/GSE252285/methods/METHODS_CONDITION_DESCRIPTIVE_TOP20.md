# Resolved-condition descriptive Top20 methods

The existing 36,111-cell robust-QC Seurat object was read only. Fixed RNA_snn_res.0.6 clusters and broad_annotation_no_doublet_sensitivity labels were retained without rerunning or modifying QC, PCA, clustering or annotation.

Conditions were assigned directly from observed barcode suffix using the user-resolved mapping GEO_-3=Control, GEO_-1=HHL and GEO_-2=NIHL. The object cell-name suffix, barcode_suffix metadata, library_id and original library_map_template.tsv were required to form a complete one-to-one audit. The legacy object condition field remained UNKNOWN and was not overwritten.

For every fixed Top20 gene, expression was summarized for each resolved condition x existing broad cell type and condition x fixed cluster. n_expressing was defined as RNA counts > 0. pct_expressing was 100*n_expressing/n_total_cells. Mean and median expression used the Seurat RNA data layer and included zero-expression cells. Cell-type fraction was cell count divided by the total cells in that condition.

Groups with fewer than 20 cells were explicitly marked LOW_CELL_COUNT. All comparisons are descriptive. No condition-level differential expression, p-values, significance tests, differential-abundance inference or cell-as-replicate claims were produced. Bars start at zero; DotPlot area encodes detection percentage and colour encodes mean expression using common scales across compared condition facets.
