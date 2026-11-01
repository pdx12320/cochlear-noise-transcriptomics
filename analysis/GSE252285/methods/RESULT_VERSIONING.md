# Selected result status and annotation guardrails

This repository contains only the 36,111-cell robust-QC/no-doublet/t-SNE branch. Older complete 31,273-cell objects and result bundles are excluded. The minimal barcode/label prior required by the adopted annotation is supplied under `inputs/`, with its source hash and exact-mapping validation.

## Stored annotation

The selected object has 28 fixed resolution-0.6 clusters and the provisional field `broad_annotation_no_doublet_sensitivity`. Its counts are Immune_myeloid 27,676; UNRESOLVED 4,726; Fibrocyte_mesenchymal 2,997; Endothelial_vascular 275; Glial 264; Supporting_epithelial 149; Hair_sensory 24; Neuronal 0.

The marker audit later provided more informative recommendations for clusters 6/9/13/19, but it did not modify the object:

- clusters 6 and 9 show a complete erythroid program and no canonical neuronal program;
- cluster 13 shows a hematopoietic-progenitor program with endothelial-associated conflict;
- cluster 19 shows a fibrocyte/OMC/perivascular program;
- SGN is not recovered.

The revised table is an annotation recommendation, not an RDS writeback. Any analysis using the stored broad field—including the resolved-condition Top20 tables—must retain the stored `UNRESOLVED` group.

## Condition mapping

The suffix mapping is verified as `GEO_-3/LIB3=Control`, `GEO_-1/LIB1=HHL`, `GEO_-2/LIB2=NIHL`. It was applied after clustering and annotation, only for descriptive summaries. No inferential condition test was performed.

## Doublet status

No doublet filter is applied in this branch. The branch reconstructs exactly 36,111 robust-QC barcodes from QC metadata, checks that no scDblFinder field enters the new object, and records `doublet_filter_applied=NO`.
