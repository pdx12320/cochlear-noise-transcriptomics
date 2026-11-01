# Inputs for the selected GSE252285 workflow

`annotation_reference.tsv` is the minimal reference required to reproduce the adopted broad annotations: 31,273 unique barcodes and their reference broad label. It was extracted read-only from the earlier annotation object. That complete object and its result bundle are not needed by the current entry point. `annotation_reference_provenance.json` records the exact source-object hash, extracted-table hash, field mapping and validation results.

The prior labels enter the annotation decision. Within each current resolution-0.6 cluster, the existing code computes the reference majority label, shared-cell count and label purity, then applies its unchanged canonical-marker and competing-lineage rules. The adapter changes input storage only. It preserves the adopted 28-cluster annotations, including the four unresolved clusters.

## Complete required-input inventory

All run-directory paths below are relative to a new `GSE252285_PROJECT_ROOT`.

| Stage | Required input | Source |
|---|---|---|
| Upstream input audit (`stage1/01_input_qc.R`) | `GSE252285_matrix.mtx.gz`, `GSE252285_features.tsv.gz`, `GSE252285_barcodes.tsv.gz` | Public GSE252285 processed Matrix Market files; set `GSE_MATRIX`, `GSE_FEATURES`, `GSE_BARCODES`, or place under `01_input/` |
| Upstream QC (`stage1/02_seurat_qc.R`) | `02_processed/01_counts_matrix.rds`, feature map and `00_manifest/library_map_template.tsv` | Produced by stage 01; the template retains the original suffix/library mapping and `UNKNOWN` condition |
| Current QC/clustering/t-SNE (`13a`) | `06_objects/01_raw_seurat.rds`, `04_tables/QC_CELL_METADATA.tsv.gz`, `04_tables/QC_THRESHOLDS.tsv` | Produced by stage 02; 13a reads the exact robust-QC pass set and ignores doublet fields |
| Current markers/annotation/Top20 (`13b`) | `06_objects/no_doublet_robust_qc_tsne_sensitivity/01_no_doublet_clustered_tsne.rds` | Produced by 13a |
| Current annotation prior (`13b`) | This repository's `inputs/annotation_reference.tsv` | Supplied; override location with `ANNOTATION_REFERENCE_TSV` |
| Optional marker checkpoint (`13b`) | Branch `04_tables/ALL_CLUSTER_MARKERS.tsv.gz` | Reused only if already present; otherwise markers are computed by the unchanged code |
| Resolution audit and figure QA | Selected `02_no_doublet_annotated_tsne.rds`; current resolution/annotation tables | Produced by 13b and `resolution_stability_audit.R` |
| Selected-branch validation | Current object, current Top20/composition/annotation/resolution tables and figures | Produced by preceding selected-branch stages |
| Resolved condition summaries (`17`) | Selected `02_no_doublet_annotated_tsne.rds` and original `00_manifest/library_map_template.tsv` | Selected object is also a release asset; template is produced by stage 01. Verified external suffix-to-condition mapping is encoded and checked in the script |
| Condition finalization (`18`) | Current condition-analysis tables, object hashes and validation records | Produced by stage 17 |
| Portable UCell programs | Selected release object; `mh.all.v2026.1.Mm.symbols.gmt`, `m2.cp.reactome.v2026.1.Mm.symbols.gmt` | Release asset and supplied `complete_results/` files |
| Portable matched-null analysis | Same object/GMT files plus GENCODE vM23 GTF | Supplied object/GMT; GTF source/hash are documented in `portable_scripts/README.md` |

A new raw-input reconstruction needs only upstream stages 01 and 02 before entering the current 13a/13b workflow. The historical `stage1/run_all.sh` also invokes stages 03–06, which are unnecessary for this selected branch. Stage 02 computes its recorded doublet audit; the current 36,111-cell branch does not use those classes to remove cells. Use the documented R/package versions; this adaptation did not rerun QC, clustering or marker estimation.

## Running the selected branch

From `analysis/GSE252285/`, after preparing the three raw/QC prerequisites:

```bash
export GSE252285_PROJECT_ROOT=/path/to/new/GSE252285_run
export PROJECT_DIR="$GSE252285_PROJECT_ROOT"
export R_BIN=/path/to/Rscript
export SC_THREADS=4
bash run_all.sh
```

The wrapper locates `inputs/annotation_reference.tsv` relative to its script directory and includes it in before/after input checksums. It no longer requires `06_annotated_final_audit_v3.rds`.

Old/new comparison exports are disabled by default (`RUN_LEGACY_COMPARISON=false`). Optional comparison mode requires an additional `frozen_cluster` column in a user-supplied reference table; that column is not required for annotation and is not included in the minimal release input. The label purity and shared-cell calculations needed by the adopted annotation always run.

## Separate descriptive annotation audits

The optional four-cluster audit scripts are not prerequisites for 13a/13b, condition summaries or UCell scoring. Their later recommendations were not written into the selected object:

- `14_four_cluster_annotation_evidence.R` additionally needs the separately reconstructed GSE217727 P7 reference object and `REFERENCE_GSE217727_PROVENANCE.tsv`.
- `14b_finalize_four_cluster_annotation_audit.R` consumes that audit's marker, lineage, reference-mapping and technical tables, plus its before/after checksums.
- `15_targeted_erythroid_cl13_audit.R` additionally reads `04_tables/DOUBLET_RESULTS.tsv.gz` and `QC_THRESHOLDS.tsv` for descriptive auditing. These audit-only doublet values do not change selected cell inclusion or labels.

Those auxiliary reference/input dependencies are stated explicitly; the archived audit outputs are not evidence that every auxiliary audit was rerun during publication preparation.

## Validation of this adaptation

The original annotation decision expression was evaluated with current cluster memberships and existing lineage-score records after reconstructing its prior mapping from the two-column reference. All 28 cluster labels, confidence values, reasons and warnings matched the adopted table exactly; the reference majority fractions and additional-cell fractions matched within 1e-12. The 13a file, marker/lineage computations, annotation-decision expression and object-label assignment remain unchanged. This is metadata/interface validation, not a new full scientific analysis.
