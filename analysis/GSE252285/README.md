# GSE252285: expression localization in 36,111 cells

This module contains the selected robust-QC single-cell branch, descriptive condition summaries and matched-gene localization analyses used by the manuscript.

| Item | Selected setting |
|---|---|
| Source | [GSE252285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252285) |
| Public matrix | 32,285 features × 39,092 barcodes |
| Retained cells | 36,111; no doublet-based exclusions in this branch |
| Normalization | LogNormalize, scale factor 10,000 |
| Variable features | 2,000, VST |
| PCA and neighbors | 50 PCs calculated; PCs 1–30, k = 20 |
| Selected clustering | Louvain resolution 0.6, seed 42; 28 clusters |
| Display | t-SNE, perplexity 30 |

## Objects and metadata

The exact selected object is available through the root [release](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.0.0). Follow the [download and metadata guide](../../docs/RELEASE_OBJECTS.md).

Use `seurat_clusters_no_doublet` for the selected 28-cluster partition and `broad_annotation_no_doublet_sensitivity` for manuscript labels. The generic `seurat_clusters` column records a 37-cluster resolution audit. The generic `condition` column is unassigned; apply the external crosswalk below.

| Library | Barcode suffix | Condition | Cells |
|---|---|---|---:|
| LIB3 | -3 | Control | 10,273 |
| LIB1 | -1 | HHL | 16,631 |
| LIB2 | -2 | NIHL | 9,207 |

`library_provenance/` supplies accession/barcode-overlap evidence. These libraries support descriptive condition profiles. Cells are not treated as independent biological replicates for condition inference.

## Reproduce analyses from the selected object

Use the [portable entry points](portable_scripts/README.md) with the documented Seurat/UCell versions and frozen GMT files. The full 10,000-iteration localization results are in `complete_results/`. The portable copies change file interfaces and integrity guards; their scientific bodies match the executed scripts. Full analysis was not rerun when adapting these interfaces.

`ucell_programs/` and `matched_null/` retain executed source code and parameters. `scripts/condition_descriptive/` contains the fixed twenty-gene summary workflow. `tables/`, `figures/`, `methods/` and `logs/` contain selected current supporting records.

## Preprocessing and annotation provenance

The preprocessing sources are in `scripts/stage1/` and `scripts/no_doublet_sensitivity/`. Broad annotation combines reference-label concordance with canonical-marker support and conflict checks. It is not a de novo marker-only assignment. The [required input inventory](inputs/README.md) includes the supplied two-column annotation reference and its source hash. The source scripts document the reference-majority, purity and shared-cell rules. The adapter was tested against all 28 adopted cluster decisions without changing the scientific expressions. The later four-cluster marker audit did not overwrite the labels used in the manuscript.

Raw-matrix reconstruction requires its declared input/reference records and dependencies; the portable current-object analyses above are the direct release entry points. Current publication files use portable paths. Existing Git history is preserved; superseded result exports and old/new comparison outputs are excluded from this release.
