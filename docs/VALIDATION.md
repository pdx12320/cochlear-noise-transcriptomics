# Validation record

## Current manuscript

The Scientific Reports manuscript contains 18 pages, Supplementary Information 9 pages, and the cover letter 1 page. All 28 pages were rendered and visually inspected. The manuscript has 36 references, 34 DOI links and 59 citation instances, in sequential numerical order. The abstract has 146 words and the title 12 words. All eight supplementary figures and Table S1 are cited. Current document hashes and checks are recorded in `manuscript/BUILD_VALIDATION.json`.

All 50 panel letters are bold lowercase. The figure audit verified identical source tables, unchanged SVG geometry and styles, and identical pixels outside the letter boxes. Bibliographic records and scientific analysis outputs are unchanged. Reference presentation, document structure, figure-panel case and supplementary figure numbering were adjusted.

The standalone LaTeX package was verified by CRC and per-file hashes, then compiled after fresh extraction. The archive hash is recorded in `manifests/RELEASE_FILES_V1_2_2.tsv`.

## Execution scope

| Component | Execution for this release | Checks |
|---|---|---|
| GSE160639 normalization and differential expression | New fit from public counts in an isolated R 4.6.1 / DESeq2 1.52.0 environment | Exact input hashes; 12 samples; 19,673 fixed genes; three Wald contrasts and one overall-dose LRT |
| GSE160639 targeted enrichment | New calculations from the new Wald statistics | Nine fixed pathways across three contrasts; frozen MSigDB 2026.1.Mm membership; one global BH correction across 27 tests |
| GSE168041 single-cell and pseudobulk analyses | Adopted saved analysis retrieved | Two annotated AnnData objects; 13 sample-level differential models; 97 supporting files matched their analysis-directory hashes |
| GSE168041 offline enrichment interface | Input staging and R parsing executed | Frozen 82,288-row membership export; scientific analysis body unchanged from the archived script |
| GSE252285 single-cell analysis and matched-null tests | Adopted saved analysis retrieved | 36,111-cell object; selected 28-cluster field; seven stored broad labels; condition crosswalk and saved null-test records |
| GSE252285 portable UCell and matched-null interfaces | R parsing and code comparison executed | Scientific bodies unchanged; explicit paths and input-hash checks added |

The full single-cell preprocessing, pseudobulk fitting, UCell and matched-null workflows were not repeated during archive preparation. Their delivered values are the adopted saved analyses. The executable entry points and required reference inputs are supplied for subsequent reproduction.

## Current bulk checks

All 45 retrieved GSE160639 raw-count-analysis files matched their source hashes. An independent calculation reproduced the median-ratio size factors (maximum absolute difference 5.33 × 10⁻¹⁵), normalized counts (3.03 × 10⁻⁹), base means, Wald/LRT P values and the four BH-adjusted gene families (maximum absolute difference 5.15 × 10⁻¹⁵). This audit checks the arithmetic and provenance of the executed fit.

| Contrast or test | Genes evaluated | Genes with BH-adjusted P < 0.05 |
|---|---:|---:|
| 94 dB vs 70 dB | 19,673 | 5,777 |
| 105 dB vs 70 dB | 19,673 | 889 |
| 105 dB vs 94 dB | 19,673 | 4,177 |
| Overall-dose LRT | 19,673 | 6,398 |

Targeted enrichment uses 17,441 unique ranked gene symbols per contrast. Seventeen of the 27 pathway comparisons meet global FDR < 0.05: eight for 94 dB vs 70 dB, nine for 105 dB vs 70 dB, and zero for 105 dB vs 94 dB. Every result, including the full 27-test family and leading-edge memberships, is supplied.

## Single-cell objects and supporting files

The annotated objects were read and inspected: GSE168041 SGN contains 31,556 cells, lateral wall contains 33,149 cells, and GSE252285 contains 36,111 cells. The original object bytes were preserved. All three SHA256 values matched the source analysis objects.

The GSE252285 selected clustering field is `seurat_clusters_no_doublet` (28 clusters at resolution 0.6). The generic `seurat_clusters` field records the separate 37-cluster resolution audit. Condition summaries use the archived GEO/library crosswalk: Control 10,273 cells, HHL 16,631 cells and NIHL 9,207 cells.

The 2.66-GB Seurat object was split into numbered release parts and reassembled locally. The reconstructed SHA256 matched the original. Sixteen objects, including the 13 pseudobulk model objects, were inspected for embedded credential and private-path strings. Publication copies of two supporting tables replace machine-specific directory fields with relative paths; scientific values remain unchanged. These transformations are listed in `manifests/RELEASE_VERIFICATION.json`.

## Portable workflow and object checks

The minimal GSE252285 annotation-reference table contains 31,273 barcode/label rows. Reading this table through the portable adapter reproduced all 28 cluster labels, confidences, decision reasons, warnings and cell counts. The marker/lineage calculations and annotation-decision expressions were preserved; preprocessing stage 13a was unchanged. Its detailed equivalence and mock-command checks are under `analysis/GSE252285/inputs/`.

The public DESeq2 model copies preserve every assay, row/column metadata field and model result. Design and dispersion-function environment references were reduced to remove captured machine-local arguments. Dispersion predictions across the fitted gene base means and all three Wald/one LRT result extractions were identical to the original objects. Source/public hashes and field-level checks accompany the release bundle.

## Archive integrity

`manifests/REPOSITORY_FILES.tsv` records the path, byte count and SHA256 of every versioned payload file other than the manifest itself. `python3 scripts/verify_archive.py` checks these bytes and dataset invariants. Release-object byte counts and hashes are listed separately in `manifests/RELEASE_OBJECTS.tsv`.

The asset retrieval and reassembly helpers were tested with local fixtures for successful reconstruction, repeated execution, preservation of corrupt existing files, partial-download cleanup and path-traversal rejection. They use Python's standard library.

The source repositories, software versions and statistical parameters are recorded with each analysis.


## Sensitivity analysis validation

All 48 GSE168041 omission fits completed. Audits checked retained sample identities, fixed gene masks, size factors, Wald P values, gene-level BH, rank ordering, gene-set overlap, leading edges and the original family-specific BH adjustment. All 56 focused pathway estimates retained their direction. The primary four-versus-four comparisons remain unchanged.

The GSE252285 sensitivity changes only the derived label of cluster 19 (319 cells). Original-label evaluation reproduced all 120,000 archived matched-null statistics with maximum absolute error 5.55e-17. Independent empirical-P and per-mode six-test BH calculations agreed within 1.1e-19. The original eleven localization assignments and all twelve significance classifications were retained. Source-object hashes are unchanged.

The revised manuscript uses six main figures and eight supplementary figures, with 36 references and 34 distinct DOI links. LaTeX logs, numerical source comparisons, PDF link checks and rendered-page inspection are recorded in `manuscript/BUILD_VALIDATION.json` and `manuscript/STATISTICAL_VALIDATION.json`. Scripted archive checks validate file integrity and declared structure; they do not rerun upstream models.
