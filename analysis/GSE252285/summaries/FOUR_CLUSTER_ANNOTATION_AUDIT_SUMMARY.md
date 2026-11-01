# GSE252285 clusters 6/9/13/19 annotation audit

STATUS=SUCCESS

This focused audit is read-only with respect to the 36,111-cell Seurat object. Resolution 0.6, all 28 clusters, PCA and t-SNE were held fixed. No QC, PCA, UMAP, t-SNE or clustering was rerun, and no annotation metadata was written back.

## Adjudication

| cluster | cells | candidate label | confidence | key evidence | reference evidence |
|---:|---:|---|---|---|---|
| 6 | 2359 | UNRESOLVED (erythroid-like hypothesis) | UNRESOLVED | 9/9 erythroid genes; Car1/Ces2g high but canonical neuronal support 0/17 | 91.6% rejected; P7 reference lacks erythroid |
| 9 | 1210 | UNRESOLVED (erythroid-like hypothesis) | UNRESOLVED | 7/9 erythroid genes; no canonical neuronal support | closed-world Immune call conflicts with erythroid program; reference lacks erythroid |
| 13 | 838 | Immune_hematopoietic (progenitor-like; not myeloid-resolved) | MODERATE | 15/17 hematopoietic-progenitor genes; endothelial/hemogenic conflict | 100% accepted broad Immune; median score 0.843, margin 0.687 |
| 19 | 319 | Fibrocyte_mesenchymal (Type I OMC-like/perivascular) | HIGH broad | 8/13 fibrocyte genes; mixed vascular/pericyte/chondrocyte programs | 99.7% TypeI_OMCs; median score 0.966, margin 0.932 |

SGN_RECOVERED=NO
AUTHOR_LABEL_CONFLICT=YES
UNRESOLVED_REMAINING=6,9

## SGN decision

Cluster 6 expresses the paper's SGN-associated Car1 (72.3%) and Ces2g (74.9%), but its dominant program is erythroid (Klf1/Gypa/Gata1/Slc4a1/Alas2) and the canonical neuronal program has 0/17 supported genes. Cluster 19 has four partial neuronal-associated genes, but lacks the pan-neuronal Tubb3/Snap25/Rbfox3/Nefh/Prph core and maps to Type I OMCs. SGN is therefore not recovered.

## Technical warnings

- No audited cluster crossed the prespecified library/QC-driver flag (library fraction >=0.90, total-variation >=0.50, or absolute robust median shift >=3 MAD).
- Library composition is imbalanced: LIB1 is only 6.2% of cluster 6 and 2.1% of cluster 9; cluster 19 is 73.4% LIB3. The mapping is now resolved, but independent biological replication remains unresolved, so these shifts are descriptive and not inferential.
- GSE217727 is a P7 whole-cochlea closed-world reference and lacks an erythroid class. Its Immune prediction for cluster 9 is explicitly treated as a reference-space failure, not lineage proof.
- Cluster 13 has a strong progenitor program and broad reference agreement, but its endothelial/hemogenic signature prevents a mature myeloid label.
- FindMarkers p-values are descriptive because clustering and testing use the same data. Labels use effect size, pct.1, pct.2, pct.1-pct.2, multi-gene coherence and reference evidence.
- Exact author cluster/cell metadata were not found in the project; author marker comparisons are not claimed as cell-level author label transfer.

## Reference

GSE217727, P7 C57BL/6J whole mouse cochlea, reconstructed 15-subtype reference; Seurat anchor transfer with score <0.50 or top-two margin <0.10 rejected. Reference PCA was rebuilt; target PCA was not rerun.

## Key files

- `04_tables/.../FOUR_CLUSTER_FINAL_ANNOTATION_AUDIT.tsv`
- `04_tables/.../FOUR_CLUSTER_TOP50_ANNOTATION_MARKERS.tsv`
- `04_tables/.../FOUR_CLUSTER_LINEAGE_PROGRAM_SCORES.tsv`
- `04_tables/.../FOUR_CLUSTER_REFERENCE_MAPPING_SUMMARY.tsv`
- `04_tables/.../FOUR_CLUSTER_TECHNICAL_AUDIT.tsv`
- `05_figures/.../FOUR_CLUSTER_LINEAGE_MARKER_DOTPLOT.pdf/png`
- `05_figures/.../FOUR_CLUSTER_SGN_AUTHOR_VS_CANONICAL_DOTPLOT.pdf/png`
- `05_figures/.../FOUR_CLUSTER_REFERENCE_MAPPING_HEATMAP.pdf/png`
