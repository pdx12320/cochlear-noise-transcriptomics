# Cluster 19 annotation sensitivity: expression and cached UCell summaries

This derived sensitivity reassigns only the 319 cells in selected resolution-0.6 cluster 19 from `UNRESOLVED` to `Fibrocyte_mesenchymal`. The original Seurat object, primary annotation, all other cells, expression values and cached UCell scores are unchanged. Fibrocyte/mesenchymal cells increase from 2,997 to 3,316 and unresolved cells decrease from 4,726 to 4,407. Cluster 19 contains 234 Control, 25 HHL and 60 NIHL cells.

The later marker/reference audit supports the broad mesenchymal assignment. Cluster 19 expresses Coch (99.69%), Otor (82.13%), Pou3f4 (75.86%), Col1a2 (83.07%) and Pdgfra (74.61%), with low Pecam1 (1.88%), Kdr (9.09%) and undetected Vwf. Its earlier unresolved status reflected a conflict between a fibrocyte reference prior and the original canonical-lineage score, rather than absence of fibrocyte evidence. The later reference audit mapped 99.69% of these cells to Type-I OMC-like cells. The sensitivity uses only the supported broad label and does not force a subtype.

All 11 originally reported genes retain their highest mean-normalized-expression broad category in each of the three libraries after reassignment. Cebpb additionally satisfies this criterion in the sensitivity because its NIHL maximum shifts from myeloid to fibrocyte/mesenchymal. The original main-analysis set remains the same 11 genes; the sensitivity is reported separately. Across the full twenty-gene panel, pooled Fos switches from fibrocyte/mesenchymal to endothelial and pooled Jun switches from endothelial to fibrocyte/mesenchymal. Detection-fraction rankings show further changes, including pooled Errfi1 shifting to supporting epithelium. These do not change the mean-expression definition used for the original eleven genes.

Apold1 remains highest among the tested broad categories by both mean normalized expression and detection fraction in pooled cells and in each library. Its pooled endothelial mean is 1.3375 and detection fraction 66.91%, while the expanded fibrocyte/mesenchymal group changes from 0.2416 to 0.3709 and from 15.82% to 21.59%. Apold1 is nevertheless strongly expressed within cluster 19 itself (mean 1.5854; 242/319 cells, 75.86%). The supported conclusion is robust ranking at the broad-category level; these data do not support endothelial-specific expression.

The cached six-program UCell scores are byte-identical to those in the user-supplied original integrated archive. None was rescored for this descriptive regrouping. All 24 highest-median assignments (six programs × pooled plus three library summaries) remain unchanged. The pooled assignments are: endothelial for MAPK/AP1, VEGF/angiogenesis, TNFA/NFKB and IL6/JAK/STAT3; unresolved for MYC and E2F. This is a descriptive result. Matched-null significance must be evaluated separately with the revised grouping because both the fibrocyte median and the comparator statistic can change.

Numerical checks reproduced 560 original expression-summary rows (maximum mean error 1.1e-14; detection-fraction error 1.0e-13) and 42 original UCell broad-summary rows (maximum median/quartile error 1.0e-15). No condition-level inferential test was added.

## Tables and plotting inputs

- `results/ALL20_TOP_ASSIGNMENTS_ORIGINAL_VS_C19.tsv`: all 160 mean/detection × pooled/library × gene comparisons, including unchanged assignments.
- `results/ALL20_MEAN_LOCALIZATION_STABILITY.tsv`: original main-eleven retention and newly stable genes for all twenty candidates.
- `results/TOP20_EXPRESSION_BOTH_ANNOTATIONS.tsv`: group sizes, numbers expressing, detection fractions, raw-count means, normalized means and ranks.
- `results/CLUSTER19_MARKER_EXPRESSION.tsv`: cluster 19 markers by library and pooled.
- `results/UCELL_SUMMARY_BOTH_ANNOTATIONS.tsv`: all six programs, all broad categories, pooled and by library.
- `results/CLUSTER19_ORIGINAL_ANNOTATION_EVIDENCE.tsv` and `CLUSTER19_LATER_AUDIT_EVIDENCE.tsv`: exact existing annotation evidence records.
- `figure_inputs/CELL_TSNE_ANNOTATION_EXPRESSION_UCELL.tsv.gz`: 36,111 cells with fixed t-SNE coordinates, cluster IDs, both label variants, resolved condition, 20 genes plus 8 lineage markers, and 6 cached UCell scores. Expression columns use the unchanged RNA data layer.
