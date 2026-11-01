# Animal-deletion sensitivity analysis

The complete four-versus-four-animal analysis remains the primary result. Forty-eight leave-one-animal-out fits evaluated six tissue–cell-type models, each retaining three animals from one condition and four from the other. The seven displayed pathway/model combinations cover TNF–NF-kB in all six models and MAPK-family signaling in lateral-wall endothelial cells. All individual animal identities, directions, NES and FDR values are supplied in `LOO_FOCUSED_ALL_56.tsv`.

| Tissue | Cell type | Pathway | Full NES | Full FDR | Same direction | LOO FDR < 0.05 | LOO NES range | LOO FDR range |
|---|---|---|---:|---:|---:|---:|---:|---:|
| SGN | Fibrocyte / mesenchymal | TNF–NF-kB | 1.463 | 0.0305 | 8/8 | 7/8 | 1.305 to 1.731 | 0.000326 to 0.136 |
| SGN | Endothelial | TNF–NF-kB | -1.854 | 7.62e-05 | 8/8 | 8/8 | -2.312 to -1.491 | 3.49e-10 to 0.0178 |
| SGN | Immune / myeloid | TNF–NF-kB | -2.375 | 1.38e-12 | 8/8 | 8/8 | -3.132 to -1.977 | 8.59e-32 to 1.44e-06 |
| Lateral_wall | Fibrocyte / mesenchymal | TNF–NF-kB | 1.648 | 0.00371 | 8/8 | 6/8 | 1.317 to 2.007 | 3.94e-06 to 0.169 |
| Lateral_wall | Endothelial | TNF–NF-kB | -2.813 | 4.27e-20 | 8/8 | 8/8 | -3.051 to -2.649 | 9.23e-29 to 3.67e-16 |
| Lateral_wall | Endothelial | MAPK-family signaling | 1.918 | 8.37e-07 | 8/8 | 8/8 | 1.807 to 1.925 | 5.5e-07 to 2.08e-05 |
| Lateral_wall | Immune / myeloid | TNF–NF-kB | -3.213 | 1.78e-32 | 8/8 | 8/8 | -3.338 to -2.927 | 3.12e-38 to 2.32e-24 |

## Interpretation

The eight deletion estimates for each model overlap in seven of eight source animals and must not be treated as eight independent replications. Direction stability and changes in adjusted P values are reported separately. Removing one animal changes information, normalization and dispersion estimates; crossing FDR 0.05 alone is not a direction reversal. No post hoc animal exclusion or replacement of the main analysis was made. These separate cell-type fits do not test a condition-by-cell-type interaction.

The original full-model fitted gene mask was retained. Normalization and dispersion were re-estimated for each seven-animal dataset using the original PyDESeq2 design and settings. Enrichment retained the frozen original Hallmark-plus-Reactome memberships and the original eligible testing family for each model. BH family sizes stayed fixed, including unavailable members where applicable. The LOO FDR values are family-specific sensitivity summaries; there was no pooled confirmatory FDR across deletion iterations.

Per-animal expression summaries use the original full-fit normalized pseudobulk counts. Representative-gene values and all displayed pathway-member counts are available with sample identities. The pathway expression profile is the unsigned equal-weight mean of gene-standardized log1p normalized counts, using sample SD across the original eight animals. It is a descriptive score, not NES, a cell-level test or a power calculation.

The finite sample size and broad cell annotations remain limitations. This exercise evaluates influence of observed animals under the existing annotation and modeling choices; it does not establish biological replication in the separate pooled-library localization dataset.

## Observed stability and warnings

All 56 focused deletion estimates retained the full-analysis direction. Fibrocyte/mesenchymal TNF–NF-kB remained below FDR 0.05 in seven of eight SGN deletions and six of eight lateral-wall deletions. The three exceptions were SGN_noise_1 (GSM5124295; NES 1.305459, FDR 0.135819), Lateral_wall_naive_2 (GSM5124300; NES 1.317430, FDR 0.169287) and Lateral_wall_noise_3 (GSM5124305; NES 1.357393, FDR 0.117546). The four endothelial/immune TNF models and lateral-wall endothelial MAPK model remained below FDR 0.05 in all eight deletions.

No PyDESeq2 warnings were recorded. fgsea reported tied ranking statistics in every deletion fit; the original lexical secondary ordering was retained, and the warnings and tie counts are supplied. These warnings are not evidence that bias occurred. All original pathway families remained testable. At most five retained genes became all zero after a deletion; they remain in the complete DE tables, and undefined statistics are omitted from the corresponding rankings as in the original workflow.
