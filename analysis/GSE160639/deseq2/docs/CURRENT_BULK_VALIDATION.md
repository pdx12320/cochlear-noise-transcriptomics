# Current bulk analysis: independent statistical validation

The current chain starts from public GSE160639 raw counts, fits DESeq2 1.52.0 under R 4.6.1 with newly estimated median-ratio normalization, and passes the three current Wald-statistic tables to the fixed nine-pathway by three-contrast fgsea analysis. The locked 12 libraries (four per dose), 19,673 Ensembl IDs, annotations and contrast directions were preserved. No scientific threshold, inclusion label or fixed pathway membership was changed during this validation.

## Gene-level results

| Contrast | Genes with FDR < 0.05 | Positive | Negative | Nonmissing adjusted P values |
|---|---:|---:|---:|---:|
| 94 dB vs 70 dB | 5,777 | 3,104 | 2,673 | 18,882 |
| 105 dB vs 70 dB | 889 | 511 | 378 | 16,596 |
| 105 dB vs 94 dB | 4,177 | 1,908 | 2,269 | 15,836 |
| Overall dose LRT | 6,398 | — | — | 17,738 |

Each complete table contains 19,673 rows. The 28 missing nominal P values per model and the additional independently filtered genes remain explicit; no missing adjusted P value was treated as significant. BH was independently reproduced over each model's nonmissing adjusted-P family, with maximum absolute discrepancy 5.15e-15. Wald statistics equal log2 fold changes divided by their standard errors; two-sided normal P values and the LRT chi-square P values (two degrees of freedom) were independently checked.

The size-factor calculation was reproduced from the raw fixed-gene matrix: 19,558 genes were positive in every library, the maximum factor discrepancy was 5.33e-15, and raw counts divided by these factors matched the supplied normalized counts (maximum absolute rounding discrepancy 3.03e-9). Each baseMean matched the mean normalized count. All 45 retrieved run-file hashes matched the recorded source manifest. The fit completed in 37.927 seconds without model warnings.

## Targeted pathway results

All 27 pathway–contrast combinations were present. Global BH adjustment across the 27 nominal P values was independently reproduced (maximum discrepancy 1.03e-15); the separate nine-test within-contrast columns were also checked. Current global-FDR-significant counts are 8/9 for 94 vs 70 dB, 9/9 for 105 vs 70 dB, and 0/9 for 105 vs 94 dB.

The direct 105 vs 94 dB Hallmark reactive oxygen species pathway has NES −1.48760473644908, nominal P 0.0396711946617207, and global FDR 0.059506791992581. It is not significant at the specified FDR threshold. Its direction alone does not establish a direct dose difference.

All three input table hashes match the newly completed DESeq2 outputs. Each ranking independently reconstructs 17,441 unique symbols under the locked duplicate-selection rule. All pathway overlap counts and leading-edge members agree with the frozen 2026.1.Mm memberships. All 22 retrieved enrichment-file hashes matched. Enrichment completed in 5.042 seconds without warnings. See the two machine-readable validation JSON files for all checked fields.

## Published fitted objects

The Release object copies remove captured machine-local execution environments from the design formula and dispersion function. Every other slot, all assays, row/column metadata and extra attributes are identical to the original execution objects. Predictions of the parametric dispersion function at every fitted gene baseMean are identical. Re-extraction of all three Wald comparisons and the LRT produced identical result tables and result metadata. Source/public hashes and exact transformations are recorded in `PUBLIC_OBJECT_VALIDATION.tsv` and `PUBLIC_OBJECT_NOTE.txt`.

These checks establish computational consistency of the current chain. They do not address additional biological replication, sample influence or cell-type mixture. Bulk pathway context remains distinct from a cell-type-specific treatment effect; the separate single-cell study's replication and descriptive-localization limits are unchanged.
