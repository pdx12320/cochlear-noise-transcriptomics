# Full matched-null annotation sensitivity

The original 10,000-iteration analysis was repeated with the same frozen gene sets, expression/detection matching, gene-length sensitivity, seed, UCell scoring, localization statistic and six-test BH family. Only cluster 19 entered the fibrocyte/mesenchymal group in the derived annotation. All 12 significance classifications were unchanged: five of six tests were significant in each matching mode, while MAPK/AP1-to-endothelial localization remained nonsignificant.

For the fibrocyte TNFA/NFKB comparison, the observed localization delta increased from 0.0379826157 to 0.0393930206, with q = 0.0001999800 in both matching modes. The observed deltas of the other five pairs were unchanged. The primary MAPK/AP1 q value changed slightly from 0.0723927607 to 0.0724927507; its gene-length-matched q remained 0.0816918308. All other q values were unchanged. Thus, classification robustness should not be described as exact equality of every statistic.

Paired evaluation under the original labels reproduced all 120,000 archived null values with a maximum absolute error of 5.55e-17, and the six observed deltas with error below 5.21e-17. Independent Python recalculation of empirical P values and six-test BH adjustment agreed within 1.1e-19. The full run completed in 1,069.85 seconds with one compute thread and a peak resident memory of 7.24 GiB. The selected original object was read-only; its post-run SHA256 matched the initial and released checksum.

The comparison table is `null_final/C19_VS_ORIGINAL_MATCHED_NULL.tsv`. Full null distributions and program/mode checkpoints are retained separately. Figure S6 shows all tests, both matching modes, both annotation variants, and the nonsignificant MAPK results.
