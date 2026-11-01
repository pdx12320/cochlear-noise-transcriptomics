# METHODS — matched-gene empirical null for GSE252285 UCell localization

- Fixed input object: `/path/to/cochlear_analysis/GSE252285_stage1_20260814/06_objects/no_doublet_robust_qc_tsne_sensitivity/02_no_doublet_annotated_tsne.rds` (read only; 32,285 genes x 36,111 cells).
- No QC, normalization, clustering, annotation, or condition-level analysis was run.
- Tested only the six prespecified program/cell-type localizations; MYC and E2F were excluded.
- Gene universe: genes in the RNA counts layer detected in at least one of 36,111 cells.
- Effective signature: official frozen MSigDB 2026.1.Mm signature genes present in that detectable universe.
- Observed statistic: target cell-type median UCell minus the median of the other five resolved broad cell-type medians.
- UNRESOLVED was excluded from the localization comparator; all six resolved annotation groups were otherwise retained.
- Primary matching: exact strata from rank-based deciles of global mean count and detection rate (10 x 10 grid).
- Random sets were sampled without replacement within strata and exactly matched the detectable signature size.
- Primary iterations: 10000; one-sided empirical P=(1+#null>=observed)/(B+1).
- FDR: Benjamini-Hochberg across the six prespecified tests, separately for primary and sensitivity analyses.
- Gene-length sensitivity: the same expression/detection strata plus gene-span-length tertiles; missing lengths retained in an explicit LMISSING stratum.
- Gene length source: GENCODE mouse vM23 comprehensive GTF (GRCm38/mm10); gene span=end-start+1.
- Seed: 20260823. UCell maxRank=1500; average ties; rankings computed once and reused.
- Batched null scoring used sparse matrix multiplication implementing the exact positive-signature UCell u_stat formula; maximum absolute validation error versus the official API was 5.551115e-17.
- Control/HHL/NIHL labels were not used; no condition-level P values were calculated.

## Software

- R: 4.4.3
- UCell: 2.10.1
- Seurat: 5.5.1
- SeuratObject: 5.4.0
- Matrix: 1.7.6
- matrixStats: 1.5.0
- data.table: 1.18.4
- ggplot2: 4.0.3
