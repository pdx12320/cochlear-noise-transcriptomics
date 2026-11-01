# METHODS — GSE168041 sample-aware cell-type differential-state validation

## Scope and frozen order

The analysis used 16 independent official processed 10X count matrices from `GSE168041_RAW.tar`: SGN/cochlea neurons and cochlea lateral wall were analyzed separately, each with four Naive and four Noise biological samples. The immutable input tar SHA256 is `48acf3a156264d0eacd197eccd21b24fcde6ccf2279a6d8b3f3592318dc85704`. Phases were executed in the prespecified order: input audit, QC, per-sample doublets, separate clustering/annotation, replicate gate, pseudobulk, DE, GSEA, prior-evidence integration, and finally candidate extraction.

## Input and raw counts

All 16/16 GSMs had matching matrix/barcode/feature dimensions and 27,998 Ensembl gene IDs. Barcodes were prefixed with GSM. Raw nonnegative integer counts were retained in `layers["counts"]`; normalization changed the analysis matrix only. Pseudobulk was reconstructed exclusively by summing this layer. Per-file SHA256 values and provenance are in `01_input_audit/`.

## QC and doublets

QC was audited by sample before freezing thresholds: 500–7,500 detected genes, mitochondrial fraction <=30%, hemoglobin fraction <=20%; no upper UMI threshold. Ribosomal fraction and genes-per-UMI were descriptive. Exact rules are in `02_qc/QC_RULES.md`. scDblFinder was run separately for every GSM with sample-specific seed `20260825 + sample index`, default data-driven doublet calling, and serial execution. Versions: R 4.4.3, scDblFinder 1.23.4, SingleCellExperiment 1.28.0, BiocParallel 1.40.0, Matrix 1.7-6.

## Scanpy/AnnData workflow

Python 3.12.13, Scanpy 1.12.1, AnnData 0.12.16, NumPy 2.2.6, pandas 2.3.2, SciPy 1.16.3. Each compartment used normalize_total target 10,000, log1p, 3,000 HVGs (`seurat`, GSM used only as HVG batch key), scaling capped at 10 on the HVG analysis copy, PCA 50 components, 15-neighbor graph using 40 PCs, UMAP min_dist 0.5, and Leiden resolution 0.6 with igraph. Seed: 20260825. No batch correction was applied. Markers used Wilcoxon tests for annotation support, not downstream inference.

## Broad annotation and freeze

Cluster labels used current cluster markers and canonical cochlear markers while excluding Apold1, Errfi1, Fos/Jun, Rela/Stat-family genes and current pathway hypotheses as deciding evidence. Conflicted rare clusters remained UNRESOLVED. Evidence is in `06_annotation/S1_ANNOTATION_EVIDENCE.tsv`; annotation freeze hash is `28ad95e6596e164dec4807d294ea0ff1e52618881ef53c8dbf26723bac5d6864`.

## Replicate gate and pseudobulk

Before DE results were viewed, formal inference required at least three samples per condition with >=20 cells each, no sample contributing >50% of a condition's cells, a resolved annotation, and positive integer library sizes. All other types were `DESCRIPTIVE_ONLY`. For every sample x broad cell type, raw integer counts were summed. Integer status, sample identity, and aggregate traceability are recorded in `08_pseudobulk/S1_PSEUDOBULK_AUDIT.tsv`.

## Differential expression

PyDESeq2 0.5.4 was run independently for each of 13 gate-passing broad types with `design = ~ condition`, contrast Noise vs Naive (positive effect = Noise), Naive reference, gene total-count prefilter >=10, alpha 0.05, Cook's filtering and independent filtering enabled, no LFC shrinkage, 8 CPUs, seed 20260825. Complete gene tables were retained.

## Preranked GSEA

All genes with non-missing symbol and finite PyDESeq2 Wald statistic were ranked descending; no DEG filter was applied. Duplicate symbols used the frozen GSE160639 rule: maximum absolute stat, then maximum signed stat, then lexical Ensembl ID. Mouse MSigDB 2026.1.Mm was obtained with msigdbr 26.1.1: Hallmark (`MH`, 50 source sets) plus Reactome (`M2/CP:REACTOME`, 1,333 source sets). fgsea 1.36.2 `fgseaMultilevel` used weight 1, minSize 10, maxSize 500, eps 0, single process, seed 20260825; BH FDR was computed across the combined eligible Hallmark+Reactome universe within each cell type. Positive NES = Noise; negative NES = Naive. Frozen TERM2GENE SHA256 is `98dcdf8081f93b433396b8cbc1eaa6ad0e34ce2afc84ae82fbda99df2e641caf`.

## Cross-dataset validation and candidates

The six prespecified GSE252285 localization tests were imported from the frozen primary matched-null result (10,000 iterations, seed 20260823) without recomputation. Replicated = matched-null PASS plus positive Noise NES with GSEA FDR <=0.05. Positive but nonsignificant GSEA = directionally supported; significant negative NES = opposite; unavailable eligible term = not testable. MAPK/AP-1 remained secondary because its prior localization failed matched-null. Only after annotation, DE, and GSEA freezes were Apold1 in lateral-wall Endothelial and Errfi1 in lateral-wall Fibrocyte/mesenchymal extracted. Single-cell detection summaries are descriptive; DE inference remains sample-level pseudobulk.

## Reproducibility

All executable scripts are frozen in `14_methods/scripts/`; parameters and full environment records are in `14_methods/`; logs are in `15_logs/`. Output checksums and the release manifest are produced in `16_release/`. No SRA/Cell Ranger rebuild, scVI, TF inference, cell-cell communication, trajectory, or additional dataset analysis was performed.
