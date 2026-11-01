# DECISION LOG

- Input provenance: PASS, 16/16 GSMs unambiguously matched; no FASTQ reconstruction.
- QC freeze: 500–7,500 genes, mt <=30%, hemoglobin <=20%; frozen before candidate/pathway inspection.
- Doublets: scDblFinder independently per sample; no pooled calling.
- Compartment separation: SGN and lateral wall processed independently.
- HVG flavor: Scanpy `seurat`; `seurat_v3` was not used because its optional dependency was unavailable in the isolated environment. This choice preceded biological results and did not alter raw counts.
- Batch correction: none. Sample-cluster AMI was low (SGN 0.104; lateral wall 0.128), although local same-sample-neighbor enrichment was 3.24x and 3.38x. Sample-dominant rare/conflicted clusters were retained as UNRESOLVED, rather than corrected or forced into labels.
- Annotation: candidate and hypothesis genes were excluded as deciding evidence; annotation frozen at `28ad95e6596e164dec4807d294ea0ff1e52618881ef53c8dbf26723bac5d6864`.
- Gate: >=3 samples/condition with >=20 cells, <=50% sample dominance, resolved label, positive library.
- DE: independent sample-level pseudobulk models; cells were never treated as replicates.
- GSEA: full finite Wald-stat ranking, frozen Mouse MSigDB 2026.1.Mm Hallmark+Reactome, no result-driven gene-set changes.
- Direct Reactome AP-1 activation in lateral-wall Endothelial did not meet the eligible minSize after intersection; MAPK cascade is reported as secondary pathway support, not direct AP-1 replication.
- VEGF Endothelial shows positive nonsignificant Hallmark angiogenesis but negative nonsignificant Reactome VEGF signaling; therefore classified only as directionally supported, not replicated.
- Candidate genes were viewed only after Phase 8; neither candidate changed any earlier decision.
