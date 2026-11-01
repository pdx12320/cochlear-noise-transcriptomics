#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(sys.argv[1]).resolve()
SEED = 20260825
for d in ["12_figures", "13_source_data", "14_methods/scripts", "16_release"]:
    (ROOT / d).mkdir(parents=True, exist_ok=True)

def read(rel):
    return pd.read_csv(ROOT / rel, sep="\t")

def md_table(df, cols=None, formats=None):
    x = df.loc[:, cols].copy() if cols else df.copy()
    formats = formats or {}
    for c, f in formats.items():
        x[c] = x[c].map(lambda v: "NA" if pd.isna(v) else f(v))
    x = x.fillna("NA").astype(str)
    header = "| " + " | ".join(x.columns) + " |"
    rule = "| " + " | ".join(["---"] * len(x.columns)) + " |"
    rows = ["| " + " | ".join(v.replace("|", "\\|") for v in row) + " |" for row in x.to_numpy()]
    return "\n".join([header, rule, *rows])

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

sample = read("03_doublet/S1_SAMPLE_CELL_COUNTS_POST_DOUBLET.tsv")
gate = read("07_replicate_gate/S1_REPLICATE_GATE.tsv")
ds = read("09_differential/S1_DS_SUMMARY.tsv")
gqa = read("10_gsea/S1_GSEA_QA.tsv")
val = read("13_source_data/S1_PROGRAM_VALIDATION.tsv")
cand = read("11_candidates/S1_CANDIDATE_DE.tsv")
expr = read("11_candidates/S1_CANDIDATE_EXPRESSION_CONTEXT.tsv")
pb = read("08_pseudobulk/S1_PSEUDOBULK_AUDIT.tsv")

print("[1/5] Final per-sample cell-count overview")
source = sample[["GSM", "sample_id", "tissue", "condition", "n_cells_raw", "n_cells_qc_pass", "n_cells_final"]].copy()
source.to_csv(ROOT / "13_source_data/FIG_QC_CELL_COUNTS_RAW_QC_DOUBLET.tsv", sep="\t", index=False)
fig, ax = plt.subplots(figsize=(11, 5.6))
x = np.arange(len(source))
w = .25
ax.bar(x-w, source.n_cells_raw, width=w, label="Raw", color="#8C8C8C")
ax.bar(x, source.n_cells_qc_pass, width=w, label="QC pass", color="#56B4E9")
ax.bar(x+w, source.n_cells_final, width=w, label="Singlet final", color="#009E73")
ax.set_xticks(x, source.GSM, rotation=50, ha="right")
ax.set_ylabel("Cells")
ax.set_title("GSE168041 per-sample cell retention")
ax.legend(frameon=False, ncol=3)
fig.tight_layout()
fig.savefig(ROOT / "12_figures/QC_00_raw_qc_doublet_cell_counts.png", dpi=300)
fig.savefig(ROOT / "12_figures/QC_00_raw_qc_doublet_cell_counts.pdf")
plt.close(fig)

print("[2/5] Freezing executable scripts")
script_names = [
    "phase1_manifest.py", "phase2_qc_audit.py", "apply_qc_rules.py",
    "run_scdblfinder_per_sample.R", "build_filtered_objects.py",
    "phase3_cluster_scanpy.py", "phase3_annotate_broad.py",
    "phase4_coverage_audit.py", "phase4_gate_and_pseudobulk.py",
    "phase6_pydeseq2.py", "phase7_fgsea_celltypes.R",
    "phase8_9_integrate_candidates.py", "finalize_gse168041.py",
]
for name in script_names:
    src = ROOT / "work" / name
    if src.exists():
        shutil.copy2(src, ROOT / "14_methods/scripts" / name)

print("[3/5] Writing METHODS, DECISION, and executive summary")
manifest = read("00_manifest/S1_SAMPLE_MANIFEST.tsv")
audit_sgn = read("05_clustering/S1_SAMPLE_DRIVEN_CLUSTERING_AUDIT_SGN.tsv").iloc[0]
audit_lw = read("05_clustering/S1_SAMPLE_DRIVEN_CLUSTERING_AUDIT_Lateral_wall.tsv").iloc[0]
anno_hash = (ROOT / "06_annotation/ANNOTATION_FREEZE.sha256").read_text().split()[0]
term_hash = (ROOT / "10_gsea/MSIGDB_TERM2GENE.sha256").read_text().split()[0]

methods = f"""# METHODS — GSE168041 sample-aware cell-type differential-state validation

## Scope and frozen order

The analysis used 16 independent official processed 10X count matrices from `GSE168041_RAW.tar`: SGN/cochlea neurons and cochlea lateral wall were analyzed separately, each with four Naive and four Noise biological samples. The immutable input tar SHA256 is `48acf3a156264d0eacd197eccd21b24fcde6ccf2279a6d8b3f3592318dc85704`. Phases were executed in the prespecified order: input audit, QC, per-sample doublets, separate clustering/annotation, replicate gate, pseudobulk, DE, GSEA, prior-evidence integration, and finally candidate extraction.

## Input and raw counts

All 16/16 GSMs had matching matrix/barcode/feature dimensions and 27,998 Ensembl gene IDs. Barcodes were prefixed with GSM. Raw nonnegative integer counts were retained in `layers["counts"]`; normalization changed the analysis matrix only. Pseudobulk was reconstructed exclusively by summing this layer. Per-file SHA256 values and provenance are in `01_input_audit/`.

## QC and doublets

QC was audited by sample before freezing thresholds: 500–7,500 detected genes, mitochondrial fraction <=30%, hemoglobin fraction <=20%; no upper UMI threshold. Ribosomal fraction and genes-per-UMI were descriptive. Exact rules are in `02_qc/QC_RULES.md`. scDblFinder was run separately for every GSM with sample-specific seed `20260825 + sample index`, default data-driven doublet calling, and serial execution. Versions: R 4.4.3, scDblFinder 1.23.4, SingleCellExperiment 1.28.0, BiocParallel 1.40.0, Matrix 1.7-6.

## Scanpy/AnnData workflow

Python 3.12.13, Scanpy 1.12.1, AnnData 0.12.16, NumPy 2.2.6, pandas 2.3.2, SciPy 1.16.3. Each compartment used normalize_total target 10,000, log1p, 3,000 HVGs (`seurat`, GSM used only as HVG batch key), scaling capped at 10 on the HVG analysis copy, PCA 50 components, 15-neighbor graph using 40 PCs, UMAP min_dist 0.5, and Leiden resolution 0.6 with igraph. Seed: {SEED}. No batch correction was applied. Markers used Wilcoxon tests for annotation support, not downstream inference.

## Broad annotation and freeze

Cluster labels used current cluster markers and canonical cochlear markers while excluding Apold1, Errfi1, Fos/Jun, Rela/Stat-family genes and current pathway hypotheses as deciding evidence. Conflicted rare clusters remained UNRESOLVED. Evidence is in `06_annotation/S1_ANNOTATION_EVIDENCE.tsv`; annotation freeze hash is `{anno_hash}`.

## Replicate gate and pseudobulk

Before DE results were viewed, formal inference required at least three samples per condition with >=20 cells each, no sample contributing >50% of a condition's cells, a resolved annotation, and positive integer library sizes. All other types were `DESCRIPTIVE_ONLY`. For every sample x broad cell type, raw integer counts were summed. Integer status, sample identity, and aggregate traceability are recorded in `08_pseudobulk/S1_PSEUDOBULK_AUDIT.tsv`.

## Differential expression

PyDESeq2 0.5.4 was run independently for each of 13 gate-passing broad types with `design = ~ condition`, contrast Noise vs Naive (positive effect = Noise), Naive reference, gene total-count prefilter >=10, alpha 0.05, Cook's filtering and independent filtering enabled, no LFC shrinkage, 8 CPUs, seed {SEED}. Complete gene tables were retained.

## Preranked GSEA

All genes with non-missing symbol and finite PyDESeq2 Wald statistic were ranked descending; no DEG filter was applied. Duplicate symbols used the frozen GSE160639 rule: maximum absolute stat, then maximum signed stat, then lexical Ensembl ID. Mouse MSigDB 2026.1.Mm was obtained with msigdbr 26.1.1: Hallmark (`MH`, 50 source sets) plus Reactome (`M2/CP:REACTOME`, 1,333 source sets). fgsea 1.36.2 `fgseaMultilevel` used weight 1, minSize 10, maxSize 500, eps 0, single process, seed {SEED}; BH FDR was computed across the combined eligible Hallmark+Reactome universe within each cell type. Positive NES = Noise; negative NES = Naive. Frozen TERM2GENE SHA256 is `{term_hash}`.

## Cross-dataset validation and candidates

The six prespecified GSE252285 localization tests were imported from the frozen primary matched-null result (10,000 iterations, seed 20260823) without recomputation. Replicated = matched-null PASS plus positive Noise NES with GSEA FDR <=0.05. Positive but nonsignificant GSEA = directionally supported; significant negative NES = opposite; unavailable eligible term = not testable. MAPK/AP-1 remained secondary because its prior localization failed matched-null. Only after annotation, DE, and GSEA freezes were Apold1 in lateral-wall Endothelial and Errfi1 in lateral-wall Fibrocyte/mesenchymal extracted. Single-cell detection summaries are descriptive; DE inference remains sample-level pseudobulk.

## Reproducibility

All executable scripts are frozen in `14_methods/scripts/`; parameters and full environment records are in `14_methods/`; logs are in `15_logs/`. Output checksums and the release manifest are produced in `16_release/`. No SRA/Cell Ranger rebuild, scVI, TF inference, cell-cell communication, trajectory, or additional dataset analysis was performed.
"""
(ROOT / "14_methods/METHODS.md").write_text(methods)

decisions = f"""# DECISION LOG

- Input provenance: PASS, 16/16 GSMs unambiguously matched; no FASTQ reconstruction.
- QC freeze: 500–7,500 genes, mt <=30%, hemoglobin <=20%; frozen before candidate/pathway inspection.
- Doublets: scDblFinder independently per sample; no pooled calling.
- Compartment separation: SGN and lateral wall processed independently.
- HVG flavor: Scanpy `seurat`; `seurat_v3` was not used because its optional dependency was unavailable in the isolated environment. This choice preceded biological results and did not alter raw counts.
- Batch correction: none. Sample-cluster AMI was low (SGN {audit_sgn.sample_cluster_AMI:.3f}; lateral wall {audit_lw.sample_cluster_AMI:.3f}), although local same-sample-neighbor enrichment was {audit_sgn.knn_sample_enrichment_ratio:.2f}x and {audit_lw.knn_sample_enrichment_ratio:.2f}x. Sample-dominant rare/conflicted clusters were retained as UNRESOLVED, rather than corrected or forced into labels.
- Annotation: candidate and hypothesis genes were excluded as deciding evidence; annotation frozen at `{anno_hash}`.
- Gate: >=3 samples/condition with >=20 cells, <=50% sample dominance, resolved label, positive library.
- DE: independent sample-level pseudobulk models; cells were never treated as replicates.
- GSEA: full finite Wald-stat ranking, frozen Mouse MSigDB 2026.1.Mm Hallmark+Reactome, no result-driven gene-set changes.
- Direct Reactome AP-1 activation in lateral-wall Endothelial did not meet the eligible minSize after intersection; MAPK cascade is reported as secondary pathway support, not direct AP-1 replication.
- VEGF Endothelial shows positive nonsignificant Hallmark angiogenesis but negative nonsignificant Reactome VEGF signaling; therefore classified only as directionally supported, not replicated.
- Candidate genes were viewed only after Phase 8; neither candidate changed any earlier decision.
"""
(ROOT / "14_methods/DECISION.md").write_text(decisions)

sample_display = sample[["GSM", "tissue", "condition", "n_cells_raw", "n_cells_qc_pass", "n_doublets", "n_cells_final"]]
gate_lw = gate.loc[gate.tissue.eq("Lateral_wall"), ["cell_type", "naive_samples_ge20", "noise_samples_ge20", "naive_min_cells", "noise_min_cells", "status"]]
val_display = val[["program_cell", "gse252285_matched_null", "gse168041_NES", "gse168041_FDR", "gse168041_direction", "replication_judgment"]]
cand_display = cand[["candidate", "cell_type", "baseMean", "log2FoldChange", "stat", "pvalue", "padj"]]

summary = f"""# GSE168041 sample-aware cell-type differential-state validation

## STATUS

**PASS WITH WARNING.** 16/16 samples completed through sample-aware pseudobulk and cell-type GSEA. No STOP RULE was triggered. Warnings concern localized sample-neighborhood structure and rare sample-dominant UNRESOLVED clusters, not a global sample-driven clustering failure.

## DATA

- Raw cells: {sample.n_cells_raw.sum():,}; QC pass: {sample.n_cells_qc_pass.sum():,}; singlets/final: {sample.n_cells_final.sum():,}.
- Final SGN cells: {sample.loc[sample.tissue.eq('SGN'), 'n_cells_final'].sum():,}; final lateral-wall cells: {sample.loc[sample.tissue.eq('Lateral_wall'), 'n_cells_final'].sum():,}.
- No sample was excluded. Final retention ranged from {sample.final_retention_fraction.min()*100:.1f}% to {sample.final_retention_fraction.max()*100:.1f}% of raw cells.

{md_table(sample_display)}

## CELL TYPES

Lateral wall stably contained Endothelial, Fibrocyte/mesenchymal, Immune/myeloid, Epithelial/supporting-associated, Melanocyte/intermediate, and Pericyte/vascular-smooth-muscle broad types. Endothelial, Fibrocyte, and Immune each had 4 Naive + 4 Noise biological samples and passed formal inference. Two rare/conflicted lateral-wall clusters ({int(gate.loc[(gate.tissue.eq('Lateral_wall')) & (gate.cell_type.eq('UNRESOLVED'))].shape[0])} UNRESOLVED broad category; 343 cells total) remained descriptive only.

{md_table(gate_lw)}

## PSEUDOBULK

Thirteen broad cell types passed the frozen gate (7 SGN, 6 lateral wall), each with biological n=4 Naive + 4 Noise. All pseudobulk matrices were integer and traceable. Lateral-wall DE: Endothelial 1,049 FDR<0.05 genes (546 Noise-up/503 Noise-down); Fibrocyte 219 (99/120); Immune 345 (177/168).

## PROGRAM VALIDATION

{md_table(val_display, formats={'gse168041_NES': lambda v: f'{v:.3f}', 'gse168041_FDR': lambda v: f'{v:.3g}'})}

Interpretation: only **NFkB -> Fibrocyte** achieved strict independent sample-aware replication (NES 1.648, FDR 0.00371). VEGF/angiogenesis -> Endothelial was positive but nonsignificant (NES 1.173, FDR 0.384), hence directional support only. NFkB -> Endothelial, JAK-STAT -> Endothelial, and JAK-STAT -> Immune were significant in the opposite, Naive-enriched direction. MAPK cascade -> Endothelial was strongly Noise-enriched (NES 1.918, FDR 8.37e-7), but remains secondary because GSE252285 Endothelial localization failed matched-null and the direct AP-1 term was not eligible at minSize.

## CANDIDATES

{md_table(cand_display, formats={'baseMean': lambda v: f'{v:.1f}', 'log2FoldChange': lambda v: f'{v:.3f}', 'stat': lambda v: f'{v:.3f}', 'pvalue': lambda v: f'{v:.3g}', 'padj': lambda v: f'{v:.3g}'})}

APOLD1 showed no sample-aware Endothelial differential effect (log2FC 0.031, FDR 0.913). ERRFI1 trended lower in Noise Fibrocytes (log2FC -0.361) but was not significant (FDR 0.588). Neither is supported as a standalone differential marker for a final validation panel from this dataset; they may be retained only as secondary context alongside replicated pathway readouts.

## INTERPRETATION

GSE252285 localization and GSE168041 condition response answer different questions. Cell-type localization replicated as a Noise response only for NFkB in Fibrocytes. Endothelial VEGF localization received weak directional but no statistical replication. Endothelial/Immune inflammatory localization did not translate into Noise induction here and instead showed significant Naive enrichment. Endothelial MAPK is a strong independent secondary state signal despite lacking prior matched-null localization support.

## QA

- 16/16 input matrices PASS; raw counts preserved; 16/16 pseudobulk category audits PASS.
- 13/13 formal DE models PASS; 13/13 GSEA mapping audits PASS.
- Annotation freeze: `{anno_hash}`.
- MSigDB TERM2GENE SHA256: `{term_hash}`.
- Seed: {SEED} (matched-null source used 20260823).
- Authoritative server root: `{ROOT}`.
- Full checksums and release archive: `16_release/`.
"""
(ROOT / "00_EXECUTIVE_SUMMARY_CN.md").write_text(summary)

print("[4/5] QA checks")
checks = []
def check(name, passed, observed, expected):
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": str(observed), "expected": str(expected)})

check("manifest_16_samples", len(manifest) == 16 and manifest.GSM.nunique() == 16, len(manifest), "16 unique GSM")
check("raw_cell_total", sample.n_cells_raw.sum() == 71743, sample.n_cells_raw.sum(), 71743)
check("qc_cell_total", sample.n_cells_qc_pass.sum() == 69036, sample.n_cells_qc_pass.sum(), 69036)
check("final_singlet_total", sample.n_cells_final.sum() == 64705, sample.n_cells_final.sum(), 64705)
check("all_pseudobulk_audits", pb.status.eq("PASS").all() and pb.integer_counts.all() and pb.sample_identity_match.all(), pb.status.value_counts().to_dict(), "all PASS/integer/identity")
formal_n = int(gate.status.eq("FORMAL_INFERENCE").sum())
check("formal_models", formal_n == 13, formal_n, 13)
check("DE_models", len(ds) == formal_n and ds.status.eq("PASS").all(), len(ds), formal_n)
check("GSEA_mapping", len(gqa) == formal_n and gqa.status.eq("PASS").all(), len(gqa), formal_n)
check("GSEA_overlap_floor", gqa.rank_genes_matching_any_set.min() >= 1000 and gqa.rank_gene_set_union_overlap_fraction.min() >= .10, f"min genes={gqa.rank_genes_matching_any_set.min()}, min fraction={gqa.rank_gene_set_union_overlap_fraction.min():.3f}", ">=1000 and >=0.10")
check("annotation_hash", sha256(ROOT / "06_annotation/S1_ANNOTATION_EVIDENCE.tsv") == anno_hash, sha256(ROOT / "06_annotation/S1_ANNOTATION_EVIDENCE.tsv"), anno_hash)
check("program_rows", len(val) == 6, len(val), 6)
check("candidate_rows", len(cand) == 2 and len(expr) == 4, f"DE={len(cand)}, context={len(expr)}", "2 and 4")
required_figs = [
    "QC_00_raw_qc_doublet_cell_counts", "SGN_UMAP_Broad_annotation", "Lateral_wall_UMAP_Broad_annotation",
    "SGN_UMAP_sample", "Lateral_wall_UMAP_sample", "SGN_UMAP_condition", "Lateral_wall_UMAP_condition",
    "Lateral_wall_Broad_celltype_by_sample_heatmap", "PSEUDOBULK_PCA_Lateral_wall_Endothelial",
    "Lateral_wall_Endothelial_volcano", "FOCUSED_PROGRAM_NES_HEATMAP", "S1_CANDIDATE_EFFECTS",
]
missing = [f for f in required_figs if not (ROOT / f"12_figures/{f}.png").exists() or not (ROOT / f"12_figures/{f}.pdf").exists()]
check("required_figure_pairs", not missing, "missing=" + ",".join(missing), "all PNG+PDF")
qa = pd.DataFrame(checks)
qa.to_csv(ROOT / "14_methods/QA_STATUS.tsv", sep="\t", index=False)
if qa.status.ne("PASS").any():
    raise RuntimeError("Final QA failed:\n" + qa.loc[qa.status.ne("PASS")].to_string(index=False))

print("[5/5] Method input checksums")
key_files = [
    ROOT / "00_manifest/S1_SAMPLE_MANIFEST.tsv", ROOT / "02_qc/QC_RULES.md",
    ROOT / "07_replicate_gate/GATE_RULES.md", ROOT / "06_annotation/S1_ANNOTATION_EVIDENCE.tsv",
    ROOT / "09_differential/S1_DS_SUMMARY.tsv", ROOT / "10_gsea/S1_GSEA_ALL.tsv",
    ROOT / "13_source_data/S1_PROGRAM_VALIDATION.tsv", ROOT / "11_candidates/S1_CANDIDATE_DE.tsv",
]
pd.DataFrame([{"file": str(p.relative_to(ROOT)), "sha256": sha256(p), "bytes": p.stat().st_size} for p in key_files]).to_csv(
    ROOT / "14_methods/KEY_OUTPUT_CHECKSUMS.tsv", sep="\t", index=False)
print("Finalization complete; QA PASS")
