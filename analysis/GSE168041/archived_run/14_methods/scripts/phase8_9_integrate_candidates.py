#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(sys.argv[1]).resolve()
NULL_ROOT = Path("/path/to/cochlear_analysis/GSE252285_UCELL_MATCHED_NULL_20260823/outputs")
SEED = 20260825
np.random.seed(SEED)

for d in ["11_candidates", "12_figures", "13_source_data", "14_methods"]:
    (ROOT / d).mkdir(parents=True, exist_ok=True)

print("[1/4] Integrating frozen GSE252285 matched-null evidence")
null_file = NULL_ROOT / "MATCHED_NULL_RESULTS.tsv"
null_all = pd.read_csv(null_file, sep="\t")
null = null_all.loc[null_all["match_mode"].eq("primary")].copy()
null.to_csv(ROOT / "13_source_data/GSE252285_MATCHED_NULL_PRIMARY.tsv", sep="\t", index=False)
shutil.copy2(NULL_ROOT / "RESULT_SUMMARY.md", ROOT / "13_source_data/GSE252285_MATCHED_NULL_RESULT_SUMMARY.md")

focus = pd.read_csv(ROOT / "10_gsea/S1_GSEA_FOCUSED_PROGRAMS.tsv", sep="\t")
focus = focus.loc[focus["tissue"].eq("Lateral_wall")].copy()

specs = [
    ("NFkB -> Endothelial", "TNFA_NFKB__Endothelial_vascular", "Endothelial", "Hallmark::HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PRIMARY"),
    ("VEGF -> Endothelial", "VEGF_ANGIOGENESIS__Endothelial_vascular", "Endothelial", "Hallmark::HALLMARK_ANGIOGENESIS", "PRIMARY"),
    ("JAK-STAT -> Endothelial", "IL6_JAK_STAT3__Endothelial_vascular", "Endothelial", "Hallmark::HALLMARK_IL6_JAK_STAT3_SIGNALING", "PRIMARY"),
    ("NFkB -> Fibrocyte", "TNFA_NFKB__Fibrocyte_mesenchymal", "Fibrocyte / mesenchymal", "Hallmark::HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PRIMARY"),
    ("JAK-STAT -> Immune", "IL6_JAK_STAT3__Immune_myeloid", "Immune / myeloid", "Hallmark::HALLMARK_IL6_JAK_STAT3_SIGNALING", "PRIMARY"),
    ("MAPK/AP-1 -> Endothelial", "MAPK_AP1__Endothelial_vascular", "Endothelial", "Reactome::REACTOME_MAPK_FAMILY_SIGNALING_CASCADES", "SECONDARY"),
]

rows = []
for label, test_id, cell_type, pathway_id, role in specs:
    n = null.loc[null["test_id"].eq(test_id)]
    g = focus.loc[focus["cell_type"].eq(cell_type) & focus["pathway_id"].eq(pathway_id)]
    if len(n) != 1:
        raise RuntimeError(f"Expected one matched-null row for {test_id}, found {len(n)}")
    if len(g) > 1:
        raise RuntimeError(f"Multiple GSEA rows for {cell_type} / {pathway_id}")
    nr = n.iloc[0]
    matched_pass = bool(nr["significant_FDR05"])
    if len(g) == 0:
        nes = pval = fdr = np.nan
        direction = "NOT_TESTABLE"
        judgment = "not testable"
    else:
        gr = g.iloc[0]
        nes, pval, fdr = float(gr.NES), float(gr.P), float(gr.FDR)
        direction = "Noise" if nes > 0 else "Naive"
        if role == "SECONDARY" and not matched_pass:
            judgment = "secondary pathway support; localization not replicated" if (nes > 0 and fdr <= 0.05) else "unsupported"
        elif matched_pass and nes > 0 and fdr <= 0.05:
            judgment = "replicated"
        elif matched_pass and nes > 0:
            judgment = "directionally supported"
        elif matched_pass and nes < 0 and fdr <= 0.05:
            judgment = "opposite"
        else:
            judgment = "unsupported"
    rows.append({
        "program_cell": label,
        "hypothesis_role": role,
        "gse252285_empirical_P": nr["empirical_P"],
        "gse252285_empirical_FDR": nr["empirical_FDR"],
        "gse252285_matched_null": "PASS" if matched_pass else "FAIL",
        "gse168041_pathway_id": pathway_id,
        "gse168041_NES": nes,
        "gse168041_P": pval,
        "gse168041_FDR": fdr,
        "gse168041_direction": direction,
        "replication_judgment": judgment,
    })

validation = pd.DataFrame(rows)
validation.to_csv(ROOT / "13_source_data/S1_PROGRAM_VALIDATION.tsv", sep="\t", index=False)

# Preserve the prespecified corroborative VEGF and direct AP-1 results, including absence under minSize.
corroborative_specs = [
    ("Endothelial", "Reactome::REACTOME_SIGNALING_BY_VEGF"),
    ("Endothelial", "Reactome::REACTOME_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS"),
]
cor = []
for ct, pid in corroborative_specs:
    hit = focus.loc[focus["cell_type"].eq(ct) & focus["pathway_id"].eq(pid)]
    if len(hit):
        z = hit.iloc[0]
        cor.append({"cell_type": ct, "pathway_id": pid, "status": "TESTED", "NES": z.NES, "P": z.P, "FDR": z.FDR, "size": z["size"]})
    else:
        cor.append({"cell_type": ct, "pathway_id": pid, "status": "NOT_TESTED_MIN_SIZE_OR_NO_OVERLAP", "NES": np.nan, "P": np.nan, "FDR": np.nan, "size": np.nan})
pd.DataFrame(cor).to_csv(ROOT / "13_source_data/S1_PROGRAM_CORROBORATIVE.tsv", sep="\t", index=False)

print("[2/4] Extracting candidates only after annotation, DE, and GSEA freeze")
candidate_specs = [
    ("APOLD1", "Apold1", "Endothelial", ROOT / "09_differential/Lateral_wall/S1_DS_ALL_Endothelial.tsv"),
    ("ERRFI1", "Errfi1", "Fibrocyte / mesenchymal", ROOT / "09_differential/Lateral_wall/S1_DS_ALL_Fibrocyte_mesenchymal.tsv"),
]
de_rows = []
for display, symbol, ct, f in candidate_specs:
    de = pd.read_csv(f, sep="\t")
    hit = de.loc[de["gene_symbol"].eq(symbol)].copy()
    if len(hit) != 1:
        raise RuntimeError(f"Expected exactly one DE row for {symbol} in {ct}; found {len(hit)}")
    h = hit.iloc[0]
    de_rows.append({
        "candidate": display, "gene_symbol": symbol, "tissue": "Lateral_wall", "cell_type": ct,
        "gene_id": h.gene_id, "baseMean": h.baseMean, "log2FoldChange": h.log2FoldChange,
        "lfcSE": h.lfcSE, "stat": h.stat, "pvalue": h.pvalue, "padj": h.padj,
        "contrast": "Noise_vs_Naive", "positive_direction": "Noise",
    })
cand_de = pd.DataFrame(de_rows)
cand_de.to_csv(ROOT / "11_candidates/S1_CANDIDATE_DE.tsv", sep="\t", index=False)

adata = ad.read_h5ad(ROOT / "06_annotation/Lateral_wall_annotated.h5ad")
if "counts" not in adata.layers:
    raise RuntimeError("Raw counts layer missing")
counts = adata.layers["counts"]
expr_rows = []
sample_rows = []
for display, symbol, ct, _ in candidate_specs:
    idx = np.flatnonzero(adata.var["gene_symbol"].astype(str).to_numpy() == symbol)
    if len(idx) != 1:
        raise RuntimeError(f"Expected exactly one var row for {symbol}; found {len(idx)}")
    j = int(idx[0])
    for condition in ["Naive", "Noise"]:
        mask = ((adata.obs["broad_annotation"].astype(str) == ct) &
                (adata.obs["condition"].astype(str) == condition)).to_numpy()
        v = counts[mask, j]
        v = v.toarray().ravel() if sparse.issparse(v) else np.asarray(v).ravel()
        expr_rows.append({
            "candidate": display, "gene_symbol": symbol, "cell_type": ct, "condition": condition,
            "n_cells": int(mask.sum()), "detected_cells": int((v > 0).sum()),
            "detection_fraction": float((v > 0).mean()), "mean_raw_counts_per_cell": float(v.mean()),
            "median_raw_counts_per_cell": float(np.median(v)),
        })
    for gsm in sorted(adata.obs["GSM"].astype(str).unique()):
        mask = ((adata.obs["broad_annotation"].astype(str) == ct) &
                (adata.obs["GSM"].astype(str) == gsm)).to_numpy()
        v = counts[mask, j]
        v = v.toarray().ravel() if sparse.issparse(v) else np.asarray(v).ravel()
        sample_rows.append({
            "candidate": display, "gene_symbol": symbol, "cell_type": ct, "GSM": gsm,
            "condition": str(adata.obs.loc[mask, "condition"].iloc[0]) if mask.sum() else "NA",
            "n_cells": int(mask.sum()), "detected_cells": int((v > 0).sum()),
            "detection_fraction": float((v > 0).mean()) if len(v) else np.nan,
            "mean_raw_counts_per_cell": float(v.mean()) if len(v) else np.nan,
            "sum_raw_counts": int(v.sum()),
        })

expr = pd.DataFrame(expr_rows)
sample_expr = pd.DataFrame(sample_rows)
expr.to_csv(ROOT / "11_candidates/S1_CANDIDATE_EXPRESSION_CONTEXT.tsv", sep="\t", index=False)
sample_expr.to_csv(ROOT / "11_candidates/S1_CANDIDATE_SAMPLE_CONTEXT.tsv", sep="\t", index=False)

print("[3/4] Candidate effect figure")
plot = cand_de.copy()
plot["label"] = plot["candidate"] + "\n" + plot["cell_type"]
fig, ax = plt.subplots(figsize=(6.4, 4.1))
x = np.arange(len(plot))
colors = np.where(plot["log2FoldChange"] >= 0, "#B2182B", "#2166AC")
for i, r in plot.iterrows():
    ax.errorbar(i, r["log2FoldChange"], yerr=1.96 * r["lfcSE"], fmt="none",
                ecolor=colors[i], elinewidth=1.8, capsize=4)
    ax.scatter(i, r["log2FoldChange"], c=colors[i], s=75, edgecolor="black", linewidth=.5, zorder=3)
ax.axhline(0, color="0.35", lw=.8)
ax.set_xticks(x, plot["label"])
ax.set_ylabel("Pseudobulk log2 fold change (Noise vs Naive)")
ax.set_title("Prespecified candidate effects")
for i, r in plot.iterrows():
    q = "NA" if pd.isna(r.padj) else f"{r.padj:.2g}"
    ax.text(i, r.log2FoldChange + (0.10 if r.log2FoldChange >= 0 else -0.10), f"FDR={q}", ha="center", va="bottom" if r.log2FoldChange >= 0 else "top", fontsize=8)
fig.tight_layout()
fig.savefig(ROOT / "12_figures/S1_CANDIDATE_EFFECTS.png", dpi=300)
fig.savefig(ROOT / "12_figures/S1_CANDIDATE_EFFECTS.pdf")
plt.close(fig)
plot.to_csv(ROOT / "13_source_data/S1_CANDIDATE_EFFECTS_PLOT_SOURCE.tsv", sep="\t", index=False)

print("[4/4] Provenance and checksums")
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

prov = pd.DataFrame([
    {"source": str(null_file), "sha256": sha256(null_file), "role": "frozen GSE252285 primary matched-null results"},
    {"source": str(ROOT / "10_gsea/S1_GSEA_FOCUSED_PROGRAMS.tsv"), "sha256": sha256(ROOT / "10_gsea/S1_GSEA_FOCUSED_PROGRAMS.tsv"), "role": "GSE168041 frozen focused GSEA"},
    {"source": str(ROOT / "06_annotation/ANNOTATION_FREEZE.sha256"), "sha256": sha256(ROOT / "06_annotation/ANNOTATION_FREEZE.sha256"), "role": "annotation freeze record"},
])
prov.to_csv(ROOT / "13_source_data/S1_INTEGRATION_PROVENANCE.tsv", sep="\t", index=False)
print("Phase 8-9 complete")
