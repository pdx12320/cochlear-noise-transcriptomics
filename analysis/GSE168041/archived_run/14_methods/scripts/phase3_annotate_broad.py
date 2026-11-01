#!/usr/bin/env python3
import argparse
import hashlib
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


MAPPINGS = {
    "SGN": {
        "0": "Epithelial / supporting-associated", "1": "Epithelial / supporting-associated",
        "2": "Glial / Schwann-like", "3": "Neuronal / SGN", "4": "Neuronal / SGN",
        "5": "Neuronal / SGN", "6": "Endothelial", "7": "UNRESOLVED",
        "8": "Epithelial / supporting-associated", "9": "Fibrocyte / mesenchymal",
        "10": "Fibrocyte / mesenchymal", "11": "Fibrocyte / mesenchymal",
        "12": "Fibrocyte / mesenchymal", "13": "Epithelial / supporting-associated",
        "14": "Pericyte / vascular smooth muscle", "15": "Immune / myeloid",
        "16": "Epithelial / supporting-associated", "17": "Fibrocyte / mesenchymal",
        "18": "Fibrocyte / mesenchymal", "19": "Epithelial / supporting-associated",
        "20": "Glial / Schwann-like", "21": "Epithelial / supporting-associated",
        "22": "Epithelial / supporting-associated", "23": "Immune / myeloid",
        "24": "Glial / Schwann-like", "25": "Epithelial / supporting-associated",
        "26": "Melanocyte / intermediate",
    },
    "Lateral_wall": {
        "0": "Melanocyte / intermediate", "1": "Endothelial", "2": "Melanocyte / intermediate",
        "3": "Melanocyte / intermediate", "4": "Epithelial / supporting-associated",
        "5": "Epithelial / supporting-associated", "6": "Epithelial / supporting-associated",
        "7": "Epithelial / supporting-associated", "8": "Pericyte / vascular smooth muscle",
        "9": "Fibrocyte / mesenchymal", "10": "Fibrocyte / mesenchymal",
        "11": "Fibrocyte / mesenchymal", "12": "Fibrocyte / mesenchymal",
        "13": "Epithelial / supporting-associated", "14": "Epithelial / supporting-associated",
        "15": "Immune / myeloid", "16": "Melanocyte / intermediate",
        "17": "Fibrocyte / mesenchymal", "18": "Immune / myeloid",
        "19": "Epithelial / supporting-associated", "20": "UNRESOLVED", "21": "UNRESOLVED",
    },
}

EVIDENCE = {
    "SGN": {
        "0": ("Vmo1;Rnase1;Lrp2;Gpx3", "mixed epithelial identity", "Moderate"),
        "1": ("Epcam;Cldn9;Cldn3;Tspan8", "none", "High"),
        "2": ("Pmp22;Mbp;Mpz;Cnp;Plp1", "none", "High"),
        "3": ("Nefh;Nefl;Nefm;Elavl4;Map2", "none", "High"),
        "4": ("Uchl1;Gng3;Stmn3;Eno2;Cplx1", "condition-skewed neuronal state", "Moderate"),
        "5": ("Snap25;Pvalb;Rab3a;Nefl", "condition-skewed neuronal state", "High"),
        "6": ("Flt1;Cldn5;Egfl7;Klf2;Esam;Sox17", "none", "High"),
        "7": ("Snap25;Ank2;Kcnc3", "S100a8/S100a9 and hemoglobin; 97% GSM5124294", "UNRESOLVED"),
        "8": ("Gjb2;Gjb6;Sox9;Rbp1", "glial-like support markers", "Moderate"),
        "9": ("Dcn;Lum;Coch;Cxcl12", "none", "High"),
        "10": ("Igfbp4;Efemp1;Col8a1;Aebp1", "none", "High"),
        "11": ("Mgp;Rbp1;Efemp1;Ptn", "Cldn11/Gjb6 support-like markers", "Moderate"),
        "12": ("Emilin2;Coch;Otor;Igfbp2", "Emcn without core endothelial panel", "Moderate"),
        "13": ("Tppp3;Wfdc2;Ccdc153", "rare ciliated epithelial population", "Moderate"),
        "14": ("Acta2;Tagln;Myh11;Mylk", "none", "High"),
        "15": ("C1qa;C1qb;C1qc;Lyz2;Csf1r", "none", "High"),
        "16": ("Krt8;Krt18;Epcam;Cldn4", "stress-response markers", "High"),
        "17": ("Col1a1;Col1a2;Lum;Bglap;Pth1r", "osteogenic mesenchymal subtype", "High"),
        "18": ("Dcn;Cldn11;Rbp1;Igfbp4", "support-like markers", "Moderate"),
        "19": ("Otof;Espn;Myo6;Bdnf", "sensory epithelial specialization", "High"),
        "20": ("Plp1;Fabp7;S100a6;Ndrg2", "none", "High"),
        "21": ("Krt19;Epcam;Cldn4;Wfdc2", "none", "High"),
        "22": ("Epcam;Cldn9;Tspan8;S100a1", "none", "High"),
        "23": ("Cd52;Tyrobp;Fcer1g;Coro1a;Cd53", "mixed immune lineage", "High"),
        "24": ("Plp1;Mbp;Mog;Mag;Cnp", "rare myelinating glia", "High"),
        "25": ("Kcne1;Slc12a2;Lrp2;Atp1b2", "rare marginal-like epithelium", "High"),
        "26": ("Dct;Slc45a2;Gpnmb;Cited1", "none", "High"),
    },
    "Lateral_wall": {
        "0": ("Dct;Mlana;Tyr;Slc45a2;Gsta1", "strongly Naive-skewed but multi-sample", "High"),
        "1": ("Flt1;Egfl7;Pecam1;Cd34;Lmo2", "none", "High"),
        "2": ("Kcnj10;Dct;Tyr;Ednrb;Gpnmb", "strongly Noise-skewed but multi-sample", "High"),
        "3": ("Kcnj10;Ednrb;Plp1;Gpnmb", "Noise-only state; melanocyte/glial overlap", "Moderate"),
        "4": ("Kcne1;Slc12a2;Lrp2;Atp1b2", "none", "High"),
        "5": ("Krt8;Krt18;Epcam;Cldn4;Cldn7", "none", "High"),
        "6": ("Fbxo2;S100a1;Kcnj16;Gata3", "supporting epithelial subtype", "Moderate"),
        "7": ("Kcne1;Dnase1;Igf1;Slc12a2", "elevated mitochondrial transcripts", "High"),
        "8": ("Rgs5;Ndufa4l2;Myl9;Cox4i2;Higd1b", "none", "High"),
        "9": ("Dcn;Coch;Car3;Slc4a11;Tnfrsf11b", "none", "High"),
        "10": ("Coch;Igfbp4;Bgn;Efemp1;Ogn", "none", "High"),
        "11": ("Col1a1;Col1a2;Bglap;Pth1r;Dcn", "rare osteogenic mesenchymal subtype", "High"),
        "12": ("Dcn;Ptn;Mdk;Igfbp4;Efemp1", "none", "High"),
        "13": ("Cldn11;Gjb2;Gjb6;Clu", "basal/supporting epithelium", "High"),
        "14": ("Epcam;Slc26a4;Kcnj16;Fbxo2", "none", "High"),
        "15": ("C1qa;C1qb;C1qc;Lyz2;Cd68;Trem2", "none", "High"),
        "16": ("Mlana;Gsta1;Gsta4;S100b", "weaker lineage marker coherence", "Moderate"),
        "17": ("Emilin2;Coch;Otor;Igfbp2;Efemp1", "Emcn without core endothelial panel", "Moderate"),
        "18": ("Coro1a;Cd52;Rac2;Cd53;Gmfg", "rare mixed immune lineage", "High"),
        "19": ("Vmo1;Rnase1;Lrp2;Wfdc2", "mixed epithelial subtype", "Moderate"),
        "20": ("Wfdc2;Bpifa1;Krt8;Krt18;Krt19", "98% GSM5124303", "UNRESOLVED"),
        "21": ("Foxj1;Tppp3;Ccdc153;Rsph1", "49 cells; 100% GSM5124303", "UNRESOLVED"),
    },
}

MARKER_PANEL = {
    "Endothelial": ["Pecam1", "Cdh5", "Flt1", "Kdr", "Esam", "Emcn", "Egfl7", "Vwf", "Klf2"],
    "Fibrocyte_mesenchymal": ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Coch", "Pdgfra", "Igfbp4"],
    "Pericyte_vascular": ["Rgs5", "Pdgfrb", "Cspg4", "Acta2", "Tagln", "Myh11", "Cox4i2"],
    "Immune_myeloid": ["Ptprc", "Lyz2", "C1qa", "C1qb", "C1qc", "Tyrobp", "Aif1", "Csf1r", "Cd68"],
    "Glial_Schwann": ["Sox10", "S100b", "Plp1", "Mpz", "Mbp", "Cnp", "Fabp7"],
    "Neuronal_SGN": ["Snap25", "Syt1", "Tubb3", "Nefh", "Nefm", "Elavl4", "Map2"],
    "Epithelial_supporting": ["Epcam", "Krt8", "Krt18", "Krt19", "Gjb2", "Gjb6", "Kcne1", "Slc12a2"],
    "Melanocyte_intermediate": ["Dct", "Tyr", "Mlana", "Pmel", "Mitf", "Slc45a2", "Kcnj10", "Ednrb"],
}


def save_pair(fig, stem):
    fig.savefig(stem + ".png", dpi=240, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    plt.close(fig)


def plot_umap(adata, tissue, stem):
    values = adata.obs["broad_annotation"].astype(str)
    levels = sorted(values.unique())
    cmap = plt.get_cmap("tab10")
    colors = {level: cmap(i % 10) for i, level in enumerate(levels)}
    coords = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for level in levels:
        mask = values.to_numpy() == level
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2.2, alpha=0.65, linewidths=0,
                   rasterized=True, color=colors[level], label=level)
    ax.set(title=f"{tissue}: frozen Broad annotation", xlabel="UMAP1", ylabel="UMAP2")
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, markerscale=3, fontsize=8)
    fig.tight_layout(); save_pair(fig, stem)


def marker_panel_table(adata, tissue):
    rows = []
    symbol_to_indices = {}
    symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    for family, genes in MARKER_PANEL.items():
        for gene in genes:
            idx = np.flatnonzero(symbols == gene)
            symbol_to_indices[gene] = idx
            for cluster in sorted(adata.obs["leiden"].astype(str).unique(), key=int):
                mask = adata.obs["leiden"].astype(str).to_numpy() == cluster
                if len(idx) == 0:
                    mean = np.nan; pct = np.nan
                else:
                    block = adata.X[mask, :][:, idx]
                    mean = float(block.mean())
                    pct = float((block > 0).sum() / (block.shape[0] * block.shape[1]))
                rows.append({"tissue": tissue, "cluster": cluster, "marker_family": family,
                             "gene_symbol": gene, "n_gene_ids": len(idx), "mean_log_normalized": mean,
                             "detection_fraction": pct})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    evidence_rows = []
    panel_parts = []
    for tissue in ["SGN", "Lateral_wall"]:
        comp = pd.read_csv(os.path.join(root, "05_clustering", f"S1_CLUSTER_COMPOSITION_AUDIT_{tissue}.tsv"), sep="\t")
        comp["cluster"] = comp["cluster"].astype(str)
        for cluster, label in MAPPINGS[tissue].items():
            positive, conflicting, confidence = EVIDENCE[tissue][cluster]
            row = comp.loc[comp.cluster == cluster].iloc[0]
            evidence_rows.append({
                "tissue": tissue, "cluster": cluster, "positive_markers": positive,
                "conflicting_markers_or_evidence": conflicting, "final_broad_annotation": label,
                "confidence": confidence, "reference_support": "canonical cochlear markers + within-dataset cluster markers; no external automated reference",
                "n_cells": int(row.n_cells), "n_samples_present": int(row.n_samples_present),
                "dominant_sample": row.dominant_sample, "dominant_sample_fraction": row.dominant_sample_fraction,
                "normalized_sample_entropy": row.normalized_sample_entropy,
            })
    evidence = pd.DataFrame(evidence_rows).sort_values(["tissue", "cluster"], key=lambda x: x.map(lambda y: int(y) if str(y).isdigit() else y))
    evidence_path = os.path.join(root, "06_annotation", "S1_ANNOTATION_EVIDENCE.tsv")
    os.makedirs(os.path.dirname(evidence_path), exist_ok=True)
    evidence.to_csv(evidence_path, sep="\t", index=False)
    digest = hashlib.sha256(open(evidence_path, "rb").read()).hexdigest()
    with open(os.path.join(root, "06_annotation", "ANNOTATION_FREEZE.sha256"), "w") as handle:
        handle.write(f"{digest}  S1_ANNOTATION_EVIDENCE.tsv\n")

    for tissue in ["SGN", "Lateral_wall"]:
        path = os.path.join(root, "05_clustering", f"{tissue}_clustered.h5ad")
        adata = ad.read_h5ad(path)
        counts_sum = int(adata.layers["counts"].sum())
        mapping = MAPPINGS[tissue]
        labels = adata.obs["leiden"].astype(str).map(mapping)
        if labels.isna().any() or set(mapping) != set(adata.obs["leiden"].astype(str).unique()):
            raise RuntimeError(f"Incomplete cluster mapping for {tissue}")
        adata.obs["broad_annotation"] = pd.Categorical(labels)
        adata.uns["annotation_freeze_sha256"] = digest
        adata.uns["annotation_prohibited_genes_used"] = False
        panel_parts.append(marker_panel_table(adata, tissue))
        if int(adata.layers["counts"].sum()) != counts_sum or adata.layers["counts"].dtype.kind not in "iu":
            raise RuntimeError(f"Counts changed during annotation: {tissue}")
        output = os.path.join(root, "06_annotation", f"{tissue}_annotated.h5ad")
        adata.write_h5ad(output, compression="gzip")
        plot_umap(adata, tissue, os.path.join(root, "12_figures", f"{tissue}_UMAP_Broad_annotation"))

    panel = pd.concat(panel_parts, ignore_index=True)
    panel.to_csv(os.path.join(root, "06_annotation", "S1_ANNOTATION_MARKER_PANEL.tsv"), sep="\t", index=False)
    detected = panel.dropna(subset=["mean_log_normalized"]).copy()
    means = detected.pivot_table(index=["tissue", "cluster"], columns="gene_symbol", values="mean_log_normalized")
    z = (means - means.mean(axis=0)) / means.std(axis=0, ddof=0).replace(0, 1)
    for tissue in ["SGN", "Lateral_wall"]:
        block = z.loc[tissue]
        fig, ax = plt.subplots(figsize=(16, max(5, 0.32 * len(block))))
        sns.heatmap(block, cmap="vlag", center=0, ax=ax, cbar_kws={"label": "Cluster mean z-score"})
        ax.set(title=f"{tissue}: canonical marker evidence", xlabel="Canonical marker", ylabel="Leiden cluster")
        fig.tight_layout(); save_pair(fig, os.path.join(root, "12_figures", f"{tissue}_canonical_marker_heatmap"))
        block.to_csv(os.path.join(root, "13_source_data", f"FIG_{tissue}_CANONICAL_MARKER_ZSCORES.tsv"), sep="\t")
    print(f"Annotation frozen: {digest}")


if __name__ == "__main__":
    main()
