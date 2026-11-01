#!/usr/bin/env python3
import argparse
import os
import sys
from importlib.metadata import version

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from sklearn.metrics import adjusted_mutual_info_score


SEED = 20260825
N_HVG = 3000
N_PCS = 50
NEIGHBOR_PCS = 40
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 0.6


def save_pair(fig, stem):
    fig.savefig(stem + ".png", dpi=240, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    plt.close(fig)


def categorical_umap(adata, column, title, stem):
    coords = adata.obsm["X_umap"]
    values = adata.obs[column].astype(str)
    levels = sorted(values.unique())
    cmap = plt.get_cmap("tab20")
    colors = {level: cmap(i % 20) for i, level in enumerate(levels)}
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    for level in levels:
        mask = values.to_numpy() == level
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2.2, alpha=0.65, linewidths=0,
                   rasterized=True, color=colors[level], label=level)
    ax.set(title=title, xlabel="UMAP1", ylabel="UMAP2")
    ax.set_xticks([]); ax.set_yticks([])
    ncol = 2 if len(levels) > 10 else 1
    ax.legend(title=column, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
              markerscale=3, fontsize=8, ncol=ncol)
    fig.tight_layout()
    save_pair(fig, stem)


def knn_sample_mixing(adata):
    distances = adata.obsp["distances"].tocsr()
    samples = adata.obs["GSM"].astype(str).to_numpy()
    same_fractions = np.zeros(adata.n_obs, dtype=float)
    for i in range(adata.n_obs):
        idx = distances.indices[distances.indptr[i]:distances.indptr[i + 1]]
        same_fractions[i] = np.mean(samples[idx] == samples[i]) if len(idx) else np.nan
    proportions = pd.Series(samples).value_counts(normalize=True)
    expected = float((proportions ** 2).sum())
    observed = float(np.nanmean(same_fractions))
    return same_fractions, expected, observed


def normalized_entropy(counts):
    values = counts.to_numpy(dtype=float)
    total = values.sum()
    if total <= 0:
        return np.nan
    p = values[values > 0] / total
    return float(-(p * np.log(p)).sum() / np.log(len(values))) if len(values) > 1 else 0.0


def run_compartment(root, tissue):
    input_path = os.path.join(root, "04_objects", f"{tissue}_counts_singlets.h5ad")
    adata = ad.read_h5ad(input_path)
    counts = adata.layers["counts"]
    if not sparse.issparse(counts) or counts.dtype.kind not in "iu" or (counts.data < 0).any():
        raise RuntimeError(f"Raw counts invariant failed before clustering: {tissue}")
    count_sum_before = int(counts.sum())
    if not adata.obs_names.is_unique:
        raise RuntimeError(f"Non-unique barcodes: {tissue}")

    adata.X = counts.astype(np.float32).copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat", n_top_genes=N_HVG, batch_key="GSM", inplace=True
    )
    if int(adata.var["highly_variable"].sum()) != N_HVG:
        raise RuntimeError(f"Unexpected HVG count: {tissue}")
    analysis = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(analysis, max_value=10)
    sc.tl.pca(analysis, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(analysis, n_neighbors=N_NEIGHBORS, n_pcs=NEIGHBOR_PCS, metric="euclidean", random_state=SEED)
    sc.tl.umap(analysis, random_state=SEED, min_dist=0.5)
    sc.tl.leiden(
        analysis, resolution=LEIDEN_RESOLUTION, random_state=SEED, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )
    adata.obs["leiden"] = analysis.obs["leiden"].astype(str).astype("category")
    adata.obsm["X_pca"] = analysis.obsm["X_pca"].copy()
    adata.obsm["X_umap"] = analysis.obsm["X_umap"].copy()
    adata.obsp["distances"] = analysis.obsp["distances"].copy()
    adata.obsp["connectivities"] = analysis.obsp["connectivities"].copy()
    adata.uns["neighbors"] = analysis.uns["neighbors"].copy()
    adata.uns["pca"] = analysis.uns["pca"].copy()
    adata.uns["umap"] = analysis.uns["umap"].copy()
    adata.uns["clustering_parameters"] = {
        "normalization_target_sum": 10000.0,
        "hvg_flavor": "seurat", "n_top_genes": N_HVG, "hvg_batch_key": "GSM",
        "scale_max_value": 10.0, "pca_n_comps": N_PCS, "neighbors_n_neighbors": N_NEIGHBORS,
        "neighbors_n_pcs": NEIGHBOR_PCS, "leiden_resolution": LEIDEN_RESOLUTION,
        "random_seed": SEED, "batch_correction": "none",
    }

    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon", use_raw=False, pts=True, key_added="cluster_markers"
    )
    markers = sc.get.rank_genes_groups_df(adata, group=None, key="cluster_markers")
    markers = markers.rename(columns={"group": "cluster", "names": "gene_id"})
    markers["gene_symbol"] = adata.var.loc[markers.gene_id, "gene_symbol"].to_numpy()
    markers.insert(0, "tissue", tissue)
    markers = markers.sort_values(["cluster", "scores"], ascending=[True, False])
    markers = markers[markers.scores > 0].groupby("cluster", observed=True).head(300).reset_index(drop=True)
    marker_path = os.path.join(root, "05_clustering", f"S1_CLUSTER_MARKERS_{tissue}.tsv")
    markers.to_csv(marker_path, sep="\t", index=False)

    cluster_sample = pd.crosstab(adata.obs["leiden"], adata.obs["GSM"])
    cluster_condition = pd.crosstab(adata.obs["leiden"], adata.obs["condition"])
    cluster_sample.to_csv(os.path.join(root, "05_clustering", f"S1_CLUSTER_BY_SAMPLE_{tissue}.tsv"), sep="\t")
    cluster_condition.to_csv(os.path.join(root, "05_clustering", f"S1_CLUSTER_BY_CONDITION_{tissue}.tsv"), sep="\t")

    composition_rows = []
    for cluster in cluster_sample.index:
        row = cluster_sample.loc[cluster]
        cond = cluster_condition.loc[cluster]
        composition_rows.append({
            "tissue": tissue, "cluster": cluster, "n_cells": int(row.sum()),
            "n_samples_present": int((row > 0).sum()), "n_samples_ge_1pct_cluster": int((row / row.sum() >= 0.01).sum()),
            "dominant_sample": str(row.idxmax()), "dominant_sample_fraction": float(row.max() / row.sum()),
            "normalized_sample_entropy": normalized_entropy(row),
            "naive_fraction": float(cond.get("Naive", 0) / cond.sum()),
            "noise_fraction": float(cond.get("Noise", 0) / cond.sum()),
        })
    composition = pd.DataFrame(composition_rows)
    composition.to_csv(os.path.join(root, "05_clustering", f"S1_CLUSTER_COMPOSITION_AUDIT_{tissue}.tsv"), sep="\t", index=False)

    same_fraction, expected_same, observed_same = knn_sample_mixing(adata)
    adata.obs["knn_same_sample_fraction"] = same_fraction
    sample_ami = adjusted_mutual_info_score(adata.obs["leiden"].astype(str), adata.obs["GSM"].astype(str))
    condition_ami = adjusted_mutual_info_score(adata.obs["leiden"].astype(str), adata.obs["condition"].astype(str))
    audit = pd.DataFrame([{
        "tissue": tissue, "n_cells": adata.n_obs, "n_genes": adata.n_vars,
        "n_clusters": adata.obs["leiden"].nunique(), "sample_cluster_AMI": sample_ami,
        "condition_cluster_AMI": condition_ami, "knn_same_sample_observed": observed_same,
        "knn_same_sample_expected_random": expected_same,
        "knn_sample_enrichment_ratio": observed_same / expected_same,
        "median_cluster_dominant_sample_fraction": composition.dominant_sample_fraction.median(),
        "max_cluster_dominant_sample_fraction": composition.dominant_sample_fraction.max(),
        "min_cluster_sample_entropy": composition.normalized_sample_entropy.min(),
        "batch_correction": "none",
    }])
    audit.to_csv(os.path.join(root, "05_clustering", f"S1_SAMPLE_DRIVEN_CLUSTERING_AUDIT_{tissue}.tsv"), sep="\t", index=False)

    fig_dir = os.path.join(root, "12_figures")
    categorical_umap(adata, "leiden", f"{tissue}: Leiden clusters", os.path.join(fig_dir, f"{tissue}_UMAP_Leiden"))
    categorical_umap(adata, "GSM", f"{tissue}: sample identity", os.path.join(fig_dir, f"{tissue}_UMAP_sample"))
    categorical_umap(adata, "condition", f"{tissue}: condition", os.path.join(fig_dir, f"{tissue}_UMAP_condition"))

    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(cluster_sample))))
    sample_fraction = cluster_sample.div(cluster_sample.sum(axis=1), axis=0)
    sns.heatmap(sample_fraction, cmap="Blues", ax=ax, cbar_kws={"label": "Fraction within cluster"})
    ax.set(title=f"{tissue}: cluster composition by sample", xlabel="GSM", ylabel="Leiden cluster")
    fig.tight_layout()
    save_pair(fig, os.path.join(fig_dir, f"{tissue}_cluster_by_sample_heatmap"))
    sample_fraction.to_csv(os.path.join(root, "13_source_data", f"FIG_{tissue}_CLUSTER_SAMPLE_FRACTIONS.tsv"), sep="\t")

    fig, ax = plt.subplots(figsize=(6, 4))
    variance = analysis.uns["pca"]["variance_ratio"]
    ax.plot(np.arange(1, len(variance) + 1), variance, marker="o", markersize=2, linewidth=1)
    ax.axvline(NEIGHBOR_PCS, color="firebrick", linestyle="--", linewidth=1)
    ax.set(xlabel="PC", ylabel="Explained variance ratio", title=f"{tissue}: PCA variance")
    fig.tight_layout()
    save_pair(fig, os.path.join(fig_dir, f"{tissue}_PCA_variance"))

    source = adata.obs[["sample_id", "GSM", "tissue", "condition", "replicate", "leiden", "knn_same_sample_fraction"]].copy()
    source.insert(0, "barcode", source.index)
    source["UMAP1"] = adata.obsm["X_umap"][:, 0]
    source["UMAP2"] = adata.obsm["X_umap"][:, 1]
    source.to_csv(os.path.join(root, "13_source_data", f"S1_UMAP_SOURCE_{tissue}.tsv.gz"), sep="\t", index=False, compression="gzip")

    if int(adata.layers["counts"].sum()) != count_sum_before or adata.layers["counts"].dtype.kind not in "iu":
        raise RuntimeError(f"Counts layer changed during clustering: {tissue}")
    output_path = os.path.join(root, "05_clustering", f"{tissue}_clustered.h5ad")
    adata.write_h5ad(output_path, compression="gzip")
    print(audit.to_string(index=False), flush=True)
    return marker_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    os.makedirs(os.path.join(root, "05_clustering"), exist_ok=True)
    paths = [run_compartment(root, tissue) for tissue in ["SGN", "Lateral_wall"]]
    combined = pd.concat([pd.read_csv(x, sep="\t") for x in paths], ignore_index=True)
    combined.to_csv(os.path.join(root, "05_clustering", "S1_CLUSTER_MARKERS.tsv"), sep="\t", index=False)
    params = pd.DataFrame({
        "parameter": ["seed", "normalization_target_sum", "hvg_flavor", "n_top_genes", "hvg_batch_key",
                      "scale_max_value", "pca_n_comps", "neighbors_k", "neighbors_n_pcs", "umap_min_dist",
                      "leiden_resolution", "leiden_flavor", "batch_correction", "marker_method"],
        "value": [SEED, 10000, "seurat", N_HVG, "GSM", 10, N_PCS, N_NEIGHBORS, NEIGHBOR_PCS, 0.5,
                  LEIDEN_RESOLUTION, "igraph", "none", "wilcoxon"],
    })
    params.to_csv(os.path.join(root, "14_methods", "SCANPY_PARAMETERS.tsv"), sep="\t", index=False)
    print("Phase 3 clustering complete", flush=True)


if __name__ == "__main__":
    main()
