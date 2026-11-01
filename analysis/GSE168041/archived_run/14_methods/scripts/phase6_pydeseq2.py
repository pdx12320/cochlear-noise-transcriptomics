#!/usr/bin/env python3
import argparse
import os
import re
import sys
from importlib.metadata import version

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA


SEED = 20260825
MIN_TOTAL_COUNT = 10
ALPHA = 0.05
N_CPUS = 8


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def save_pair(fig, stem):
    fig.savefig(stem + ".png", dpi=240, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    plt.close(fig)


def plot_pb_pca(normed_counts, metadata, tissue, cell_type, root):
    transformed = np.log1p(normed_counts.astype(float))
    variances = transformed.var(axis=0)
    keep = np.argsort(variances)[-min(2000, transformed.shape[1]):]
    coords = PCA(n_components=2, random_state=SEED).fit_transform(transformed[:, keep])
    source = metadata.copy()
    source["PC1"] = coords[:, 0]
    source["PC2"] = coords[:, 1]
    safe = safe_name(cell_type)
    source.to_csv(os.path.join(root, "13_source_data", f"PSEUDOBULK_PCA_{tissue}_{safe}.tsv"), sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    colors = {"Naive": "#4C78A8", "Noise": "#E45756"}
    for condition, group in source.groupby("condition"):
        ax.scatter(group.PC1, group.PC2, s=60, color=colors[str(condition)], label=str(condition), edgecolor="white", linewidth=0.6)
        for row in group.itertuples(index=False):
            ax.text(row.PC1, row.PC2, str(row.GSM)[-4:], fontsize=7, ha="left", va="bottom")
    ax.set(title=f"{tissue} — {cell_type}\npseudobulk PCA", xlabel="PC1", ylabel="PC2")
    ax.legend(frameon=False)
    fig.tight_layout(); save_pair(fig, os.path.join(root, "12_figures", f"PSEUDOBULK_PCA_{tissue}_{safe}"))


def run_model(root, gate_row):
    tissue = gate_row.tissue
    cell_type = gate_row.cell_type
    safe = safe_name(cell_type)
    pb_dir = os.path.join(root, "08_pseudobulk", tissue, safe)
    counts = pd.read_csv(os.path.join(pb_dir, "counts_samples_by_gene.tsv.gz"), sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(pb_dir, "metadata.tsv"), sep="\t", index_col="sample_id")
    gene_meta = pd.read_csv(os.path.join(pb_dir, "gene_metadata.tsv"), sep="\t").set_index("gene_id")
    if list(counts.index) != list(meta.index):
        raise RuntimeError(f"Sample order mismatch: {tissue} {cell_type}")
    if counts.shape[0] != 8 or meta.condition.value_counts().to_dict() != {"Naive": 4, "Noise": 4}:
        raise RuntimeError(f"Biological replicate mismatch: {tissue} {cell_type}")
    array = counts.to_numpy()
    if not np.issubdtype(array.dtype, np.integer) or (array < 0).any():
        raise RuntimeError(f"Counts invalid: {tissue} {cell_type}")
    keep_genes = counts.sum(axis=0) >= MIN_TOTAL_COUNT
    filtered = counts.loc[:, keep_genes].astype(np.int64)
    metadata = meta.copy()
    metadata["condition"] = pd.Categorical(metadata["condition"], categories=["Naive", "Noise"], ordered=True)
    inference = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(
        counts=filtered,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
        low_memory=False,
    )
    dds.deseq2()
    ds = DeseqStats(
        dds,
        contrast=["condition", "Noise", "Naive"],
        alpha=ALPHA,
        cooks_filter=True,
        independent_filter=True,
        inference=inference,
        quiet=False,
    )
    ds.summary()
    result = ds.results_df.copy()
    result.index.name = "gene_id"
    result.insert(0, "gene_symbol", gene_meta.loc[result.index, "gene_symbol"].to_numpy())
    result.insert(0, "gene_id", result.index)
    required = ["gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    result = result[required]
    out_dir = os.path.join(root, "09_differential", tissue)
    os.makedirs(out_dir, exist_ok=True)
    result_path = os.path.join(out_dir, f"S1_DS_ALL_{safe}.tsv")
    result.to_csv(result_path, sep="\t", index=False)
    portable = dds.to_picklable_anndata()
    for key, value in list(portable.uns.items()):
        if isinstance(value, pd.Series):
            portable.uns[key] = value.rename("value").to_frame()
    portable.write_h5ad(os.path.join(out_dir, f"S1_DDS_OBJECT_{safe}.h5ad"), compression="gzip")

    size_factors = pd.DataFrame({
        "sample_id": dds.obs_names, "GSM": meta.loc[dds.obs_names, "GSM"].to_numpy(),
        "condition": meta.loc[dds.obs_names, "condition"].to_numpy(),
        "size_factor": dds.obs["size_factors"].to_numpy(),
        "n_cells": meta.loc[dds.obs_names, "n_cells"].to_numpy(),
        "raw_library_size": meta.loc[dds.obs_names, "raw_library_size"].to_numpy(),
    })
    size_factors.to_csv(os.path.join(out_dir, f"S1_DDS_QA_{safe}.tsv"), sep="\t", index=False)
    normed_layer = dds.layers["normed_counts"]
    normed = normed_layer.toarray() if hasattr(normed_layer, "toarray") else np.asarray(normed_layer)
    plot_pb_pca(normed, meta.loc[dds.obs_names].reset_index(), tissue, cell_type, root)

    if tissue == "Lateral_wall" and cell_type == "Endothelial":
        plot = result.copy()
        plot["minus_log10_padj"] = -np.log10(plot.padj.clip(lower=np.nextafter(0, 1)))
        sig = plot.padj < ALPHA
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(plot.loc[~sig, "log2FoldChange"], plot.loc[~sig, "minus_log10_padj"], s=5, alpha=0.3, color="#9E9E9E", rasterized=True)
        ax.scatter(plot.loc[sig, "log2FoldChange"], plot.loc[sig, "minus_log10_padj"], s=7, alpha=0.65, color="#D95F02", rasterized=True)
        ax.axvline(0, color="black", linewidth=0.7); ax.axhline(-np.log10(ALPHA), color="black", linewidth=0.7, linestyle="--")
        ax.set(xlabel="log2 fold change (Noise vs Naive)", ylabel="-log10(FDR)", title="Lateral-wall Endothelial pseudobulk DE")
        fig.tight_layout(); save_pair(fig, os.path.join(root, "12_figures", "Lateral_wall_Endothelial_volcano"))
        plot.to_csv(os.path.join(root, "13_source_data", "FIG_Lateral_wall_Endothelial_volcano.tsv"), sep="\t", index=False)

    sig = result[result.padj.notna() & (result.padj < ALPHA)]
    return {
        "tissue": tissue, "cell_type": cell_type,
        "n_naive_samples": int((meta.condition == "Naive").sum()),
        "n_noise_samples": int((meta.condition == "Noise").sum()),
        "naive_cells": int(meta.loc[meta.condition == "Naive", "n_cells"].sum()),
        "noise_cells": int(meta.loc[meta.condition == "Noise", "n_cells"].sum()),
        "genes_input": counts.shape[1], "genes_tested": len(result),
        "genes_filtered_total_count_lt10": int((~keep_genes).sum()),
        "fdr_lt_0_05": len(sig), "fdr_up_noise": int((sig.log2FoldChange > 0).sum()),
        "fdr_down_noise": int((sig.log2FoldChange < 0).sum()),
        "finite_wald_stat": int(np.isfinite(result.stat).sum()), "status": "PASS",
        "result_file": result_path,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    np.random.seed(SEED)
    gate = pd.read_csv(os.path.join(root, "07_replicate_gate", "S1_REPLICATE_GATE.tsv"), sep="\t")
    eligible = gate[gate.status == "FORMAL_INFERENCE"]
    rows = []
    for row in eligible.itertuples(index=False):
        print(f"PyDESeq2: {row.tissue} / {row.cell_type}", flush=True)
        rows.append(run_model(root, row))
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(root, "09_differential", "S1_DS_SUMMARY.tsv"), sep="\t", index=False)
    params = pd.DataFrame({
        "parameter": ["Python", "PyDESeq2", "design", "contrast", "reference", "min_total_count",
                      "alpha", "refit_cooks", "cooks_filter", "independent_filter", "n_cpus", "seed"],
        "value": [sys.version.split()[0], version("pydeseq2"), "~condition", "Noise_vs_Naive", "Naive",
                  MIN_TOTAL_COUNT, ALPHA, True, True, True, N_CPUS, SEED],
    })
    params.to_csv(os.path.join(root, "14_methods", "PYDESEQ2_PARAMETERS.tsv"), sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)
    print("Phase 6 complete", flush=True)


if __name__ == "__main__":
    main()
