#!/usr/bin/env python3
import argparse
import gzip
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


SEED = 20260825
np.random.seed(SEED)


def read_gzip_lines(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def load_sample(row):
    matrix = sc.read_mtx(row.matrix_file).X.T.tocsr()
    if matrix.dtype.kind not in "iu":
        if matrix.data.size and not np.allclose(matrix.data, np.rint(matrix.data)):
            raise ValueError(f"Non-integer counts in {row.GSM}")
        matrix.data = np.rint(matrix.data).astype(np.int32)
        matrix = matrix.astype(np.int32)
    else:
        matrix = matrix.astype(np.int32)
    barcodes = read_gzip_lines(row.barcode_file)
    features = [x.split("\t") for x in read_gzip_lines(row.feature_file)]
    gene_ids = [x[0] for x in features]
    symbols = [x[1] if len(x) > 1 else "" for x in features]
    if matrix.shape != (len(barcodes), len(gene_ids)):
        raise ValueError(f"Dimension mismatch for {row.GSM}: {matrix.shape}")
    obs_names = pd.Index([f"{row.GSM}_{x}" for x in barcodes], name="barcode")
    obs = pd.DataFrame(index=obs_names)
    obs["original_barcode"] = barcodes
    for col in ["sample_id", "GSM", "tissue", "condition", "replicate"]:
        obs[col] = getattr(row, col)
    var = pd.DataFrame(index=pd.Index(gene_ids, name="gene_id"))
    var["gene_symbol"] = symbols
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.layers["counts"] = matrix.copy()
    symbol = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = symbol.str.match(r"(?i)^mt-")
    adata.var["hb"] = symbol.str.match(r"(?i)^(hba|hbb|hbz)")
    adata.var["ribo"] = symbol.str.match(r"(?i)^rp[sl][0-9]")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "hb", "ribo"], percent_top=None, log1p=False, inplace=True
    )
    adata.obs["genes_per_umi"] = adata.obs["n_genes_by_counts"] / adata.obs["total_counts"].clip(lower=1)
    counts = adata.layers["counts"]
    assert sparse.isspmatrix_csr(counts)
    assert counts.dtype.kind in "iu"
    assert counts.data.size == 0 or (counts.data >= 0).all()
    return adata


def quantile(series, q):
    return float(series.quantile(q))


def save_pair(fig, stem):
    fig.savefig(stem + ".png", dpi=240, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    manifest = pd.read_csv(os.path.join(root, "00_manifest", "S1_SAMPLE_MANIFEST.tsv"), sep="\t")
    metrics_parts = []
    summary_rows = []

    for row in manifest.itertuples(index=False):
        print(f"Loading {row.GSM}", flush=True)
        adata = load_sample(row)
        obs = adata.obs.copy()
        keep_cols = [
            "original_barcode", "sample_id", "GSM", "tissue", "condition", "replicate",
            "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb",
            "pct_counts_ribo", "genes_per_umi",
        ]
        obs = obs[keep_cols]
        obs.insert(0, "barcode", obs.index)
        metrics_parts.append(obs.reset_index(drop=True))
        summary_rows.append({
            "sample_id": row.sample_id,
            "GSM": row.GSM,
            "tissue": row.tissue,
            "condition": row.condition,
            "replicate": row.replicate,
            "n_cells_raw": adata.n_obs,
            "n_genes_raw": adata.n_vars,
            "aggregate_counts": int(adata.obs["total_counts"].sum()),
            "median_total_counts": float(adata.obs["total_counts"].median()),
            "p01_total_counts": quantile(adata.obs["total_counts"], 0.01),
            "p99_total_counts": quantile(adata.obs["total_counts"], 0.99),
            "median_n_genes": float(adata.obs["n_genes_by_counts"].median()),
            "p01_n_genes": quantile(adata.obs["n_genes_by_counts"], 0.01),
            "p99_n_genes": quantile(adata.obs["n_genes_by_counts"], 0.99),
            "median_pct_mt": float(adata.obs["pct_counts_mt"].median()),
            "p95_pct_mt": quantile(adata.obs["pct_counts_mt"], 0.95),
            "p99_pct_mt": quantile(adata.obs["pct_counts_mt"], 0.99),
            "median_pct_hb": float(adata.obs["pct_counts_hb"].median()),
            "p95_pct_hb": quantile(adata.obs["pct_counts_hb"], 0.95),
            "median_pct_ribo": float(adata.obs["pct_counts_ribo"].median()),
            "median_genes_per_umi": float(adata.obs["genes_per_umi"].median()),
        })

    metrics = pd.concat(metrics_parts, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    qc_dir = os.path.join(root, "02_qc")
    fig_dir = os.path.join(root, "12_figures")
    source_dir = os.path.join(root, "13_source_data")
    os.makedirs(qc_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(source_dir, exist_ok=True)
    metrics.to_csv(os.path.join(qc_dir, "S1_QC_METRICS_UNFILTERED.tsv.gz"), sep="\t", index=False, compression="gzip")
    summary.to_csv(os.path.join(qc_dir, "S1_SAMPLE_QC_SUMMARY_RAW.tsv"), sep="\t", index=False)
    summary.to_csv(os.path.join(source_dir, "FIG_QC_SAMPLE_SUMMARY.tsv"), sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="notebook")
    order = manifest["GSM"].tolist()
    palette = {"SGN": "#4C78A8", "Lateral_wall": "#F58518"}

    fig, ax = plt.subplots(figsize=(11, 4.5))
    tmp = summary.set_index("GSM").loc[order].reset_index()
    sns.barplot(data=tmp, x="GSM", y="n_cells_raw", hue="tissue", palette=palette, dodge=False, ax=ax)
    ax.tick_params(axis="x", rotation=60)
    ax.set(xlabel=None, ylabel="Raw cells", title="GSE168041 raw cells per sample")
    ax.legend(title=None, frameon=False)
    save_pair(fig, os.path.join(fig_dir, "QC_01_raw_cell_counts"))

    plot_specs = [
        ("n_genes_by_counts", "Detected genes per cell", "QC_02_nFeature_distribution", True),
        ("total_counts", "UMIs per cell", "QC_03_nCount_distribution", True),
        ("pct_counts_mt", "Mitochondrial counts (%)", "QC_04_mito_distribution", False),
    ]
    for metric, ylabel, filename, logy in plot_specs:
        fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
        for ax, tissue in zip(axes, ["SGN", "Lateral_wall"]):
            d = metrics[metrics.tissue == tissue]
            sns.violinplot(data=d, x="GSM", y=metric, order=[x for x in order if x in set(d.GSM)],
                           inner="quartile", cut=0, linewidth=0.7, color=palette[tissue], ax=ax)
            if logy:
                ax.set_yscale("log")
            ax.set(title=tissue, xlabel=None, ylabel=ylabel)
            ax.tick_params(axis="x", rotation=45)
        fig.tight_layout()
        save_pair(fig, os.path.join(fig_dir, filename))

    sampled = pd.concat(
        [x.sample(min(len(x), 1500), random_state=SEED) for _, x in metrics.groupby("GSM", observed=True)],
        ignore_index=True,
    )
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, tissue in zip(axes, ["SGN", "Lateral_wall"]):
        d = sampled[sampled.tissue == tissue]
        sns.scatterplot(data=d, x="total_counts", y="n_genes_by_counts", hue="GSM", s=7, alpha=0.45,
                        linewidth=0, ax=ax, legend=False)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(tissue)
        ax.set_xlabel("UMIs per cell"); ax.set_ylabel("Detected genes per cell")
    fig.tight_layout()
    save_pair(fig, os.path.join(fig_dir, "QC_05_nFeature_vs_nCount"))

    heat_cols = ["n_cells_raw", "median_total_counts", "median_n_genes", "median_pct_mt",
                 "median_pct_hb", "median_pct_ribo", "median_genes_per_umi"]
    heat = summary.set_index("GSM")[heat_cols]
    heat_z = (heat - heat.mean(axis=0)) / heat.std(axis=0, ddof=0).replace(0, 1)
    heat_z.to_csv(os.path.join(source_dir, "FIG_QC_HEATMAP_ZSCORES.tsv"), sep="\t")
    fig, ax = plt.subplots(figsize=(9, 7))
    sns.heatmap(heat_z, cmap="vlag", center=0, linewidths=0.25, ax=ax, cbar_kws={"label": "Across-sample z-score"})
    ax.set(title="Sample-level QC audit", xlabel=None, ylabel=None)
    fig.tight_layout()
    save_pair(fig, os.path.join(fig_dir, "QC_06_sample_metrics_heatmap"))

    versions = pd.DataFrame({"software": ["Python", "scanpy", "anndata", "numpy", "pandas", "scipy"],
                             "version": [sys.version.split()[0], version("scanpy"), version("anndata"),
                                         version("numpy"), version("pandas"), version("scipy")]})
    versions.to_csv(os.path.join(root, "14_methods", "PYTHON_CORE_VERSIONS.tsv"), sep="\t", index=False)
    print(f"QC audit complete: {len(metrics)} cells", flush=True)


if __name__ == "__main__":
    main()
