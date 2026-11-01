#!/usr/bin/env python3
import argparse
import gzip
import os
from glob import glob

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def read_gzip_lines(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def read_sample_counts(row):
    matrix = sc.read_mtx(row.matrix_file).X.T.tocsr()
    if matrix.dtype.kind not in "iu":
        if matrix.data.size and not np.allclose(matrix.data, np.rint(matrix.data)):
            raise ValueError(f"Non-integer counts in {row.GSM}")
        matrix = matrix.astype(np.int32)
    else:
        matrix = matrix.astype(np.int32)
    barcodes = read_gzip_lines(row.barcode_file)
    features = [x.split("\t") for x in read_gzip_lines(row.feature_file)]
    gene_ids = [x[0] for x in features]
    symbols = [x[1] if len(x) > 1 else "" for x in features]
    obs_names = pd.Index([f"{row.GSM}_{x}" for x in barcodes], name="barcode")
    obs = pd.DataFrame(index=obs_names)
    obs["original_barcode"] = barcodes
    for col in ["sample_id", "GSM", "tissue", "condition", "replicate"]:
        obs[col] = getattr(row, col)
    var = pd.DataFrame({"gene_symbol": symbols}, index=pd.Index(gene_ids, name="gene_id"))
    return ad.AnnData(X=matrix, obs=obs, var=var)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    manifest = pd.read_csv(os.path.join(root, "00_manifest", "S1_SAMPLE_MANIFEST.tsv"), sep="\t")
    qc = pd.read_csv(os.path.join(root, "02_qc", "S1_QC_CELLS.tsv"), sep="\t")
    qc["keep"] = qc["keep"].astype(str).str.lower().eq("true")

    doublet_files = sorted(glob(os.path.join(root, "03_doublet", "per_sample", "GSM*_scDblFinder.tsv")))
    if len(doublet_files) != 16:
        raise RuntimeError(f"Expected 16 doublet files, found {len(doublet_files)}")
    evaluated = pd.concat([pd.read_csv(x, sep="\t") for x in doublet_files], ignore_index=True)
    evaluated["keep"] = evaluated["keep"].astype(str).str.lower().eq("true")
    if evaluated.barcode.duplicated().any() or len(evaluated) != int(qc.keep.sum()):
        raise RuntimeError("Doublet barcode completeness/uniqueness failed")
    evaluated = evaluated.rename(columns={"keep": "doublet_keep"})
    evaluated["drop_reason"] = np.where(evaluated.doublet_keep, "PASS", "scDblFinder_doublet")
    evaluated.to_csv(os.path.join(root, "03_doublet", "S1_DOUBLET.tsv"), sep="\t", index=False)

    merged = qc.merge(
        evaluated[["barcode", "doublet_score", "doublet_class", "doublet_keep", "seed"]],
        on="barcode", how="left", validate="one_to_one",
    )
    merged["doublet_class"] = merged["doublet_class"].fillna("not_evaluated_qc_fail")
    merged["doublet_keep"] = merged["doublet_keep"].fillna(False).astype(bool)
    merged["final_keep"] = merged["keep"] & merged["doublet_keep"]
    merged["final_reason"] = np.select(
        [~merged["keep"], merged["keep"] & ~merged["doublet_keep"]],
        ["QC:" + merged["reason"], "scDblFinder_doublet"],
        default="PASS",
    )
    merged.to_csv(os.path.join(root, "03_doublet", "S1_CELL_FILTER_STATUS.tsv.gz"), sep="\t", index=False, compression="gzip")

    summary = merged.groupby(["sample_id", "GSM", "tissue", "condition", "replicate"], observed=True).agg(
        n_cells_raw=("barcode", "size"),
        n_cells_qc_pass=("keep", "sum"),
        n_doublets=("doublet_class", lambda x: int((x == "doublet").sum())),
        n_cells_final=("final_keep", "sum"),
    ).reset_index()
    summary["doublet_rate_among_qc_pass"] = summary.n_doublets / summary.n_cells_qc_pass
    summary["final_retention_fraction"] = summary.n_cells_final / summary.n_cells_raw
    summary.to_csv(os.path.join(root, "03_doublet", "S1_SAMPLE_CELL_COUNTS_POST_DOUBLET.tsv"), sep="\t", index=False)

    final_ids = set(merged.loc[merged.final_keep, "barcode"])
    qc_by_barcode = merged.set_index("barcode")
    compartments = {"SGN": [], "Lateral_wall": []}
    sample_audit = []
    reference_gene_ids = None
    reference_symbols = None
    for row in manifest.itertuples(index=False):
        print(f"Building filtered object: {row.GSM}", flush=True)
        a = read_sample_counts(row)
        if reference_gene_ids is None:
            reference_gene_ids = a.var_names.copy()
            reference_symbols = a.var["gene_symbol"].copy()
        else:
            if not a.var_names.equals(reference_gene_ids) or not a.var["gene_symbol"].equals(reference_symbols):
                raise RuntimeError(f"Feature order/content differs in {row.GSM}")
        keep_mask = a.obs_names.isin(final_ids)
        b = a[keep_mask, :].copy()
        meta = qc_by_barcode.loc[b.obs_names]
        for col in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb",
                    "pct_counts_ribo", "genes_per_umi", "doublet_score", "doublet_class"]:
            b.obs[col] = meta[col].to_numpy()
        b.layers["counts"] = b.X.copy()
        original_sum = int(a.X[keep_mask, :].sum())
        layer_sum = int(b.layers["counts"].sum())
        if original_sum != layer_sum:
            raise RuntimeError(f"Count traceability failed for {row.GSM}")
        sample_audit.append({
            "GSM": row.GSM, "tissue": row.tissue, "n_final": b.n_obs, "n_genes": b.n_vars,
            "original_subset_count_sum": original_sum, "counts_layer_sum": layer_sum,
            "integer_counts": b.layers["counts"].dtype.kind in "iu",
            "barcode_unique": b.obs_names.is_unique, "status": "PASS",
        })
        compartments[row.tissue].append(b)

    os.makedirs(os.path.join(root, "04_objects"), exist_ok=True)
    combined_audit = []
    for tissue, objects in compartments.items():
        combined = ad.concat(objects, axis=0, join="inner", merge="same", index_unique=None)
        if combined.n_vars != len(reference_gene_ids):
            raise RuntimeError(f"Gene intersection unexpectedly changed in {tissue}")
        combined.strings_to_categoricals()
        counts = combined.layers["counts"]
        if not sparse.issparse(counts) or counts.dtype.kind not in "iu":
            raise RuntimeError(f"Counts-layer type failure in {tissue}")
        if not combined.obs_names.is_unique or combined.obs[["sample_id", "GSM", "condition", "tissue", "replicate"]].isna().any().any():
            raise RuntimeError(f"Metadata/barcode failure in {tissue}")
        combined.uns["counts_provenance"] = "Official GSE168041 Matrix Market integer counts; QC and per-GSM scDblFinder singlets only"
        combined.uns["random_seed"] = 20260825
        path = os.path.join(root, "04_objects", f"{tissue}_counts_singlets.h5ad")
        combined.write_h5ad(path, compression="gzip")
        combined_audit.append({
            "object": os.path.basename(path), "tissue": tissue, "n_cells": combined.n_obs,
            "n_genes": combined.n_vars, "counts_sum": int(counts.sum()), "integer_counts": True,
            "barcode_unique": True, "metadata_complete": True, "status": "PASS",
        })
        print(f"Wrote {path}: {combined.n_obs} x {combined.n_vars}", flush=True)

    pd.DataFrame(sample_audit).to_csv(os.path.join(root, "04_objects", "S1_COUNTS_LAYER_AUDIT_BY_SAMPLE.tsv"), sep="\t", index=False)
    pd.DataFrame(combined_audit).to_csv(os.path.join(root, "04_objects", "S1_COUNTS_LAYER_AUDIT_OBJECTS.tsv"), sep="\t", index=False)
    if int(summary.n_cells_final.sum()) != 64705:
        raise RuntimeError("Unexpected final cell total")
    print("Phase 2 objects PASS: 64705 singlets", flush=True)


if __name__ == "__main__":
    main()
