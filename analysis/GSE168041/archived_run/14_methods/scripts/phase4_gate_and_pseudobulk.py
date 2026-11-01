#!/usr/bin/env python3
import argparse
import os
import re

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


MIN_CELLS = 20
MIN_REPLICATES = 3
MAX_DOMINANCE = 0.50


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def save_pair(fig, stem):
    fig.savefig(stem + ".png", dpi=240, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    gate_rows = []
    audit_rows = []
    os.makedirs(os.path.join(root, "08_pseudobulk"), exist_ok=True)

    for tissue in ["SGN", "Lateral_wall"]:
        adata = ad.read_h5ad(os.path.join(root, "06_annotation", f"{tissue}_annotated.h5ad"))
        counts = adata.layers["counts"].tocsr()
        if counts.dtype.kind not in "iu" or (counts.data < 0).any():
            raise RuntimeError(f"Non-integer/negative counts: {tissue}")
        total_input_sum = int(counts.sum())
        sample_meta = adata.obs[["sample_id", "GSM", "condition", "replicate"]].drop_duplicates().copy()
        sample_meta["condition"] = sample_meta["condition"].astype(str)
        sample_meta["condition_order"] = sample_meta.condition.map({"Naive": 0, "Noise": 1})
        sample_meta = sample_meta.sort_values(["condition_order", "replicate"]).drop(columns="condition_order")
        order = sample_meta.GSM.tolist()
        cell_types = sorted(adata.obs["broad_annotation"].astype(str).unique())
        reconstructed_sum = 0

        cell_count_matrix = pd.crosstab(adata.obs["broad_annotation"], adata.obs["GSM"]).reindex(index=cell_types, columns=order, fill_value=0)
        fig, ax = plt.subplots(figsize=(10, max(4, 0.48 * len(cell_types))))
        sns.heatmap(np.log10(cell_count_matrix + 1), cmap="viridis", annot=cell_count_matrix, fmt="d", ax=ax,
                    cbar_kws={"label": "log10(cells + 1)"})
        ax.set(title=f"{tissue}: Broad cell type × sample counts", xlabel="GSM", ylabel="Broad cell type")
        fig.tight_layout(); save_pair(fig, os.path.join(root, "12_figures", f"{tissue}_Broad_celltype_by_sample_heatmap"))
        cell_count_matrix.to_csv(os.path.join(root, "13_source_data", f"FIG_{tissue}_BROAD_CELLTYPE_SAMPLE_COUNTS.tsv"), sep="\t")

        for cell_type in cell_types:
            pb_rows = []
            pb_meta_rows = []
            for sm in sample_meta.itertuples(index=False):
                mask = ((adata.obs["GSM"].astype(str).to_numpy() == sm.GSM) &
                        (adata.obs["broad_annotation"].astype(str).to_numpy() == cell_type))
                vector = np.asarray(counts[mask, :].sum(axis=0)).ravel().astype(np.int64)
                if (vector < 0).any():
                    raise RuntimeError("Negative pseudobulk count")
                n_cells = int(mask.sum())
                library_size = int(vector.sum())
                reconstructed_sum += library_size
                pb_rows.append(vector)
                pb_meta_rows.append({
                    "sample_id": sm.sample_id, "GSM": sm.GSM, "condition": sm.condition,
                    "replicate": int(sm.replicate), "tissue": tissue, "cell_type": cell_type,
                    "n_cells": n_cells, "raw_library_size": library_size,
                })
            pb = np.vstack(pb_rows)
            meta = pd.DataFrame(pb_meta_rows).set_index("sample_id", drop=False)
            pb_df = pd.DataFrame(pb, index=meta.index, columns=adata.var_names)
            if not np.issubdtype(pb_df.values.dtype, np.integer):
                raise RuntimeError("Pseudobulk matrix not integer")

            condition_stats = {}
            for condition in ["Naive", "Noise"]:
                x = meta[meta.condition == condition]
                eligible = x.n_cells >= MIN_CELLS
                condition_stats[condition] = {
                    "n_ge_min": int(eligible.sum()),
                    "min_cells": int(x.n_cells.min()),
                    "total_cells": int(x.n_cells.sum()),
                    "dominance": float(x.n_cells.max() / x.n_cells.sum()) if x.n_cells.sum() else np.nan,
                    "min_library": int(x.raw_library_size.min()),
                }
            reasons = []
            if cell_type == "UNRESOLVED": reasons.append("unresolved_annotation")
            for condition in ["Naive", "Noise"]:
                st = condition_stats[condition]
                if st["n_ge_min"] < MIN_REPLICATES: reasons.append(f"{condition.lower()}_replicates_lt_{MIN_REPLICATES}")
                if st["dominance"] > MAX_DOMINANCE: reasons.append(f"{condition.lower()}_sample_dominance_gt_{MAX_DOMINANCE}")
                if st["min_library"] <= 0: reasons.append(f"{condition.lower()}_nonpositive_library")
            status = "FORMAL_INFERENCE" if not reasons else "DESCRIPTIVE_ONLY"
            gate_rows.append({
                "tissue": tissue, "cell_type": cell_type,
                "naive_samples_ge20": condition_stats["Naive"]["n_ge_min"],
                "noise_samples_ge20": condition_stats["Noise"]["n_ge_min"],
                "naive_min_cells": condition_stats["Naive"]["min_cells"],
                "noise_min_cells": condition_stats["Noise"]["min_cells"],
                "naive_dominant_sample_fraction": condition_stats["Naive"]["dominance"],
                "noise_dominant_sample_fraction": condition_stats["Noise"]["dominance"],
                "naive_min_library_size": condition_stats["Naive"]["min_library"],
                "noise_min_library_size": condition_stats["Noise"]["min_library"],
                "status": status, "reason": "PASS" if not reasons else ";".join(reasons),
            })

            out_dir = os.path.join(root, "08_pseudobulk", tissue, safe_name(cell_type))
            os.makedirs(out_dir, exist_ok=True)
            pb_df.to_csv(os.path.join(out_dir, "counts_samples_by_gene.tsv.gz"), sep="\t", compression="gzip", index_label="sample_id")
            meta.to_csv(os.path.join(out_dir, "metadata.tsv"), sep="\t", index=False)
            gene_meta = adata.var[["gene_symbol"]].copy()
            gene_meta.insert(0, "gene_id", gene_meta.index)
            gene_meta.to_csv(os.path.join(out_dir, "gene_metadata.tsv"), sep="\t", index=False)
            audit_rows.append({
                "tissue": tissue, "cell_type": cell_type, "directory": out_dir,
                "n_samples": pb_df.shape[0], "n_genes": pb_df.shape[1],
                "aggregate_count_sum": int(pb_df.to_numpy().sum()), "integer_counts": True,
                "sample_identity_match": list(pb_df.index) == list(meta.index), "status": "PASS",
            })

        if reconstructed_sum != total_input_sum:
            raise RuntimeError(f"Pseudobulk conservation failed for {tissue}: {reconstructed_sum} != {total_input_sum}")

    gate = pd.DataFrame(gate_rows)
    gate.to_csv(os.path.join(root, "07_replicate_gate", "S1_REPLICATE_GATE.tsv"), sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(os.path.join(root, "08_pseudobulk", "S1_PSEUDOBULK_AUDIT.tsv"), sep="\t", index=False)
    print(gate.to_string(index=False))
    print("Phase 4/5 pseudobulk conservation PASS", flush=True)


if __name__ == "__main__":
    main()
