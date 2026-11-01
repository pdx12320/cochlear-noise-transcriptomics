#!/usr/bin/env python3
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    os.makedirs(os.path.join(root, "07_replicate_gate"), exist_ok=True)
    all_rows = []
    for tissue in ["SGN", "Lateral_wall"]:
        adata = ad.read_h5ad(os.path.join(root, "06_annotation", f"{tissue}_annotated.h5ad"))
        samples = adata.obs[["GSM", "condition", "replicate"]].drop_duplicates().copy()
        samples["condition_order"] = samples["condition"].astype(str).map({"Naive": 0, "Noise": 1})
        samples = samples.sort_values(["condition_order", "replicate"]).drop(columns="condition_order")
        order = samples.GSM.tolist()
        cell_counts = pd.crosstab(adata.obs["broad_annotation"], adata.obs["GSM"]).reindex(columns=order, fill_value=0)
        cell_counts.index.name = "cell_type"
        cell_counts.to_csv(os.path.join(root, "07_replicate_gate", f"S1_CELLTYPE_SAMPLE_COUNTS_{tissue}.tsv"), sep="\t")

        group_keys = adata.obs[["GSM", "broad_annotation"]].astype(str).agg("||".join, axis=1)
        library_rows = []
        for key in sorted(group_keys.unique()):
            gsm, cell_type = key.split("||", 1)
            mask = group_keys.to_numpy() == key
            library_rows.append({"tissue": tissue, "cell_type": cell_type, "GSM": gsm,
                                 "n_cells": int(mask.sum()), "raw_library_size": int(adata.layers["counts"][mask, :].sum())})
        library = pd.DataFrame(library_rows)
        library.to_csv(os.path.join(root, "07_replicate_gate", f"S1_CELLTYPE_SAMPLE_LIBRARY_SIZES_{tissue}.tsv"), sep="\t", index=False)

        meta = adata.obs[["GSM", "condition", "replicate"]].drop_duplicates().set_index("GSM")
        for cell_type, row in cell_counts.iterrows():
            rec = {"tissue": tissue, "cell_type": str(cell_type), "total_cells": int(row.sum())}
            for condition in ["Naive", "Noise"]:
                gsm_ids = meta.index[meta.condition.astype(str) == condition]
                values = row.reindex(gsm_ids, fill_value=0)
                rec[f"{condition.lower()}_samples_any"] = int((values > 0).sum())
                rec[f"{condition.lower()}_min_cells"] = int(values.min())
                rec[f"{condition.lower()}_median_cells"] = float(values.median())
                rec[f"{condition.lower()}_max_cells"] = int(values.max())
                rec[f"{condition.lower()}_dominant_sample_fraction"] = float(values.max() / values.sum()) if values.sum() else np.nan
            all_rows.append(rec)
    summary = pd.DataFrame(all_rows)
    summary.to_csv(os.path.join(root, "07_replicate_gate", "S1_CELLTYPE_COVERAGE_PREFREEZE.tsv"), sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
