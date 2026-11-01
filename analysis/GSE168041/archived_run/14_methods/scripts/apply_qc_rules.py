#!/usr/bin/env python3
import argparse
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    path = os.path.join(root, "02_qc", "S1_QC_METRICS_UNFILTERED.tsv.gz")
    data = pd.read_csv(path, sep="\t")
    reasons = []
    for row in data.itertuples(index=False):
        failed = []
        if row.n_genes_by_counts < 500:
            failed.append("n_genes_lt_500")
        if row.n_genes_by_counts > 7500:
            failed.append("n_genes_gt_7500")
        if row.pct_counts_mt > 30:
            failed.append("pct_mt_gt_30")
        if row.pct_counts_hb > 20:
            failed.append("pct_hb_gt_20")
        reasons.append("PASS" if not failed else ";".join(failed))
    data["keep"] = [x == "PASS" for x in reasons]
    data["reason"] = reasons
    out_cols = [
        "barcode", "original_barcode", "sample_id", "GSM", "tissue", "condition", "replicate",
        "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo",
        "genes_per_umi", "keep", "reason",
    ]
    data[out_cols].to_csv(os.path.join(root, "02_qc", "S1_QC_CELLS.tsv"), sep="\t", index=False)
    summary = data.groupby(["sample_id", "GSM", "tissue", "condition", "replicate"], observed=True).agg(
        n_cells_raw=("barcode", "size"), n_cells_qc_pass=("keep", "sum")
    ).reset_index()
    summary["n_cells_qc_drop"] = summary.n_cells_raw - summary.n_cells_qc_pass
    summary["qc_retention_fraction"] = summary.n_cells_qc_pass / summary.n_cells_raw
    summary.to_csv(os.path.join(root, "02_qc", "S1_QC_SAMPLE_COUNTS.tsv"), sep="\t", index=False)
    reason_counts = data.loc[~data.keep].groupby(["GSM", "reason"], observed=True).size().rename("n_cells").reset_index()
    reason_counts.to_csv(os.path.join(root, "02_qc", "S1_QC_DROP_REASONS.tsv"), sep="\t", index=False)
    assert len(data) == 71743
    assert int(data.keep.sum()) == 69036
    print(f"QC frozen: {int(data.keep.sum())}/{len(data)} cells pass")


if __name__ == "__main__":
    main()
