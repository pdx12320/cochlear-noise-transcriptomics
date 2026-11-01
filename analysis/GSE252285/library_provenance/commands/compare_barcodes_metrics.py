#!/usr/bin/env python3
"""Compute the required 3x3 barcode overlap and metrics comparison."""

import csv
import gzip
import json
import os
import re


PROJECT = os.environ.get(
    "PROJECT", "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
COMPARISON = os.path.join(PROJECT, "comparison")
LIBRARIES = ["A1_Control", "B1_HHL", "C1_NIHL"]
CONDITIONS = {
    "A1_Control": "Control",
    "B1_HHL": "HHL",
    "C1_NIHL": "NIHL",
}
SUFFIXES = ["-1", "-2", "-3"]
OFFICIAL_COUNTS = {"-1": 17720, "-2": 10028, "-3": 11344}


def read_barcodes(path):
    with gzip.open(path, "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


def strip_suffix(barcode):
    stripped = re.sub(r"-\d+$", "", barcode)
    if stripped == barcode:
        raise RuntimeError("barcode has no GEM suffix: {}".format(barcode))
    return stripped


def base_set(barcodes, label):
    stripped = [strip_suffix(barcode) for barcode in barcodes]
    result = set(stripped)
    if len(result) != len(stripped):
        raise RuntimeError("duplicate suffix-stripped barcode in {}".format(label))
    return result


def parse_number(value):
    text = value.strip().replace(",", "")
    if text.endswith("%"):
        return float(text[:-1]) / 100.0
    try:
        return int(text)
    except ValueError:
        try:
            return float(text)
        except ValueError:
            return value


def main():
    os.makedirs(COMPARISON, exist_ok=True)
    official_barcodes = read_barcodes(
        os.path.join(PROJECT, "official", "GSE252285_barcodes.tsv.gz")
    )
    official_sets = {}
    for suffix in SUFFIXES:
        group = [barcode for barcode in official_barcodes if barcode.endswith(suffix)]
        if len(group) != OFFICIAL_COUNTS[suffix]:
            raise RuntimeError("official suffix count mismatch for {}".format(suffix))
        official_sets[suffix] = base_set(group, "GEO_{}".format(suffix))

    reconstruction_sets = {}
    metrics_by_library = {}
    for library in LIBRARIES:
        barcode_path = os.path.join(
            PROJECT,
            "cellranger",
            library,
            "outs",
            "filtered_feature_bc_matrix",
            "barcodes.tsv.gz",
        )
        reconstruction_sets[library] = base_set(read_barcodes(barcode_path), library)
        metrics_path = os.path.join(
            PROJECT, "cellranger", library, "outs", "metrics_summary.csv"
        )
        with open(metrics_path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        if len(rows) != 1:
            raise RuntimeError("unexpected metrics_summary rows for {}".format(library))
        metrics_by_library[library] = {
            key: parse_number(value) for key, value in rows[0].items()
        }

    overlap_rows = []
    overlap_lookup = {}
    for library in LIBRARIES:
        reconstructed = reconstruction_sets[library]
        for suffix in SUFFIXES:
            official = official_sets[suffix]
            intersection = len(reconstructed & official)
            union = len(reconstructed | official)
            row = {
                "library": library,
                "condition": CONDITIONS[library],
                "geo_suffix": suffix,
                "reconstruction_barcodes": len(reconstructed),
                "official_barcodes": len(official),
                "intersection": intersection,
                "union": union,
                "jaccard": intersection / float(union),
                "official_recall": intersection / float(len(official)),
                "reconstruction_recall": intersection / float(len(reconstructed)),
            }
            overlap_rows.append(row)
            overlap_lookup[(library, suffix)] = row

    long_path = os.path.join(COMPARISON, "barcode_overlap_3x3_long.tsv")
    long_fields = [
        "library",
        "condition",
        "geo_suffix",
        "reconstruction_barcodes",
        "official_barcodes",
        "intersection",
        "union",
        "jaccard",
        "official_recall",
        "reconstruction_recall",
    ]
    with open(long_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=long_fields)
        writer.writeheader()
        for row in overlap_rows:
            rendered = dict(row)
            for key in ("jaccard", "official_recall", "reconstruction_recall"):
                rendered[key] = "{:.8f}".format(row[key])
            writer.writerow(rendered)

    wide_path = os.path.join(COMPARISON, "barcode_overlap_3x3_matrix.tsv")
    with open(wide_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["library"] + ["GEO_{}".format(s) for s in SUFFIXES])
        for library in LIBRARIES:
            values = []
            for suffix in SUFFIXES:
                row = overlap_lookup[(library, suffix)]
                values.append(
                    "I={};J={:.6f};Orec={:.6f};Rrec={:.6f}".format(
                        row["intersection"],
                        row["jaccard"],
                        row["official_recall"],
                        row["reconstruction_recall"],
                    )
                )
            writer.writerow([library] + values)

    row_best = {
        library: max(
            SUFFIXES,
            key=lambda suffix: overlap_lookup[(library, suffix)]["intersection"],
        )
        for library in LIBRARIES
    }
    col_best = {
        suffix: max(
            LIBRARIES,
            key=lambda library: overlap_lookup[(library, suffix)]["intersection"],
        )
        for suffix in SUFFIXES
    }
    resolved = len(set(row_best.values())) == 3
    evidence = {}
    for library, suffix in row_best.items():
        matched = overlap_lookup[(library, suffix)]
        off_row = max(
            overlap_lookup[(library, other)]["intersection"]
            for other in SUFFIXES
            if other != suffix
        )
        off_col = max(
            overlap_lookup[(other, suffix)]["intersection"]
            for other in LIBRARIES
            if other != library
        )
        row_ratio = matched["intersection"] / float(max(off_row, 1))
        col_ratio = matched["intersection"] / float(max(off_col, 1))
        mutual_best = col_best[suffix] == library
        strong = (
            mutual_best
            and matched["intersection"] >= 1000
            and matched["jaccard"] >= 0.20
            and matched["official_recall"] >= 0.30
            and matched["reconstruction_recall"] >= 0.30
            and row_ratio >= 20.0
            and col_ratio >= 20.0
        )
        evidence[library] = {
            "suffix": suffix,
            "mutual_best": mutual_best,
            "matched_intersection": matched["intersection"],
            "matched_jaccard": matched["jaccard"],
            "matched_official_recall": matched["official_recall"],
            "matched_reconstruction_recall": matched["reconstruction_recall"],
            "row_intersection_ratio_to_next_best": row_ratio,
            "column_intersection_ratio_to_next_best": col_ratio,
            "passes_strong_unique_rule": strong,
        }
        resolved = resolved and strong

    suffix_to_condition = (
        {suffix: CONDITIONS[col_best[suffix]] for suffix in SUFFIXES} if resolved else {}
    )
    mapping = {
        "CONDITION_MAPPING_RESOLVED": "YES" if resolved else "NO",
        "mapping_confidence": "HIGH" if resolved else "UNRESOLVED",
        "suffix_to_condition": suffix_to_condition,
        "strong_unique_rule": {
            "one_to_one_mutual_best": True,
            "minimum_intersection": 1000,
            "minimum_jaccard": 0.20,
            "minimum_each_recall": 0.30,
            "minimum_row_and_column_intersection_ratio": 20.0,
        },
        "evidence": evidence,
    }
    with open(os.path.join(COMPARISON, "condition_mapping.json"), "w") as handle:
        json.dump(mapping, handle, indent=2, sort_keys=True)
        handle.write("\n")

    official_summary = {}
    with open(os.path.join(COMPARISON, "official_matrix_summary.tsv"), newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            official_summary[row["suffix"]] = row

    selected_metrics = [
        "Estimated Number of Cells",
        "Mean Reads per Cell",
        "Median Genes per Cell",
        "Median UMI Counts per Cell",
        "Fraction Reads in Cells",
        "Valid Barcodes",
        "Sequencing Saturation",
        "Reads Mapped Confidently to Transcriptome",
    ]
    comparison_rows = []
    for library in LIBRARIES:
        suffix = row_best[library] if resolved else "UNRESOLVED"
        metrics = metrics_by_library[library]
        estimated = metrics.get("Estimated Number of Cells")
        official_count = OFFICIAL_COUNTS[suffix] if resolved else "NA"
        comparison_rows.append(
            {
                "library": library,
                "condition": CONDITIONS[library],
                "matched_geo_suffix": suffix,
                "reconstructed_filtered_barcodes": len(reconstruction_sets[library]),
                "cellranger_estimated_cells": estimated,
                "official_processed_cells": official_count,
                "estimated_minus_official": estimated - official_count if resolved else "NA",
                "reconstruction_median_genes": metrics.get("Median Genes per Cell", "NA"),
                "official_matrix_median_genes": (
                    official_summary[suffix]["median_genes_per_cell"] if resolved else "NA"
                ),
                "reconstruction_median_umi": metrics.get("Median UMI Counts per Cell", "NA"),
                "official_matrix_median_umi": (
                    official_summary[suffix]["median_umi_per_cell"] if resolved else "NA"
                ),
                "reconstruction_fraction_reads_in_cells": metrics.get(
                    "Fraction Reads in Cells", "NA"
                ),
                "official_fraction_reads_in_cells": "NOT_AVAILABLE_FROM_GEO_MATRIX",
            }
        )
    comparison_fields = list(comparison_rows[0].keys())
    with open(os.path.join(COMPARISON, "metrics_comparison.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=comparison_fields)
        writer.writeheader()
        writer.writerows(comparison_rows)

    metrics_fields = ["library", "condition"] + selected_metrics
    with open(os.path.join(COMPARISON, "reconstruction_metrics.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=metrics_fields)
        writer.writeheader()
        for library in LIBRARIES:
            row = {"library": library, "condition": CONDITIONS[library]}
            row.update(
                {name: metrics_by_library[library].get(name, "NA") for name in selected_metrics}
            )
            writer.writerow(row)

    total_estimated = sum(
        metrics_by_library[library]["Estimated Number of Cells"] for library in LIBRARIES
    )
    total_difference = {
        "reconstructed_total_estimated_cells": total_estimated,
        "official_total_processed_cells": 39092,
        "difference": total_estimated - 39092,
        "percent_difference_vs_official": (total_estimated - 39092) / 39092.0 * 100.0,
    }
    with open(os.path.join(COMPARISON, "total_cell_difference.json"), "w") as handle:
        json.dump(total_difference, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(mapping, indent=2, sort_keys=True))
    print(json.dumps(total_difference, indent=2, sort_keys=True))
    print("BARCODE_AND_METRICS_COMPARISON_COMPLETE")


if __name__ == "__main__":
    main()
