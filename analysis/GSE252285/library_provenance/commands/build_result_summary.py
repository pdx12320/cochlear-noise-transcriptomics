#!/usr/bin/env python3
"""Build the concise final result requested by the user."""

import csv
import json
import os


PROJECT = os.environ.get(
    "PROJECT", "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
COMPARISON = os.path.join(PROJECT, "comparison")


def load_json(name):
    with open(os.path.join(COMPARISON, name)) as handle:
        return json.load(handle)


def load_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    provenance = load_tsv(
        os.path.join(PROJECT, "metadata", "gsm_srr_library_provenance.tsv")
    )
    manifest = load_tsv(os.path.join(PROJECT, "metadata", "fastq_manifest.tsv"))
    metrics = load_tsv(os.path.join(COMPARISON, "metrics_comparison.tsv"))
    overlaps = load_tsv(os.path.join(COMPARISON, "barcode_overlap_3x3_long.tsv"))
    mapping = load_json("condition_mapping.json")
    difference = load_json("total_cell_difference.json")
    reference = load_json("reference_feature_audit.json")

    if mapping["CONDITION_MAPPING_RESOLVED"] != "YES":
        raise RuntimeError("mapping is unresolved; SUCCESS cannot be emitted")
    if len(manifest) != 12 or any(row["checksum_status"] != "PASS" for row in manifest):
        raise RuntimeError("FASTQ manifest is incomplete or has failed checksums")

    required = []
    for library in ("A1_Control", "B1_HHL", "C1_NIHL"):
        base = os.path.join(PROJECT, "cellranger", library, "outs")
        required.extend(
            [
                os.path.join(base, "web_summary.html"),
                os.path.join(base, "metrics_summary.csv"),
                os.path.join(base, "raw_feature_bc_matrix", "barcodes.tsv.gz"),
                os.path.join(base, "raw_feature_bc_matrix", "features.tsv.gz"),
                os.path.join(base, "raw_feature_bc_matrix", "matrix.mtx.gz"),
                os.path.join(base, "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
                os.path.join(base, "filtered_feature_bc_matrix", "features.tsv.gz"),
                os.path.join(base, "filtered_feature_bc_matrix", "matrix.mtx.gz"),
            ]
        )
    missing = [path for path in required if not os.path.isfile(path) or not os.path.getsize(path)]
    if missing:
        raise RuntimeError("required Cell Ranger outputs missing: {}".format(missing))

    by_gsm = {}
    for row in provenance:
        by_gsm.setdefault(row["gsm"], []).append(row)
    total_bytes = sum(int(row["actual_bytes"]) for row in manifest)

    lines = ["SUCCESS", "", "GSM-SRR-library mapping:"]
    for gsm in ("GSM7998628", "GSM7998629", "GSM7998630"):
        rows = sorted(by_gsm[gsm], key=lambda item: item["srr"])
        lines.append(
            "{} -> {} -> {} ({}) [{}]".format(
                gsm,
                ",".join(row["srr"] for row in rows),
                rows[0]["library"],
                rows[0]["condition"],
                "one 10x library; {} technical partition(s); shared {} / {}".format(
                    len(rows), rows[0]["srx"], rows[0]["biosample"]
                ),
            )
        )

    lines.extend(
        [
            "",
            "FASTQ: {} files; {:.3f} GiB ({:.3f} GB); ENA MD5 {}/{} PASS".format(
                len(manifest),
                total_bytes / (1024.0 ** 3),
                total_bytes / 1e9,
                sum(row["checksum_status"] == "PASS" for row in manifest),
                len(manifest),
            ),
            "Cell Ranger: 7.1.0",
            "Reference: {}; official/reference feature-pair set exact_match={} (intersection={}, official_only={}, reference_only={})".format(
                reference["selected_reference"],
                reference["exact_pair_set_match"],
                reference["pair_intersection"],
                reference["official_only"],
                reference["candidate_only"],
            ),
            "",
            "Estimated cells:",
        ]
    )
    for row in metrics:
        lines.append(
            "{}: {} (matched GEO {} official {}; difference {:+d})".format(
                row["library"],
                row["cellranger_estimated_cells"],
                row["matched_geo_suffix"],
                row["official_processed_cells"],
                int(row["estimated_minus_official"]),
            )
        )

    lookup = {(row["library"], row["geo_suffix"]): row for row in overlaps}
    lines.extend(
        [
            "",
            "3x3 barcode overlap (I/Jaccard/official_recall/reconstruction_recall):",
            "library\tGEO_-1\tGEO_-2\tGEO_-3",
        ]
    )
    for library in ("A1_Control", "B1_HHL", "C1_NIHL"):
        cells = []
        for suffix in ("-1", "-2", "-3"):
            row = lookup[(library, suffix)]
            cells.append(
                "{}/{:.6f}/{:.6f}/{:.6f}".format(
                    int(row["intersection"]),
                    float(row["jaccard"]),
                    float(row["official_recall"]),
                    float(row["reconstruction_recall"]),
                )
            )
        lines.append("{}\t{}".format(library, "\t".join(cells)))

    lines.extend(["", "Suffix -> condition mapping:"])
    for suffix in ("-1", "-2", "-3"):
        lines.append("{} = {}".format(suffix, mapping["suffix_to_condition"][suffix]))
    lines.append("CONDITION_MAPPING_RESOLVED=YES")
    ratio_summary = "; ".join(
        "{} row_ratio={:.1f}x col_ratio={:.1f}x".format(
            library,
            mapping["evidence"][library]["row_intersection_ratio_to_next_best"],
            mapping["evidence"][library]["column_intersection_ratio_to_next_best"],
        )
        for library in ("A1_Control", "B1_HHL", "C1_NIHL")
    )
    lines.append("Mapping confidence: HIGH ({})".format(ratio_summary))
    lines.extend(
        [
            "",
            "Official 39092-cell difference: reconstructed total {} cells; difference {:+d} ({:+.3f}%).".format(
                difference["reconstructed_total_estimated_cells"],
                difference["difference"],
                difference["percent_difference_vs_official"],
            ),
            "Result directory: {}".format(PROJECT),
        ]
    )

    output = os.path.join(PROJECT, "RESULT_SUMMARY.txt")
    with open(output, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
