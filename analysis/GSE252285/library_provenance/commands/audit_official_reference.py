#!/usr/bin/env python3
"""Audit the official GEO matrix and candidate 10x mouse references."""

import csv
import gzip
import json
import os
import re


PROJECT = os.environ.get(
    "PROJECT", "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
OFFICIAL = os.path.join(PROJECT, "official")
COMPARISON = os.path.join(PROJECT, "comparison")
REFERENCE_ROOT = os.path.join(PROJECT, "reference")
REFERENCE_NAMES = ["refdata-gex-GRCm39-2024-A", "refdata-gex-mm10-2020-A"]


def read_gzip_lines(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n\r") for line in handle]


def parse_gtf_attributes(text):
    values = {}
    for part in text.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        key, value = part.split(" ", 1)
        values[key] = value.strip().strip('"')
    return values


def reference_gene_pairs(gtf_path):
    pairs = []
    seen = set()
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            pair = (attrs.get("gene_id", ""), attrs.get("gene_name", ""))
            if pair not in seen:
                seen.add(pair)
                pairs.append(pair)
    return pairs


def main():
    os.makedirs(COMPARISON, exist_ok=True)
    barcode_path = os.path.join(OFFICIAL, "GSE252285_barcodes.tsv.gz")
    feature_path = os.path.join(OFFICIAL, "GSE252285_features.tsv.gz")
    matrix_path = os.path.join(OFFICIAL, "GSE252285_matrix.mtx.gz")

    barcodes = read_gzip_lines(barcode_path)
    feature_rows = [line.split("\t") for line in read_gzip_lines(feature_path)]
    feature_pairs = [(row[0], row[1]) for row in feature_rows]
    if any(len(row) < 3 or row[2] != "Gene Expression" for row in feature_rows):
        raise RuntimeError("official feature schema is not 10x Gene Expression")

    suffix_indices = {"-1": [], "-2": [], "-3": []}
    for index, barcode in enumerate(barcodes):
        match = re.search(r"(-[123])$", barcode)
        if not match:
            raise RuntimeError("unexpected official barcode suffix: {}".format(barcode))
        suffix_indices[match.group(1)].append(index)

    expected_suffix_counts = {"-1": 17720, "-2": 10028, "-3": 11344}
    observed_suffix_counts = {key: len(value) for key, value in suffix_indices.items()}
    if observed_suffix_counts != expected_suffix_counts:
        raise RuntimeError(
            "official suffix counts differ: {}".format(observed_suffix_counts)
        )

    with gzip.open(matrix_path, "rt") as handle:
        dimensions = None
        official_matrix_metadata = {}
        for line in handle:
            if line.startswith("%metadata_json:"):
                official_matrix_metadata = json.loads(line.split(":", 1)[1].strip())
                continue
            if line.startswith("%"):
                continue
            dimensions = tuple(int(value) for value in line.split())
            break
    if dimensions != (len(feature_rows), len(barcodes), 89879619):
        raise RuntimeError("matrix dimensions disagree with validated inputs")

    summary_path = os.path.join(COMPARISON, "official_matrix_summary.tsv")
    if not os.path.isfile(summary_path):
        raise RuntimeError("official matrix summary must be generated before reference audit")
    with open(summary_path, newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    summary_counts = {
        row["suffix"]: int(row["cell_count"])
        for row in summary_rows
        if row["suffix"] in expected_suffix_counts
    }
    if summary_counts != expected_suffix_counts:
        raise RuntimeError("R matrix summary suffix counts failed validation")

    official_set = set(feature_pairs)
    candidates = {}
    for reference_name in REFERENCE_NAMES:
        reference_path = os.path.join(REFERENCE_ROOT, reference_name)
        gtf_paths = [
            os.path.join(reference_path, "genes", "genes.gtf.gz"),
            os.path.join(reference_path, "genes", "genes.gtf"),
        ]
        gtf_path = next((path for path in gtf_paths if os.path.isfile(path)), None)
        if gtf_path is None:
            continue
        candidate_pairs = reference_gene_pairs(gtf_path)
        candidate_set = set(candidate_pairs)
        candidate_audit = {
            "candidate_reference_gene_count": len(candidate_pairs),
            "candidate_unique_gene_pairs": len(candidate_set),
            "pair_intersection": len(official_set & candidate_set),
            "official_only": len(official_set - candidate_set),
            "candidate_only": len(candidate_set - official_set),
            "exact_pair_set_match": official_set == candidate_set,
            "exact_pair_order_match": feature_pairs == candidate_pairs,
        }
        for metadata_name in ("reference.json", "version"):
            metadata_path = os.path.join(reference_path, metadata_name)
            if os.path.isfile(metadata_path):
                with open(metadata_path) as handle:
                    candidate_audit[metadata_name] = handle.read().strip()
        candidates[reference_name] = candidate_audit

    exact_candidates = [
        name for name, audit in candidates.items() if audit["exact_pair_set_match"]
    ]
    if len(exact_candidates) != 1:
        raise RuntimeError(
            "reference is not uniquely identified by feature pairs: {}".format(
                exact_candidates
            )
        )
    selected_reference = exact_candidates[0]
    selected = candidates[selected_reference]
    reference_audit = {
        "official_feature_count": len(feature_pairs),
        "official_unique_feature_pairs": len(official_set),
        "selected_reference": selected_reference,
        "candidate_reference_gene_count": selected["candidate_reference_gene_count"],
        "candidate_unique_gene_pairs": selected["candidate_unique_gene_pairs"],
        "pair_intersection": selected["pair_intersection"],
        "official_only": selected["official_only"],
        "candidate_only": selected["candidate_only"],
        "exact_pair_set_match": selected["exact_pair_set_match"],
        "exact_pair_order_match": selected["exact_pair_order_match"],
        "candidates": candidates,
        "matrix_rows": dimensions[0],
        "matrix_columns": dimensions[1],
        "matrix_nonzero_entries": dimensions[2],
        "official_matrix_metadata": official_matrix_metadata,
        "suffix_counts": observed_suffix_counts,
    }
    audit_path = os.path.join(COMPARISON, "reference_feature_audit.json")
    with open(audit_path, "w") as handle:
        json.dump(reference_audit, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(reference_audit, indent=2, sort_keys=True))
    print("OFFICIAL_MATRIX_AUDIT_COMPLETE {}".format(summary_path))


if __name__ == "__main__":
    main()
