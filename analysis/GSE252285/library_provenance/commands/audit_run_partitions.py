#!/usr/bin/env python3
"""Use 10x barcode evidence to test whether SRRs are technical partitions."""

import collections
import csv
import gzip
import json
import os


PROJECT = os.environ.get(
    "PROJECT", "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
MANIFEST = os.path.join(PROJECT, "metadata", "fastq_manifest.tsv")
WHITELIST = os.path.join(PROJECT, "metadata", "3M-february-2018.txt.gz")
OUTPUT = os.path.join(PROJECT, "qc", "run_partition_barcode_overlap.tsv")
STATUS = os.path.join(PROJECT, "qc", "run_partition_status.json")
SAMPLE_RECORDS = int(os.environ.get("PARTITION_AUDIT_RECORDS", "2000000"))
TOP_BARCODES = int(os.environ.get("PARTITION_AUDIT_TOP_BARCODES", "15000"))


def load_whitelist():
    with gzip.open(WHITELIST, "rt") as handle:
        return {line.strip() for line in handle if line.strip()}


def top_barcode_set(path, whitelist):
    counts = collections.Counter()
    sampled = 0
    exact = 0
    with gzip.open(path, "rt") as handle:
        while sampled < SAMPLE_RECORDS:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline()
            if not header.startswith("@") or not plus.startswith("+") or not quality:
                raise RuntimeError("malformed FASTQ: {}".format(path))
            barcode = sequence[:16]
            sampled += 1
            if barcode in whitelist:
                counts[barcode] += 1
                exact += 1
    top = {barcode for barcode, count in counts.most_common(TOP_BARCODES)}
    return top, sampled, exact


def main():
    whitelist = load_whitelist()
    with open(MANIFEST, newline="") as handle:
        all_rows = list(csv.DictReader(handle, delimiter="\t"))
    r1_rows = sorted(
        [row for row in all_rows if row["read"] == "R1"],
        key=lambda row: row["srr"],
    )
    if len(r1_rows) != 6:
        raise RuntimeError("expected six R1 files")

    sets = {}
    sample_stats = {}
    srx = {}
    library = {}
    for row in r1_rows:
        print("sampling barcode evidence from {}".format(row["archive_path"]), flush=True)
        top, sampled, exact = top_barcode_set(row["archive_path"], whitelist)
        sets[row["srr"]] = top
        srx[row["srr"]] = row["srx"]
        library[row["srr"]] = row["library"]
        sample_stats[row["srr"]] = {
            "sampled_records": sampled,
            "exact_whitelist_records": exact,
            "exact_whitelist_fraction": exact / float(sampled),
            "top_barcode_count": len(top),
        }

    pair_rows = []
    within_hhl = []
    between_srx = []
    srrs = sorted(sets)
    for left_index, left in enumerate(srrs):
        for right in srrs[left_index + 1 :]:
            intersection = len(sets[left] & sets[right])
            union = len(sets[left] | sets[right])
            same_srx = srx[left] == srx[right]
            row = {
                "srr_1": left,
                "srr_2": right,
                "library_1": library[left],
                "library_2": library[right],
                "srx_1": srx[left],
                "srx_2": srx[right],
                "same_srx": "YES" if same_srx else "NO",
                "top_barcode_intersection": intersection,
                "top_barcode_union": union,
                "jaccard": intersection / float(union),
            }
            pair_rows.append(row)
            if same_srx and srx[left] == "SRX23063846":
                within_hhl.append(row["jaccard"])
            elif not same_srx:
                between_srx.append(row["jaccard"])

    fields = list(pair_rows[0].keys())
    with open(OUTPUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in pair_rows:
            rendered = dict(row)
            rendered["jaccard"] = "{:.8f}".format(row["jaccard"])
            writer.writerow(rendered)

    min_hhl = min(within_hhl)
    max_between = max(between_srx)
    unique_ratio = min_hhl / float(max(max_between, 1e-12))
    confirmed = min_hhl >= 0.50 and unique_ratio >= 20.0
    result = {
        "HHL_TECHNICAL_PARTITIONS_CONFIRMED": "YES" if confirmed else "NO",
        "method": "top {} exact-v3-whitelist barcodes from first {} R1 records per SRR".format(
            TOP_BARCODES, SAMPLE_RECORDS
        ),
        "hhl_srx": "SRX23063846",
        "hhl_runs": [
            "SRR27387602",
            "SRR27387603",
            "SRR27387604",
            "SRR27387605",
        ],
        "minimum_within_hhl_pairwise_jaccard": min_hhl,
        "maximum_between_srx_pairwise_jaccard": max_between,
        "minimum_within_to_maximum_between_ratio": unique_ratio,
        "strong_rule": "min within-HHL Jaccard >= 0.50 and >=20x max between-SRX Jaccard",
        "sample_stats": sample_stats,
    }
    with open(STATUS, "w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not confirmed:
        raise RuntimeError("HHL SRRs are not strongly supported as technical partitions")


if __name__ == "__main__":
    main()
