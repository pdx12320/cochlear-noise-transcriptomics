#!/usr/bin/env python3
"""Sample every raw FASTQ to audit pairing, lengths, and v3 barcode structure."""

import csv
import gzip
import os
import statistics
import sys


PROJECT = os.environ.get(
    "PROJECT", "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
MANIFEST = os.path.join(PROJECT, "metadata", "fastq_manifest.tsv")
OUTPUT = os.path.join(PROJECT, "qc", "read_structure_audit.tsv")
SAMPLE_RECORDS = int(os.environ.get("READ_STRUCTURE_SAMPLE_RECORDS", "100000"))


def load_whitelist(path):
    if not path:
        return None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        return {line.strip() for line in handle if line.strip()}


def audit_fastq(path, read_name, whitelist):
    lengths = []
    records = 0
    exact = 0
    within_one = 0
    n_barcodes = 0
    bases = "ACGT"
    with gzip.open(path, "rt") as handle:
        while records < SAMPLE_RECORDS:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not header.startswith("@") or not plus.startswith("+"):
                raise RuntimeError("malformed FASTQ record in {}".format(path))
            if len(sequence) != len(quality):
                raise RuntimeError("sequence/quality length mismatch in {}".format(path))
            lengths.append(len(sequence))
            records += 1
            if read_name == "R1" and whitelist is not None and len(sequence) >= 16:
                barcode = sequence[:16]
                n_barcodes += 1
                if barcode in whitelist:
                    exact += 1
                    within_one += 1
                else:
                    found = False
                    for index, original in enumerate(barcode):
                        for base in bases:
                            if base != original and barcode[:index] + base + barcode[index + 1 :] in whitelist:
                                found = True
                                break
                        if found:
                            break
                    if found:
                        within_one += 1
    if not records:
        raise RuntimeError("no records sampled from {}".format(path))
    return {
        "sampled_records": records,
        "min_read_length": min(lengths),
        "median_read_length": statistics.median(lengths),
        "max_read_length": max(lengths),
        "v3_barcode_whitelist_exact_fraction": (
            "{:.8f}".format(exact / float(n_barcodes)) if n_barcodes else "NA"
        ),
        "v3_barcode_whitelist_within_1bp_fraction": (
            "{:.8f}".format(within_one / float(n_barcodes)) if n_barcodes else "NA"
        ),
    }


def main():
    whitelist_path = sys.argv[1] if len(sys.argv) > 1 else ""
    whitelist = load_whitelist(whitelist_path)
    with open(MANIFEST, newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for item in manifest:
        result = audit_fastq(item["archive_path"], item["read"], whitelist)
        result.update(
            {
                "library": item["library"],
                "srr": item["srr"],
                "read": item["read"],
                "fastq": item["archive_path"],
            }
        )
        rows.append(result)
        print("audited {}".format(item["archive_path"]), flush=True)

    fields = [
        "library",
        "srr",
        "read",
        "fastq",
        "sampled_records",
        "min_read_length",
        "median_read_length",
        "max_read_length",
        "v3_barcode_whitelist_exact_fraction",
        "v3_barcode_whitelist_within_1bp_fraction",
    ]
    with open(OUTPUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print("READ_STRUCTURE_AUDIT_COMPLETE {}".format(OUTPUT))


if __name__ == "__main__":
    main()
