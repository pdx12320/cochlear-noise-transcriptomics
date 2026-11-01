#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import os
import re
import shutil
import sys
import tarfile
from collections import Counter, defaultdict


EXPECTED = {}
for gsm_num, tissue, condition, replicate in [
    *[(n, "SGN", "Naive", n - 5124290) for n in range(5124291, 5124295)],
    *[(n, "SGN", "Noise", n - 5124294) for n in range(5124295, 5124299)],
    *[(n, "Lateral_wall", "Naive", n - 5124298) for n in range(5124299, 5124303)],
    *[(n, "Lateral_wall", "Noise", n - 5124302) for n in range(5124303, 5124307)],
]:
    gsm = "GSM{}".format(gsm_num)
    EXPECTED[gsm] = {
        "sample_id": "{}_{}_{}".format(tissue, condition.lower(), replicate),
        "GSM": gsm,
        "tissue": tissue,
        "condition": condition,
        "replicate": replicate,
    }


def sha256_file(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def gzip_line_count(path):
    count = 0
    with gzip.open(path, "rt") as handle:
        for _ in handle:
            count += 1
    return count


def matrix_header(path):
    with gzip.open(path, "rt") as handle:
        banner = handle.readline().strip()
        if not banner.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError("Unexpected Matrix Market banner in {}: {}".format(path, banner))
        for line in handle:
            if not line.startswith("%"):
                fields = line.strip().split()
                if len(fields) != 3:
                    raise ValueError("Malformed Matrix Market dimensions in {}".format(path))
                return tuple(int(x) for x in fields)
    raise ValueError("Missing Matrix Market dimensions in {}".format(path))


def audit_features(path):
    ids = []
    symbols = []
    widths = Counter()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            widths[len(fields)] += 1
            ids.append(fields[0] if fields else "")
            symbols.append(fields[1] if len(fields) > 1 else "")
    return {
        "feature_rows": len(ids),
        "feature_field_widths": ",".join("{}:{}".format(k, widths[k]) for k in sorted(widths)),
        "duplicate_gene_ids": len(ids) - len(set(ids)),
        "empty_gene_ids": sum(not x for x in ids),
        "duplicate_gene_symbols": len(symbols) - len(set(symbols)),
        "empty_gene_symbols": sum(not x for x in symbols),
    }


def audit_barcodes(path):
    values = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            values.append(line.rstrip("\n"))
    return {
        "barcode_rows": len(values),
        "duplicate_barcodes": len(values) - len(set(values)),
        "empty_barcodes": sum(not x for x in values),
    }


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-tar", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()

    tar_path = os.path.abspath(args.input_tar)
    root = os.path.abspath(args.output_root)
    extract_dir = os.path.join(root, "work", "input_matrices")
    manifest_dir = os.path.join(root, "00_manifest")
    audit_dir = os.path.join(root, "01_input_audit")
    os.makedirs(extract_dir, exist_ok=False)

    pattern = re.compile(r"^(GSM\d+)_([^/]+)_(barcodes\.tsv\.gz|genes\.tsv\.gz|features\.tsv\.gz|matrix\.mtx\.gz)$")
    members_by_gsm = defaultdict(dict)
    extracted = []

    with tarfile.open(tar_path, "r") as archive:
        members = archive.getmembers()
        if len(members) != 48:
            raise RuntimeError("Expected exactly 48 tar members, found {}".format(len(members)))
        for member in members:
            if not member.isfile() or os.path.basename(member.name) != member.name:
                raise RuntimeError("Unsafe or non-file tar member: {}".format(member.name))
            match = pattern.match(member.name)
            if not match:
                raise RuntimeError("Unexpected filename: {}".format(member.name))
            gsm, descriptor, kind = match.groups()
            if gsm not in EXPECTED:
                raise RuntimeError("Unexpected GSM: {}".format(gsm))
            role = "barcode" if kind.startswith("barcodes") else "feature" if kind.startswith(("genes", "features")) else "matrix"
            if role in members_by_gsm[gsm]:
                raise RuntimeError("Duplicate {} file for {}".format(role, gsm))
            expected = EXPECTED[gsm]
            expected_tissue_token = "MouseCochleaNeurons" if expected["tissue"] == "SGN" else "MouseCochleaStria"
            expected_condition_token = "{}{}".format(expected["condition"].lower(), expected["replicate"])
            if descriptor != "{}_{}".format(expected_tissue_token, expected_condition_token):
                raise RuntimeError("Filename metadata conflict for {}: {}".format(gsm, descriptor))
            destination = os.path.join(extract_dir, member.name)
            source = archive.extractfile(member)
            if source is None:
                raise RuntimeError("Cannot read tar member {}".format(member.name))
            with source, open(destination, "wb") as out:
                shutil.copyfileobj(source, out)
            members_by_gsm[gsm][role] = destination
            extracted.append((gsm, role, destination, member.size))

    if set(members_by_gsm) != set(EXPECTED):
        raise RuntimeError("GSM set mismatch")
    if any(set(files) != {"matrix", "barcode", "feature"} for files in members_by_gsm.values()):
        raise RuntimeError("One or more GSMs lack matrix/barcode/feature triplets")

    file_rows = []
    manifest_rows = []
    feature_rows = []
    for gsm in sorted(EXPECTED, key=lambda x: int(x[3:])):
        meta = EXPECTED[gsm]
        files = members_by_gsm[gsm]
        hashes = {role: sha256_file(path) for role, path in files.items()}
        n_genes, n_cells, nnz = matrix_header(files["matrix"])
        barcode_audit = audit_barcodes(files["barcode"])
        feature_audit = audit_features(files["feature"])
        status = "PASS"
        reasons = []
        if barcode_audit["barcode_rows"] != n_cells:
            status = "FAIL"; reasons.append("barcode_count_mismatch")
        if feature_audit["feature_rows"] != n_genes:
            status = "FAIL"; reasons.append("feature_count_mismatch")
        if barcode_audit["duplicate_barcodes"] or barcode_audit["empty_barcodes"]:
            status = "FAIL"; reasons.append("invalid_barcodes")
        if feature_audit["duplicate_gene_ids"] or feature_audit["empty_gene_ids"]:
            status = "FAIL"; reasons.append("invalid_gene_ids")
        bundle = hashlib.sha256("".join(hashes[x] for x in ("matrix", "barcode", "feature")).encode()).hexdigest()
        manifest_rows.append({
            **meta,
            "matrix_file": files["matrix"],
            "barcode_file": files["barcode"],
            "feature_file": files["feature"],
            "n_genes_raw": n_genes,
            "n_cells_raw": n_cells,
            "checksum": bundle,
            "provenance_status": status if not reasons else status + ":" + ",".join(reasons),
        })
        feature_rows.append({
            "GSM": gsm,
            "n_genes_matrix": n_genes,
            "n_cells_matrix": n_cells,
            "nnz_matrix": nnz,
            **barcode_audit,
            **feature_audit,
            "status": status,
        })
        for role in ("matrix", "barcode", "feature"):
            path = files[role]
            file_rows.append({
                "GSM": gsm,
                "role": role,
                "file": path,
                "bytes": os.path.getsize(path),
                "sha256": hashes[role],
            })

    tar_sha = sha256_file(tar_path)
    write_tsv(os.path.join(manifest_dir, "S1_SAMPLE_MANIFEST.tsv"), manifest_rows, [
        "sample_id", "GSM", "tissue", "condition", "replicate", "matrix_file", "barcode_file",
        "feature_file", "n_genes_raw", "n_cells_raw", "checksum", "provenance_status",
    ])
    write_tsv(os.path.join(audit_dir, "S1_INPUT_FILES_SHA256.tsv"), file_rows, [
        "GSM", "role", "file", "bytes", "sha256",
    ])
    write_tsv(os.path.join(audit_dir, "S1_MATRIX_DIMENSION_AUDIT.tsv"), feature_rows, [
        "GSM", "n_genes_matrix", "n_cells_matrix", "nnz_matrix", "barcode_rows", "duplicate_barcodes",
        "empty_barcodes", "feature_rows", "feature_field_widths", "duplicate_gene_ids", "empty_gene_ids",
        "duplicate_gene_symbols", "empty_gene_symbols", "status",
    ])
    qa_rows = [
        {"check": "input_tar_sha256", "status": "PASS", "detail": tar_sha},
        {"check": "expected_gsm", "status": "PASS" if len(manifest_rows) == 16 else "FAIL", "detail": "{}/16".format(len(manifest_rows))},
        {"check": "expected_files", "status": "PASS" if len(file_rows) == 48 else "FAIL", "detail": "{}/48".format(len(file_rows))},
        {"check": "sample_provenance", "status": "PASS" if all(x["provenance_status"] == "PASS" for x in manifest_rows) else "FAIL", "detail": "matrix/barcode/feature and filename metadata concordance"},
    ]
    write_tsv(os.path.join(audit_dir, "PHASE1_QA.tsv"), qa_rows, ["check", "status", "detail"])
    if any(row["status"] == "FAIL" for row in qa_rows):
        return 2
    print("Phase 1 PASS: 16/16 GSM, 48/48 files, tar SHA256 {}".format(tar_sha))
    return 0


if __name__ == "__main__":
    sys.exit(main())
