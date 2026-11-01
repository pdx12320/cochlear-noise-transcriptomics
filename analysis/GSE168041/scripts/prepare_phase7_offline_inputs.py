#!/usr/bin/env python3
"""Stage exact archived GSE168041 DE tables for the offline phase7 script.

Usage: python prepare_phase7_offline_inputs.py SOURCE_GSE168041 NEW_STAGED_ROOT
No statistical analysis is performed. Existing output paths are never reused.
"""
from pathlib import Path
import csv
import hashlib
import re
import shutil
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
source = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2]).resolve()
if target.exists():
    raise SystemExit("Use a new staged root; refusing to overwrite an existing path")

summary_path = source / "09_differential/S1_DS_SUMMARY.tsv"
with summary_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    columns = list(reader.fieldnames or [])
    rows = list(reader)
if len(rows) != 13 or any(row["status"] != "PASS" for row in rows):
    raise SystemExit("Expected exactly thirteen eligible archived models")

planned = []
for row in rows:
    stem = re.sub(r"[^A-Za-z0-9]+", "_", row["cell_type"]).strip("_")
    relative = Path("differential") / row["tissue"] / f"S1_DS_ALL_{stem}.tsv"
    path = source / relative
    if not path.is_file():
        raise SystemExit(f"Missing required differential table: {relative}")
    with path.open() as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if not {"gene_id", "gene_symbol", "stat"}.issubset(header):
        raise SystemExit(f"Required ranking columns are absent: {relative}")
    planned.append((row, relative, path))
if len({str(relative) for _, relative, _ in planned}) != 13:
    raise SystemExit("Differential table mapping is not one-to-one")

target.mkdir(parents=True)
for name in ["09_differential", "12_figures", "13_source_data", "14_methods"]:
    (target / name).mkdir()
manifest = []
for row, relative, path in planned:
    copied = target / relative
    copied.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, copied)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if hashlib.sha256(copied.read_bytes()).hexdigest() != digest:
        raise RuntimeError("Input-copy hash mismatch")
    row["result_file"] = str(copied)
    manifest.append({"relative_input": str(relative), "sha256": digest})
if "result_file" not in columns:
    columns.append("result_file")
with (target / "09_differential/S1_DS_SUMMARY.tsv").open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
with (target / "14_methods/STAGED_INPUT_SHA256.tsv").open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=["relative_input", "sha256"], delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest)
print(f"Staged thirteen unchanged differential tables in {target}; no enrichment run")
