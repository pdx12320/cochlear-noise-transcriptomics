#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

required=(README.md run_all.sh scripts tables figures methods logs environment manifests summaries)
for path in "${required[@]}"; do
  test -e "$path" || { echo "MISSING_REQUIRED=$path"; exit 1; }
done

for script in $(find scripts -type f -name '*.sh' -print); do
  bash -n "$script"
done

PYTHONDONTWRITEBYTECODE=1 python3 - <<'PY'
import csv
import hashlib
import os
import re
from collections import Counter
from pathlib import Path

root = Path.cwd()
manifest = root / "manifests" / "SHA256SUMS.tsv"
if not manifest.is_file():
    raise SystemExit("MISSING_SHA256_MANIFEST")

with manifest.open(encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
for row in rows:
    path = root / row["relative_path"]
    if not path.is_file():
        raise SystemExit(f"MANIFEST_MISSING_FILE={row['relative_path']}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != row["sha256"] or path.stat().st_size != int(row["bytes"]):
        raise SystemExit(f"MANIFEST_MISMATCH={row['relative_path']}")

forbidden_suffixes = (".rds", ".h5", ".h5ad", ".hdf5", ".loom", ".bam", ".cram", ".tar.gz")
for path in root.rglob("*"):
    if not path.is_file() or ".git" in path.parts or ".local" in path.parts:
        continue
    lower = path.name.lower()
    if lower.endswith(forbidden_suffixes):
        raise SystemExit(f"FORBIDDEN_FILE={path.relative_to(root)}")
    if path.stat().st_size >= 100_000_000:
        raise SystemExit(f"GITHUB_OVERSIZE_FILE={path.relative_to(root)}")

excluded_old = (
    "tables/stage1_annotation_v2", "tables/final_annotation_v3",
    "figures/stage1_annotation_v2", "figures/final_closure",
    "figures/focused_endothelial_tsne", "scripts/annotation_v2",
    "scripts/final_annotation_v3", "scripts/final_closure",
)
if any((root / path).exists() for path in excluded_old):
    raise SystemExit("SUPERSEDED_RESULT_DIRECTORY_PRESENT")

with (root / "tables/no_doublet_sensitivity/CELL_FLOW.tsv").open() as handle:
    flow = {row["step"]: int(row["n_cells"]) for row in csv.DictReader(handle, delimiter="\t")}
if flow != {"INITIAL_CELLS": 39092, "ROBUST_QC_RETAINED": 36111,
            "DOUBLETS_REMOVED": 0, "FINAL_CELLS": 36111}:
    raise SystemExit(f"CELL_FLOW_INVARIANT_FAILED={flow}")

with (root / "tables/no_doublet_sensitivity/CLUSTER_RESOLUTION_AUDIT.tsv").open() as handle:
    resolution = {float(row["resolution"]): int(row["number_of_clusters"])
                  for row in csv.DictReader(handle, delimiter="\t")}
if resolution != {0.4: 27, 0.6: 28, 0.8: 33, 1.0: 37}:
    raise SystemExit(f"RESOLUTION_INVARIANT_FAILED={resolution}")

with (root / "tables/no_doublet_sensitivity/CLUSTER_SIZES_RES0.6.tsv").open() as handle:
    cluster_sizes = list(csv.DictReader(handle, delimiter="\t"))
if len(cluster_sizes) != 28 or sum(int(row["cell_number"]) for row in cluster_sizes) != 36111:
    raise SystemExit("CLUSTER_SIZE_INVARIANT_FAILED")

with (root / "tables/no_doublet_sensitivity/BROAD_CELLTYPE_COUNTS.tsv").open() as handle:
    broad_counts = {row["broad_cell_type"]: int(row["cell_number"])
                    for row in csv.DictReader(handle, delimiter="\t")}
expected_broad = {
    "Immune_myeloid": 27676, "UNRESOLVED": 4726,
    "Fibrocyte_mesenchymal": 2997, "Endothelial_vascular": 275,
    "Glial": 264, "Supporting_epithelial": 149, "Hair_sensory": 24,
}
if broad_counts != expected_broad or sum(broad_counts.values()) != 36111:
    raise SystemExit(f"BROAD_COUNT_INVARIANT_FAILED={broad_counts}")

with (root / "tables/no_doublet_sensitivity/ALL_CLUSTERS_MARKER_CLASSIFICATION_REVISED.tsv").open() as handle:
    revised = list(csv.DictReader(handle, delimiter="\t"))
if len(revised) != 28 or {row["cluster"] for row in revised} != {str(i) for i in range(28)}:
    raise SystemExit("REVISED_MARKER_CLASSIFICATION_INVARIANT_FAILED")

with (root / "tables/condition_descriptive/CONDITION_TOTAL_CELLS.tsv").open() as handle:
    condition_rows = list(csv.DictReader(handle, delimiter="\t"))
condition_counts = {row["condition"]: int(row["total_cells"]) for row in condition_rows}
if condition_counts != {"Control": 10273, "HHL": 16631, "NIHL": 9207}:
    raise SystemExit("CONDITION_COUNT_INVARIANT_FAILED")

with (root / "manifests/source/GEO_CONDITION_MAPPING_AUDIT.tsv").open() as handle:
    mapping = list(csv.DictReader(handle, delimiter="\t"))
observed = {(row["geo_library_key"], row["library_id"]): row["condition_resolved"] for row in mapping}
expected_mapping = {("GEO_-3", "LIB3"): "Control", ("GEO_-1", "LIB1"): "HHL", ("GEO_-2", "LIB2"): "NIHL"}
if observed != expected_mapping or any(row["CONDITION_MAPPING_RESOLVED"] != "YES" for row in mapping):
    raise SystemExit("CONDITION_MAPPING_INVARIANT_FAILED")

with (root / "tables/condition_descriptive/TOP20_FEATURE_AUDIT_CONDITION.tsv").open() as handle:
    top20 = list(csv.DictReader(handle, delimiter="\t"))
if len(top20) != 20 or any(row["status"] != "FOUND" for row in top20):
    raise SystemExit("TOP20_PRESENCE_INVARIANT_FAILED")

figure_files = [path for path in (root / "figures").rglob("*") if path.suffix.lower() in {".pdf", ".png"}]
if not figure_files or any(path.stat().st_size == 0 for path in figure_files):
    raise SystemExit("FIGURE_NONEMPTY_CHECK_FAILED")

print(f"SHA256_VERIFIED={len(rows)}")
print(f"SELECTED_CELLS={flow['FINAL_CELLS']}")
print(f"CLUSTERS_RES0.6={resolution[0.6]}")
print(f"CONDITION_CELLS={sum(condition_counts.values())}")
print(f"FIGURES_NONEMPTY={len(figure_files)}")
PY

if rg -n --hidden -g '!.git/**' -g '!.local/**' -g '!scripts/utilities/verify_repository.sh' -g '!*.pdf' -g '!*.png' -g '!*.gz' \
  '/Users/[^/]+/|/data/[a-z][a-z0-9_-]+/|10\.[0-9]+\.[0-9]+\.[0-9]+|github_pat_[A-Za-z0-9_]+|ghp_[A-Za-z0-9]+|BEGIN (RSA|OPENSSH|EC|DSA) PRIVATE KEY' .; then
  echo "PRIVATE_PATH_OR_SECRET_PATTERN_FOUND"
  exit 1
fi

echo "REPOSITORY_VERIFICATION=PASS"
