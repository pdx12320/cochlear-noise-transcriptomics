#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${GSE252285_PROJECT_ROOT:-.}}"
SCRIPT_DIR="${SCRIPT_DIR:-${PROJECT_DIR}/03_scripts/no_doublet_robust_qc_tsne_sensitivity}"
R_BIN="${R_BIN:-${PROJECT_DIR}/work/r_env/bin/Rscript}"
ANNOTATION_REFERENCE_TSV="${ANNOTATION_REFERENCE_TSV:-${SCRIPT_DIR}/../../inputs/annotation_reference.tsv}"
RUN_LEGACY_COMPARISON="${RUN_LEGACY_COMPARISON:-false}"
BRANCH="no_doublet_robust_qc_tsne_sensitivity"
LOG_DIR="${PROJECT_DIR}/07_logs/${BRANCH}"
MANIFEST_DIR="${PROJECT_DIR}/00_manifest/${BRANCH}"
mkdir -p "${LOG_DIR}" "${MANIFEST_DIR}"

SOURCE_FILES=(
  "${PROJECT_DIR}/06_objects/01_raw_seurat.rds"
  "${ANNOTATION_REFERENCE_TSV}"
  "${PROJECT_DIR}/04_tables/QC_CELL_METADATA.tsv.gz"
  "${PROJECT_DIR}/04_tables/QC_THRESHOLDS.tsv"
)

sha256sum "${SOURCE_FILES[@]}" > "${MANIFEST_DIR}/SOURCE_SHA256_BEFORE.tsv"

export PROJECT_DIR SCRIPT_DIR ANNOTATION_REFERENCE_TSV RUN_LEGACY_COMPARISON SC_THREADS="${SC_THREADS:-4}"
"${R_BIN}" "${SCRIPT_DIR}/13a_no_doublet_cluster_tsne.R" 2>&1 | tee "${LOG_DIR}/13a_cluster_tsne.log"
"${R_BIN}" "${SCRIPT_DIR}/13b_no_doublet_markers_annotation_top20.R" 2>&1 | tee "${LOG_DIR}/13b_markers_annotation_top20.log"
"${R_BIN}" "${SCRIPT_DIR}/resolution_stability_audit.R" 2>&1 | tee "${LOG_DIR}/resolution_stability_audit.log"
"${R_BIN}" "${SCRIPT_DIR}/figure_qa_regenerate.R" 2>&1 | tee "${LOG_DIR}/figure_qa_regenerate.log"
"${R_BIN}" "${SCRIPT_DIR}/validate_no_doublet_branch.R" 2>&1 | tee "${LOG_DIR}/validation.log"

sha256sum "${SOURCE_FILES[@]}" > "${MANIFEST_DIR}/SOURCE_SHA256_AFTER.tsv"
cmp "${MANIFEST_DIR}/SOURCE_SHA256_BEFORE.tsv" "${MANIFEST_DIR}/SOURCE_SHA256_AFTER.tsv"

find \
  "${PROJECT_DIR}/00_manifest/${BRANCH}" \
  "${SCRIPT_DIR}" \
  "${PROJECT_DIR}/04_tables/${BRANCH}" \
  "${PROJECT_DIR}/05_figures/${BRANCH}" \
  "${PROJECT_DIR}/06_objects/${BRANCH}" \
  "${PROJECT_DIR}/07_logs/${BRANCH}" \
  "${PROJECT_DIR}/08_methods/${BRANCH}" \
  "${PROJECT_DIR}/09_qc/${BRANCH}" \
  "${PROJECT_DIR}/10_release/${BRANCH}" \
  -type f ! -name 'BUNDLE_MANIFEST_SHA256.tsv' -print0 | sort -z | xargs -0 sha256sum \
  > "${MANIFEST_DIR}/BUNDLE_MANIFEST_SHA256.tsv"

echo "SUCCESS no-doublet robust-QC t-SNE sensitivity branch"
